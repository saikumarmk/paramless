using Statistics
using SpecialFunctions
using KahanSummation
include("paramless.jl")

function gaussian_j(q::Float64, q_threshold::Float64, ϵ::Float64)::Float64

    #edge cases:
    if q == q_threshold && ϵ == 0.0
        return 0.5
    end
    result = 1.0 - ((1 + erf((q_threshold - q) / ϵ / sqrt(2))) / 2)
    if isnan(result)
        throw(error("Error computing gaussian_j: q = $q, q_threshold = $q_threshold, ϵ = $ϵ"))
    end
    return Float64(result)
end


function payoff_scientist(q::Float64, q_threshold::Float64, c::Float64, ϵ::Float64, q_bar::Float64, version::Int)::Float64
    #= 
    q = quality of scientist
    q_threshold = used to determine probability of acceptance 
    ϵ = noise generated by the journal 
    q_bar = average quality of the paper, has a role in the benefit  
    =#
    p = gaussian_j(q, q_threshold, ϵ)

    if version == 3
        return max(p * q_bar - (1.0 - p) * c, 0.0)
    elseif version == 4
        return max(p * q_bar - (1.0 - p) * q * c, 0.0)
    else
        throw(ErrorException("Invalid Version, version ∈ {3,4}"))
    end
end

"""
    compute_scientist_payoff_curve(domain, submission_curve, q_threshold, c, ϵ, q_bar)

For a submission curve over [0,1], determine the collective payoff.
"""
function compute_scientist_payoff_curve(domain::Array{Float64}, submission_curve::Array{Float64}, q_threshold::Float64, c::Float64, ϵ::Float64, q_bar::Float64, version::Int)::Array{Float64}

    return [submission_curve[idx] * payoff_scientist(domain[idx], q_threshold, 
                                                        c, ϵ, 
                                                        q_bar, version) for idx = 1:length(domain)]
end


"""

Compute the fitness of a submission set and its mutant.

"""
function science_fitness(resident_submission::Array{Float64}, mutant_submission::Array{Float64}; kwargs...)
    if haskey(kwargs, :q_threshold)
        q_threshold = kwargs[:q_threshold]
    end

    if haskey(kwargs, :ϵ)
        ϵ = kwargs[:ϵ]
    end

    if haskey(kwargs, :c)
        c = kwargs[:c]
    end

    if haskey(kwargs, :domain)
        domain = kwargs[:domain]
    end

    if haskey(kwargs, :q_bar)
        q_bar = kwargs[:q_bar]
    end

    if haskey(kwargs, :version)
        version = kwargs[:version]
    else
        version = 3
    end

    fitness_resident = sum(compute_scientist_payoff_curve(domain, resident_submission, q_threshold, c, ϵ, q_bar, version))
    fitness_mutant = sum(compute_scientist_payoff_curve(domain, mutant_submission, q_threshold, c, ϵ, q_bar, version))

    return fitness_resident, fitness_mutant
end

"""
    binary_mutation(vector::Array; kwargs...)
Binary mutation in which a random submission is flipped.

"""
function binary_mutation(vector::Array; kwargs...)

    mutant = deepcopy(vector)

    index = rand(1:length(vector))

    mutant[index] = 1 - mutant[index]
    return mutant

end

"""
    binary_mutations(vector::Array; kwargs...)

Experimental mutation function that flips multiple scientists.

"""
function binary_mutations(vector::Array; kwargs...)
    mutant = deepcopy(vector)

    indices = rand(1:length(vector), 5)

    for idx = indices
        mutant[idx] = 1 - mutant[idx]
    end
    return mutant
end


"""
    average_quality_accepted(q_threshold::Float64, ϵ::Float64, bit_vector::Array{Float64}, domain::Array{Float64})::Float64

Computes q_bar for a submission curve. This is given by the discretisation of the integral of the centroid.    
"""
function average_quality_accepted(q_threshold::Float64, ϵ::Float64, bit_vector::Array{Float64}, domain::Array{Float64})::Float64
    #=
    Produces an area curve, then computes noise 
    (q1*e1+q2*e2+..)/(e1+e2+..)
    Computes q_bar of new batch
    =#
    numerator_vector = Float64[]
    denominator_vector = Float64[]
    for idx = 1:length(domain) 
        noise = gaussian_j(domain[idx], q_threshold, ϵ) 
        push!(numerator_vector, noise * bit_vector[idx] * domain[idx])
        push!(denominator_vector, noise * bit_vector[idx])
    end

    numerator = sum_kbn(numerator_vector) 
    denominator = sum_kbn(denominator_vector)

    if denominator == 0.0
        return 0.0
    else
        result = numerator / denominator
    end
    return result
end

"""

Computes the acceptance rate.

"""
function acceptance_rate(q_threshold::Float64, ϵ::Float64, submission_set::Array{Float64}, domain::Array{Float64})::Float64
    numerator_vector = Float64[]
    for idx = 1:length(submission_set) 
        push!(numerator_vector, submission_set[idx] * gaussian_j(domain[idx], q_threshold, ϵ))
    end
    numerator = sum_kbn(numerator_vector)
    denominator = sum(submission_set) # the rough length of accepted 
    if denominator == 0.0
        return 1.0
    else
        return numerator / denominator
    end
end

"""

Computes rejection rate.

"""
function rejection_rate(q_threshold::Float64, ϵ::Float64, submission_set::Array{Float64}, domain::Array{Float64})::Float64
    return 1.0 - acceptance_rate(q_threshold, ϵ, submission_set, domain)
end

"""

Computes π_journal = benefit-cost, with cost = 1/(1+ϵ)^k and benefit either q_bar or rejection rate.

cost not
"""
function payoff_journal(q_threshold::Float64, ϵ::Float64, k::Float64, submission::Array, domain::Array, quality::Bool)::Float64

    if quality
        benefit = average_quality_accepted(q_threshold, ϵ, submission, domain) # Has a version that depends on a vector, submission, probability, domain
    else
        benefit = rejection_rate(q_threshold, ϵ, submission, domain)
    end

    cost = 1 / (1 + ϵ)^k
    result = benefit - cost

    return result

end

""" 

Checks fitness by recomputing submission curve under new parameters.

"""
function journal_fitness(resident::Float64, mutant::Float64; kwargs...)

    if haskey(kwargs, :mutate_qt)
        mutate_qt = kwargs[:mutate_qt]
    end

    if haskey(kwargs, :ϵ)
        ϵ = kwargs[:ϵ]
    end

    if haskey(kwargs, :q_threshold)
        q_threshold = kwargs[:q_threshold]
    end

    if haskey(kwargs, :k)
        k = kwargs[:k]
    end

    if haskey(kwargs, :submission)
        submission = kwargs[:submission]
    end

    if haskey(kwargs, :domain)
        domain = kwargs[:domain]
    end

    if haskey(kwargs, :quality)
        quality = kwargs[:quality]
    else
        quality = true
    end

    if mutate_qt
        fitness_resident = payoff_journal(resident, ϵ, k, submission, domain, quality)
        fitness_mutant = payoff_journal(mutant, ϵ, k, submission, domain, quality)
    else
        fitness_resident = payoff_journal(q_threshold, resident, k, submission, domain, quality)
        fitness_mutant = payoff_journal(q_threshold, mutant, k, submission, domain, quality)
    end

    return fitness_resident, fitness_mutant

end

"""

Mutates a scalar quantity for the journal. 

"""

function journal_mutation(param::Float64; kwargs...)::Float64
    is_inside = false
    attempt = 0
    mutant = param

    if haskey(kwargs, :mutation_epsilon_journal)
        mutation_epsilon = kwargs[:mutation_epsilon_journal]
    end
    if haskey(kwargs, :lower_bound)
        lower_bound = kwargs[:lower_bound]
    else
        lower_bound = nothing
    end

    if haskey(kwargs, :upper_bound)
        upper_bound = kwargs[:upper_bound]
    else
        upper_bound = nothing
    end

    while !is_inside
        mutant = _attempt_journal_mutation(param, mutation_epsilon) # Modify their DNA
        is_inside = _is_within_bounds_scalar(mutant, lower_bound, upper_bound) # Check if everything is contained
        attempt += 1

        if attempt > MAX_ITER
            throw(ErrorException("Attempted too many mutations without producing anythin within bounds"))
        end
    end
    return mutant
end


function _attempt_journal_mutation(param::Float64, mutation_epsilon::Float64)::Float64
    return param + rand([-1, 1], 1)[1] * mutation_epsilon
end

""" 

Evolves a submission set in conjunction with a journal with specified parameters.
Need to specify:
k ∈ (0, ∞)
c ∈ [0,1]
mutation_epsilon ∈ (0,1)
lower_bound = 0.0
q_bar = 0.4

"""
function co_evolve(domain::Array{Float64}, resident_submission::Array{Float64}, resident_ϵ::Float64, resident_q_threshold::Float64,
    iterations, atol::Float64, seed = 0, snapshot_iter = 1e3; kwargs...)

    Random.seed!(seed)


    time_series_submission = [[0, deepcopy(resident_submission)]]
    time_series_ϵ = [[0, resident_ϵ]]
    time_series_q_threshold = [[0, resident_q_threshold]]
    previous_submission = zeros(size(resident_submission))
    previous_q_threshold = 0
    previous_ϵ = 0


    for step = 1:iterations

        q_bar = average_quality_accepted(resident_q_threshold, resident_ϵ, resident_submission, domain)

        # Compute and mutate science fitness at the same time,
        previous_submission = resident_submission
        previous_q_threshold = resident_q_threshold
        previous_ϵ = resident_ϵ
        # is q_bar passed to science? This is a big issue
        resident_submission, invasion = evolution_step(resident_submission, science_fitness,
            gaussian_mutation, atol;
            kwargs..., domain=domain, q_threshold=previous_q_threshold, ϵ=previous_ϵ, q_bar=q_bar)



        # Compute for ϵ
        resident_ϵ, invasion = evolution_step(resident_ϵ, journal_fitness,
            journal_mutation, atol;
            kwargs..., submission=previous_submission,
            q_threshold=previous_q_threshold, mutate_qt=false, domain=domain)

        # compute for threshold
        resident_q_threshold, invasion = evolution_step(resident_q_threshold, journal_fitness,
            journal_mutation, atol;
            kwargs..., submission=previous_submission,
            ϵ=previous_ϵ, mutate_qt=true, domain=domain)

        if step % snapshot_iter == 0
            push!(time_series_submission, [step, resident_submission])
            push!(time_series_ϵ, [step, resident_ϵ])
            push!(time_series_q_threshold, [step, resident_q_threshold])
        end


    end

    return Dict("resident_submission" => resident_submission,
        "resident_ϵ" => resident_ϵ,
        "resident_q_threshold" => resident_q_threshold,
        "time_series_submission" => time_series_submission,
        "time_series_ϵ" => time_series_ϵ,
        "time_series_q_threshold" => time_series_q_threshold
    )

end



function main()

    domain = Array(range(0.0,1.0, length=100));
    resident_sub = 0.01*ones(size(domain)); # 0.0, not 0, everyone submits with pr 0
    resident_eps = 0.0
    resident_qt = 0.0

    result_dict = co_evolve(domain, resident_sub, resident_eps, resident_qt,
    1e5, 1e-9, 0, 1e3; k=8.0, 
    c= 0.6, mutation_epsilon=0.001,mutation_epsilon_journal=0.001, 
    lower_bound=0.0, upper_bound=1.0, 
    width=0.01, quality=false, version=3);
    println(result_dict)
end

