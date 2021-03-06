using Statistics
using SpecialFunctions
using KahanSummation
using NumericalIntegration
include("paramless.jl")

"""
    payoff_scientist(q::Float64, c::Float64, q_bar::Float64, confusion::Array{Float64}, index::Int, version::Int)::Float64

Computes the payoff of a scientist given an evolving submsision function.
"""
function payoff_scientist(q::Float64, c::Float64, q_bar::Float64, confusion::Array{Float64}, index::Int, version::Int)::Float64
    p = confusion[index]

    if version == 3
        return max(p * q_bar - (1.0 - p) * c,0)
    elseif version == 4
        return max(p * q_bar - (1.0 - p) * q * c, 0)
    else
        throw(ErrorException("Invalid Version, version ∈ {3,4}"))
    end
end

"""
    compute_scientist_payoff_curve(domain, submission_curve, q_threshold, c, ϵ, q_bar)

For a submission curve over [0,1], determine the collective payoff.
"""
function compute_scientist_payoff_curve(domain::Array{Float64}, submission_curve::Array{Float64}, c::Float64, confusion::Array{Float64}, q_bar::Float64, version::Int)::Array{Float64}

    return [submission_curve[idx] * payoff_scientist(domain[idx], c, q_bar, confusion, idx, version) for idx = 1:length(domain)]
end


"""
    science_fitness(resident_submission::Array{Float64}, mutant_submission::Array{Float64}; kwargs...)

Compute the fitness of a submission set and its mutant.

"""
function science_fitness(resident_submission::Array{Float64}, mutant_submission::Array{Float64}; kwargs...)
    if haskey(kwargs, :confusion)
        confusion = kwargs[:confusion]
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

    fitness_resident = sum(compute_scientist_payoff_curve(domain, resident_submission, c, confusion, q_bar, version))
    fitness_mutant = sum(compute_scientist_payoff_curve(domain, mutant_submission, c, confusion, q_bar, version))

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


    average_quality_accepted(confusion::Array{Float64}, bit_vector::Array{Float64}, domain::Array{Float64})::Float64

Computes q_bar for a submission curve. This is given by the discretisation of the integral of the centroid.    
"""
function average_quality_accepted(confusion::Array{Float64}, bit_vector::Array{Float64}, domain::Array{Float64})::Float64
    #=
    Produces an area curve, then computes noise 
    (q1*e1+q2*e2+..)/(e1+e2+..)
    Computes q_bar of new batch
    =#
    numerator_vector = Float64[]
    denominator_vector = Float64[]
    for idx = 1:length(domain) 
        noise = confusion[idx] 
        push!(numerator_vector, noise * bit_vector[idx] * domain[idx])
        push!(denominator_vector, noise * bit_vector[idx])
    end

    numerator = sum_kbn(numerator_vector) # sum_kbn avoids numerical error
    denominator = sum_kbn(denominator_vector)

    if denominator == 0.0
        return 0.0
    else
        result = numerator / denominator
    end
    return result
end

"""
    acceptance_rate(confusion::Array{Float64}, submission_set::Array{Float64})::Float64

Computes the acceptance rate.

"""
function acceptance_rate(confusion::Array{Float64}, submission_set::Array{Float64})::Float64
    numerator_vector = Float64[]
    for idx = 1:length(submission_set) 
        push!(numerator_vector, submission_set[idx] * confusion[idx])
    end
    numerator = sum_kbn(numerator_vector)
    denominator = sum(submission_set) 
    if denominator == 0.0
        return 1.0
    else
        return numerator / denominator
    end
end

"""
    rejection_rate(confusion::Array{Float64}, submission_set::Array{Float64})::Float64

Computes rejection rate.

"""
function rejection_rate(confusion::Array{Float64}, submission_set::Array{Float64})::Float64
    return 1.0 - acceptance_rate(confusion, submission_set)
end

"""
    payoff_journal(confusion::Array, k::Float64, submission::Array, domain::Array, quality::Bool)::Float64

Computes π_journal = benefit-cost, with cost = 1/(1+ϵ)^k and benefit either q_bar or rejection rate.


"""
function payoff_journal(confusion::Array, k::Float64, submission::Array, domain::Array, quality::Bool)::Float64

    if quality
        benefit = average_quality_accepted(confusion, submission, domain) 
    else
        benefit = rejection_rate(confusion, submission)
    end

    cost = k * (1 - integrate(domain, confusion)) # Compute area of noise function, penalise strict curves
    result = benefit - cost

    return result

end

""" 

    journal_fitness(resident::Array{Float64}, mutant::Array{Float64}; kwargs...)

Measures journal fitness of resident and mutant.
"""
function journal_fitness(resident::Array{Float64}, mutant::Array{Float64}; kwargs...)
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

    fitness_resident = payoff_journal(resident, k, submission, domain, quality)
    fitness_mutant = payoff_journal(mutant, k, submission, domain, quality)

    return fitness_resident, fitness_mutant

end

"""

    journal_mutation(confusion::Array{Float64}; kwargs...)

Evolves a confusion curve.
"""

function journal_mutation(confusion::Array{Float64}; kwargs...)

    mutant = deepcopy(confusion)
    if haskey(kwargs, :mutation_epsilon_journal)
        mutation_epsilon = kwargs[:mutation_epsilon_journal]
    end
    if haskey(kwargs, :domain)
        domain = kwargs[:domain]
    end
    if haskey(kwargs, :journal_width)
        width = kwargs[:width]
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


    is_inside = false
    attempt = 0

    while !is_inside
        mutant = _attempt_gaussian_mutation(confusion, mutation_epsilon, domain, width) #, lower_bound, upper_bound)
        is_inside = _is_within_bounds(mutant, lower_bound, upper_bound)
        attempt += 1
        if attempt > MAX_ITER
            throw(ErrorException("Attempted too many mutations without producing anythin within bounds"))
        end
    end
    return mutant
end


""" 

    co_evolve(domain::Array{Float64}, resident_submission::Array{Float64},
resident_confusion::Array{Float64}, iterations=1e5,
atol::Float64=1e-8, seed::Int=0, snapshot_iter=1e3; kwargs...)

Evolves a submission set in conjunction with a journal with specified parameters.
Need to specify:
k ∈ [0, 1]
c ∈ [0,1]


"""
function co_evolve(domain::Array{Float64}, resident_submission::Array{Float64},
    resident_confusion::Array{Float64}, iterations=1e5,
    atol::Float64=1e-8, seed::Int=0, snapshot_iter=1e3; kwargs...)

    if haskey(kwargs, :mutation_epsilon_journal)
        mutation_epsilon_journal = kwargs[:mutation_epsilon_journal]
    end
    if haskey(kwargs, :journal_width)
        journal_width = kwargs[:journal_width]
    end


    Random.seed!(seed)

    time_series_submission = [[0, deepcopy(resident_submission)]]
    time_series_confusion = [[0, resident_confusion]]

    previous_submission = zeros(size(resident_submission))
    previous_confusion = zeros(size(resident_confusion))

    time_series_q_bar = []

    for step = 1:iterations

        q_bar = average_quality_accepted(resident_confusion, resident_submission, domain)
        push!(time_series_q_bar, q_bar)

        # Compute and mutate science fitness at the same time,
        previous_submission = resident_submission
        previous_confusion = resident_confusion

        resident_submission, invasion = evolution_step(resident_submission, science_fitness,
            gaussian_mutation, atol;
            kwargs..., domain=domain, confusion=previous_confusion, q_bar=q_bar)

        # Compute for confusion
        resident_confusion, invasion = evolution_step(resident_confusion, journal_fitness,
            journal_mutation, atol;
            kwargs..., submission=previous_submission, domain=domain, mutation_epsilon_journal=mutation_epsilon_journal, journal_width=journal_width)

        if step % snapshot_iter == 0
            push!(time_series_submission, [step, resident_submission])
            push!(time_series_confusion, [step, resident_confusion])
        end


    end

    return Dict{String,Any}("resident_submission" => resident_submission,
        "resident_confusion" => resident_confusion,
        "time_series_submission" => time_series_submission,
        "time_series_confusion" => time_series_confusion
    )



end



function main()

    domain = Array(range(0.0, 1.0, length=100))
    resident_sub = 0.01 * ones(size(domain))
    resident_confuse = 0.01 * zeros(size(domain))

    result_dict = co_evolve(domain, resident_sub, resident_confuse,
        1e5, 1e-9, 0, 1e3; k=8.0,
        c=0.6, mutation_epsilon=0.001, mutation_epsilon_journal=0.001,
        lower_bound=0.0, upper_bound=1.0,
        width=0.01, journal_width=0.01, quality=true, version=3)

    println(result_dict)

end
