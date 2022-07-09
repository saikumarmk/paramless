using Statistics
using SpecialFunctions
using KahanSummation


function gaussian_j(q, q_threshold,ϵ)

    #edge cases:
    if q == q_threshold && ϵ==0.0
        return 0.5
    end
    result = 1.0 - ((1 + erf((q_threshold - q) / ϵ  / sqrt(2))) / 2)
    if isnan(result)    
        throw(error("Error computing gaussian_j: q = $q, q_threshold = $q_threshold, ϵ = $ϵ"))
    end
    return Float64(result)
end


function payoff_scientist(q, q_threshold, c, ϵ, q_bar, version = 3)
    #= 
    q = quality of scientist
    q_threshold = used to determine probability of acceptance 
    ϵ = noise generated by the journal 
    q_bar = average quality of the paper, has a role in the benefit  
    =#
    p = gaussian_j(q, q_threshold, ϵ)

    if version == 3
        return max(p * q_bar - (1.0 - p) * c ,0)
    elseif version == 4
        return max(p * q_bar - (1.0-p)*q*c, 0)
    else
        throw(ErrorException("Invalid Version, version ∈ {3,4}"))
end 

"""
    compute_scientist_payoff_curve(domain, submission_curve, q_threshold, c, ϵ, q_bar)

For a submission curve over [0,1], determine the collective payoff.
"""
function compute_scientist_payoff_curve(domain, submission_curve, q_threshold, c, ϵ, q_bar)
    return [payoff_scientist(domain[idx], q_threshold, c, ϵ, q_bar) 
            for idx = 1:length(domain) 
            if submission_curve[idx]==1]
end 


# Avoid kwargs?
function science_fitness(resident_submission, mutant_submission; kwargs...)


    fitness_resident = sum(compute_scientist_payoff_curve(domain, resident_submission, q_threshold, c, ϵ, q_bar))
    fitness_mutant = sum(compute_scientist_payoff_curve(domain, mutant_submission, q_threshold, c, ϵ, q_bar))

    return fitness_resident, fitness_mutant
end 

"""

Binary mutation in which a random submission is flipped.

"""
function binary_mutation(vector::BitArray; kwargs...)

    mutant = deepcopy(vector)

    index = rand(1:length(vector))

    mutant[index] = 1 - mutant[index]
    return mutant 

end 

function prepare_for_integral(bit_vector)
    #= 
    If bit is on, then add it for integration 

    =#
    ans = []
    @assert length(bit_vector) == length(all_q)
    for (index, bit) in enumerate(bit_vector)
        @assert bit ==1 || bit==0
        if bit ==1 
            push!(ans,all_q[index])
        end
    end
return ans
end

function average_quality_accepted(q_threshold, ϵ, bit_vector)
    #=
    Produces an area curve, then computes noise 
    (q1*e1+q2*e2+..)/(e1+e2+..)
    Computes q_bar of new batch
    =#
    submission_set = prepare_for_integral(bit_vector)
    numerator_vector = Float64[]
    denominator_vector = Float64[]
    for q in submission_set
        noise = gaussian_j(q, q_threshold, ϵ)
        push!(numerator_vector, noise*q)
        push!(denominator_vector, noise)
    end

    numerator = sum_kbn(numerator_vector) # sum_kbn avoids numerical error
    denominator = sum_kbn(denominator_vector)

    if denominator == 0.0 
        return 0.0
    else
        result = numerator/denominator
    end
    return result
end

function acceptance_rate(q_threshold, ϵ, submission_set)
    #=
    Computes number of accepted 

    =#
    numerator_vector = Float64[]
    for (i, mask_value) in enumerate(submission_set)
        if mask_value == 1
            q = all_q[i] # Check the quality associated, need to change 
            push!(numerator_vector,  gaussian_j(q, q_threshold, ϵ))
        end
    end
    numerator = sum_kbn(numerator_vector)
    denominator = sum(submission_set) # the rough length of accepted 
    if denominator == 0.0
        # we define that if nothing is sibmitted the acceptance rate is 1
        return 1.0 
    else
        return numerator/denominator
    end
end


function rejection_rate(q_threshold, ϵ, submission_set)
    return 1.0 - acceptance_rate(q_threshold, ϵ, submission_set)
end

function payoff_journal(q_threshold::Real, c::Real, ϵ::Real, k::Real, submission::BitArray, quality::Bool = true)

    if quality 
        benefit = average_quality_accepted(q_threshold, ϵ, submission)
    else 
        benefit = rejection_rate(q_threshold, ϵ, submission) 
    end 

    cost = 1/(1+ϵ)^k
    result = benefit - cost 

    return result

end 

""" 

Checks fitness by recomputing submission curve under new parameters.

"""
function fitness_journal(resident, mutant; kwargs...)

    if haskey(kwargs, :mutate_qt)
        mutate_qt = kwargs[:mutate_qt] 
    end 

    if haskey(kwargs, :ϵ)
        ϵ = kwargs[:ϵ]
    end # Parse at fitness level
    
    if mutate_qt 
        fitness_resident = payoff_journal(q_threshold, c, ϵ)
        fitness_mutant = payoff_journal(q_threshold, c, ϵ)
    else 
        fitness_resident = payoff_journal(q_threshold, c, ϵ)
        fitness_mutant = payoff_journal(q_threshold, c, ϵ)
    end 

    return fitness_resident, fitness_mutant

end 

""" 



:param domain: Discretisation of [0,1]


"""
function co_evolve(domain, resident_submission, resident_ϵ, resident_q_threshold,
                                                     iterations, atol, seed=0; kwargs...)

    Random.seed!(seed)

    time_series_submission = [deepcopy(resident_submission)]
    time_series_ϵ = [resident_ϵ]
    time_series_q_threshold = [resident_q_threshold]


    for step = 1:iterations
        # Compute and mutate science fitness
        resident_submission, invasion = evolution_step(resident_submission, science_fitness, 
                                    binary_mutation, atol; 
                                    kwargs..., domain=domain)

        if invasion 
            push!(time_series_submission, [step, resident_submission]) 
        end 

        # Compute for ϵ
        resident_ϵ, invasion = evolution_step()

        if invasion 
            push!(time_series_ϵ, [step, resident_ϵ])
        end 
        
        # compute for threshold
        resident_q_threshold, invasion = evolution_step()

        if invasion 
            push!(time_series_q_threshold, [step, resident_q_threshold])
        end 


    end 

    return resident_submission, resident_ϵ, resident_q_threshold, 
            time_series_resident, time_series_ϵ, time_series_q_threshold


end 



