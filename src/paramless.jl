using Random

DEF_ATOL = 1e-8

MAX_ITER = 1e6

"""
    _is_within_bounds(vector, lower_bound, upper_bound)

Checks if all elements of an array are contained within [`lower_bound`, `upper_bound`].

"""
function _is_within_bounds(vector::Array, lower_bound, upper_bound)

    if lower_bound !== nothing && minimum(vector) < lower_bound
        return false
    end
    if upper_bound !== nothing && maximum(vector) > upper_bound
        return false
    end
    return true
end

""" 

    _is_within_bounds_scalar(param, lower_bound, upper_bound) 

Checks if a parameter satisfies lower_bound <= param <= upper_bound

"""
function _is_within_bounds_scalar(param, lower_bound, upper_bound)
    if lower_bound !== nothing && param < lower_bound
        return false
    end
    if upper_bound !== nothing && param > upper_bound
        return false
    end
    return true

end 

"""

Binary mutation in which a random submission is flipped.

"""
function binary_mutation(vector::BitArray)

    mutant = deepcopy(vector)

    index = rand(1:length(vector))

    mutant[index] = 1 - mutant[index]
    return mutant 

end 

"""

Performs a point mutation by shifting a point up/down.

"""
function _attempt_point_mutation(vector::Array, mutation_epsilon::Float64)::Array
    mutant = deepcopy(vector)
    position = rand(1:length(vector), 1)

    mutant[position] += rand([-1, 1], 1) * mutation_epsilon

    return mutant
end

"""

Mutates the line until it lies within bounds.

"""
function point_mutation(vector::Array; kwargs...)::Array
    is_inside = false
    attempt = 0
    mutant = deepcopy(vector)

    if haskey(kwargs, :mutation_epsilon)
        mutation_epsilon = kwargs[:mutation_epsilon]
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
        mutant = _attempt_point_mutation(vector, mutation_epsilon) # Modify their DNA
        is_inside = _is_within_bounds(mutant, lower_bound, upper_bound) # Check if everything is contained
        attempt += 1

        if attempt > MAX_ITER
            throw(ErrorException("Attempted too many mutations without producing anythin within bounds"))
        end
    end
    return mutant
end

"""

Performs point mutation on multiple points.

"""
function point_mutation_distribution(vector::Array; kwargs...)::Array
    # moves up a random point, and down a random point at the same time
    if haskey(kwargs, :mutation_epsilon)
        mutation_epsilon = kwargs[:mutation_epsilon]
    end

    mutant = deepcopy(vector)
    pos_up, pos_down = rand(1:length(vector), 2)
    adjusted_epsilon = mutation_epsilon
    # taking care that no negative values are allowed
    if mutant[pos_down] - mutation_epsilon < 0
        adjusted_epsilon = mutant[pos_down]
    end
    mutant[pos_up] = mutant[pos_up] + adjusted_epsilon
    mutant[pos_down] = mutant[pos_down] - adjusted_epsilon
    return mutant
end

function one_norm_distance(u, v)
    return sum(abs.(u - v))
end

function distance_fitness_function(resident, mutant; kwargs...)

    if haskey(kwargs, :target_surface)
        target_surface = kwargs[:target_surface]
    end

    fitness_res = 1.0 / one_norm_distance(resident, target_surface)
    fitness_mut = 1.0 / one_norm_distance(mutant, target_surface)
    return fitness_res, fitness_mut
end

""" 

Gaussian smoothing function.

"""
function _gaussian_mutation_helper(x, mutation_epsilon::Real, loc::Real, width::Real)
    return mutation_epsilon .* (ℯ .^ (-(x .- loc) .^ 2 ./ width))
end

function _attempt_gaussian_mutation(vector::Array, mutation_epsilon::Float64, domain, width::Float64)::Array
    location_idx = rand(1:length(vector), 1)[1]
    location_val = domain[location_idx]
    mutant = deepcopy(vector)

    perturbation = _gaussian_mutation_helper(domain, mutation_epsilon, location_val, width * rand(Float64, 1)[1])
    mutant .+= rand([-1, 1], 1)[1] .* perturbation

    return mutant
end


function gaussian_mutation(vector::Array; kwargs...)::Array
    #=
    Should generalise this, abstract away method of mutation
    , mutation_epsilon, domain, width, lower_bound=nothing, upper_bound=nothing
    =#
    # Parse args 
    mutant = deepcopy(vector)
    if haskey(kwargs, :mutation_epsilon)
        mutation_epsilon = kwargs[:mutation_epsilon]
    end
    if haskey(kwargs, :domain)
        domain = kwargs[:domain]
    end
    if haskey(kwargs, :width)
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
        mutant = _attempt_gaussian_mutation(vector, mutation_epsilon, domain, width) #, lower_bound, upper_bound)
        is_inside = _is_within_bounds(mutant, lower_bound, upper_bound)
        attempt += 1
        if attempt > MAX_ITER
            throw(ErrorException("Attempted too many mutations without producing anythin within bounds"))
        end
    end
    return mutant
end

function gaussian_mutation_distribution(vector::Array; kwargs...)::Array
    #= moves up a random point, and down a random point at the same time
    - Produces a smoothing effect around each of the points like a bump
    =#
    if haskey(kwargs, :mutation_epsilon)
        mutation_epsilon = kwargs[:mutation_epsilon]
    end
    if haskey(kwargs, :domain)
        domain = kwargs[:domain]
    end
    if haskey(kwargs, :width)
        width = kwargs[:width]
    end



    mutant = deepcopy(vector)
    pos_up, pos_down = rand(1:length(vector), 2)
    min_val = minimum(vector)

    if mutant[pos_down] - mutation_epsilon < 0
        adjusted_epsilon = min_val
    else
        adjusted_epsilon = mutation_epsilon
    end
    width = rand(Float64, 1)[1] * width
    perturbation_up = _gaussian_mutation_helper(
        domain, adjusted_epsilon, pos_up, width)
    perturbation_down = _gaussian_mutation_helper(
        domain, adjusted_epsilon, pos_down, width)

    mutant .+= perturbation_up
    mutant .-= perturbation_down

    return mutant
end




function evolution_step(resident_surface, fitness_function, mutation_function, atol; kwargs...)
    """
    One generation iteration
    """
    invasion = false
    mutant = mutation_function(resident_surface; kwargs...)
    fitness_resident, fitness_mutant = fitness_function(
        resident_surface, mutant; kwargs...)
    if fitness_resident < fitness_mutant && abs(fitness_resident - fitness_mutant) > atol
        resident_surface = deepcopy(mutant)
        invasion = true
    end

    return resident_surface, invasion
end

# return_time_series=false
function evolve(initial_surface, fitness_function, mutation_function, iterations, atol=DEF_ATOL, seed=0; kwargs...)
    """
    Evolve 
    Returns last resident , plus time series data if required
    """
    Random.seed!(seed)

    last_entry_time = 0
    resident = deepcopy(initial_surface)
    seq = 0

    time_series = []

    previous_resident = zeros(size(initial_surface))

    for step = 1:iterations
        previous_resident = deepcopy(resident)
        resident, invasion = evolution_step(
            resident, fitness_function, mutation_function, atol; kwargs...)

        if (invasion) # Update when a mutant gene takes over
            push!(time_series, [step - last_entry_time, previous_resident])
            last_entry_time = step
        end
    end
    return resident, time_series
end