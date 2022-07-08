using NumericalIntegration


function energy(surface, domain)
    integrate(domain, surface, SimpsonEven())
end

function prob_seed(a, x, k)
    return ℯ .^ (-x ./ k(a))
end

function k(a)
    return 100 * (2 + sin(2 * π * ((a - 1)^2 - 0.25)))
end

function c(e, α)
    return 2.0 / (1 + ℯ^(-α * e))
end


function birth(x_mutant, x_resident, domain, α, β, k_function)

    integrand = [x_mutant[idx] * c((x_mutant[idx] - x_resident[idx]) / (x_mutant[idx]^β), α) * prob_seed(domain[idx], x_resident[idx], k_function)
                 for idx = 1:length(domain)]
    return integrate(domain, integrand, SimpsonEven())
end

function flowering_fitness(x_resident, x_mutant; kwargs...)
    if haskey(kwargs, :domain)
        domain = kwargs[:domain]
    end
    if haskey(kwargs, :alpha)
        alpha = kwargs[:alpha]
    end
    if haskey(kwargs, :beta)
        beta = kwargs[:beta]
    end
    if haskey(kwargs, :k_function)
        k_function = kwargs[:k_function]
    end

    fitness_resident = birth(x_resident, x_resident, domain, alpha, beta, k_function)
    fitness_mutant = birth(x_mutant, x_resident, domain, alpha, beta, k_function)
    return fitness_resident, fitness_mutant
end