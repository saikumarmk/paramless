using NumericalIntegration



function energy(surface, domain)
    return integrate(domain, surface, SimpsonEven())
end

function resource_density(a)
    return 4.0*(1.0-a)
end 

function total_intake(surface, domain)
    
    integrand = [resource_density(domain[idx])*(surface[idx]/(surface[idx]+domain[idx])) 
                            for idx = 1:length(surface)]
    return integrate(domain, integrand)
end

function metabolic_fitness(resident, mutant; kwargs...)
    if haskey(kwargs, :c)
        c = kwargs[:c]
    end 
    if haskey(kwargs, :domain)
        domain = kwargs[:domain]
    end
    
    fitness_resident = total_intake(resident, domain) - c*energy(resident, domain)
    fitness_mutant = total_intake(mutant, domain) - c*energy(mutant, domain)
    return fitness_resident, fitness_mutant
end