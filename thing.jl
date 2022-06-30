
include("paramless.jl")

Random.seed!(10)


#=
x = Array(range(0.01, 1.0, length = 100))
parent = zeros(size(x))
mutant = point_mutation(parent; mutation_epsilon= 0.2, lower_bound=0.0, upper_bound=0.5)
println(minimum(mutant))
=#
function target_function(x::Float64)
    return x^2
end 

x = Array(range(-1.0,1.0,length=100));
target = target_function.(x);
initial_surface = zeros(size(x));

iterations = 1e4 
mutation_epsilon = 0.01
width = 0.05

ans_soft, series_soft = evolve(initial_surface, distance_fitness_function, gaussian_mutation, iterations, 1e-12, 777; target_surface=target, mutation_epsilon=mutation_epsilon, domain=x, width=width)