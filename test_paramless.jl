
using Test
include("paramless.jl")

@testset "testPointMutationWithBounds" begin
    x = Array(range(0.01, 1.0, length = 100))
    parent = zeros(size(x))
    mutant = point_mutation(parent; mutation_epsilon= 0.2, lower_bound=0.0, upper_bound=0.5)
    @test maximum(mutant) <= 0.5 
    @test minimum(mutant) >= 0.0

end


@testset "testWithinBounds" begin 
    z = zeros(10)
    @test _within_bounds(z, -1.0, 1.0)
    @test _within_bounds(z, nothing,nothing)

end 