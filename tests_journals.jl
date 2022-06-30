include("journals.jl")

using Test

@testset "sortednearest" begin
    @test  searchsortednearest([1, 2, 3], 1) == 1
    @test  searchsortednearest([1, 3, 4], 1.5) == 1
    @test  searchsortednearest([1, 3, 4], 5) == 3
end

@testset "mask_submit_above_threshold" begin
    resample(11)
    # println(all_q)
    @test mask_submit_above_threshold(0.1) == [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test mask_submit_above_threshold(0.2) == [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test mask_submit_above_threshold(1.0) == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    @test mask_submit_above_threshold(0.0) == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    resample()
end

@testset "mask_submit_between(low_index, high_index)" begin
    resample(11)
    @test mask_submit_between(1, 11) == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test mask_submit_between(2, 11) == [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    @test mask_submit_between(2, 10) == [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
    @test mask_submit_between(3, 9) == [0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0]
end

@testset "complement_set(set_of_masks)" begin
    @test complement_set([[0, 1, 1],[1, 1, 1], [0, 0, 0]]) ==  [[1, 0, 0],[0, 0, 0], [1, 1, 1]]
end


@testset "all_above_threshold_sets" begin
resample(5)
result = [[0, 1, 1, 1, 1],
 [0, 0, 1, 1, 1],
 [0, 0, 0, 1, 1],
 [0, 0, 0, 0, 1]]
@test all_above_threshold_sets() == result
resample()
end

@testset "all_convex_sets" begin
resample(5)    
result = [[1, 1, 0, 0, 0], [1, 1, 1, 0, 0], [1, 1, 1, 1, 0], [1, 1, 1, 1, 1], [0, 1, 1, 0, 0],
[0, 1, 1, 1, 0], [0, 1, 1, 1, 1], [0, 0, 1, 1, 0], [0, 0, 1, 1, 1], [0, 0, 0, 1, 1]]
@test all_convex_sets() == result
resample()
end

@testset "gaussian_j" begin
@test gaussian_j(0, 0, 0.01) == 0.5
@test gaussian_j(0, 3, 0.01) == 0.0
@test gaussian_j(3, 0, 0.01) == 1.0
end


@testset "payoff" begin
@test payoff_journal(0.5, 1, 0.1; version = 1, quality = false, rate = true, k = 0.1) <= 0 #-0.9210430274501143
end


@testset "acceptance_rate" begin
    q_threshold = 1.0
    ϵ = 1.0
    version = 2
    c = 0.5
    submission_set = [q for q in all_q if payoff_scientist(q, q_threshold, c, ϵ, version=version) ≥ 0]
    #print(acceptance_rate(q_threshold, ϵ, submission_set))
    @test 0 < acceptance_rate(q_threshold, ϵ, submission_set) < 1
    
end

####################################
## Kevin's sanity checks begin here
####################################

@testset "Test 0.1" begin
    k = 4.0
    for quality in [true, false]
        rate = ~quality
        results = []
        for dupla in Iterators.product(all_q, all_q)
            c, qt = dupla
            # ACHTUNG  This one has edge problems for qt=1, check what is q bar and r for qt=1, qt=0.99


            # TASK: For test 0.1, the issue is with the gaussian definition, for ϵ =0 should approach -1/2 when qt goes to 0, and 0 when qt goes to 1 (check the edges)

            @test payoff_journal(qt, c, 0.0; version=0, quality=quality, rate=rate, k=k) <= 0  
        end
    end   
end


@testset "Test 0.2" begin
    k = 4.0
    rango = 0.01:0.01:0.99
    for quality in [true, false]
        rate = ~quality
        for dupla in Iterators.product(rango, rango)
            c, epsilon = dupla
            results = []
            for qt in rango
                #quality_threshold, c, epsilon
                push!(results, payoff_journal(qt, c, epsilon; version=0, quality=quality, rate=rate, k=k))
            end
            # Max is always at largest qt for all available values
            @test(argmax(results)== length(results))
        end    
    end
end

@testset "Test 1.1" begin
    k = 4.0
    rango = 0.01:0.01:0.99
    for quality in [true, false]
        rate = ~quality
        results = []
        for dupla in Iterators.product(rango, rango)
            c, qt = dupla
            @test payoff_journal(qt, c, 0.0; version=1, quality=quality, rate=rate, k=k) <= 0  
        end
    end   
end


@testset "Test 1.2" begin
    # quality 
    @test all(equilibrium(0; version=1, quality=true, rate=false,  k=4.0)["mask_submission"] .== 1)
    # rate
    @test all(equilibrium(0; version=1, quality=false, rate=true,  k=4.0)["mask_submission"] .== 1)
    c = 0
    ϵ = 0
    for quality in [true, false]
        rate = ~quality
        for q_threshold in all_q
            for q in all_q
                if q < q_threshold
                    @test payoff_scientist(q, q_threshold, c, ϵ, version=1) == 0.0
                else 
                    @test payoff_scientist(q, q_threshold, c, ϵ, version=1) != 0.0
                end
            end
        end

    end
end


# ACHTUNG This one is failing at q= 1
@testset "Test 1.3" begin
    # c = 1
    @test all(equilibrium(1.0; version=1, quality=true, rate=false,  k=4.0)["mask_submission"][1:end-1]  .== 0)
    @test all(equilibrium(1.0; version=1, quality=false, rate=true,  k=4.0)["mask_submission"][1:end-1]   .== 0)
    # c > 1
    @test all(equilibrium(1.0; version=1, quality=true, rate=false,  k=4.0)["mask_submission"][1:end-1]  .== 0)
    @test all(equilibrium(1.0; version=1, quality=false, rate=true,  k=4.0)["mask_submission"][1:end-1]   .== 0)
  
end


function is_step(array)
    index_submit = findfirst(isequal(1),array)
    if isnothing(index_submit)
         return true
    end
    return all(array[1:index_submit-1].==0)  && all(array[index_submit:end].==1) 
end



@testset "Test 1.4" begin
    k = 4.0
    rango = 0.01:0.01:0.99
    for quality in [true, false] # ACHTUNG fails for True
        rate = ~quality
        results = []
        for c in rango
            @test is_step(equilibrium(c; version=1, quality=quality, rate=rate,  k=k)["mask_submission"])
        end
    end   
  
end


@testset "Test 2.1" begin
    k = 4.0
    rango = 0.01:0.01:0.99
    for quality in [true, false]
        rate = ~quality
        results = []
        for dupla in Iterators.product(rango, rango)
            c, qt = dupla
            @test payoff_journal(qt, c, 0.0; version=2, quality=quality, rate=rate, k=k) <= 0  
        end
    end   
end



@testset "Test 2.2" begin
    # quality 
    @test all(equilibrium(0; version=2, quality=true, rate=false,  k=4.0)["mask_submission"] .== 1)
    # rate
    @test all(equilibrium(0; version=2, quality=false, rate=true,  k=4.0)["mask_submission"] .== 1)
    c = 0
    ϵ = 0
    for quality in [true, false]
        rate = ~quality
        for q_threshold in all_q
            for q in all_q
                if q < q_threshold
                    @test payoff_scientist(q, q_threshold, c, ϵ, version=2) == 0.0
                else 
                    @test payoff_scientist(q, q_threshold, c, ϵ, version=2) != 0.0
                end
            end
        end

    end
  
end

@testset "Test 2.3" begin

    # First part
    k = 4.0
    for quality in [false, true] # ACHTUNG fails for True
        rate = ~quality
        results = []
        for c in all_q
            @test equilibrium(c; version=2, quality=true, rate=false,  k=4.0)["mask_submission"][1] == 1
        end
    end   


    #second part
    rango = 0.1:0.1:1.0
    ϵ = 0.0
    for q_threshold in rango
        for c in all_q
            @test payoff_scientist(0.0, q_threshold, c, ϵ; version=2) == 0.0 
        end
    end
   
  
 end