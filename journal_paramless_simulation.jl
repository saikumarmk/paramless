include("src/journal_gauss_free.jl")
using JSON

if length(ARGS) > 1

    version = parse(Int, ARGS[1])
    N = parse(Int, ARGS[2])
    c = parse(Float64, ARGS[3])
    k = parse(Float64, ARGS[4])
    quality = parse(Bool, ARGS[5])
    iterations = Int(parse(Float64, ARGS[6]))
    snapshot_iter = Int(parse(Float64, ARGS[7]))

    # Can remove from here 
    if length(ARGS) == 7
        seed = 0
        sci_epsilon = journal_epsilon = 1e-3
        width = journal_width = 1e-2
    elseif length(ARGS) == 8
        seed = parse(Int, ARGS[8])
        sci_epsilon = journal_epsilon = 1e-3
        width = journal_width = 1e-2
    else
        sci_epsilon = parse(Float64, ARGS[9])
        journal_epsilon = parse(Float64, ARGS[10])
        width = parse(Float64, ARGS[11])
        journal_width = parse(Float64, ARGS[12])
    end 

    domain = Array(range(0.0, 1.0, length=N))
    resident_sub = 0.01 * ones(size(domain))
    resident_confuse = 0.5 * ones(size(domain))

    settings = Dict{String,Any}("version" => version, "N" => N, 
                    "c" => c, "k" => k,
                    "quality"=>quality, "iterations"=>iterations,
                    "snapshot_iter" => snapshot_iter, "seed"=> seed,
                    "sci_epsilon"=>sci_epsilon,"journal_epsilon"=>journal_epsilon,
                    "width" => width, "journal_width" => journal_width)

    result_dict = co_evolve(domain, resident_sub,
        resident_confuse, iterations, 
        1e-9, seed, 
        snapshot_iter; k=k,
        c=c, mutation_epsilon=sci_epsilon,
        mutation_epsilon_journal=journal_epsilon,
        lower_bound=0.0, upper_bound=1.0,
        width=width, journal_width = journal_width, 
        quality=quality, version=version)

    filename = string("journal_pl_version_", version, "_c_", c, "_k_", k, "_quality_", quality,
        "_iterations_", ARGS[6], "_snapshot_iter_", ARGS[7], "_seed_", seed, ".json")

    open(filename, "w") do file
        JSON.print(file, merge!(Dict{String,Any}("settings"=>settings),result_dict), 4)
    end
end