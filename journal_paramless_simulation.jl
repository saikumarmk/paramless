include("src/journal_paramless.jl")
using JSON

if length(ARGS) > 1

    version = parse(Int, ARGS[1])
    N = parse(Int, ARGS[2])
    c = parse(Float64, ARGS[3])
    k = parse(Float64, ARGS[4])
    quality = parse(Bool, ARGS[5])
    iterations = Int(parse(Float64, ARGS[6]))
    snapshot_iter = Int(parse(Float64, ARGS[7]))
    sci_epsilon = parse(Float64, ARGS[8])
    journal_epsilon = parse(Float64, ARGS[9])
    width = parse(Float64, ARGS[10])

    domain = Array(range(0.0, 1.0, length=N))
    resident_sub = 0.01 * ones(size(domain))
    resident_eps = 0.0
    resident_qt = 0.0


    result_dict = co_evolve(domain, resident_sub,
        resident_eps, resident_qt,
        iterations, 1e-9,
        0; k=k,
        c=c, mutation_epsilon=sci_epsilon,
        mutation_epsilon_journal=journal_epsilon,
        lower_bound=0.0, upper_bound=1.0,
        width=width, snapshot_iter=snapshot_iter)

    filename = string("journal_pl_version_", version, "_c_", c, "_k_", k, "_quality_", quality,
        "_iterations_", ARGS[6], "_snapshot_iter_", ARGS[7], "_sci_epsilon_",
        sci_epsilon, "_journal_epsilon_",
        journal_epsilon, "_width_", width, ".json")

    open(filename, "w") do file
        JSON.print(file, result_dict, 4)
    end
end