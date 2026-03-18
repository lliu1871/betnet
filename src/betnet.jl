module betnet

    using Distributions
    using BenchmarkTools
    using Base.Threads: @spawn, @threads
    using Random
    using StatsBase
    using FreqTables  
      
    const ERROR = 1e-5

    export
        #data and tree definition    
        TransNet,
        TransNet_N,
        Tborg
    include("betnet1.1.jl")
end # module
