module betnet
    using Distributions
    using BenchmarkTools
    using Base.Threads: @spawn, @threads
    using Random
    using StatsBase
    using FreqTables 
    using DataFrames
    using CSV
    using LinearAlgebra
    using Statistics
    using Base.Threads

    const NO_ERROR = 0
    const DEBUG = 1  # set to 1 to turn on debugging mode

    mutable struct modelParameters
        theta::Float64
        theta_lb::Float64 
        theta_ub::Float64
        theta_window::Float64
        theta_priorMean::Float64
        muRate::Float64
        muRate_lb::Float64
        muRate_ub::Float64
        muRate_window::Float64
        muRate_priorMean::Float64
        infectionRate::Float64
        infectionRate_lb::Float64
        infectionRate_ub::Float64
        infectionRate_window::Float64
        infectionRate_priorMean::Float64
        removalRate::Float64
        latent_priorMean::Float64
        Net_infID::Vector{Int64}
        Child_Vec::Vector{Int64}
        InfectionPeriod::Matrix{Float64}
        InfTime::Vector{Float64}            
        ContactProb::Matrix{Float64}
        genomeSize::Int64
        
        function modelParameters(n::Int64)
            theta = 1.0e-6
            theta_lb = 1.0e-10
            theta_ub = 1.0
            theta_window = 1.0e-7
            theta_priorMean = 1.0e-6
            muRate = 1.0e-6
            muRate_lb = 1.0e-10
            muRate_ub = 1.0
            muRate_window = 1.0e-7
            muRate_priorMean = 1.0e-6
            infectionRate = 2.0
            infectionRate_lb = 1.0e-10
            infectionRate_ub = 10.0
            infectionRate_window = 0.1
            infectionRate_priorMean = 2.0
            removalRate = 0.1 
            latent_priorMean = 0.5
            Net_infID = zeros(Int64, n)
            Child_Vec = zeros(Int64, n)
            InfectionPeriod = zeros(Float64, n, n)
            InfTime = zeros(Float64, n)
            ContactProb = zeros(Float64, n, n)
            genomeSize = 4411532
            new(theta, theta_lb,theta_ub, theta_window, theta_priorMean, muRate,muRate_lb,muRate_ub,muRate_window, muRate_priorMean,infectionRate, infectionRate_lb,infectionRate_ub,infectionRate_window,infectionRate_priorMean,removalRate,latent_priorMean,Net_infID,Child_Vec, InfectionPeriod, InfTime, ContactProb, genomeSize)
        end      
    end

    mutable struct mcmcParameter
        numIter::Int64
        burnIn::Int64
        thin::Int64
        adaptInterval::Int64
        adaptFactor::Float64
        adaptCount::Int64
        acceptCount::Int64
        proposalSD::Float64
        logLikelihood::Float64
        logPrior::Float64
        function mcmcParameter()
            numIter = 100000
            burnIn = 10000
            thin = 100
            adaptInterval = 1000
            adaptFactor = 1.1
            adaptCount = 0
            acceptCount = 0
            proposalSD = 1.0
            logLikelihood = 0.0
            logPrior = 0.0
            new(numIter, burnIn, thin, adaptInterval, adaptFactor, adaptCount, acceptCount, proposalSD, logLikelihood, logPrior)
        end
    end

    mutable struct caseData 
        numCases::Int64
        caseID::Vector{String}
        tempData::Matrix{Float64}
        SNPData::Matrix{Float64}
        contactData::Matrix{Float64}

        function caseData(x::Int64)
            numCases = x
            caseID = string.(1:x)
            tempData = zeros(Float64, x, 3)
            SNPData = zeros(Float64, x, x)
            contactData = zeros(Float64, x, x)
            new(numCases, caseID, tempData, SNPData, contactData)
        end
    end

    export transNetworkInference

    include("betnet2.0.jl")
end # module
