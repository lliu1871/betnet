function Initialization(; genomeSize::Int64, itr_MCMC::Int64, burn_in::Int64, subsample::Int64)
    # initialize mcmcAlgorithm
    global mcmc = mcmcParameter()
    mcmc.numIter = itr_MCMC
    mcmc.burnIn = burn_in
    mcmc.thin = subsample
 
    # initialize the parameters
    global parameters = modelParameters(data.numCases)
 
    # genonme size 4411532 for M. tuberculosis
    parameters.genomeSize = genomeSize

    # theta 
    parameters.theta_lb = 0.0 # lower bound of theta
    parameters.theta_ub = 1E-3 # upper bound of theta
    parameters.theta_window = (parameters.theta_ub - parameters.theta_lb)/100
    parameters.theta= rand(Uniform(parameters.theta_lb, parameters.theta_ub))
    parameters.theta_priorMean = 1E-5
   
    # mutation rate
    parameters.muRate_lb = 0.0 # lower bound of mu
    parameters.muRate_ub = 1E-3 # upper bound of mu
    parameters.muRate_window = (parameters.muRate_ub - parameters.muRate_lb)/100
    parameters.muRate = rand(Uniform(parameters.muRate_lb, parameters.muRate_ub))
    parameters.muRate_priorMean = 1E-5
    
    # infection rate
    parameters.infectionRate_lb = 0.0 # lower bound of inf
    parameters.infectionRate_ub = 20.0 # upper bound of inf
    parameters.infectionRate_window = (parameters.infectionRate_ub - parameters.infectionRate_lb)/100
    parameters.infectionRate = rand(Uniform(parameters.infectionRate_lb, parameters.infectionRate_ub))
    parameters.infectionRate_priorMean = 1.0

    # latent period prior mean
    parameters.latent_priorMean = 0.05 # mean of latent period ~ Chisq(0.5)
    
    # initial infectors with minimum SNP difference
    parameters.Net_infID[1] = 1 # the first case is the index case
    for j in 2:data.numCases
        parameters.Net_infID[j] = findmin(data.SNPData[j,1:(j-1)])[2] # initialize the infector of case j as 0
    end

    if DEBUG == 1
        if (minimum(parameters.Net_infID) < 1) || (maximum(parameters.Net_infID) > data.numCases)
            println("Infectors ID: ", parameters.Net_infID)
            error("Error: infectors!")
        end
    end

    # count how many children each patient has of the given transmission tree
    for j in 1:data.numCases
        parameters.Child_Vec[j] = count(i->(i==j), parameters.Net_infID[2:end])
    end

    if DEBUG == 1
        if sum(parameters.Child_Vec) != (data.numCases - 1)
            println("children vector ", parameters.Child_Vec)
            error("Error: children vector!")
        end
    end

    # the period in which i can infect j
    for j in 2:data.numCases
        for i in 1:(j-1) 
            if data.tempData[j,1] > data.tempData[i,2] #onset time of j > removal time of i
                parameters.InfectionPeriod[i,j] = data.tempData[i,2] - data.tempData[i,1] 
            else 
                parameters.InfectionPeriod[i,j] = data.tempData[j,1] - data.tempData[i,1]
            end
        end
    end

    if DEBUG == 1
        if minimum(parameters.InfectionPeriod) < 0
            println("Infection Period: ", parameters.InfectionPeriod)
            error("Error: negative latent periods!")
        end
    end
    
    # infection time
    for j in 2:data.numCases
        infector = parameters.Net_infID[j]
        parameters.InfTime[j] = data.tempData[infector,1] + parameters.InfectionPeriod[infector,j] * rand()
    end
        
    if DEBUG == 1
        if (minimum(parameters.InfTime) < 0) || (maximum(parameters.InfTime - data.tempData[:,1]) > 0)
            println("Infection Time: ", parameters.InfTime)
            error("Error: infection times!")
        end
    end

    # estimate removal rate using the average of observed removal times
    parameters.removalRate = 1/mean(data.tempData[:,2] - data.tempData[:,1])
    if DEBUG == 1
        if parameters.removalRate <= 0
            println("Removal Rate: ", parameters.removalRate)
            error("Error: removal rate!")
        end
    end
    
    # contact probability
    if maximum(data.contactData) > 0.0
        parameters.ContactProb = data.contactData
        # normalize the contact probability matrix
        for j in 2:data.numCases
            s = sum(data.contactData[j,1:(j-1)])
            for i in 1:(j-1)
                parameters.ContactProb[j,i] = parameters.ContactProb[i,j] = data.contactData[i,j]/s
            end
        end
    else
        # if no contact data, assume all cases are equally likely to infect case j
        for j in 2:data.numCases
            for i in 1:(j-1)
                parameters.ContactProb[j,i] = parameters.ContactProb[i,j] = 1.0/float(j-1)
            end
        end
    end

    if DEBUG == 1
        if (minimum(parameters.ContactProb) < 0.0) || (maximum(parameters.ContactProb) > 1.0)
            println("Contact Probability: ", parameters.ContactProb)
            error("Error: contact probabilities!")
        end
    end
    
    return NO_ERROR
end

function Loglikelihood(theta::Float64,mu::Float64,Net_infID::Vector{Int64},genomeSize::Int64,InfTime::Vector{Float64})
    loglikelihood = 0.0
    for j in 2:data.numCases
        i = Int64(Net_infID[j]) # infector of j
        snp = data.SNPData[i,j]
        time = data.tempData[j,2] + data.tempData[i,2] - 2 * InfTime[j]
        p = 3/4-3/4*exp(-mu * time)/(2 * theta + 1)
        loglikelihood += snp*log(p) + (genomeSize-snp)*log(1-p)                       
    end
    return loglikelihood
end

function LoglikelihoodSingleEvent(individual::Int64,infector::Int64,InfTime::Float64,theta::Float64,mu::Float64,genomeSize::Int64)
    j = individual
    i = infector # individual i infected individual j in the last loop
    snp = data.SNPData[i,j]
    time = data.tempData[j,2] + data.tempData[i,2] - 2 * InfTime
    p = 3/4-3/4*exp(-mu * time)/(2 * theta + 1)
    loglikelihood = snp*log(p) + (genomeSize-snp)*log(1-p)                       
    return loglikelihood
end

function moveTheta()   
    ## update theta
    theta_old = parameters.theta
    theta_new = rand(Uniform(max(parameters.theta_lb, theta_old - parameters.theta_window), min(parameters.theta_ub,theta_old + parameters.theta_window)))

    # loglikelihood ratio new/old
    loglike_old = mcmc.logLikelihood
    loglike_new = Loglikelihood(theta_new, parameters.muRate, parameters.Net_infID, parameters.genomeSize,parameters.InfTime)
    diff_loglikelihood = loglike_new - loglike_old
    
    # prior ratio new/old
    diff_logprior = (theta_old - theta_new)/parameters.theta_priorMean
 
    # proposal ratio old/new
    diff_logproposal = log(1/(min(parameters.theta_ub,theta_old + parameters.theta_window)-max(parameters.theta_lb, theta_old - parameters.theta_window))) - 
                     log(1/(min(parameters.theta_ub,theta_new + parameters.theta_window)-max(parameters.theta_lb, theta_new - parameters.theta_window)))
    
    # accept new theta
    hastingRatio = diff_loglikelihood + diff_logprior + diff_logproposal
    if rand() < exp(hastingRatio)
        # accept new theta
        parameters.theta = theta_new
        mcmc.logLikelihood += diff_loglikelihood
        mcmc.logPrior += diff_logprior
    end

    return NO_ERROR
end

function moveMutationRate()
    ## update muRate
    muRate_old = parameters.muRate
    muRate_new = rand(Uniform(max(parameters.muRate_lb, muRate_old - parameters.muRate_window), min(parameters.muRate_ub,muRate_old + parameters.muRate_window)))

    # loglikelihood ratio new/old
    loglike_old = mcmc.logLikelihood
    loglike_new = Loglikelihood(parameters.theta, muRate_new, parameters.Net_infID, parameters.genomeSize,parameters.InfTime)
    diff_loglikelihood = loglike_new - loglike_old

    # prior ratio new/old
    diff_logprior = (muRate_old - muRate_new)/parameters.muRate_priorMean
                    
    # proposal ratio old/new
    diff_logproposal = log(1/(min(parameters.muRate_ub,muRate_old + parameters.muRate_window)-max(parameters.muRate_lb, muRate_old - parameters.muRate_window))) - 
                     log(1/(min(parameters.muRate_ub,muRate_new + parameters.muRate_window)-max(parameters.muRate_lb, muRate_new - parameters.muRate_window)))
    
    # hastings ratio
    hastingRatio = diff_loglikelihood + diff_logprior + diff_logproposal

    # accept new muRate
    if rand() < exp(hastingRatio) 
        parameters.muRate = muRate_new
        mcmc.logLikelihood += diff_loglikelihood
        mcmc.logPrior += diff_logprior
    end

    return NO_ERROR
end

function moveInfectionRate()
    # update infectionRate
    infectionRate_old = parameters.infectionRate
    infectionRate_new = rand(Uniform(max(parameters.infectionRate_lb, infectionRate_old - parameters.infectionRate_window), min(parameters.infectionRate_ub,infectionRate_old + parameters.infectionRate_window)))

    # logprior ratio new/old
    logprior_phi_old = 0.0
    for i in 1: (data.numCases-1)
        logprior_phi_old += log(pdf(Geometric(parameters.removalRate/(infectionRate_old+parameters.removalRate)),parameters.Child_Vec[i]))
    end

    logprior_phi_new = 0.0
    for i in 1: (data.numCases-1)
        logprior_phi_new += log(pdf(Geometric(parameters.removalRate/(infectionRate_new+parameters.removalRate)),parameters.Child_Vec[i]))
    end

    diff_logprior = logprior_phi_new - logprior_phi_old 
    diff_logprior += (infectionRate_old - infectionRate_new)/parameters.infectionRate_priorMean
    
    # logproposal ratio old/new
    diff_logproposal = log(1/(min(parameters.infectionRate_ub,infectionRate_old + parameters.infectionRate_window)-max(parameters.infectionRate_lb, infectionRate_old - parameters.infectionRate_window))) - 
                     log(1/(min(parameters.infectionRate_ub,infectionRate_new + parameters.infectionRate_window)-max(parameters.infectionRate_lb, infectionRate_new - parameters.infectionRate_window)))
    
    # hastings ratio
    hastingRatio = diff_logprior + diff_logproposal

    # accept new infection rate
    if rand() < exp(hastingRatio)
        parameters.infectionRate = infectionRate_new
        mcmc.logPrior += diff_logprior
    end
    
    return NO_ERROR
end

#Update transmission tree
function moveInfector()     
    # select one individual to update its infector
    # no need to update the first two individuals
    probs = (3:data.numCases)/sum(3:data.numCases)
    if DEBUG == 1
        if length(probs) != (data.numCases - 2)
            error("Error: probs length!")
        end
        if sum(probs) != 1.0
            error("Error: probs sum!")
        end
    end
    individual_update = sample(3:data.numCases, Weights(probs)) 
    
    # randomly select a new infector
    infector_old = Int64(parameters.Net_infID[individual_update])
    infector_new = rand(1:(individual_update-1))

    # update InfTime for the new infector
    InfTime_old = parameters.InfTime[individual_update]
    InfTime_new = data.tempData[infector_new,1] + parameters.InfectionPeriod[infector_new,individual_update] * rand()

    # loglikelihood ratio new/old
    loglike_old = LoglikelihoodSingleEvent(individual_update, infector_old, InfTime_old, parameters.theta, parameters.muRate, parameters.genomeSize)
    loglike_new = LoglikelihoodSingleEvent(individual_update, infector_new, InfTime_new, parameters.theta, parameters.muRate, parameters.genomeSize)
    diff_loglikelihood = loglike_new - loglike_old
 
    # logprior ratio new/old
    if infector_old == infector_new
        latent_time_old = data.tempData[individual_update,1] - InfTime_old
        latent_time_new = data.tempData[individual_update,1] - InfTime_new 
        if DEBUG == 1
            if latent_time_new < 0 || latent_time_old < 0
                error("Error: negative latent time!")
            end
        end
        # check if latent time is 0, if yes, use a small value to calculate logprior to avoid log(0)
        if latent_time_new == 0.0
            diff_logprior = log(pdf(Chisq(parameters.latent_priorMean), 1.0e-10)) - log(pdf(Chisq(parameters.latent_priorMean), latent_time_old))
        elseif latent_time_old == 0.0
             diff_logprior = log(pdf(Chisq(parameters.latent_priorMean), latent_time_new)) - log(pdf(Chisq(parameters.latent_priorMean), 1.0e-10))
        else
            diff_logprior = log(pdf(Chisq(parameters.latent_priorMean), latent_time_new)) - log(pdf(Chisq(parameters.latent_priorMean), latent_time_old))
        end #diff_logprior = log(pdf(Chisq(parameters.latent_priorMean), latent_time_new)) - log(pdf(Chisq(parameters.latent_priorMean), latent_time_old))
    else
        latent_time_old = data.tempData[individual_update,1] - InfTime_old
        latent_time_new = data.tempData[individual_update,1] - InfTime_new 
        if DEBUG == 1
            if latent_time_new < 0 || latent_time_old < 0
                error("Error: negative latent time!")
            end
        end
        # check if latent time is 0, if yes, use a small value to calculate logprior to avoid log(0)
        if latent_time_new == 0.0
            diff_logprior = log(pdf(Chisq(parameters.latent_priorMean), 1.0e-10)) - log(pdf(Chisq(parameters.latent_priorMean), latent_time_old))
        elseif latent_time_old == 0.0
             diff_logprior = log(pdf(Chisq(parameters.latent_priorMean), latent_time_new)) - log(pdf(Chisq(parameters.latent_priorMean), 1.0e-10))
        else
            diff_logprior = log(pdf(Chisq(parameters.latent_priorMean), latent_time_new)) - log(pdf(Chisq(parameters.latent_priorMean), latent_time_old))
        end
        #diff_logprior = log(pdf(Chisq(parameters.latent_priorMean), latent_time_new)) - log(pdf(Chisq(parameters.latent_priorMean), latent_time_old))
        diff_logprior += log(parameters.ContactProb[individual_update,infector_new]) - log(parameters.ContactProb[individual_update,infector_old])   
        diff_logprior += log(pdf(Geometric(parameters.removalRate/(parameters.infectionRate+parameters.removalRate)),parameters.Child_Vec[infector_new]+1)) +
                     log(pdf(Geometric(parameters.removalRate/(parameters.infectionRate+parameters.removalRate)),parameters.Child_Vec[infector_old]-1)) -
                     log(pdf(Geometric(parameters.removalRate/(parameters.infectionRate+parameters.removalRate)),parameters.Child_Vec[infector_old])) -
                     log(pdf(Geometric(parameters.removalRate/(parameters.infectionRate+parameters.removalRate)),parameters.Child_Vec[infector_new]))  
    end

    # hastings ratio
    hastingRatio = diff_loglikelihood + diff_logprior

    # accept new infector
    if rand() < exp(hastingRatio)
        # update loglikelihood, logprior
        mcmc.logLikelihood += diff_loglikelihood
        mcmc.logPrior += diff_logprior
        # update Net_infID, InfTime, Child_Vec
        parameters.Net_infID[individual_update] = infector_new
        parameters.InfTime[individual_update] = InfTime_new
        parameters.Child_Vec[infector_old] -= 1
        parameters.Child_Vec[infector_new] += 1           
        if DEBUG == 1
            if sum(parameters.Child_Vec) != (data.numCases - 1)
                println("children vector ", parameters.Child_Vec)
                error("Error: children vector!")
            end

            for j in 1:data.numCases
                x = count(i->(i==j), parameters.Net_infID[2:end])
                if x != parameters.Child_Vec[j]
                    println("children vector ", parameters.Child_Vec)
                    error("Error: children vector!")
                end
            end
        end    
    end
                     
    return NO_ERROR
end

function mcmcAlgorithm(outputfile::String)
    # initial loglikelihood
    mcmc.logLikelihood = Loglikelihood(parameters.theta, parameters.muRate, parameters.Net_infID, parameters.genomeSize,parameters.InfTime)
    
    # initial logprior
    mcmc.logPrior = log(pdf(Exponential(parameters.theta_priorMean), parameters.theta)) + log(pdf(Exponential(parameters.muRate_priorMean), parameters.muRate)) + log(pdf(Exponential(parameters.infectionRate_priorMean), parameters.infectionRate))
    for k in 1:(data.numCases-1)
        mcmc.logPrior += log(pdf(Geometric(4.0/(parameters.infectionRate+4.0)),count(i->(i== k), parameters.Net_infID[2:end])))
    end

    for j in 2:data.numCases
        infID = parameters.Net_infID[j]
        latent_time = data.tempData[j,1] - parameters.InfTime[j]        
        if DEBUG == 1
            if latent_time < 0
                error("Error: negative latent time!")
            end
        end
        # if latent_time == 0, we use a small value 1.0e-10 to calculate the logprior to avoid log(0)
        if(latent_time == 0.0)
            mcmc.logPrior += log(pdf(Chisq(parameters.latent_priorMean), 1.0e-10)) - log(cdf(Chisq(parameters.latent_priorMean), 1.0e-10)) + log(parameters.ContactProb[infID,j])
        else
            mcmc.logPrior += log(pdf(Chisq(parameters.latent_priorMean), latent_time)) - log(cdf(Chisq(parameters.latent_priorMean), parameters.InfectionPeriod[infID,j])) + log(parameters.ContactProb[infID,j])
        end      
        #mcmc.logPrior += log(pdf(Chisq(parameters.latent_priorMean), latent_time)) - log(cdf(Chisq(parameters.latent_priorMean), parameters.InfectionPeriod[infID,j])) + log(parameters.ContactProb[infID,j])
    end
    
    for round in 1:mcmc.numIter       
        r_random = rand()    
        if r_random <= 0.1
            moveTheta()
        elseif r_random >0.1 && r_random <=0.2
            moveMutationRate()    
        elseif r_random >0.2 && r_random <=0.3
            moveInfectionRate()
        else
            moveInfector()
        end  

        if round == 1 || (round > mcmc.burnIn && round%mcmc.thin == 0)
            println("...... Iteration ",round,": Loglikelihood ",trunc(mcmc.logLikelihood,digits=2),", theta ",parameters.theta,", mu ",parameters.muRate,", inf_rate ",parameters.infectionRate)
            CSV.write(outputfile,DataFrame(reshape(append!([string(round)],map(string,[mcmc.logLikelihood,mcmc.logPrior,parameters.theta,parameters.muRate,parameters.infectionRate]),data.caseID[parameters.Net_infID],string.(parameters.InfTime)), 1, :), :auto),append=true)
        end
    end

    return NO_ERROR
end

function transNetworkInference(; tempfile::String, SNPfile::String, Contactfile::String, genomeSize::Int64, itr_MCMC::Int64, burn_in::Int64, subsample::Int64, outputfile::String="")
    # read temporal Data
    x = DataFrame(CSV.File(tempfile))    
    global data = caseData(size(x,1))  #create an instance of caseData
    data.caseID = string.(x[:,1])
    data.tempData[:,1:2] = float(Array(x[:,2:3]))
    
    if (data.numCases != length(unique(data.caseID)))
        error("duplicate ID labels!")
    end
    if issorted(data.tempData[:,1]) == false
        error("Error: temporal data not sorted by infectious date!")
    end
    if minimum(data.tempData[:,2] - data.tempData[:,1]) < 0
        error("Error: removal time earlier than infectious time!")
    end
    println("...... Reading temporal data of ",data.numCases," cases ........ completed!")

    # SNP Difference Data
    x = Array(DataFrame(CSV.File(SNPfile)))
    data.SNPData = float(x)
    if (size(data.SNPData,1) != data.numCases) || (size(data.SNPData,2) != data.numCases)
        error("Error: SNP difference matrix size incorrect!")
    end
    if issymmetric(data.SNPData) == false
        error("Error: SNP difference matrix not symmetric!")
    end
    println("...... Reading SNP difference data ............... completed!")

    # network information
    if(!isempty(Contactfile))
        data.contactData = float(Array(DataFrame(CSV.File(Contactfile))))
        println("...... Reading contact data ..................... completed!")
    else
        println("...... No contact data provided!")
    end

    Initialization(genomeSize=genomeSize, itr_MCMC=itr_MCMC, burn_in=burn_in, subsample=subsample)
    println("...... Initialization ............................ completed!\n")

    # output File
    outputfile = isempty(outputfile) ? "parameter_bestnet.csv" : outputfile
    CSV.write(outputfile,DataFrame(reshape(append!(["iteration","logLikelihood","logPrior","theta","mutation_rate","infection_rate"],data.caseID, data.caseID .* "_InfTime"), 1, :), :auto),writeheader = false)

    # run MCMC
    mcmcAlgorithm(outputfile)

    return "The MCMC is completed!"
end


