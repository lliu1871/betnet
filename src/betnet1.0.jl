#Update theta (effective population size)
function thetaUpdate(theta::Float64, JumpingWidth::Float64, theta_lb::Float64=0.0)
     
    theta_new = rand(Uniform(max(theta_lb, theta-JumpingWidth), theta+JumpingWidth))
    length = theta+JumpingWidth-max(theta_lb, theta-JumpingWidth)
    return(theta_new=theta_new,length=length)
end

#Update tb (latent period)
function tbUpdate(tb::Float64, tb_lb::Float64, tb_ub::Float64, JumpingWidth::Float64)

    tb_new = rand(Uniform(max(tb_lb, tb-JumpingWidth), min(tb+JumpingWidth,tb_ub)))
    length = tb+JumpingWidth-max(tb_lb, tb-JumpingWidth)
    return(tb_new=tb_new,length=length)
end

#Update mu (mutation rate)
function muUpdate(mu::Float64, mu_lb::Float64, JumpingWidth::Float64)

    mu_new = rand(Uniform(max(mu_lb, mu-JumpingWidth), mu+JumpingWidth))
    length = mu+JumpingWidth-max(mu_lb, mu-JumpingWidth)
    return(mu_new=mu_new,length=length)
end

#Update alpha (infection rate)
function infrateUpdate(infrate::Float64,inf_lb::Float64, JumpingWidth::Float64)
    
    inf_new = rand(Uniform(max(inf_lb, infrate-JumpingWidth), infrate+JumpingWidth))
    length = infrate+JumpingWidth-max(inf_lb, infrate-JumpingWidth)
    return(inf_new=inf_new,length=length)
end



#Initialize the parameters

function Initialization(Pevent::Matrix{Float64},n::Int64,input::Matrix{Float64}, D::Matrix{Float64}, mu::Float64, N::Float64, inf_rate::Float64,
    Jump_theta::Float64,Jump_mu::Float64,Jump_inf::Float64,
    itr_MCMC::Int64=100000,burn_in::Int64=5000)

    # initialize the parameters
    num_out = Int64(floor((itr_MCMC - burn_in)/100)) #vector/list length of each output
    Net_infIDout = zeros(Int64, num_out, n) #output for the most likely individual who infected individual i (i is the sequence/row number)
    Param_out = zeros(Float64, num_out, 4) #output to store loglikelihood & theta
    Tb_out = zeros(Float64, n, n, num_out) #output of tb (in matrix) corresponding to each plausible transmission event.
    ProbNet_event_out = zeros(Float64, num_out, n) #output of ProbNet_event (in every 250 iterations)

    ProbNet_event = zeros(Float64, n) #the transmission probability of the most likely transmission event that each individual got infected.
    Net_infID = zeros(Int64, n) #the most likely individual who infected individual i (i is the sequence/row number)
    TRATIA_Vec= zeros(Float64, n) #vector to store the infectious period of individuals
    Child_Vec= zeros(Int64, n)
    Pro_Phi_Vec = zeros(Float64, n)
    TAB_Mtx= zeros(Float64, n, n) #matrix to store whether the symptom onset time time difference between each pairs of individuals.
    TRATIB_Mtx = zeros(Float64, n, n)
    Tb_sim = zeros(Float64, n, n) #matrix to store the tb values (latent period) of whom got infected in the corresponding transmission event.
    Tb_Proposal_Dis = zeros(Float64, n, n) #matrix to store the probability of proposal distribution of tb (truncated chi-square) corresponding 
    # to each plausible transmission event.
    T1_Mtx = zeros(Float64, n, n) #matrix to store t1 corresponding to each plausible transmission event
    Jump_tb_Mtx = zeros(Float64, n, n)
    Tb_lb = zeros(Float64, n, n)
    Prob_prod = zeros(Float64, n, n)
    TRA_TIB_TIA_min =  zeros(Float64, n, n)
    scal = 2/28 #2/28 is within 180 - 270 days (0.5 - 0.75 years). In other words, the expectation of latent period lies within 6 - 9 months
    dum_term=zeros(Int64,1)
   
    Pro_Tb = zeros(Float64, n, n) #matrix to store the pdf of given tb corresponding to each plausible transmission event.
    #t_star_sim_all = Matrix{Vector{Float64}}(undef,n,n)

    P_Mtx=zeros(Float64, n, n) #save the probability different site nucleotide given the same lineage
    Likeli_Mtx=zeros(Float64, n, n) #save the probability different site nucleotide given the same lineage
    Log_Likeli_Mtx=zeros(Float64, n, n)



    # 1)simulate initial value of theta 2)calculate the pdf of selected theta 3) calculate the range 

    theta_lb = 1E-6 # lower bound of theta
    theta_ub = 5E-6 # upper bound of theta
    theta= rand(Uniform(theta_lb, theta_ub))
    Pro_theta = pdf(Exponential(1), theta) #pdf of given theta
    theta_length= theta+Jump_theta-max(0, theta-Jump_theta)

     
    # 1)simulate initial value of mu 2)calculate the pdf of selected mu 3) calculate the range
    mu_lb = 5E-7 # lower bound of mu
    mu_ub = 1E-6# upper bound of mu
    mu=rand(Uniform(mu_lb, mu_ub))
    Pro_mu = pdf(Exponential(1/10), mu) 
    mu_length=mu+Jump_mu-max(0, mu-Jump_mu)

    # calculate the Likelihood at starting point. Calculate and store Fact_D & InfTime_diff at the same time.
    # Consider individual j was infected by individual i.
    Tb_sim[1,1] = 0.57 # mean of expected latent period.
    Net_infID[1] = 1
    TRATIA_Vec[1]=input[1,3]-input[1,2]

    for j in 2:n
        TIB = input[j,2] #symptom on set time of individual j
        TRB = input[j,3] #removal time of indiviudal j
        SNP = 1000

        TRATIA_Vec[j] = TRB - TIB

        for i in 1:(j-1) 
            TIA = input[i,2] #symptom on set time of individual i
            #Fact_D[i,j] = factorialSNP(D[i,j],N) #Calculate the coefficients of I1 and I2 based on SNP data
            TRA = input[i,3] #removal time of individual i
            TAB_Mtx[i,j] = TIB-TIA   #time period from sympton onset time of ind i to symptom onset time of ind j
            TRATIB_Mtx[i,j] = TIB-TRA # time period from the removal time of ind i to symptom onset time of ind j
            
            println("i ",i)
            println("j ",j)
            println("TAB_Mtx[i,j]",TAB_Mtx[i,j])
            println("TRATIB_Mtx[i,j]",TRATIB_Mtx[i,j])
            if TRATIB_Mtx[i,j] <= 2 && TAB_Mtx[i,j] > 0   #1.The gap from removel time of i to infection time of j should less than 2   2.infection time of j is later than infection time of i
                #calculate the maximum latent period that individual B could have
                tAB = TAB_Mtx[i,j] #the upper limit of tb and t_star   #from infection time of i to infection time of j
                TRA_TIB_TIA_min[i,j] = min(TRA-TIA, tAB)     #min of length infectious period or gap between two onset time
                
                #fact = Fact_D[i,j]
                d = D[i,j]

                if d < SNP
                    SNP = d
                    Net_infID[j] = i
                end


                Tb_lb[i,j] = max(0,TIB-TRA)   
                tb_l = Tb_lb[i,j]#the lower limit of tb
                Jump_tb_Mtx[i,j] = (tAB-tb_l)/5.0     

                Tb_sim[i,j] = tb_l+(tAB-tb_l)*rand()   #generated from Uniform(tb lower bound, tb upper bound)
                tb = Tb_sim[i,j]
                
                Tb_Proposal_Dis[i,j] = cdf(Chisq(8), tAB/scal)-cdf(Chisq(8), tb_l/scal)
                Pro_Tb[i,j] = pdf(Chisq(8), tb/scal)    

                T1_Mtx[i,j] = TRA+TRB+2*tb-2*TIB   #t1 is t1+t2 in pdf
                t1 = T1_Mtx[i,j]

                #compute the transition probability when two sequence from same lineage 
                P_Mtx[i,j]=3/4-3/4*exp(-mu*t1)/(2*theta+1)
                p = P_Mtx[i,j]
                
                Likeli_Mtx[i,j]= p^d *(1-p)^(N-d)
                likeli = Likeli_Mtx[i,j]

                Log_Likeli_Mtx[i,j]=log(likeli)
                logLike = Log_Likeli_Mtx[i,j]

                Prob_prod[i,j] = max(-1E10,logLike+Pevent[i,j]+log(Pro_Tb[i,j])-log(Tb_Proposal_Dis[i,j]))              
            end        
           
        end

        println("j ",j)
        println("Net_infID[j] ",Net_infID[j])
        ProbNet_event[j] = Prob_prod[Net_infID[j],j]

    end


    #estimate Beta
    beta=1/ mean(TRATIA_Vec[1:end-1])
    println("beta",beta)
  
    #count how many children each patient has of the given transmission tree
    for k in 1:n
        Child_Vec[k] = count(i->(i== k), Net_infID[2:end])
    end

    for z in 1: (n-1)
        Pro_Phi_Vec[z] = log(pdf(Geometric(beta/(inf_rate+beta)),Child_Vec[z]))
    end
    
    Pro_Phi= sum(Pro_Phi_Vec)

    println("Pro_Phi_Vec ",Pro_Phi_Vec)
    println("Pro_Phi ",Pro_Phi)
    

    # 1)simulate initial value of infection rate 2)calculate the pdf of selected infection rate 3) calculate the range

    inf_lb = 2 # lower bound of inf
    inf_ub = 6 # upper bound of inf
    inf_rate=rand(Uniform(inf_lb, inf_ub))
    alpha_rate = 4.0
    Pro_inf =  pdf(Exponential(alpha_rate), inf_rate)
    
    inf_length=inf_rate+Jump_inf-max(0, inf_rate-Jump_inf)

    results=(theta=theta,mu=mu,Net_infID=Net_infID,ProbNet_event=ProbNet_event,TAB_Mtx=TAB_Mtx,T1_Mtx=T1_Mtx,
    TRATIA_Vec=TRATIA_Vec,beta=beta,Child_Vec=Child_Vec,Pro_Phi_Vec=Pro_Phi_Vec,Pro_Phi=Pro_Phi,
    Pro_Tb=Pro_Tb,Tb_Proposal_Dis=Tb_Proposal_Dis,Pro_theta=Pro_theta,Pro_mu=Pro_mu,Pro_inf=Pro_inf,TRATIB_Mtx=TRATIB_Mtx,
    Tb_sim=Tb_sim,Tb_lb=Tb_lb,Jump_tb_Mtx=Jump_tb_Mtx,scal=scal,dum_term=dum_term,TRA_TIB_TIA_min=TRA_TIB_TIA_min,Prob_prod=Prob_prod,
    Net_infIDout=Net_infIDout,Tb_out =Tb_out,ProbNet_event_out=ProbNet_event_out,
    Param_out=Param_out,P_Mtx=P_Mtx,Likeli_Mtx=Likeli_Mtx,
    Log_Likeli_Mtx=Log_Likeli_Mtx,theta_length=theta_length,mu_length=mu_length,inf_length=inf_length)

    return results

end


# Update theta
function thetaUpdate_itr(n::Int64,Pro_theta::Float64,theta::Float64,
    Net_infID::Vector{Int64},ProbNet_event::Vector{Float64},Prob_prod::Matrix{Float64},
    TAB_Mtx::Matrix{Float64}, TRATIB_Mtx::Matrix{Float64},T1_Mtx::Matrix{Float64},
    Pro_Tb::Matrix{Float64},Tb_Proposal_Dis::Matrix{Float64},
    D::Matrix{Float64}, mu::Float64, N::Float64, Jump_theta::Float64,Tb_sim::Matrix{Float64},
    P_Mtx::Matrix{Float64},Likeli_Mtx::Matrix{Float64},
    Log_Likeli_Mtx::Matrix{Float64},
    Child_Vec::Vector{Int64},Pro_Phi_Vec::Vector{Float64},
    Pro_Phi::Float64,beta::Float64,inf_rate::Float64,TRATIA_Vec::Vector{Float64},
    theta_length::Float64,
    dum_term::Vector{Int64},
    Pevent=0)
     
    
    ## update theta
    sim=thetaUpdate(theta, Jump_theta, 0.0)
    theta_new = sim.theta_new
    theta_length_new=sim.length

    Pro_theta_new = pdf(Exponential(1), theta_new)
 
    Prob_prod_new = zeros(Float64, n, n)
    ProbNet_event_new = zeros(Float64, n)
    

    P_Mtx_new = zeros(Float64, n,n)
    Likeli_Mtx_new = zeros(Float64, n,n)
    Log_Likeli_Mtx_new = zeros(Float64, n,n)


    #loglike_new = zeros(Float64, n)

    for j in 2:n
        #println("infect",j)
        i = Int64(Net_infID[j]) # individual i infected individual j in the last loop
        
        if ProbNet_event[j] != 0

            #tAB = TAB_Mtx[i,j]
            t1 = T1_Mtx[i,j]
            d = D[i,j] 
            


            #compute the transition probability when two sequence from same lineage 
            P_Mtx_new[i,j]=3/4-3/4*exp(-mu*t1)/(2*theta_new+1)
            p = P_Mtx_new[i,j]
            
            Likeli_Mtx_new[i,j]= p^d *(1-p)^(N-d)
            likeli = Likeli_Mtx_new[i,j]

            Log_Likeli_Mtx_new[i,j]= log(likeli)
            logLike = Log_Likeli_Mtx_new[i,j]


           #loglike_new[j]=logLike
           ProbNet_event_new[j] = max(-1E10,logLike+Pevent[i,j]+log(Pro_Tb[i,j])-log(Tb_Proposal_Dis[i,j]))

           Prob_prod_new[i,j] = ProbNet_event_new[j]

          
        end
    end
    
    #log ratio of the proposal density in two directions 
    log_ratio_pd=log((1/theta_length)/(1/theta_length_new))

    H = min(1, exp(sum(ProbNet_event_new)+log(Pro_theta_new)-sum(ProbNet_event)-log(Pro_theta)+log_ratio_pd))
    
    println("H",H)


    if H > rand()
        println("accept")

        theta = theta_new
        theta_length=theta_length
        ProbNet_event = ProbNet_event_new 
        Prob_prod = Prob_prod_new
        Pro_theta = Pro_theta_new 
        P_Mtx=P_Mtx_new
        Likeli_Mtx=Likeli_Mtx_new
        Log_Likeli_Mtx = Log_Likeli_Mtx_new



        for j = 2:n
           for q in 1:(j-1)
              if q != Int64(Net_infID[j])
                   if TRATIB_Mtx[q,j]<=2 && TAB_Mtx[q,j] > 0

                      #tAB = TAB_Mtx[q,j]
                      t1 = T1_Mtx[q,j]
                      
                      d = D[q,j]
                      
                       P_Mtx[q,j] = 3/4-3/4*exp(-mu*t1)/(2*theta+1)
                       p= P_Mtx[q,j]
                       Likeli_Mtx[q,j]= p^d *(1-p)^(N-d)
                       likeli = Likeli_Mtx[q,j]
                       
                       Log_Likeli_Mtx[q,j] = log(likeli)
                       logLike = Log_Likeli_Mtx[q,j]
                       Prob_prod[q,j] = max(-1E10,logLike+Pevent[q,j]+log(Pro_Tb[q,j])-log(Tb_Proposal_Dis[q,j]))
                   end 
              end    
           end
        end 
    end

    result=(theta = theta,Prob_prod = Prob_prod,ProbNet_event = ProbNet_event,Pro_theta = Pro_theta,
    P_Mtx=P_Mtx,Likeli_Mtx=Likeli_Mtx,Log_Likeli_Mtx=Log_Likeli_Mtx,theta_length=theta_length,dum_term=dum_term)
    return result
end


#Update mutation rate
function muUpdate_itr(n::Int64,Net_infID::Vector{Int64},ProbNet_event::Vector{Float64},Prob_prod::Matrix{Float64},
    TAB_Mtx::Matrix{Float64}, TRATIB_Mtx::Matrix{Float64},T1_Mtx::Matrix{Float64},
    Pro_Tb::Matrix{Float64},Tb_Proposal_Dis::Matrix{Float64},
    D::Matrix{Float64}, theta::Float64,mu::Float64, N::Float64,Jump_mu::Float64,
    P_Mtx::Matrix{Float64},Likeli_Mtx::Matrix{Float64},Log_Likeli_Mtx::Matrix{Float64},
    Child_Vec::Vector{Int64},Pro_Phi_Vec::Vector{Float64},
    Pro_Phi::Float64,beta::Float64,inf_rate::Float64,TRATIA_Vec::Vector{Float64},
    Pro_mu::Float64,mu_length::Float64,
    dum_term::Vector{Int64},Pevent=0)
    

    ## update mu
    sim=muUpdate(mu, 0.0,Jump_mu)
    mu_new = sim.mu_new
    mu_length_new=sim.length
    Pro_mu_new = pdf(Exponential(1/10), mu_new)
    #Pro_mu_new = pdf(Exponential(1), mu_new)

    ProbNet_event_new = zeros(Float64, n)
    Prob_prod_new = zeros(Float64, n, n)
    
    P_Mtx_new = zeros(Float64, n,n)
    Likeli_Mtx_new = zeros(Float64, n,n)
    Log_Likeli_Mtx_new = zeros(Float64, n, n)

    ProbNet_event_old =zeros(Float64, n)

    for j in 2:n

        i = Int64(Net_infID[j]) 

        if ProbNet_event[j] != 0

            t1 = T1_Mtx[i,j]

            d = D[i,j] 

            #compute the transition probability when two sequence from same lineage 
            P_Mtx_new[i,j]=3/4-3/4*exp(-mu_new*t1)/(2*theta+1)
            p = P_Mtx_new[i,j]
            
            Likeli_Mtx_new[i,j]= p^d *(1-p)^(N-d)
            likeli = Likeli_Mtx_new[i,j]

            Log_Likeli_Mtx_new[i,j] = log(likeli)
            logLike_new = Log_Likeli_Mtx_new[i,j]

            ProbNet_event_new[j] = max(-1E10,logLike_new+Pevent[i,j]+log(Pro_Tb[i,j])-log(Tb_Proposal_Dis[i,j]))

            Prob_prod_new[i,j] = ProbNet_event_new[j]

            ProbNet_event_old[j] = max(-1E10,Log_Likeli_Mtx[i,j]+Pevent[i,j]+log(Pro_Tb[i,j])-log(Tb_Proposal_Dis[i,j]))
           
        end
    end

    #log ratio of the proposal density in two directions 
    log_ratio_pd=log((1/mu_length)/(1/mu_length_new))

    H = min(1, exp(sum(ProbNet_event_new)+log(Pro_mu_new)-sum(ProbNet_event)-log(Pro_mu)+log_ratio_pd))
    println("H ",H)



    if H > rand()
        
        println("accept")

        mu = mu_new
        mu_length=mu_length_new
        Pro_mu=Pro_mu_new
        ProbNet_event = ProbNet_event_new 
        Prob_prod = Prob_prod_new

        P_Mtx = P_Mtx_new
        Likeli_Mtx = Likeli_Mtx_new

        #Sum_temp_Mtx = Sum_temp_Mtx_new
        Log_Likeli_Mtx = Log_Likeli_Mtx_new

        
        for j = 2:n
           for q in 1:(j-1)
              if q != Int64(Net_infID[j])
                   if TRATIB_Mtx[q,j]<=2 && TAB_Mtx[q,j] > 0

                      #update P_sl_Mtx & Likeli_sl_Mtx
                      t1 = T1_Mtx[q,j]

                      d = D[q,j]

                      P_Mtx[q,j]=3/4-3/4*exp(-mu*t1)/(2*theta+1)
                      p = P_Mtx[q,j]
            
                      Likeli_Mtx[q,j]= p^d *(1-p)^(N-d)
                      likeli = Likeli_Mtx[q,j]

                      Log_Likeli_Mtx[q,j] = log(likeli)
                      logLike_new = Log_Likeli_Mtx[q,j]

                      Prob_prod[q,j] = max(-1E10,logLike_new+Pevent[q,j]+log(Pro_Tb[q,j])-log(Tb_Proposal_Dis[q,j]))


                   end 
              end    
           end
        end 
    end
    result=(mu = mu,Prob_prod = Prob_prod, ProbNet_event = ProbNet_event,P_Mtx = P_Mtx,
            Likeli_Mtx = Likeli_Mtx,Log_Likeli_Mtx=Log_Likeli_Mtx,mu_length=mu_length,Pro_mu=Pro_mu,
            dum_term=dum_term)
    return result

end


#Update infection rate
function infrateUpdate_itr(n::Int64,inf_rate::Float64,Net_infID::Vector{Int64},ProbNet_event::Vector{Float64},
    Prob_prod::Matrix{Float64},Pro_inf::Float64,beta::Float64, Child_Vec::Vector{Int64},
    TRATIA_Vec::Vector{Float64},Pro_Phi::Float64,Pro_Phi_Vec::Vector{Float64},
    Tb_sim::Matrix{Float64},TAB_Mtx::Matrix{Float64}, TRATIB_Mtx::Matrix{Float64},
    TRA_TIB_TIA_min::Matrix{Float64},
    Log_Likeli_Mtx::Matrix{Float64},
    Pro_Tb::Matrix{Float64},Tb_Proposal_Dis::Matrix{Float64},
    Jump_inf::Float64, inf_length::Float64,
    dum_term::Vector{Int64},Pevent=0)
    
    println("Pro_Phi ",Pro_Phi)


    ## update inf_rate
    sim=infrateUpdate(inf_rate,0.0,Jump_inf)
    inf_rate_new = sim.inf_new
    inf_length_new=sim.length

    Pro_inf_new = pdf(Exponential(4.0), inf_rate_new)

    Pro_Phi_Vec_new = zeros(Float64, n)
    
    #println("maxnum ",max_num)
    #println("Child_Vec ",Child_Vec)
    #println("inf_rate_new ",inf_rate_new)



    for l in 1: (n-1)
        Pro_Phi_Vec_new[l] = log(pdf(Geometric(beta/(inf_rate_new+beta)),Child_Vec[l]))
        #println("l ",l)
        #println("Child_Vec[l] ",Child_Vec[l])
        #println("Pro_Phi_Vec[l]  ",Pro_Phi_Vec[l] )
    end

    #println("Pro_Phi_Vec  ",Pro_Phi_Vec )
    Pro_Phi_new = sum(Pro_Phi_Vec_new)

    #log ratio of the proposal density in two directions 
    log_ratio_pd=log((1/inf_length)/(1/inf_length_new))
    H = min(1, exp(Pro_Phi_new+log(Pro_inf_new)-Pro_Phi-log(Pro_inf)+log_ratio_pd))


    if H > rand()
        println("accept inf")
        
        inf_rate = inf_rate_new
        inf_length=inf_length_new
        Pro_inf = Pro_inf_new
        Pro_Phi_Vec = Pro_Phi_Vec_new
        Pro_Phi = Pro_Phi_new

    end
    result=(inf_rate = inf_rate,Pro_inf = Pro_inf,Pro_Phi_Vec = Pro_Phi_Vec,Pro_Phi=Pro_Phi,inf_length=inf_length,dum_term=dum_term)
    
    return result
end


#Update transmission tree
function transUpdate_itr(n::Int64,TRATIB_Mtx::Matrix{Float64},TAB_Mtx::Matrix{Float64},Tb_sim::Matrix{Float64},
    Tb_lb::Matrix{Float64},Jump_tb_Mtx::Matrix{Float64},scal::Float64,inf_rate::Float64,
    TRA_TIB_TIA_min::Matrix{Float64},
    TRATIA_Vec::Vector{Float64},Pro_Phi::Float64,Pro_Phi_Vec::Vector{Float64},
    Child_Vec::Vector{Int64},beta::Float64,
    T1_Mtx::Matrix{Float64},D::Matrix{Float64},
    mu::Float64,N::Float64, 
    Tb_Proposal_Dis::Matrix{Float64},ProbNet_event::Vector{Float64},Prob_prod::Matrix{Float64},
    Pro_Tb::Matrix{Float64},Net_infID::Vector{Int64},Pevent::Matrix{Float64},
    P_Mtx::Matrix{Float64},Likeli_Mtx::Matrix{Float64},
    Log_Likeli_Mtx::Matrix{Float64},
    dum_term::Vector{Int64},
    theta::Float64)
    
    #println("Prob ",ProbNet_event)
    #println("Tb_lb ",Tb_lb)

   
    ## update one transmission event and its corresponding tb
    # now randomly change one transmission event
    m = rand(2:n)
    
    #println("select ID ",m) #new
    
    # randomly select one individual (before m) that infected m, making sure the individual that infected m is not the 
    # the one in the previous network.
    
    ID_pool = shuffle(1:(m-1))
    #println("ID pool ", ID_pool)
    
    for q in 1:(m-1)
        infID = ID_pool[q]
        if infID != Int64(Net_infID[m])
            if (TRATIB_Mtx[infID,m] <= 2 && TAB_Mtx[infID,m] > 0)

                Net_infID_new = zeros(Int64,n)

                ProbeventID_new = 0.0 
                tAB = TAB_Mtx[infID,m]
                tb = Tb_sim[infID,m]
                tb_l = Tb_lb[infID,m]
                #print("tb_l",tb_l)
                jump_tb = Jump_tb_Mtx[infID,m]
                
                sim= tbUpdate(tb, tb_l, tAB, jump_tb)
                tb_new=sim.tb_new
                tb_length_new=sim.length
                tb_length=tb+jump_tb-max(tb_l,tb-jump_tb)

                Pro_Tb_new = pdf(Chisq(8), tb_new/scal)

                t1_new = T1_Mtx[infID,m]-2*tb+2*tb_new   #TRA+TRB+2*tb-2*TIB
            
                d = D[infID,m]

                p = 3/4-3/4*exp(-mu*t1_new)/(2*theta+1)
                likeli = p^d *(1-p)^(N-d)
                logLike = log(likeli)

               ProbeventID_newtb = max(-1E10,logLike+Pevent[infID,m]+log(Pro_Tb_new)-log(Tb_Proposal_Dis[infID,m]))
               
               #println("ProbeventID_newtb",ProbeventID_newtb)

               ProbeventID = ProbNet_event[m]
               
               #log ratio of the proposal density in two directions 
               log_ratio_pd=log((1/tb_length)/(1/tb_length_new))
               rv2=rand()
               H2 = min(1,exp(ProbeventID_newtb-Prob_prod[infID,m]+log_ratio_pd))

               if H2 > rv2
                  #println("update tb", m)
                  #println("update ID", infID)
                  #println("Prob ", ProbeventID_newtb)
                  T1_Mtx[infID,m] = t1_new
                  Tb_sim[infID,m] = tb_new
                  Pro_Tb[infID,m] = Pro_Tb_new
                  Prob_prod[infID,m] = ProbeventID_newtb
                  P_Mtx[infID,m] = p
                  Likeli_Mtx[infID,m] = likeli
                  #T_birth_Mtx[infID,m] = t_birth
                  #E_nl_Mtx[infID,m] = E_nl
                  #Sum_temp_Mtx[infID,m] = sum_temp
                  Log_Likeli_Mtx[infID,m] = logLike
                  #t_star_sim_all[infID,m] =t_star_sim
                  ProbeventID_new = ProbeventID_newtb
                  #P_tstar_Mtx[infID,m] = P_tstar_val
                  
                  #println("update tb")
                  #println("update m",m)
                  #println("update ID",infID)
               else 
                  ProbeventID_new=Prob_prod[infID,m]
                  #println("update tb")
               end
              
              #println("before NetID ", Net_infID)
              Net_infID_new = copy(Net_infID)
              Net_infID_new[m] = infID

              Child_Vec_new = zeros(Int64,n)
              Pro_Phi_Vec_new = zeros(Float64,n)
              
              for k in 1:n
                  Child_Vec_new[k] = count(i->(i== k), Net_infID_new[2:end])
               end

               for l in 1: (n-1)
                   Pro_Phi_Vec_new[l] = log(pdf(Geometric(beta/(inf_rate+beta)),Child_Vec_new[l]))
               end
        

               Pro_Phi_new = sum(Pro_Phi_Vec_new)

              H3 = min(1, exp(ProbeventID_new-ProbeventID+Pro_Phi_new-Pro_Phi))
              #println("H3 ",H3)
              #println("H3_1 ",ProbeventID_new-ProbeventID)
              #println("old ProbeventID",ProbeventID)
              #println("H3_2 ",Pro_Phi_new-Pro_Phi)
             
              if H3>rand()   #new
                  #println("update inf ID", m)
                  #println("update inf ID to ", infID)
                  #println("Prob ", ProbeventID_new)
                  ProbNet_event[m] = ProbeventID_new
                  Net_infID[m] = infID
                  Child_Vec = Child_Vec_new
                  Pro_Phi_Vec = Pro_Phi_Vec_new
                  Pro_Phi = Pro_Phi_new
              end   
            end 
        end
    end

    println("update", m)

    #println("ProbNet ", ProbNet_event)
    #println("NetID ", Net_infID)

    results=(ProbNet_event=ProbNet_event,Net_infID=Net_infID,T1_Mtx=T1_Mtx,Tb_sim=Tb_sim,
    Pro_Tb=Pro_Tb,Prob_prod=Prob_prod,P_Mtx=P_Mtx,Likeli_Mtx=Likeli_Mtx,
    Log_Likeli_Mtx=Log_Likeli_Mtx,
    Child_Vec = Child_Vec,Pro_Phi_Vec = Pro_Phi_Vec,Pro_Phi = Pro_Phi,beta = beta,
    dum_term=dum_term)

    return results
end


#MCMC Algorithm
function TransNet(input::Matrix{Float64}, D::Matrix{Float64}, mu::Float64, N::Float64, 
    inf_rate::Float64, Jump_theta::Float64,Jump_inf::Float64,Jump_mu::Float64,
    itr_MCMC::Int64=100000,Pevent::Int64=0, burn_in::Int64=5000)
    # This function is aiming to find the most possible transmission network using temporal information: 
    # input, which is a n*3 matrix containing the individual ID, symptom on set time, and removal time of n infected individuals. 
    # D is a n*n matrix including the SNP data between each pair of infects. 
    # mu is the mutation rate (per site per day) of Tuberculosis. 
    # N is the number of pathogen base pairs. 
    # itr_t_star is the number of iterations that used to compute the MCMC integral of t_star. 
    # itr_MCMC is the number of iterationsof MCMC. 
    # burnin is the number of iterations that we discard (default is the first 5000 iterations).
    
    n = size(input, 1)
    # initialize 
    if Pevent == 0
        Pevent = zeros(Float64, n, n)
        for i in 2:n
            for j in 1:(i-1)
                Pevent[j,i] = log(1.0/float(i-1))
            end
        end
    end
    if size(input, 2) != 3
        return "Input data has incorrect format: there should be three columns containing individual ID, symptom onset time and removal time"
    elseif (n != length(unique(input[:,1])))
        return "incorrect ID labeling of infectious individuals (duplicates or incontinuous)"
    elseif size(Pevent, 1) != n | size(Pevent, 2) != n
        return "incorrect assignment of Pevent (Matrix of Prior Probability of transmission event )"
    else
        init=Initialization(Pevent,n,input,D,mu,N,inf_rate,Jump_theta,Jump_mu,Jump_inf,itr_MCMC,burn_in)
        
        theta=init.theta
        mu=init.mu
        Net_infID=init.Net_infID
        ProbNet_event=init.ProbNet_event
        TAB_Mtx=init.TAB_Mtx
        T1_Mtx=init.T1_Mtx
        TRATIA_Vec = init.TRATIA_Vec
        Child_Vec = init.Child_Vec
        beta = init.beta
        Pro_Phi_Vec = init.Pro_Phi_Vec
        Pro_Phi = init.Pro_Phi
        Pro_inf = init.Pro_inf
        Pro_Tb=init.Pro_Tb
        Tb_Proposal_Dis=init.Tb_Proposal_Dis
        Pro_theta=init.Pro_theta
        Pro_mu=init.Pro_mu
        TRATIB_Mtx=init.TRATIB_Mtx
        Tb_sim=init.Tb_sim
        Tb_lb=init.Tb_lb
        Jump_tb_Mtx=init.Jump_tb_Mtx
        scal=init.scal
        dum_term=init.dum_term
        TRA_TIB_TIA_min=init.TRA_TIB_TIA_min
        Prob_prod=init.Prob_prod
        Net_infIDout=init.Net_infIDout

        #Direct_infout = init.Direct_infout
        Tb_out =init.Tb_out 
        ProbNet_event_out=init.ProbNet_event_out
        Param_out=init.Param_out
        Tb_sim = init.Tb_sim

        #t_star_sim_all = init.t_star_sim_all
        P_Mtx = init.P_Mtx
        Likeli_Mtx = init.Likeli_Mtx
        Log_Likeli_Mtx = init.Log_Likeli_Mtx
        theta_length=init.theta_length
        mu_length=init.mu_length
        inf_length=init.inf_length
        #println("Tb_lb",Tb_lb)
        
        initial_theta=theta
        initial_mu=mu
        for k in 2:itr_MCMC 
            println("itr",k)
            #nonzero_index = (ProbNet_event .!= 0) #get the index of individual possess non-zeros transmission probability
            #Prob_net_nonzero = ProbNet_event[nonzero_index] #get the maximum non-zeros probability of transmission event of each individual.
            H = 0.0 # store hasting ratio
            
            r_random = rand()
        
            if r_random <= 0.3

                ## update theta
                println("update theta")
                theta_itr =thetaUpdate_itr(n,Pro_theta,theta,Net_infID,ProbNet_event,Prob_prod,TAB_Mtx,TRATIB_Mtx,T1_Mtx,Pro_Tb,Tb_Proposal_Dis,
                D, mu, N,Jump_theta,Tb_sim,P_Mtx,Likeli_Mtx,Log_Likeli_Mtx,Child_Vec,Pro_Phi_Vec,
                Pro_Phi,beta,inf_rate,TRATIA_Vec,theta_length,dum_term,Pevent)


                theta = theta_itr.theta
                #t_star_sim_all = theta_itr.t_star_sim_all
                Prob_prod = theta_itr.Prob_prod
                ProbNet_event = theta_itr.ProbNet_event
                Pro_theta = theta_itr.Pro_theta
                #E_nl_Mtx = theta_itr.E_nl_Mtx
                #Sum_temp_Mtx = theta_itr.Sum_temp_Mtx
                P_Mtx=theta_itr.P_Mtx
                Likeli_Mtx=theta_itr.Likeli_Mtx
                Log_Likeli_Mtx = theta_itr.Log_Likeli_Mtx
                theta_length=theta_itr.theta_length
                dum_term = theta_itr.dum_term

                println("First ProbNet ",ProbNet_event[1])

            elseif r_random >0.3 && r_random <=0.4
                println("update mu")

                mu_itr =muUpdate_itr(n,Net_infID,ProbNet_event,Prob_prod,TAB_Mtx,TRATIB_Mtx,T1_Mtx,Pro_Tb,Tb_Proposal_Dis,
                D, theta,mu, N, Jump_mu,P_Mtx,Likeli_Mtx,Log_Likeli_Mtx,Child_Vec,Pro_Phi_Vec,
                Pro_Phi,beta,inf_rate,TRATIA_Vec,Pro_mu,mu_length,dum_term,Pevent)

                mu = mu_itr.mu
                Prob_prod = mu_itr.Prob_prod
                ProbNet_event = mu_itr.ProbNet_event
                P_Mtx = mu_itr.P_Mtx
                Likeli_Mtx = mu_itr.Likeli_Mtx
                #Sum_temp_Mtx = mu_itr.Sum_temp_Mtx
                Log_Likeli_Mtx = mu_itr.Log_Likeli_Mtx
                Pro_mu=mu_itr.Pro_mu
                mu_length=mu_itr.mu_length
                dum_term = mu_itr.dum_term

                #println("First ProbNet ",ProbNet_event[1])
            
            elseif r_random >0.5 && r_random <=0.8
                println("update infrate")
                infrate_itr =infrateUpdate_itr(n,inf_rate,Net_infID,ProbNet_event,Prob_prod,Pro_inf,beta,Child_Vec,
                TRATIA_Vec,Pro_Phi,Pro_Phi_Vec,Tb_sim,TAB_Mtx,TRATIB_Mtx,TRA_TIB_TIA_min,Log_Likeli_Mtx,Pro_Tb,Tb_Proposal_Dis,Jump_inf,inf_length,dum_term,Pevent)

                inf_rate = infrate_itr.inf_rate
                Pro_Phi_Vec  =infrate_itr.Pro_Phi_Vec
                Pro_Phi =infrate_itr.Pro_Phi
                Pro_inf = infrate_itr.Pro_inf
                inf_length=infrate_itr.inf_length
                dum_term = infrate_itr.dum_term

            else
                ## update one transmission event and its corresponding tb
                # now randomly change one transmission event
                # randomly select one individual (before m) that infected m, making sure the individual that infected m is not the 
                # the one in the previous network.
                println("update tran")
                tran_itr=transUpdate_itr(n,TRATIB_Mtx,TAB_Mtx,Tb_sim,Tb_lb,Jump_tb_Mtx,scal,inf_rate,TRA_TIB_TIA_min,
                TRATIA_Vec,Pro_Phi,Pro_Phi_Vec,Child_Vec,beta,T1_Mtx,D,mu,N,Tb_Proposal_Dis,ProbNet_event,Prob_prod,
                Pro_Tb,Net_infID,Pevent,P_Mtx,Likeli_Mtx,Log_Likeli_Mtx,dum_term,theta)
                
                ProbNet_event=tran_itr.ProbNet_event
                Net_infID=tran_itr.Net_infID
                T1_Mtx=tran_itr.T1_Mtx
                Tb_sim=tran_itr.Tb_sim
                Pro_Tb=tran_itr.Pro_Tb
                Prob_prod=tran_itr.Prob_prod
                P_Mtx = tran_itr.P_Mtx
                Likeli_Mtx = tran_itr.Likeli_Mtx
                #T_birth_Mtx =tran_itr.T_birth_Mtx
                #E_nl_Mtx = tran_itr.E_nl_Mtx
                #Sum_temp_Mtx = tran_itr.Sum_temp_Mtx
                Log_Likeli_Mtx = tran_itr.Log_Likeli_Mtx
                #P_tstar_Mtx = tran_itr.P_tstar_Mtx
                Child_Vec = tran_itr.Child_Vec
                Pro_Phi_Vec = tran_itr.Pro_Phi_Vec
                Pro_Phi = tran_itr.Pro_Phi
                #max_num = tran_itr.max_num
                dum_term = tran_itr.dum_term

            end  


            if k > burn_in && k%100 == 0
                h = Int64(ceil((k-burn_in)/100))
                Net_infIDout[h,:] = Net_infID
                #Direct_infout[h,:] = DirectTran_event
                Tb_out[:,:,h] = Tb_sim
                ProbNet_event_out[h,:] = ProbNet_event
                Param_out[h,1] = sum(ProbNet_event)+log(Pro_theta)+Pro_Phi+log(Pro_inf)
                Param_out[h,2] = theta
                Param_out[h,3] = mu
                Param_out[h,4] = inf_rate

            end

        end
       
        println("intial theta ",initial_theta)
        println("intial mu ",initial_mu)
        println("beta ",beta)

        return Net_infIDout, Tb_out, ProbNet_event_out, Param_out,
               Pro_theta,Prob_prod,
               TAB_Mtx,TRATIB_Mtx,T1_Mtx,TRA_TIB_TIA_min,
               P_Mtx,Likeli_Mtx,
               Pevent,Pro_Tb,Tb_Proposal_Dis,
               Jump_mu,Jump_theta,Jump_inf,Jump_tb_Mtx,
               Tb_lb,scal
    end
end


# This is the code to store Tb
function Tborg(x, pos, infpos)
    k = size(x)[3]
    l = size(pos)[1]
    LatentP_out = zeros(Float64, k,l)

    for i in 1:k
        for j in 1:l
            LatentP_out[i,j] = x[infpos[j], pos[j], i]
        end
    end
    return LatentP_out
end


