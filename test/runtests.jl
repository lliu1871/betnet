using bestNet
using Distributions
using Test

#@testset "betnet.jl" begin
    #Load the temporal and Genomic Data
    @time transNetworkInference(tempfile="./data/sim_temporal.csv",SNPfile="./data/sim_snp.csv",Contactfile="",genomeSize=100000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_sim_50.csv")
    
    # simulation data with/without contact network
    @time transNetworkInference(tempfile="./bestNet/data/TemporalData1_100.csv",SNPfile="./betnet2.0/data/SNP1_100.csv",Contactfile="",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_nonet.csv")
    @time transNetworkInference(tempfile="./bestNet/data/TemporalData1_100.csv",SNPfile="./betnet2.0/data/SNP1_100.csv",Contactfile="./betnet2.0/data/ContactProb1_100.csv",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_net.csv")
 
#end

cd("./bestNet/data/")
#for i in 1:10
    tempfile = "sim_temporal" * string(i) * ".csv"
    snpfile = "sim_snp" * string(i) * ".csv"
    outputfile = "parameter" * string(i) * ".csv"
    snpfile = "./bestNet/data/sim_snp.csv"
    tempfile = "./bestNet/data/sim_temporal.csv"
    outputfile = "./bestNet/data/parameter.csv"
    transNetworkInference(tempfile=tempfile,SNPfile=snpfile,Contactfile="",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile=outputfile)
#end