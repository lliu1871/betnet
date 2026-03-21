

using betnet
@time transNetworkInference(tempfile="./data/TemporalData1_100.csv",SNPfile="./data/SNP1_100.csv",Contactfile="./data/ContactProb1_100.csv",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_net.csv")

