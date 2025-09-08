# Bayesian Estimation of Transmission NETworks 
#  <center>(betnet v1.1)</center>
This Julia package implements a Bayesian model for reconstructing transmission networks in infectious diseases by integrating genomic and temporal data.

## Installation
### Local installation
1. Download the package from GitHub
2. In the package mode, run "dev path-to-betnet" to add the package in Julia.

### Installation from Github
In the package mode, run "add https://github.com/lliu1871/betnet" to add the package in Julia.

### Global installation (unavailable)
In the package mode, run "add betnet" to add the package in Julia.

## How to run the program
The input datasets time_real_data.csv (temporal data) and SNP_real_data.csv (SNP distances) are available in the folder "test". To analyze the datasets, run the following Julia code

    using betnet
    using DataFrames
    using CSV
    using Test

    #Load the temporal and Genomic Data
    #Test on 69 Patients

    #Temporal Data
    testDDD = Array(DataFrame(CSV.File("./test/time_real_data.csv")))
    testDD = testDDD[:,3:5]
    testD = zeros(Float64,length(testDDD[:,1]),3)
    testD[:,2:3] = float(testDD[:,2:3])/365
    testD[:,1] = float(testDD[:,1])

    #SNP Difference Data
    DDD = Array(DataFrame(CSV.File("./test/SNP_real_data.csv")))
    DD = DDD[:,2:70]
    D = float(DD)

    #MCMC
    @time rd_ppp_111_1 = TransNet(testD, D, 1.0e-6, 4411532.0, 2.0, 5E-7,0.3,5e-7, 100000,0,0)

    inf_Net_rd_ppp_111_1 = DataFrame(rd_ppp_111_1[1],:auto)
    Param_rd_ppp_111_1 = DataFrame(rd_ppp_111_1[4],:auto)
    Tb_rd_ppp_111_1 = rd_ppp_111_1[2]

    #Here input the inferred prediction ID to get the corresponding tb values
    Latent_rd_ppp_111_1 = DataFrame(Tborg(Tb_rd_ppp_111_1,[2:1:69;], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 11, 6, 7, 2, 1, 1, 13, 1, 19, 1, 16, 1, 1, 22, 7, 7, 27, 1, 1, 5, 5, 27, 27, 27, 27, 18, 37, 36, 34, 34, 37, 1, 27, 5, 1, 1, 1, 34, 22, 46, 27, 14, 27, 49, 36, 5, 27, 48, 23, 25, 38, 35, 1, 53, 19, 27, 30, 27]), :auto)

    CSV.write("./test/inf_NetID_rd1111.csv", inf_Net_rd_ppp_111_1)
    CSV.write("./test/Para_Net_rd1111.csv", Param_rd_ppp_111_1)
    CSV.write("./test/Latent_rd1111.csv", Latent_rd_ppp_111_1)

    #Simulation Data
    testDDD = Array(DataFrame(CSV.File("./data/TemporalData1_100.csv")))
    testD = float(testDDD)

    #SNP Difference Data
    DDD = Array(DataFrame(CSV.File("./data/SNP1_100.csv")))
    D = float(DDD)

    CCC = Array(DataFrame(CSV.File("./data/ContactProb1_100.csv")))
    C = float(CCC)

    #run for model with network
    @time pppfull = TransNet_N(testD, D, C,  1000000.0,2.0, 1E-6,0.3,5e-7,100000,0)

    #uncomment the next row to run for model without network
    @time pppfull=TransNet(testD, D, 1.0e-6, 1000000.0, 2.0, 1E-6,0.3,5e-7,100000,0,0)   

    inf_Net_ppp_1 = DataFrame(pppfull[1],:auto)
    Param_ppp_1 = DataFrame(pppfull[4],:auto)
    Tb_ppp_1 = pppfull[2]

    #based on the inferred ID to get the correspond latent period
    #Fill in the list in the last bracket; currently filled with [1,2,1,3,.....]

    Latent_rd_ppp_1 = DataFrame(Tborg(Tb_ppp_1,[2:1:99;], [1, 2, 1, 3, 1, 4, 1, 1, 4, 1, 6, 4, 7, 14, 11, 12, 6, 17, 18, 7, 4, 13, 6, 3, 21, 17, 2, 18, 19, 20, 1, 30, 20, 9, 35, 21, 4, 31, 28, 24, 40, 1, 41, 4, 14, 38, 20, 18, 43, 8, 34, 36, 47, 25, 29, 27, 47, 54, 49, 22, 18, 14, 2, 62, 34, 52, 26, 53, 57, 42, 68, 19, 52, 57, 73, 57, 57, 7, 40, 21, 25, 79, 37, 42, 82, 57, 16, 75, 33, 26, 57, 48, 70, 33, 62, 14, 34, 97]), :auto)


    #output the data
    CSV.write("./test/IMPACT_inf_NetID_Test.csv", inf_Net_ppp_1)
    CSV.write("./test/IMPACT_Para_Net_Test.csv", Param_ppp_1)
    CSV.write("./test/Latent_1.csv", Latent_rd_ppp_1)

## Resources
- Support: We use GitHub for the development of the Julia package Distributions itself. For support and questions, please use the Julia Discourse forum. Also, for casual conversation and quick questions, there are the channels #helpdesk and #statistics in the official Julia chat (https://julialang.slack.com). To get an invitation, please visit https://julialang.org/slack/.



## Contributing
### Workflow with Git and GitHub
To contribute to the package, fork the repository on GitHub, clone it and make modifications on a new branch, do not commit modifications on master. Once your changes are made, push them on your fork and create the Pull Request on the main repository.

### Citation
Xu, J., Hu, H., Ellison, G., Yu, L., Whalen, C.C., Liu, L. Bayesian estimation of transmission networks for infectious diseases. J. Math. Biol. 2025, 90(29) [DOI link](https://doi.org/10.1007/s00285-025-02193-1).





