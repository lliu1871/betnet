# Bayesian Estimation of Transmission NETworks (betnet)
This Julia package is designed for implementing a Bayesian model to reconstruct transmission networks in infectious diseases. The Bayesian transmission model combines genomic and temporal data to reconstruct transmission networks for infectious diseases.


## Installation
### Local installation
1. Download the package from GitHub
2. To add the local library betnet, run "] add path-to-betnet" in Julia.

### Global installation
In the package mode, run "add betnet" to add the package in Julia.

## How to run the program
The input datasets time_real_data.csv (temporal data) and SNP_real_data.csv (SNP distances) are available in the folder "test". To analyze a real dataset, run the following code

    using betnet
    using DataFrames
    using CSV
    
    testDDD = Array(DataFrame(CSV.File("time_real_data.csv")))
    testDD = testDDD[:,3:5]
    testD = zeros(Float64,length(testDDD[:,1]),3)
    testD[:,2:3] = float(testDD[:,2:3])/365
    testD[:,1] = float(testDD[:,1])

    #SNP Difference Data
    DDD = Array(DataFrame(CSV.File("SNP_real_data.csv")))
    DD = DDD[:,2:70]
    D = float(DD)


    #MCMC
    @time rd_ppp_111_1 = TransNet(testD, D, 1.0e-6, 4411532.0, 2.0, 5E-7,0.3,5e-7, 5000,0,0)

    inf_Net_rd_ppp_111_1 = DataFrame(rd_ppp_111_1[1],:auto)
    Param_rd_ppp_111_1 = DataFrame(rd_ppp_111_1[4],:auto)
    Tb_rd_ppp_111_1 = rd_ppp_111_1[2]

    #Here input the inferred prediction ID to get the corresponding tb values
    Latent_rd_ppp_111_1 = DataFrame(Tborg(Tb_rd_ppp_111_1,[2:1:69;], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 11, 6, 7, 2, 1, 1, 13, 1, 19, 1, 16, 1, 1, 22, 7, 7, 27, 1, 1, 5, 5, 27, 27, 27, 27, 18, 37, 36, 34, 34, 37, 1, 27, 5, 1, 1, 1, 34, 22, 46, 27, 14, 27, 49, 36, 5, 27, 48, 23, 25, 38, 35, 1, 53, 19, 27, 30, 27]), :auto)

    CSV.write("inf_NetID_rd1111.csv", inf_Net_rd_ppp_111_1)
    CSV.write("Para_Net_rd1111.csv", Param_rd_ppp_111_1)
    CSV.write("Latent_rd1111.csv", Latent_rd_ppp_111_1)

## Resources
- Support: We use GitHub for the development of the Julia package Distributions itself. For support and questions, please use the Julia Discourse forum. Also, for casual conversation and quick questions, there are the channels #helpdesk and #statistics in the official Julia chat (https://julialang.slack.com). To get an invitation, please visit https://julialang.org/slack/.



## Contributing
### Workflow with Git and GitHub
To contribute to the package, fork the repository on GitHub, clone it and make modifications on a new branch, do not commit modifications on master. Once your changes are made, push them on your fork and create the Pull Request on the main repository.

### Citation
Xu, J., Hu, H., Ellison, G., Yu, L., Whalen, C.C., Liu, L.* Bayesian estimation of transmission networks for infectious diseases. J. Math. Biol. 2025, 90(29) [DOI link](https://doi.org/10.1007/s00285-025-02193-1).





