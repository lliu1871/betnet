# Bayesian Estimation of Transmission NETworks (betnet v1.1)
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
The input datasets time_real_data.csv (temporal data) and SNP_real_data.csv (SNP distances) are available in the folder "data". To analyze the datasets, run the following Julia code

    using betnet
    
    #Load the temporal and Genomic Data
    @time transNetworkInference(tempfile="./bestNet/data/time_real_data.csv",SNPfile="./betnet2.0/data/SNP_real_data.csv",Contactfile="",genomeSize=4411532, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_69.csv")
    
    # simulation data with/without contact network
    @time transNetworkInference(tempfile="./bestNet/data/TemporalData1_100.csv",SNPfile="./betnet2.0/data/SNP1_100.csv",Contactfile="",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_nonet.csv")
    
    @time transNetworkInference(tempfile="./bestNet/data/TemporalData1_100.csv",SNPfile="./betnet2.0/data/SNP1_100.csv",Contactfile="./betnet2.0/data/ContactProb1_100.csv",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_net.csv")
 

## Resources
- Support: We use GitHub for the development of the Julia package Distributions itself. For support and questions, please use the Julia Discourse forum. Also, for casual conversation and quick questions, there are the channels #helpdesk and #statistics in the official Julia chat (https://julialang.slack.com). To get an invitation, please visit https://julialang.org/slack/.


## Contributing
### Workflow with Git and GitHub
To contribute to the package, fork the repository on GitHub, clone it and make modifications on a new branch, do not commit modifications on master. Once your changes are made, push them on your fork and create the Pull Request on the main repository.

### Citation
Xu, J., Hu, H., Ellison, G., Yu, L., Whalen, C.C., Liu, L. Bayesian estimation of transmission networks for infectious diseases. J. Math. Biol. 2025, 90(29) [DOI link](https://doi.org/10.1007/s00285-025-02193-1).





