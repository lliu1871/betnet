# Bayesian Estimation of Transmission NETworks (betnet v2.0)
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

1. In Julia

```{code}
    using betnet
    
    #Load the temporal and Genomic Data
    @time transNetworkInference(tempfile="./data/time_real_data.csv",SNPfile="./data/SNP_real_data.csv",Contactfile="",genomeSize=4411532, itr_MCMC=1000000, burn_in=0, subsample=1000, outputfile="parameter_betnet_69.csv")
    
    # simulation data with/without contact network
    @time transNetworkInference(tempfile="./data/TemporalData1_100.csv",SNPfile="./data/SNP1_100.csv",Contactfile="",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_nonet.csv")
    
    @time transNetworkInference(tempfile="./data/TemporalData1_100.csv",SNPfile="./data/SNP1_100.csv",Contactfile="./data/ContactProb1_100.csv",genomeSize=1000000, itr_MCMC=1000000, burn_in=0,subsample=1000, outputfile="parameter_betnet_net.csv")
```

2. The binary betnet2.0.jl is in the latest release
```{code}
    chmod +x betnet2.0.jl
    julia betnet2.0.jl --help
    julia betnet2.0.jl -t TemporalData1_100.csv -s SNP1_100.csv -c ContactProb1_100.csv
```

## Output
In the MCMC output file, the first row contains columns called "iteration",	"logLikelihood", "logPrior", "theta", "mutation_rate", "infection_rate", and each remaining column corresponds to one individual in your dataset. Inside each column, the value = ID of the inferred infector at that MCMC iteration. Since this varies by iteration, you can compute the  posterior probability that A infected B and thereafter find the most likely infector of each case.


## Resources
- Simulation: The R package [transNetwork](https://github.com/lliu1871/transNetwork) provide a unified framework for simulating, summarizing, and visualizing infectious disease transmission networks using genomic and epidemiological data
- Support: We use GitHub for the development of the Julia package Distributions itself. For support and questions, please use the Julia Discourse forum. Also, for casual conversation and quick questions, there are the channels #helpdesk and #statistics in the official Julia chat (https://julialang.slack.com). To get an invitation, please visit https://julialang.org/slack/.

## Contributing
### Workflow with Git and GitHub
To contribute to the package, fork the repository on GitHub, clone it and make modifications on a new branch, do not commit modifications on master. Once your changes are made, push them on your fork and create the Pull Request on the main repository.

### Citation
Xu, J., Hu, H., Ellison, G., Yu, L., Whalen, C.C., Liu, L. Bayesian estimation of transmission networks for infectious diseases. J. Math. Biol. 2025, 90(29) [DOI link](https://doi.org/10.1007/s00285-025-02193-1).

Xu, J., Kim, J., Ji, P., Yu, L., Whalen, C.C., Liu, L. A Bayesian framework for the network analysis of transmission dynamics in infectious disease. J Mol Evol 2026. https://doi.org/10.1007/s00239-026-10303-w  





