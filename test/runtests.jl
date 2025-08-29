using betnet
using DataFrames
using CSV
using Test

@testset "betnet.jl" begin
    #Load the temporal and Genomic Data
    #Test on 69 Patients

    #Temporal Data
    testDDD = Array(DataFrame(CSV.File("/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/betnet/test/time_real_data.csv")))
    testDD = testDDD[:,3:5]
    testD = zeros(Float64,length(testDDD[:,1]),3)
    testD[:,2:3] = float(testDD[:,2:3])/365
    testD[:,1] = float(testDD[:,1])

    #SNP Difference Data
    DDD = Array(DataFrame(CSV.File("/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/betnet/test/SNP_real_data.csv")))
    DD = DDD[:,2:70]
    D = float(DD)


    #MCMC
    @time rd_ppp_111_1 = TransNet(testD, D, 1.0e-6, 4411532.0, 2.0, 5E-7,0.3,5e-7, 5000,0,0)

    inf_Net_rd_ppp_111_1 = DataFrame(rd_ppp_111_1[1],:auto)
    Param_rd_ppp_111_1 = DataFrame(rd_ppp_111_1[4],:auto)
    Tb_rd_ppp_111_1 = rd_ppp_111_1[2]

    #Here input the inferred prediction ID to get the corresponding tb values
    Latent_rd_ppp_111_1 = DataFrame(Tborg(Tb_rd_ppp_111_1,[2:1:69;], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 11, 6, 7, 2, 1, 1, 13, 1, 19, 1, 16, 1, 1, 22, 7, 7, 27, 1, 1, 5, 5, 27, 27, 27, 27, 18, 37, 36, 34, 34, 37, 1, 27, 5, 1, 1, 1, 34, 22, 46, 27, 14, 27, 49, 36, 5, 27, 48, 23, 25, 38, 35, 1, 53, 19, 27, 30, 27]), :auto)

    CSV.write("/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/betnet/test/inf_NetID_rd1111.csv", inf_Net_rd_ppp_111_1)
    CSV.write("/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/betnet/test/Para_Net_rd1111.csv", Param_rd_ppp_111_1)
    CSV.write("/Users/lliu/Library/CloudStorage/OneDrive-UniversityofGeorgia/Dropbox/Github/Julia/betnet/test/Latent_rd1111.csv", Latent_rd_ppp_111_1)

end
