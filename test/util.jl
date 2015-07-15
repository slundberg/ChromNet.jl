using GZip

## window_bed_file
binValues = GZip.open(window_bed_file, "data/K562_tss.bed.gz")

## save_chromnet_matrix and load_chromnet_matrix
P = 26
M = trues(P,100150)
M[:,:] = rand(P,100150) .> 0.5
letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
header = [letters[i:i] for i in 1:P]
isfile("data/tmp.ChromNet") && rm("data/tmp.ChromNet")
save_chromnet_matrix("data/tmp.ChromNet", M, header)
M2,header2 = load_chromnet_matrix("data/tmp.ChromNet")
@test all(M2 .== M)
@test all(header2 .== header)

## streaming_cov from a file
C = open(f->streaming_cov(f, chunkSize=3, quiet=true), "data/tmp.ChromNet/matrix")
@test all(abs(C .- cov(M')) .< 1e-10)

## streaming_cov from an inmemory bit array
@test all(abs(streaming_cov(M, chunkSize=3, quiet=true) .- cov(M')) .< 1e-10)

## conditional_cov
C = cov(M')
smallC = conditional_cov(C, 10, 0.0)
IC = inv(C)
@test all(abs(smallC .- inv(IC[11:end,11:end])) .< 1e-10)

## cov2cor!
C = cov(M')
@test all(abs(cov2cor!(C) .- cor(M')) .< 1e-10)
