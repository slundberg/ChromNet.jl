using ChromNet
using Base.Test
using GZip


## window_bed_file
binValues = GZip.open(window_bed_file, "data/K562_tss.bed.gz")

## save_chromnet_matrix and load_chromnet_matrix
P = 30
M = trues(P,100150)
M[:,:] = rand(P,100150) .> 0.5
header = ["A" for i in 1:P]
isfile("data/tmp.BitArray") && rm("data/tmp.BitArray")
save_chromnet_matrix("data/tmp.BitArray", M, header)
M2,header2 = load_chromnet_matrix("data/tmp.BitArray")
@test all(M2 .== M)
@test all(header2 .== header)

## streaming_cov from a file
C = open(f->streaming_cov(f, chunkSize=3, quiet=true), "data/tmp.BitArray/matrix")
@test all(abs(C .- cov(M')) .< 1e-10)

## streaming_cov from an inmemory bit array
#C = open(f->streaming_cov(f, chunkSize=3), "data/tmp.BitArray/matrix")
@test all(abs(streaming_cov(M, chunkSize=3, quiet=true) .- cov(M')) .< 1e-10)

## cov2cor!
C = cov(M')
@test all(abs(cov2cor!(C) .- cor(M')) .< 1e-10)