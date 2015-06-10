using ChromNet
using Base.Test
using GZip


## window_bed_file
binValues = GZip.open(window_bed_file, "data/K562_tss.bed.gz")

## save_chromnet_matrix and load_chromnet_matrix
M = trues(30,100150)
M[:,:] = rand(30,100150) .> 0.5
header = ["A" for i in 1:30]
isfile("data/tmp.BitArray") && rm("data/tmp.BitArray")
save_chromnet_matrix("data/tmp.BitArray", M, header)
M2,header2 = load_chromnet_matrix("data/tmp.BitArray")
@test all(M2 .== M)
@test all(header2 .== header)

## streaming_cov
C = open(f->streaming_cov(f, chunkSize=3), "data/tmp.BitArray/matrix")
all(C .== cov(M'))