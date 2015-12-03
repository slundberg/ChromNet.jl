using GZip

## window_bed_file
#binValues = GZip.open(window_bed_file, "data/K562_tss.bed.gz")

## streaming_cor from a file
data = jldopen(f->read(f, "data"), "data/small.ChromNet.jld")
C = streaming_cor("data/small.ChromNet.jld", chunkSize=3, quiet=true)
@test all(abs(C .- cor(data, vardim=2)) .< 1e-10)

## streaming_cov from an in memory bit array
# @test all(abs(streaming_cov(M, chunkSize=3, quiet=true) .- cov(M')) .< 1e-10)

# ## conditional_cov
# C = cov(M')
# smallC = conditional_cov(C, 1:10, 11:P, 0.0)
# IC = inv(C)
# @test all(abs(smallC .- inv(IC[11:end,11:end])) .< 1e-10)

## cov2cor!
C = cov(data, vardim=2)
@test all(abs(cov2cor!(C) .- cor(data, vardim=2)) .< 1e-10)
