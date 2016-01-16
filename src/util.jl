using JLD
using ProgressMeter

export window_bed_file, streaming_cor, cov2cor!, streaming_edgeimpact

# creates a BitArray with a 1 for every bin overlapped by the bed file regions
function window_bed_file(stream, contigs; binSize=1000)
    numBins = ceil(Int64, sum(contigs.sizes) / binSize)
    chrOffsets = Dict{ASCIIString,Int64}()
    for i in 1:contigs.count
        chrOffsets[contigs.names[i]] = contigs.offsets[i]
    end

    # mark all bins that are touched with 1
    binValues = falses(numBins)
    for line in eachline(stream)
        parts = split(line, '\t')
        if haskey(chrOffsets, parts[1])
            startPos = ceil(Int64, (chrOffsets[parts[1]]+parse(Int64, parts[2]))/binSize)
            endPos = ceil(Int64, (chrOffsets[parts[1]]+parse(Int64, parts[3]))/binSize)
            for i in startPos:endPos
                binValues[i] = true
            end
        end
    end

    binValues
end

"Designed to compute the correlation matrix from a ChromNet data bundle."
function streaming_cor(dataBundles; chunkSize=10000, quiet=false)

    # load data arrays
    P = 0
    datas = Any[]
    N = read(dataBundles[1], "N")
    for bundle in dataBundles
        push!(datas, readmmap(bundle["data"]))
        P += read(bundle, "P")
        @assert N == read(bundle, "N") "Data bundles have differing numbers of samples!"
    end

    # get the means
    means = zeros(P)
    lastPos = 0
    for bundle in dataBundles
        bundleP = read(bundle, "P")
        means[lastPos+1:lastPos+bundleP] = read(bundle, "means")
        lastPos = lastPos+bundleP
    end

    # loop through all the chunks
    XtX = zeros(Float64, P, P)
    chunk = Array(Float64, P, chunkSize)
    numChunks = round(Int64, ceil(N/chunkSize))
    p = Progress(N, 0.1, "Computing correlation: ", 34, :green, STDERR)
    for i in 1:chunkSize:N
        if i+chunkSize-1 > N # account for the last unevenly sized chunk
            chunk = Array(Float64, P, N-i+1)
        end
        streaming_cor_chunk!(XtX, datas, (i, min(i+chunkSize-1,N)), chunk, means)
        if !quiet update!(p, min(i+chunkSize-1,N)) end
    end

    cov2cor!(XtX)
    XtX + XtX' - eye(P)
end
function streaming_cor_chunk!(XtX, datas, range, chunk, means)
    w = size(chunk)[2]

    lastPos = 0
    for data in datas
        bundleP = size(data)[1]
        chunk[lastPos+1:lastPos+bundleP,:] = data[:,range[1]:range[1]+w-1]
        lastPos = lastPos+bundleP
    end

    for j in 1:size(chunk)[2]
        @simd for i in 1:size(chunk)[1]
            @inbounds chunk[i,j] -= means[i]
        end
    end
    BLAS.syrk!('U', 'N', 1.0, chunk, 1.0, XtX)
end

function streaming_edgeimpact(dataBundles, C, rows, cols; chunkSize=10000, quiet=false)
    IC = inv(C)

    # load data arrays
    P = 0
    datas = Any[]
    N = read(dataBundles[1], "N")
    for bundle in dataBundles
        push!(datas, readmmap(bundle["data"]))
        P += read(bundle, "P")
        @assert N == read(bundle, "N") "Data bundles have differing numbers of samples!"
    end

    # get the means and norms
    means = zeros(P)
    norms = zeros(P)
    lastPos = 0
    for bundle in dataBundles
        bundleP = read(bundle, "P")
        means[lastPos+1:lastPos+bundleP] = read(bundle, "means")
        norms[lastPos+1:lastPos+bundleP] = read(bundle, "norms")
        lastPos = lastPos+bundleP
    end

    # loop through all the chunks
    out = Float64[]
    vtmp = zeros(Float64, P)
    ctmp = zeros(Float64, P)
    outChunk = zeros(Float64, chunkSize)
    chunk = Array(Float64, P, chunkSize)
    chunkTmp = Array(Float64, P, chunkSize)
    p = Progress(N, 0.1, "Computing sample relevance: ", 34, :green, STDERR)
    for i in 1:chunkSize:N
        if i+chunkSize-1 > N # account for the last unevenly sized chunk
            outChunk = Array(Float64, N-i+1)
            chunk = Array(Float64, P, N-i+1)
            chunkTmp = Array(Float64, P, N-i+1)
        end
        streaming_edgeimpact_chunk!(outChunk, datas, (i, min(i+chunkSize-1,N)), IC, rows, cols, chunk, chunkTmp, means, norms, vtmp, ctmp)
        append!(out, outChunk)
        if !quiet update!(p, min(i+chunkSize-1,N)) end
    end
    out
end
function streaming_edgeimpact_chunk!(out, datas, range, IC, rows, cols, chunk, chunkTmp, means, norms, vtmp, ctmp)
    P,w = size(chunk)
    @assert size(IC) == (P,P)
    @assert length(means) == P
    @assert length(norms) == P
    @assert length(vtmp) == P
    @assert length(ctmp) == P
    @assert length(out) == w

    # load the data chunk and normalize it
    lastPos = 0
    for data in datas
        bundleP = size(data)[1]
        chunk[lastPos+1:lastPos+bundleP,:] = data[:,range[1]:range[1]+w-1] # converting to a float array is also important to get BLAS speed
        lastPos = lastPos+bundleP
    end
    for j in 1:size(chunk)[2]
        @simd for i in 1:P
            @inbounds chunk[i,j] -= means[i]
            @inbounds chunk[i,j] /= norms[i]
        end
    end

    # precompute the matrix-vector multiplication for all samples
    BLAS.symm!('L', 'U', 1.0, IC, chunk, 0.0, chunkTmp)

    # compute the sum of the changes in the given edges for each sample
    for ind in 1:length(out)
        copy!(ctmp, 1, chunk, P*(ind-1)+1, P)
        copy!(vtmp, 1, chunkTmp, P*(ind-1)+1, P)
        s = 1 - dot(ctmp,vtmp)
        totalChange = 0.0
        for j in cols
            @inbounds dj = sqrt(1 - ctmp[j].^2)
            for i in rows
                @inbounds di = sqrt(1 - ctmp[i].^2)
                @inbounds totalChange += di*dj*(IC[i,j] + vtmp[i]*vtmp[j]/s) - IC[i,j]
            end
        end
        @inbounds out[ind] = totalChange
    end
    out
end


# "Designed to compute the correlation matrix from a ChromNet data bundle."
# function streaming_cor(dataFile, memData=nothing; chunkSize=10000, quiet=false)
#
#     # open data file
#     f = jldopen(dataFile)
#     bundleP = read(f, "P")
#     N = read(f, "N")
#     means = read(f, "means")
#     data = readmmap(f["data"])
#     means = [read(f, "means"); mean(memData, 2)]
#     close(f)
#
#     if memData == nothing
#         memData = zeros(0,N)
#     else
#         @assert size(memData)[2] == N "Passed in-memory data size does not match bundle ($(size(memData)[2]) != $N)"
#     end
#     P = size(memData)[1] + bundleP
#
#     # loop through all the chunks
#     XtX = zeros(Float64, P, P)
#     chunk = Array(Float64, P, chunkSize)
#     numChunks = round(Int64, ceil(N/chunkSize))
#     p = Progress(N, 0.1, "Computing correlation: ", 34, :green, STDERR)
#     for i in 1:chunkSize:N
#         if i+chunkSize-1 > N # account for the last unevenly sized chunk
#             chunk = Array(Float64, P, N-i+1)
#         end
#         streaming_cor_chunk!(XtX, data, memData, (i, min(i+chunkSize-1,N)), chunk, means)
#         if !quiet update!(p, min(i+chunkSize-1,N)) end
#     end
#
#     cov2cor!(XtX)
#     XtX + XtX' - eye(P)
# end
# function streaming_cor_chunk!(XtX, data, memData, range, chunk, means)
#     memP = size(memData)[1]
#     bundleP = size(chunk)[1] - memP
#     w = size(chunk)[2]
#
#     chunk[1:bundleP,:] = data[:,range[1]:range[1]+w-1]
#     chunk[bundleP+1:end,:] = memData[:,range[1]:range[1]+w-1]
#     #copy!(chunk, 1, data, P*(range[1]-1)+1, P*w) # converting to a float array is also important to get BLAS speed
#     for j in 1:size(chunk)[2]
#         @simd for i in 1:size(chunk)[1]
#             @inbounds chunk[i,j] -= means[i]
#         end
#     end
#     BLAS.syrk!('U', 'N', 1.0, chunk, 1.0, XtX)
# end

function cov2cor!(M)
    for i in 1:size(M)[1]
        val = sqrt(M[i,i])
        if val > 0.0
            M[i,:] ./= val
            M[:,i] ./= val
        end
    end
    M
end
