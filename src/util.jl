using JLD
using ProgressMeter

export window_bed_file, streaming_cor, cov2cor!

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
function streaming_cor(dataFile; chunkSize=100000, quiet=false)

    # open data data file
    f = jldopen(dataFile)
    P = read(f, "P")
    N = read(f, "N")
    means = read(f, "means")
    data = readmmap(f["data"])
    close(f)

    # get all the regular sized chunks
    XtX = zeros(Float64, P, P)
    chunkInt32 = Array(Int32, P, chunkSize)
    chunk = Array(Float64, P, chunkSize)
    numChunks = round(Int64, ceil(N/chunkSize))
    p = Progress(N, 0.1, "Computing correlation: ", 34)
    for i in 1:chunkSize:N
        if i+chunkSize-1 > N # account for the last unevenly sized chunk
            chunk = Array(Float64, P, N-i+1)
        end
        streaming_cor_chunk!(XtX, data, (i, min(i+chunkSize-1,N)), chunk, means)
        quiet || update!(p, min(i+chunkSize-1,N))
    end

    cov2cor!(XtX)
    XtX + XtX' - eye(P)
end
function streaming_cor_chunk!(XtX, data, range, chunk, means)
    P,w = size(chunk)
    copy!(chunk, 1, data, P*(range[1]-1)+1, P*w) # converting to a float array is also important to get BLAS speed
    for j in 1:size(chunk)[2]
        @simd for i in 1:size(chunk)[1]
            @inbounds chunk[i,j] -= means[i]
        end
    end
    BLAS.syrk!('U', 'N', 1.0, chunk, 1.0, XtX)
end

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
