export window_bed_file, load_chromnet_matrix, save_chromnet_matrix, streaming_cov, cov2cor!

chrNames = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"
]

# how long we consider each chromosome (a rough upper bound)
chrLengths = [
    250000000,244000000,199000000,192000000,181000000,172000000,160000000,147000000,
    142000000,136000000,136000000,134000000,116000000,108000000,103000000,91000000,
    82000000,79000000,60000000,64000000,49000000,52000000,156000000,60000000
]

# compute offsets to embed al the chromsomes in a linear sequence
chrOffsets = Dict{ASCIIString,Int64}()
for i in 1:length(chrLengths)
    chrOffsets[chrNames[i]] = sum(chrLengths[1:i-1])
end

# see how long the genome is in bins
totalBins = int(ceil(sum(chrLengths)/1000))
 
# creates a BitArray with a 1 for every bin overlapped by the bed file regions
function window_bed_file(stream)
	global chrOffsets, totalBins

	# mark all bins that are touched with 1
	binValues = falses(totalBins)
	for line in eachline(stream)
	    parts = split(line, '\t')
	    if haskey(chrOffsets, parts[1])
	        startPos = ceil((chrOffsets[parts[1]]+int(parts[2]))/1000)
	        endPos = ceil((chrOffsets[parts[1]]+int(parts[3]))/1000)
	        for i in startPos:endPos
	            binValues[i] = 1
	        end
	    end
	end

	binValues
end

function streaming_cov(stream::IOStream; chunkSize=100000, quiet=false)
    P = read(stream, Int64)
    N = read(stream, Int64)
    XtX = zeros(Float64, P, P)
    varSums = zeros(Float64, P, 1)

    # force the chunk size to line up with 64 bit word boundaries,
    # this is important for loading the BitArray in blocks.
    # we try and keep the size close to what was requested
    chunkSize = max(int(chunkSize/64),1)*64
    
    # build XtX incrementally and also the totals of every variable.
    chunkBit = BitArray(P, chunkSize)
    chunk = Array(Float32, P, chunkSize)
    numChunks = int(ceil(N/chunkSize))
    for i in 1:numChunks-1
        read!(stream, chunkBit)
        chunk[:,:] = chunkBit
        XtX .+= A_mul_Bt(chunk,chunk) # using a float array is important to get LAPACK speed
        varSums .+= sum(chunk,2)
        if !quiet println("processed $(i*chunkSize*1000) bp...") end
    end

    # get the last unevenly sized chunk
    chunkBit = BitArray(P, N - (numChunks-1)*chunkSize)
    chunk = Array(Float32, P, N - (numChunks-1)*chunkSize)
    read!(stream, chunkBit)
    chunk[:,:] = chunkBit
    XtX .+= A_mul_Bt(chunk,chunk)
    varSums .+= sum(chunk,2)
    
    # convert XtX to a covariance matrix
    XtX .-= varSums*varSums'/N
    XtX ./= (N-1)
end

function streaming_cov(bigData::BitArray{2}; chunkSize=100000, quiet=false)
    P,N = size(bigData)
    XtX = zeros(Float64, P, P)
    varSums = zeros(Float64, P, 1)

    # force the chunk size to line up with 64 bit word boundaries,
    # this is most important for loading from a file, but we also use it here.
    # we try and keep the size close to what was requested
    chunkSize = max(int(chunkSize/64),1)*64
    
    # build XtX incrementally and also the totals of every variable.
    chunk = Array(Float32, P, chunkSize)
    numChunks = int(ceil(N/chunkSize))
    for i in 1:numChunks-1
        chunk[:,:] = bigData[:,(i-1)*chunkSize+1:i*chunkSize]
        XtX .+= A_mul_Bt(chunk,chunk) # using a float array is important to get LAPACK speed
        varSums .+= sum(chunk,2)
        if !quiet println("processed $(i*chunkSize*1000) bp...") end
    end

    # get the last unevenly sized chunk
    chunk = Array(Float32, P, N - (numChunks-1)*chunkSize)
    chunk[:,:] = bigData[:,(numChunks-1)*chunkSize+1:end]
    XtX .+= A_mul_Bt(chunk,chunk)
    varSums .+= sum(chunk,2)
    
    # convert XtX to a covariance matrix
    XtX .-= varSums*varSums'/N
    XtX ./= (N-1)
end

function load_chromnet_matrix(dirName)
	f = open("$dirName/matrix")
	P = read(f, Int64)
	N = read(f, Int64)
	data = falses(P, N)
	read!(f, data)
	close(f)
	data,vec(readdlm("$dirName/header"))
end

function save_chromnet_matrix(outDir, data::BitArray{2}, header::Array)
	@assert size(data)[1] == length(header)

	mkpath(outDir)

	f = open("$outDir/matrix", "w")
	write(f, size(data)[1])
    write(f, size(data)[2])
    write(f, data)
    close(f)

    writedlm("$outDir/header", header)
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