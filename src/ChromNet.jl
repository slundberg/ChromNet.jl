module ChromNet

export window_bed_file, load_chromnet_matrix, save_chromnet_matrix

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

function cov_sparse(data::BitArray)
    N,P = size(data)
    sumData = sum(data,1);
    C = zeros(Float64, P, P)
    for i in 1:P
        C[:,i] = sum(data[data[:,i],:], 1)
        if i % 10 == 0
            println(i)
        end
    end
    C .-= transpose(sumData)*sumData/N
    C ./= (N-1);
end

function load_chromnet_matrix(dirName)
	f = open("$dirName/matrix")
	N = read(f, Int64)
	P = read(f, Int64)
	data = falses(N, P)
	read!(f, data)
	close(f)
	data,vec(readdlm("$dirName/header"))
end

function save_chromnet_matrix(outDir, data::BitArray{2}, header::Array)
	@assert size(data)[2] == length(header)

	mkpath(outDir)

	f = open("$outDir/matrix", "w")
	write(f, size(data)[1])
    write(f, size(data)[2])
    write(f, data)
    close(f)

    writedlm("$outDir/header", header)
end

end