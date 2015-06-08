module ChromNet

export window_bed_file!

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

# this fills in the columns of the passed matrix from the given stream
function window_bed_file!(matrix, columnIndex::Int64, stream)
	global chrOffsets

	# mark all bins that are touched with 1
	for line in eachline(stream)
	    parts = split(line)
	    if haskey(chrOffsets, parts[1])
	        startPos = int(ceil((chrOffsets[parts[1]]+int(parts[2]))/1000))
	        endPos = int(ceil((chrOffsets[parts[1]]+int(parts[3]))/1000))
	        for i in startPos:endPos
	            matrix[i,columnIndex] = 1
	        end
	    end
	end
end

end