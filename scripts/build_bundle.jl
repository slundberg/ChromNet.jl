using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "output_file"
        help = "Output location to save the new data bundle."
        arg_type = ASCIIString
        required = true
    "input_config"
        help = "File listing which external data files to integrate into the network. These can be"*
               "either BAM files or BED files."*
               "For each line trailing fields may be omitted and follow the tab separated format: "*
               "FILE_NAME SHORT_TITLE LONG_TITLE CELL_TYPE LAB EXPERIMENT_ID ANTIBODY_ID TREATMENTS ORGANISM LIFE_STAGE LINK"
        arg_type = ASCIIString
        required = true
    "--assembly"
        help = "Which assembly genome coordinates are in (currently GRCh38, GRCh38_alt, hg19, mm10, mm9)."
        arg_type = ASCIIString
        default = "GRCh38_alt"
    "--quiet", "-q"
        help = "Don't print anything"
        action = :store_true
end
args = parse_args(s)

using ChromNet
using JSON
using GZip
using JLD
using SamIO
using ProgressMeter

# parse the input arguments
inputConfig = args["input_config"]
outputFile = args["output_file"]
if !ismatch(r"\.ChromNet\.jld$", outputFile)
    outputFile = "$outputFile.ChromNet.jld"
end
quiet = args["quiet"]

# load the metadata
inStream = STDIN
dataFileRoot = "."
if inputConfig != "-"
    inStream = open(inputConfig)
    dataFileRoot = dirname(inputConfig)
end
metadata = parse_config(inStream, dataFileRoot)
inputConfig != "-" && close(inStream)

# extract all the data files and custom dataset ids
fileTypes = ASCIIString[]
fileNames = ASCIIString[]
ids = ASCIIString[]
for (k,v) in metadata
    if v["fileType"] == "unknown"
        quiet || println(STDERR, "Ignoring file of unknown type: ", v["dataFile"])
        continue
    end
    push!(fileTypes, v["fileType"])
    push!(fileNames, v["dataFile"])
    push!(ids, k)
end

referenceContigs = SamIO.assembly[args["assembly"]]

# build the custom data matrix
numBins = ceil(Int64, sum(referenceContigs.sizes) / 1000)
data = zeros(Float32, length(ids), numBins)
@showprogress "Binning BAM/BED files: " for i in 1:length(ids)
    if fileTypes[i] == "bed" || fileTypes[i] == "narrowPeak"
        reader = BinningMap(BedReader(fileNames[i], referenceContigs), 1000)
        while !eof(reader)
            data[i,position(reader)] = value(reader)
            advance!(reader)
        end
        close(reader)

    elseif fileTypes[i] == "bed.gz" || fileTypes[i] == "narrowPeak.gz"
        f = GZip.open(fileNames[i])
        reader = BinningMap(BedReader(f, referenceContigs), 1000)
        while !eof(reader)
            data[i,position(reader)] = value(reader)
            advance!(reader)
        end
        close(reader)

    elseif fileTypes[i] == "bam"
        reader = BinningMap(BamReader(fileNames[i], :any, referenceContigs), 1000)
        while !eof(reader)
            data[i,position(reader)] = value(reader)
            advance!(reader)
        end
        close(reader)
    end
end

# build summary statistics
!quiet && println(STDERR, "Computing means and norms...")
means = vec(mean(data, 2))
norms = Float64[norm(data[i,:] - means[i]) for i in 1:size(data)[1]]

# write the JLD file
!quiet && println(STDERR, "Writing results to $outputFile...")
jldopen(outputFile, "w") do f
    f["data"] = data
    f["means"] = vec(mean(data, 2))
    f["norms"] = Float64[norm(data[i,:] - means[i]) for i in 1:size(data)[1]]
    f["ids"] = ids
    f["P"] = length(ids)
    f["N"] = size(data)[2]
    f["metadata"] = metadata
    f["assembly"] = args["assembly"]
end
