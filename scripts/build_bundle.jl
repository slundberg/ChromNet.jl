using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "bundle_name"
        help = "Output location to save the new ChromNet data bundle."
        arg_type = ASCIIString
        required = true
    "input_config"
        help = "File listing which external data files to integrate into the network. These can be"*
               "either BAM files or BED files using hg38 coordinates."
               "For each line trailing fields may be omitted and follow the tab separated format: "*
               "FILE_NAME SHORT_TITLE LONG_TITLE CELL_TYPE LAB EXPERIMENT_ID ANTIBODY_ID TREATMENTS ORGANISM LIFE_STAGE LINK"
        arg_type = ASCIIString
        required = true
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

# parse the input arguments
inputConfig = args["input_config"]
# dataBundle = nothing
# if args["base_bundle"] != "none"
#     dataBundle = jldopen(args["base_bundle"])
# end
output = args["output"]
#saveCor = args["save-cor"]
quiet = args["quiet"]

# load the metadata
inStream = STDIN
dataFileRoot = "."
if inputConfig != "-"
    inStream = open(inputConfig)
    dataFileRoot = dirname(inputConfig)
end

metadata = dataBundle == nothing ? Dict() : JSON.parse(read(dataBundle, "metadata"))
metadata = merge(metadata, parse_config(inStream, dataFileRoot))
inputConfig != "-" && close(inStream)

# extract all the data files and custom dataset ids
bedFiles = ASCIIString[]
bamFiles = ASCIIString[]
customIds = ASCIIString[]
for (k,v) in metadata
    if haskey(v, "dataFile")
        # if v["fileType"] == "bed"
        #     push!(bedFiles, v["dataFile"])
        # else
        if v["fileType"] == "bam"
            push!(bamFiles, v["dataFile"])
        else
            quiet || println(STDERR, "Ignoring file of unknown type: ", v["dataFile"])
        end
        push!(customIds, k)
    end
end

# load the main data matrix
# mainHeader = ASCIIString[]
# if dataFile != "none"
#     mainData,mainHeader = load_chromnet_matrix(dataFile)
# end
#
# # create space for the joint data matrix
# jointData = falses(length(customHeader)+length(mainHeader), ChromNet.totalBins)
# if dataFile != "none"
#     jointData[length(customHeader)+1:end,:] = mainData
# end
jointIds = [customIds; read(dataBundle, "ids")]

numBins = round(Int, ceil(sum(ReferenceContigs_hg38.sizes) / 1000))
customData = zeros(Int32, length(customIds), numBins)

# bin the bed files and add them to the joint matrix
# @showprogress "Windowing BED files: " for (i,file) in enumerate(bedFiles)
#     if ismatch(r"\.gz$", file)
#         jointData[i,:] = GZip.open(window_bed_file, file)
#     else
#         jointData[i,:] = open(window_bed_file, file)
#     end
# end

# bin the BAM files and add them to the joint matrix
@showprogress "Binning BAM files: " for (i,bamFile) in enumerate(bamFiles)
    reader = BinningMap(BamReader(bamFile, :any, ReferenceContigs_hg38), 1000)
    while !eof(reader)
        customData[i,position(reader)] = value(reader)
        advance!(reader)
    end
    close(reader)
end

# compute the joint correlation matrix
C = streaming_cov(jointData, quiet=quiet)
C .+= eye(size(C)...)*0.000001 # prevent singular matricies
if typeof(saveCor) != Void
    f = open(saveCor, "w")
    println(f, join(jointHeader, ','))
    writecsv(f, C)
    close(f)
end
cov2cor!(C)

# compute groups using hclust
groups = build_groups(C, jointHeader)

# compute groupgm matrix
G,headerG = build_groupgm(inv(C), jointHeader, groups)

# open the output stream
outStream = STDOUT
if output != "-"
    outStream = open(output)
end

# output network json data
network = build_network(G, headerG, groups, metadata, threshold=0.03, groupLinkThreshold=0.7)
println(outStream, json(network))
