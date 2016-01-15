using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "data_bundle"
        help = "Location of the ChromNet data bundle (i.e. ENCODE_buildX.ChromNet.jld)"
        arg_type = ASCIIString
        required = true
    "input_config"
        help = "File listing which external data files to integrate into the network. These are"*
               "BAM files using hg38 coordinates."
               "For each line trailing fields may be omitted and follow the tab separated format: "*
               "FILE_NAME SHORT_TITLE LONG_TITLE CELL_TYPE LAB EXPERIMENT_ID ANTIBODY_ID TREATMENTS ORGANISM LIFE_STAGE LINK"
        arg_type = ASCIIString
        required = true
    "--output", "-o"
        help = "Where to save the network JSON file (default is STDOUT)"
        default = "-"
    "--save-cor"
        help = "Save a copy of the raw data correlation matrix"
        metavar = "COR_FILE"
    "--threshold"
        help = "At what absolute value threshold do we keep links"
        arg_type = Float64
        default = 0.2
    "--group-score-threshold"
        help = "Below what group score do we discard connecting links [0.0-1.0]"
        arg_type = Float64
        default = 0.7
        metavar = "THRESHOLD"
    "--quiet", "-q"
        help = "Don't print anything"
        action = :store_true
end
args = parse_args(s)

try
    @assert Pkg.installed("SamIO") != nothing
catch
    Pkg.clone("https://github.com/slundberg/SamIO.jl.git")
end


using ChromNet
using JSON
using JLD
using SamIO

# open the data bundle and load the metadata
dataBundle = jldopen(args["data_bundle"])
inStream = STDIN
dataFileRoot = "."
if args["input_config"] != "-"
    inStream = open(args["input_config"])
    dataFileRoot = dirname(args["input_config"])
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

# build the custom data matrix
ids = [customIds; read(dataBundle, "ids")]
numBins = round(Int, ceil(sum(ReferenceContigs_hg38.sizes) / 1000))
customData = zeros(Float32, length(customIds), numBins)
# @showprogress "Windowing BED files: " for (i,file) in enumerate(bedFiles)
#     if ismatch(r"\.gz$", file)
#         jointData[i,:] = GZip.open(window_bed_file, file)
#     else
#         jointData[i,:] = open(window_bed_file, file)
#     end
# end
@showprogress "Binning BAM files: " for (i,bamFile) in enumerate(bamFiles)
    reader = BinningMap(BamReader(bamFile, :any, ReferenceContigs_hg38), 1000)
    while !eof(reader)
        customData[i,position(reader)] = value(reader)
        advance!(reader)
    end
    close(reader)
end

# compute the joint correlation matrix
C = streaming_cor(args["data_bundle"], customData, quiet=args["quiet"])
if typeof(args["save-cor"]) != Void
    f = open(args["save-cor"], "w")
    println(f, join(ids, '\t'))
    writedlm(f, C)
    close(f)
end
cov2cor!(C)

# compute groups using hclust
groups = build_groups(C, ids)

# compute groupgm matrix and create network
# we add a small amount of regularization to avoid potential singularity issues in user-provided data
G,idsG = build_groupgm(inv(C + 1e-8I), ids, groups)
network = build_network(G, groups, metadata, threshold=args["threshold"], groupScoreThreshold=args["group-score-threshold"], quiet=args["quiet"])

# send the network in JSON format to the output stream
outStream = STDOUT
if args["output_location"] != "-"
    outStream = open(args["output_location"])
end
println(outStream, json(network))
