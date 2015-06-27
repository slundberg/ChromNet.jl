using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "input_config"
        help = "File listing which bed files to integrate into the network. "*
               "For each line trailing fields may be omitted and follow the tab separated format: "*
               "BED_FILE_NAME SHORT_TITLE LONG_TITLE CELL_TYPE LAB EXPERIMENT_ID ANTIBODY_ID TREATMENTS ORGANISM LIFE_STAGE LINK"
        arg_type = ASCIIString
        required = true
    "--data","-d"
        help = "Location of the binary chromnet data matrix"
        default = dirname(Base.source_path())*"/ENCODE_build1.ChromNet"
    "--threshold"
        help = "at what absolute value threshold do we keep links"
        arg_type = Float64
        default = 0.03
    "--group_link_threshold"
        help = "below what group score do we discard connecting links [0.0-1.0]"
        arg_type = Float64
        default = 0.7
    "--quiet", "-q"
        help = "Don't print anything"
        action = :store_true
end
args = parse_args(s)
inputConfig = args["input_config"]
dataDir = args["data"]
quiet = args["quiet"]

using ChromNet
using JSON
using GZip

# load the metadata
inStream = STDIN
bedFileRoot = "."
if inputConfig != "-"
    inStream = open(inputConfig)
    bedFileRoot = dirname(inputConfig)
end

metadata = merge(JSON.parse(open(readall, "$dataDir/metadata")), parse_config(inStream, bedFileRoot))
if inputConfig != "-" close(inStream) end

# extract all the bed files and custom dataset ids
bedFiles = ASCIIString[]
customHeader = ASCIIString[]
for (k,v) in metadata
    if haskey(v, "bedFile")
        push!(bedFiles, v["bedFile"])
        push!(customHeader, k)
    end
end

# load the main data matrix
mainData,mainHeader = load_chromnet_matrix(dataDir)

# create space for the joint data matrix
jointData = falses(length(customHeader)+length(mainHeader), ChromNet.totalBins)
jointData[length(customHeader)+1:end,:] = mainData
jointHeader = [customHeader, mainHeader]

# bin the bed files and add them to the joint matrix
quiet || println(STDERR, "Windowing BED files...")
for (i,file) in enumerate(bedFiles)
    if ismatch(r"\.gz$", file)
        jointData[i,:] = GZip.open(window_bed_file, file)
    else
        jointData[i,:] = open(window_bed_file, file)
    end
end

# compute the joint correlation matrix
quiet || println(STDERR, "Computing data covariance...")
C = cov2cor!(streaming_cov(jointData, quiet=quiet))

# compute groups using hclust
groups = build_groups(C, jointHeader)

# compute groupgm matrix
G,headerG = build_groupgm(inv(C), jointHeader, groups)

# output network json data
network = build_network(G, headerG, groups, metadata, threshold=0.03, groupLinkThreshold=0.7)
println(json(network))

