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
    "--output","-o"
        help = "Where to save the network JSON file (default is STDOUT)"
        default = "-"
    "--save-cov"
        help = "Save a copy of the raw data covariance matrix"
        metavar = "COV_FILE"
    "--threshold"
        help = "At what absolute value threshold do we keep links"
        arg_type = Float64
        default = 0.03
    "--group-link-threshold"
        help = "Below what group score do we discard connecting links [0.0-1.0]"
        arg_type = Float64
        default = 0.7
        metavar = "THRESHOLD"
    "--quiet", "-q"
        help = "Don't print anything"
        action = :store_true
end
args = parse_args(s)
inputConfig = args["input_config"]
dataDir = args["data"]
output = args["output"]
saveCov = args["save-cov"]
quiet = args["quiet"]
println(saveCov)
exit()
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

metadata = dataDir == "none" ? {} : JSON.parse(open(readall, "$dataDir/metadata"))
metadata = merge(metadata, parse_config(inStream, bedFileRoot))
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
mainHeader = ASCIIString[]
if dataDir != "none"
    mainData,mainHeader = load_chromnet_matrix(dataDir)
end

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
C = streaming_cov(jointData, quiet=quiet)
if saveCov != nothing
    f = open(saveCov, "w")
    println(f, join(jointData, ','))
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
