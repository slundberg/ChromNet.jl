using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "input_config"
        help = "a config file that lists which bed files to integrate into the network"
        arg_type = ASCIIString
        nargs = '+'
    "--data", "-d"
        help = "location of the binary chromnet data matrix"
        default = "./ENCODE_build1.ChromNet"
    "--quiet", "-q"
        help = "don't print anything"
        action = :store_true
end
args = parse_args(s)
inputConfig = args["input_config"]
dataDir = args["data"]
quiet = args["quiet"]

using ChromNet

# load the main data matrix
mainData,mainHeader = load_chromnet_matrix(dataDir)

# create space for the joint data matrix
jointData = falses(length(bedFiles)+length(mainHeader), ChromNet.totalBins)
jointData[length(bedFiles)+1:end,:] = mainData'
jointHeader = [bedFiles, mainHeader]

# bin the bed files and add them to the joint matrix
for (i,file) in enumerate(bedFiles)
    jointData[i,:] = open(window_bed_file, file)
    println(i)
end

# compute the joint correlation matrix
C = cov2corr!(streaming_cov(jointData))

# compute groups using hclust

# compute groupgm matrix

# save network json data

