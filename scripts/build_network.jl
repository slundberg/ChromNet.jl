using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "data_bundle"
        help = "Location of the ChromNet data bundle (i.e. ENCODE_buildX.ChromNet.jld)"
        arg_type = ASCIIString
        required = true
    "output_location"
        help = "Where to save the network JSON file (default is STDOUT)"
        default = "-"
    "--save-cor"
        help = "Save a copy of the raw data correlation matrix"
        metavar = "COR_FILE"
    "--threshold"
        help = "At what absolute value threshold do we keep links"
        arg_type = Float64
        default = 0.03
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

using ChromNet
using JSON
using JLD

# open the data bundle and load the metadata
dataBundle = jldopen(args["data_bundle"])
metadata = read(dataBundle, "metadata")
ids = read(dataBundle, "ids")

# compute the joint correlation matrix
C = streaming_cor(args["data_bundle"], quiet=args["quiet"])
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
network = build_network(G, groups, metadata, threshold=args["threshold"], groupScoreThreshold=args["group-score-threshold"])

# send the network in JSON format to the output stream
outStream = STDOUT
if args["output_location"] != "-"
    outStream = open(args["output_location"])
end
println(outStream, json(network))
