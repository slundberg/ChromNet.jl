using ArgParse
s = ArgParseSettings()
@add_arg_table s begin
    "group1"
        help = "One end point of the edge for which to analyze context, where "*
               "| is used to separate mulitple ids. (i.e. ENCFF000ARO|ENCFF000ARK)."
        arg_type = ASCIIString
        required = true
    "group2"
        help = "The other end point of the edge for which to analyze context, where "*
               "| is used to separate mulitple ids. (i.e. ENCFF000ARO|ENCFF000ARK)."
        arg_type = ASCIIString
        required = true
    "data_bundles"
        help = "Location of the ChromNet data bundles (i.e. ENCODE_buildX.ChromNet.jld)"
        arg_type = ASCIIString
        nargs = '+'
        required = true
    "--output", "-o"
        help = "Where to save the edge context BED file"
        default = "-"
    "--load-cor"
        help = "A copy of the raw data correlation matrix"
        metavar = "COR_FILE"
    "--quiet", "-q"
        help = "Don't print anything"
        action = :store_true
end
args = parse_args(s)

using ChromNet
using JSON
using JLD
using SamIO
using ProgressMeter

# open the data bundle and load the metadata
dataBundles = [jldopen(x) for x in args["data_bundles"]]
metadata = Dict()
for bundle in dataBundles
    metadata = merge(metadata, read(bundle, "metadata"))
end
ids = vcat([read(bundle, "ids") for bundle in dataBundles]...)
numBins = round(Int, ceil(sum(ReferenceContigs_hg38.sizes) / 1000))

# find the indexes of the groups this edge is between
group1Inds = map(x->findfirst(ids, x), split(args["group1"], "|"))
group2Inds = map(x->findfirst(ids, x), split(args["group2"], "|"))
@assert all(group1Inds .> 0) "Some ids from group 1 were not found: $(group1Inds[group1Inds .== 0])"
@assert all(group2Inds .> 0) "Some ids from group 2 were not found: $(group2Inds[group2Inds .== 0])"

# compute the joint correlation matrix
C = nothing
if typeof(args["load-cor"]) != Void
    C,header = readdlm(args["load-cor"], header=true)
    @assert all(header .== ids) "Loaded correlation matrix does not match data! Perhaps a different file name ordering?"
else
    C = streaming_cor(dataBundles, quiet=args["quiet"])
end
cov2cor!(C)

# compute the edge impact
edgeScores = streaming_edgeimpact(dataBundles, C, group1Inds, group2Inds, quiet=args["quiet"])
map(close, dataBundles)

# send the data in BED format to the output stream
outStream = STDOUT
if args["output"] != "-"
    outStream = open(args["output"], "w")
end
contigIndex = 1
pos = 1
ctgs = ReferenceContigs_hg38
for i in 1:numBins
    while pos > ctgs.offsets[contigIndex] + ctgs.sizes[contigIndex]
        contigIndex += 1
    end
    println(outStream, ctgs.names[contigIndex], "\t", pos - ctgs.offsets[contigIndex], "\t", pos - ctgs.offsets[contigIndex] + 999, "\t*\t", edgeScores[i])
    pos += 1000
end
