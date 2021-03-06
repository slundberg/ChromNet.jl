export build_groups, build_groupgm, parse_config, build_network

using Clustering
using GZip

"Build GroupGM groups from a correlation matrix and a set of ids."
function build_groups(C, header)

    # force C to be perfectly symmetric (rounding errors make hclust unhappy)
    symC = deepcopy(C)
    for i in 1:size(C)[1], j in 1:size(C)[2]
        symC[j,i] = symC[i,j]
    end

    h = hclust(1 .- abs(symC), :complete)

    # build a list of all internal nodes in the hclust tree
    internalNodes = ASCIIString[]
    internalScores = Float64[]
    for i in 1:length(header)-1
        v1 = h.merge[i,1] < 0 ? header[-h.merge[i,1]] : internalNodes[h.merge[i,1]]
        v2 = h.merge[i,2] < 0 ? header[-h.merge[i,2]] : internalNodes[h.merge[i,2]]
        push!(internalNodes, "$v1 $v2") # note this keeps sub group names next to each other in larger group names
    end

    # build the results
    groups = Tuple{Float64, ASCIIString}[]
    for (i,v) in enumerate(header)
        push!(groups, (0.0, v))
    end
    for (i,v) in enumerate(internalNodes)
        push!(groups, (h.height[i], v))
    end

    groups
end

"Build a GroupGM score matrix form an inverse correlation matrix and the given groups."
function build_groupgm(I, header, groups; combiner=sum)

    # create a dense array and a map to its indexes
    G = zeros(length(groups), length(groups))
    indexMap = Dict{ASCIIString,Int64}()
    for i in 1:length(header)
        indexMap[header[i]] = i
    end

    # fill in our group matrix
    indexGroups = collect(map(g->collect(map(x->indexMap[x], split(g[2]))), groups))
    for i in 1:length(indexGroups)
        for j in i:length(indexGroups)
            G[j,i] = G[i,j] = combiner(I[indexGroups[j],indexGroups[i]])
        end
    end

    G,[g[2] for g in groups]
end

# each line is tab separated with the follow format (tab separated)
# BED_FILE_NAME SHORT_TITLE LONG_TITLE CELL_TYPE LAB EXPERIMENT_ID ANTIBODY_ID TREATMENTS ORGANISM LIFE_STAGE LINK
function parse_config(stream, fileRoot)
    metadata = Dict()
    for line in eachline(stream)
        parts = split(strip(line), '\t')
        l = length(parts)

        obj = Dict()
        if ismatch(r".*\.bam$", parts[1])
            obj["fileType"] = "bam"
        elseif ismatch(r".*\.bed$", parts[1])
            obj["fileType"] = "bed"
        elseif ismatch(r".*\.bed\.gz$", parts[1])
            obj["fileType"] = "bed.gz"
        elseif ismatch(r".*\.narrowPeak$", parts[1])
            obj["fileType"] = "narrowPeak"
        elseif ismatch(r".*\.narrowPeak\.gz$", parts[1])
            obj["fileType"] = "narrowPeak.gz"
        else
            obj["fileType"] = "unknown"
        end
        obj["dataFile"] = joinpath(fileRoot, parts[1])
        obj["name"] = l >= 2 ? parts[2] : basename(obj["dataFile"])
        obj["description"] = l >= 3 ? parts[3] : obj["name"]
        obj["cellType"] = l >= 4 ? parts[4] : "Unlabeled custom upload"
        obj["lab"] = l >= 5 ? parts[5] : "Unlabeled custom upload"
        obj["id"] = l >= 6 ? parts[6] : obj["name"]*"_"*randstring(15)
        obj["antibody"] = l >= 7 ? parts[7] : "Unlabeled custom upload"
        obj["treatments"] = l >= 8 ? parts[8] : "None"
        obj["organism"] = l >= 9 ? parts[9] : "Homo sapiens"
        obj["lifeStage"] = l >= 10 ? parts[10] : "None"
        obj["link"] = l >= 11 ? parts[11] : ""
        metadata[obj["id"]] = obj
    end
    metadata
end

"Build a network from a GroupGM matrix, groups, and metadata."
function build_network(G, groups, metadata; threshold=0.2, groupScoreThreshold=0.7, quiet=false)

    function find_parent(groups, id, ind)
        for i in ind:length(groups)
            if searchindex(groups[i][2], id) != 0
                return i
            end
        end
        0
    end

    # build the node data
    nodes = Any[]
    for (i,id) in enumerate(map(x->x[2], groups))
        if search(id, ' ') == 0
            md = metadata[id]
            push!(nodes, Dict{Any,Any}(
                "name" => get(md, "name", id),
                "group" => 0,
                "groupScore" => 1.0,
                "parent" => find_parent(groups, id, i+1)-1,
                "data" => Dict{Any,Any}(
                    "id" => id,
                    "description" => get(md, "description", get(md, "name", id)),
                    "cellType" => get(md, "cellType", "Unknown"),
                    "lab" => get(md, "lab", "Unknown"),
                    "organism" => get(md, "organism", "Unknown"),
                    "lifeStage" => get(md, "lifeStage", "Unknown"),
                    "treatments" => get(md, "treatments", "None"),
                    "antibody" => get(md, "antibody", "Unknown")
                )
            ))
        else
            push!(nodes, Dict{Any,Any}(
                "name" => "$i",
                "group" => 1,
                "groupScore" => 1-groups[i][1],
                "parent" => find_parent(groups, id, i+1)-1,
                "data" => Dict{Any,Any}(
                    "id" => "$i"
                )
            ))
        end
    end

    allNodes = Array{ASCIIString,1}[]
    for i in 1:length(groups)
        push!(allNodes, ASCIIString[])
    end
    for (i,(v,id)) in enumerate(groups)

        # if we are a leaf node then add our id to all ancestors (including ourselves)
        if !(' ' in id)
            ind = i
            while ind != 0
                push!(allNodes[ind], id)
                ind = nodes[ind]["parent"]+1
            end
        end
    end

    # build the edge array
    links = Any[]
    for i in 1:length(nodes)
        idi = nodes[i]["data"]["id"]
        for j in i+1:length(nodes)
            idj = nodes[j]["data"]["id"]

            # make sure we pass the threshold and are not between nested groups
            if (abs(G[i,j]) > threshold && searchindex(groups[i][2], groups[j][2]) == 0
                && searchindex(groups[j][2], groups[i][2]) == 0
                && groups[i][1] < 1-groupScoreThreshold && groups[j][1] < 1-groupScoreThreshold)

                d = Dict{Any,Any}(
                    "source" => i-1,
                    "target" => j-1,
                    "coeff" => -G[i,j],
                    "labels" => get(metadata, "$idi,$idj", ASCIIString[])
                )

                # find all the labels that are represented by the edges between these groups
                labels = ASCIIString[]
                for id1 in allNodes[i]
                    for id2 in allNodes[j]
                        if haskey(metadata, "$id1,$id2")
                            push!(labels, metadata["$id1,$id2"]...)
                        end
                    end
                end

                if length(labels) > 0
                    d["labels"] = unique(labels)
                end

                push!(links, d)
            end
        end
    end

    if !quiet
        println(STDERR, "Created network with $(length(links)) links and $(length(nodes)) nodes.")
    end
    Dict("nodes" => nodes, "links" => links)
end
