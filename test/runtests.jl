using ChromNet
using JLD
using Base.Test
srand(1)

# create a sample network dataset to use
jldopen("data/small.ChromNet.jld", "w") do f
    data = zeros(Int32, 4, 1000000)
    data[:,:] = rand(0:100, 4, 1000000)
    ids = ["A01", "B01", "C01", "D01"]
    f["P"] = 4
    f["N"] = 1000000
    f["ids"] = ids
    f["data"] = data
    f["means"] = mean(data, 2)
    f["norms"] = Float64[norm(data[i,:] - mean(data[i,:])) for i in 1:4]
    f["metadata"] = Dict(
        "A01" => Dict(
            "name" => "A",
            "id" => "A01",
            "description" => "A desc",
            "cellType" => "A cellType",
            "lab" => "A cellType",
            "organism" => "Homo sapiens",
            "lifeStage" => "None",
            "treatments" => "None",
            "antibody" => "Unknown"
        ),
        "B01" => Dict(
            "name" => "B",
            "id" => "B01",
            "description" => "B desc",
            "cellType" => "B cellType",
            "lab" => "B cellType",
            "organism" => "Homo sapiens",
            "lifeStage" => "None",
            "treatments" => "None",
            "antibody" => "Unknown"
        ),
        "C01" => Dict(
            "name" => "C",
            "id" => "C01",
            "description" => "C desc",
            "cellType" => "C cellType",
            "lab" => "C cellType",
            "organism" => "Homo sapiens",
            "lifeStage" => "None",
            "treatments" => "None",
            "antibody" => "Unknown"
        ),
        "D01" => Dict(
            "name" => "D",
            "id" => "D01",
            "description" => "D desc",
            "cellType" => "D cellType",
            "lab" => "D cellType",
            "organism" => "Homo sapiens",
            "lifeStage" => "None",
            "treatments" => "None",
            "antibody" => "Unknown"
        )
    )
end

include("util.jl")
include("groupgm.jl")
include("scripts.jl")
