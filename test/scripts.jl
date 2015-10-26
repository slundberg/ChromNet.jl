using JSON

## build_network.jl
M = trues(4, ChromNet.totalBins)
M[:,:] = rand(4, ChromNet.totalBins) .> 0.5
header = ["A01", "B01", "C01", "D01"]
isfile("data/small.ChromNet") && rm("data/small.ChromNet")
save_chromnet_matrix("data/small.ChromNet", M, header)
metadata = Dict(
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
open(f->println(f, json(metadata)), "data/small.ChromNet/metadata", "w")
run(pipeline(`julia ../scripts/build_network.jl data/four_bed.config --data none --quiet`, "data/small.json"))
run(pipeline(`julia ../scripts/build_network.jl data/single_bed.config --data data/small.ChromNet --quiet`, "data/small.json"))
