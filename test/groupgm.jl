using JSON

## build_groups
C = cov([M' M'])
groups = build_groups(C, [header, header])
for i in 1:4:length(groups[end][2])
    @test groups[end][2][i] == groups[end][2][i+2] # make sure duplicates are grouped together
end

## build_groupgm
C = cov([M[1,:]' M[1,:]' M[2,:]'.+0.01*M[1,:]' M[3,:]'] .+ 0.0001*randn(size(M)[2], 4))
groups = build_groups(C, header[1:4])
I = inv(C)
G,headerG = build_groupgm(I, header[1:4], groups)
@test I[1,3] + I[2,3] == G[5,3]

## parse_config
metadata = open(f->parse_config(f, "data"), "data/single_bed.config")
@test metadata["EX001"]["name"] == "K562_tss"
@test metadata["EX001"]["lab"] == "Sample lab"

## build_network
M = trues(4, 1050)
M[:,:] = rand(4, 1050) .> 0.5
header = ["A01", "B01", "C01", "D01"]
isfile("data/small.ChromNet") && rm("data/small.ChromNet")
save_chromnet_matrix("data/small.ChromNet", M, header)
metadata = {
    "A01" => {
        "name" => "A",
        "id" => "A01",
        "description" => "A desc",
        "cellType" => "A cellType",
        "lab" => "A cellType",
        "organism" => "Homo sapiens",
        "lifeStage" => "None",
        "treatments" => "None",
        "antibody" => "Unknown"
    },
    "B01" => {
        "name" => "B",
        "id" => "B01",
        "description" => "B desc",
        "cellType" => "B cellType",
        "lab" => "B cellType",
        "organism" => "Homo sapiens",
        "lifeStage" => "None",
        "treatments" => "None",
        "antibody" => "Unknown"
    },
    "C01" => {
        "name" => "C",
        "id" => "C01",
        "description" => "C desc",
        "cellType" => "C cellType",
        "lab" => "C cellType",
        "organism" => "Homo sapiens",
        "lifeStage" => "None",
        "treatments" => "None",
        "antibody" => "Unknown"
    },
    "D01" => {
        "name" => "D",
        "id" => "D01",
        "description" => "D desc",
        "cellType" => "D cellType",
        "lab" => "D cellType",
        "organism" => "Homo sapiens",
        "lifeStage" => "None",
        "treatments" => "None",
        "antibody" => "Unknown"
    }
}
C = cov(M')
groups = build_groups(C, header)
I = inv(C)
G,headerG = build_groupgm(I, header, groups)
network = build_network(G, header, groups, metadata; threshold=0.03, groupLinkThreshold=0.7)
