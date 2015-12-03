using JSON

## build_groups
f = jldopen("data/small.ChromNet.jld")
data = read(f, "data")
ids = read(f, "ids")
metadata = read(f, "metadata")
close(f)
C = cov([data' data'])
groups = build_groups(C, [ids; ids])
for i in 1:8:length(groups[end][2])
    @test groups[end][2][i] == groups[end][2][i+4] # make sure duplicates are grouped together
end

## build_groupgm
C = cov([data[1,:]' data[1,:]' data[2,:]'.+0.01*data[1,:]' data[3,:]'] .+ 0.0001*randn(size(data)[2], 4))
groups = build_groups(C, ids[1:4])
I = inv(C)
G,idsG = build_groupgm(I, ids[1:4], groups)
@test abs(I[1,3] + I[2,3] - G[5,3]) < 1e-7

## parse_config
metadataTmp = open(f->parse_config(f, "data"), "data/single_bed.config")
@test metadataTmp["EX001"]["name"] == "K562_tss"
@test metadataTmp["EX001"]["lab"] == "Sample lab"

## build_network
C = cor(data, vardim=2)
groups = build_groups(C, ids)
G,idsG = build_groupgm(inv(C), ids, groups)
network = build_network(G, groups, metadata; threshold=0.03, groupScoreThreshold=0.7)
