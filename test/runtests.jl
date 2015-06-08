using ChromNet
using Base.Test
using GZip


binValues = GZip.open(window_bed_file, "data/K562_tss.bed.gz")

data = spzeros(ChromNet.totalBins,1)
GZip.open(f->window_bed_file!(data, 1, f), "data/K562_tss.bed.gz")