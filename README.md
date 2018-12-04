
# ChromNet.jl

Integrates any BAM, BED, or narrowPeak file into the ChromNet group graphical model, the JSON network file produced can then be explored using http://chromnet.cs.washington.edu. The genomic context driving any group edge in the network can also be calculated. Details are available in our paper at http://www.genomebiology.com/2016/17/1/82.

[![Build Status](https://travis-ci.org/slundberg/ChromNet.jl.svg?branch=master)](https://travis-ci.org/slundberg/ChromNet.jl)


## Installation

- Ensure that a compatible version of [Julia](http://www.julialang.com/downloads/) is installed (v0.4 - v0.5). ChromNet is written in Julia for performance and portability reasons.
- Download the [current data package](https://drive.google.com/uc?export=download&id=0B8QcnMD1YRTXRnJFVy1BSkw0bW8) (currently at build 3) which contains the ChromNet processed version of all ENCODE ChIP-seq data. Also contained in this package is a Julia script, which when run, will generate a new ChromNet model using the ENCODE data and any user-provided BAM/BED files.

#### Upgrading
If you have already used ChromNet before, ensure you have the latest data package, and then run `Pkg.update()` in the Julia console to fetch any code updates. In a UNIX shell code updates can be fetched with:
```shell
julia -e 'Pkg.update()'
```

## Usage

Below are examples of basic usage for each command. Running each command with the `--help` option will give more detailed documentation.

### Build sorted BAM file aligned to hg38

ChromNet expects sorted BAM files aligned to hg38. Use our provided [Bowtie index](https://drive.google.com/open?id=0B8QcnMD1YRTXSzVCUEhYMFREXzQ) to ensure your data matches the ChromNet pre-built encode data bundle. You can create these using the following commands:

```shell
zcat MY_EXP.fastq.gz | bowtie2 -p 20 -x /path/bowtie/hg38 -U - | samtools view -bS - | MY_EXP.unsorted.bam;
samtools sort MY_EXP.unsorted.bam -o MY_EXP.sorted.bam -@ 10
```

### Build a custom data bundle

To build a network from custom data, a custom data bundle must be generated. To do this, decompress the downloaded data package and from inside the directory run:

```bash
julia build_bundle.jl CONFIG_FILE -o custom_bundle_name.ChromNet.jld
```

The config file lists custom BAM/BED files to be incorporated into the network. Note all BAM and BED files must be aligned to GRCh38 or the `--assembly` option must be specified (see julia build_bundle.jl --help for available reference assemblies). Each line of the config file should conform to the following TAB separated format, where trailing entries can be omitted:

```
FILE_NAME SHORT_TITLE LONG_TITLE CELL_TYPE LAB EXPERIMENT_ID ANTIBODY_ID TREATMENTS ORGANISM LIFE_STAGE LINK
```

The simplest configuration file is just a list of file names, and `build_bundle.jl` supports STDIN and STDOUT streaming. This means a one line invocation on UNIX systems is simply (where '-' denotes STDIN):

```shell
ls ~/my_bed_files/*.bam | julia build_bundle.jl custom_bundle_name -
```

### Build a custom network

Given a set of data bundles a single ChromNet network that incorperates data from all bundles can be generated using the following command:

```bash
julia build_network.jl custom_bundle_name.ChromNet.jld ENCODE_human_build3.ChromNet.jld > network.json
```

The output JSON file can then be dropped into the ChromNet interface at http://chromnet.cs.washington.edu.

### Calculate edge context

To calculate the genomic context driving a specific group edge you can use the following command:

```bash
julia compute_edge_context.jl ENCSR000DMA|ENCSR000EGM ENCSR000BGX custom_bundle_name.ChromNet.jld ENCODE_human_build3.ChromNet.jld > out.bed
```

## BED format

ChromNet only relies on the first three fields of the BED format (chrom, chromStart, and chromEnd). This means other tab delimited formats that follow the same conventions are also compatabile with `build_bundle.jl`. This includes the narrowPeak format produced by MACS2 and other peak calling software.
