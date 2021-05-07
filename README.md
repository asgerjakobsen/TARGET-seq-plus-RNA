# RNA-seq pipeline

A pipeline for pre-processing of RNA-seq data.

## Setup

1. Clone the git repository
2. Install the conda environment using:
```
nohup conda env create -f rnaseq_env.yml &
```

## Using the pipeline

### Input files

FASTQ files (or symlinks) need to be in a subdirectory called `fastq/` within the working directory where you are running the pipeline.

The pipeline can deal with two types of fastq files:
1. Single-ended reads without inline barcodes (already demultiplexed).
These need a file extension in the format: `R1_001.fastq.gz`.
In the `pipeline.yml`, specify: `SE_demultiplexed`

2. Paired-end reads where Read 1 is a cell barcode.
These need a file extension in the format: `L00[1-4]_R[1-2]_001.fq.gz`.
In the `pipeline.yml`, specify: `PE_barcoded`


### Trimming

Cutadapt trims poly(A) tails and ISPCR adapters, as well as bases with a quality score below 15.

### Demultiplexing

For paired-end reads where Read 1 is a cell barcode, cutadapt will demultiplex the cell barcode. 
This is done using the `barcodes.fasta` file, which needs to be copied into the working directory. 
Reads will be split and named according to the barcode: "HTxxx".

### Mapping

Currently mapping on a single read works.
For paired-end reads where Read 1 is a cell barcode, please specify `mapping: read2` in the `pipeline.yml`.
STAR and FeatureCounts need a GTF with gene annotations for mapping and gene assignment.
