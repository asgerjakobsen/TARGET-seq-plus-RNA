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
These need a file extension in the format: `_R[1-2]_001.fq.gz`.
In the `pipeline.yml`, specify: `PE_barcoded`

If FASTQ files are split by lane, the pipeline can merge them. In the `pipeline.yml`, specify: `lanes: split`


### Trimming

Cutadapt trims poly(A) tails and ISPCR adapters, as well as bases with a quality score below 15.

### Demultiplexing

For paired-end reads where Read 1 is a cell barcode, cutadapt will demultiplex the cell barcode. 
This is done using the `barcodes.fasta` file, which needs to be copied into the working directory. 
Reads will be split and named according to the barcode: "HTxxx".

### Mapping

The pipeline maps using STAR. STAR needs a GTF with gene annotations for mapping and gene assignment. If using ERCC spike-ins, these need to be appended to the genome fasta file before building the STAR index.
Currently mapping on a single read works.
For paired-end reads where Read 1 is a cell barcode, please specify `mapping: read2` in the `pipeline.yml`.

### Counting

Feature (gene) counting is done with either Featurecounts or STARsolo.
For featurecounts, mapping is performed first in a separate step for each cell.
For STARsolo, mapping and counting are performed in a single step for all cells. This is faster.


### Usage

For Featurecounts:
```
nohup python rnaseq_pipeline.py make full_featurecounts &
```
For STARsolo:
```
nohup python rnaseq_pipeline.py make full_starsolo &
```
