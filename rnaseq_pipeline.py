#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 17:09:47 2020

@author: jakobsen
"""

from ruffus import *
from cgatcore import pipeline as P
import sys
import os
import cgatcore.iotools as IOTools

PARAMS = P.get_parameters("pipeline.yml")


@follows(mkdir('1_fastqc'))
@transform('fastq/*q.gz', regex(r'fastq/(.*).f.*q.gz'), r'1_fastqc/\1_fastqc.zip')
def fastqc(input_file, output_file):
    statement = 'fastqc -t 2 %(input_file)s -o 1_fastqc'
    P.run(statement, job_queue=PARAMS['q'], job_threads=2)

### Merge lanes ###

## Merge fastq files where they separated by lane
@follows(mkdir('2_trimmed'))
@follows(fastqc)
@active_if(P.get_params()['input'] == "PE_barcoded")
@collate('fastq/*q.gz', regex(r'fastq/(.*?)_(.*)_(L00\d)_(R[12])_001.f.*q.gz'), r'2_trimmed/\1.merged_\4.fq.gz')
def merge_lanes(input_files, output_file):
    files = ' '.join(str(x) for x in input_files)
    statement = '''cat %(files)s
    > %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=1)
    
    
##### Trimming #####

@follows(mkdir('2_trimmed'))
## SE_demultiplexed: single-end reads with file extension in the format: R1_001.fastq.gz
@active_if(P.get_params()['input'] == "SE_demultiplexed")
@transform('fastq/*.fastq.gz', regex(r'fastq/(.*)_(.*)_R1_001.fastq.gz'), r'2_trimmed/\1/\1.trimmed.fq.gz')
def trimming_SE(input_file, output_file):
    basename = P.snip(os.path.basename(input_file),"_R1_001.fastq.gz").split("_")[0]
    statement = '''mkdir 2_trimmed/%(basename)s &&
    mkdir -p 3_mapping/%(basename)s &&
    mkdir -p 4_mapping_qc/%(basename)s &&
    cutadapt --cores=0
    --nextseq-trim=%(cutadapt_q)s
    -a "A{20}N{80};min_overlap=3" -a AGCAACTCTGCGTTGATACCACTGCTT
    %(cutadapt_options)s
    -o 2_trimmed/%(basename)s/%(basename)s.trimmed.fq.gz
    %(input_file)s
    > 2_trimmed/%(basename)s_trimming_report.txt
    '''
    P.run(statement, job_queue=PARAMS['q'], job_threads=4)


## PE_barcoded: paired-end reads with file extension in the format: L001_R1_001.fastq.gz
@follows(mkdir('2_trimmed'))
@follows(merge_lanes)
@active_if(P.get_params()['input'] == "PE_barcoded")
@subdivide('2_trimmed/*merged_R1.fq.gz', regex(r"2_trimmed/(.*).merged_R1.fq.gz"), output = r"2_trimmed/\1/\1_HT*.trimmed_R1.fq.gz")
def trimming_PE(input_file, output_files):
    input_file2 = input_file.replace("_R1.fq.gz", "_R2.fq.gz")
    outfile_prefix = P.snip(input_file,".merged_R1.fq.gz").split(".")[0]
    basename = P.snip(os.path.basename(input_file),".merged_R1.fq.gz").split("_")[0]
    #lane_suffix = input_file,r"_(.*)_L00\d_R1.fq.gz").split("_")[0]
    statement = '''mkdir 2_trimmed/%(basename)s &&
    mkdir -p 3_mapping/%(basename)s &&
    mkdir -p 4_mapping_qc/%(basename)s &&
    cutadapt --cores=0
    -g file:%(cutadapt_barcodes)s
    --nextseq-trim=%(cutadapt_q)s
    -A "A{20}N{80};min_overlap=3" -A AGCAACTCTGCGTTGATACCACTGCTT
    %(cutadapt_options)s
    -o 2_trimmed/%(basename)s/%(basename)s_{name}.trimmed_R1.fq.gz
    -p 2_trimmed/%(basename)s/%(basename)s_{name}.trimmed_R2.fq.gz
    %(input_file)s
    %(input_file2)s
    > %(outfile_prefix)s_trimming_report.txt'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=4)
    if P.get_params()["zap_files"]==1:
        IOTools.zap_file(input_file)
        IOTools.zap_file(input_file2)




##### Mapping #####

@follows(mkdir('3_mapping'))
@follows(trimming_SE)
@follows(trimming_PE)
@follows(merge_lanes)
@transform('2_trimmed/*/*.fq.gz', regex(r"2_trimmed/(.*)/(.*).trimmed(_R2)?.fq.gz"), r"3_mapping/\1/\2.bam")
def star(input_file, output_file):
    outprefix = P.snip(output_file, ".bam")
    # SE mapping
    if P.get_params()['mapping'] == "SE_demultiplexed":
        statement = '''STAR
        --runThreadN %(star_threads)s
        --genomeDir %(star_ref)s
        --readFilesIn %(input_file)s
        --readFilesCommand zcat
        --outStd SAM
        --outSAMunmapped Within
        --outFileNamePrefix %(outprefix)s_
        | samtools view -bu | samtools sort -@ %(star_threads)s -o %(output_file)s'''
        P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_memory = '30G')
        if P.get_params()["zap_files"]==1:
            IOTools.zap_file(input_file)
    
    # SE mapping of PE_barcoded reads
    elif P.get_params()['input'] == "PE_barcoded" and P.get_params()['mapping'] == "read2":
        read1 = input_file.replace("_R2.fq.gz", "_R1.fq.gz")
        #basename = P.snip(os.path.basename(input_file),".trimmed_2.fq.gz").split("_")[0]
        statement = '''STAR
        --runThreadN %(star_threads)s
        --genomeDir %(star_ref)s
        --readFilesIn %(input_file)s
        --readFilesCommand zcat
        --outStd SAM
        --outSAMunmapped Within
        --outFileNamePrefix %(outprefix)s_
        | samtools view -bu | samtools sort -@ %(star_threads)s -o %(output_file)s'''
        P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_memory = '30G')
        if P.get_params()["zap_files"]==1:
            IOTools.zap_file(input_file)
            IOTools.zap_file(read1)
            
    # PE mapping
    elif P.get_params()['input'] == "PE_barcoded" and P.get_params()['mapping'] == "both":
        read1 = input_file.replace("_R2.fq.gz", "_R1.fq.gz")
        statement = '''STAR
        --runThreadN %(star_threads)s
        --genomeDir %(star_ref)s
        --readFilesIn %(read1)s %(input_file)s
        --readFilesCommand zcat
        --outStd SAM
        --outSAMunmapped Within
        --outFileNamePrefix %(outprefix)s_
        | samtools view -bu | samtools sort -@ %(star_threads)s -o %(output_file)s'''
        P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_memory = '30G')
        if P.get_params()["zap_files"]==1:
            IOTools.zap_file(input_file)
            IOTools.zap_file(read1)
            
    else:
        print("Incorrect FASTQ input parameter")


@transform(star, suffix('.bam'), '.bam.bai')     
def bam_index(input_file, output_file):
    statement = 'samtools index %(input_file)s -@ 1 2> %(output_file)s.log'
    P.run(statement, job_queue=PARAMS['q'], job_threads=1, job_memory = '100M')


##### Mapping QC #####
    
@follows(bam_index)
@follows(mkdir('4_mapping_qc'))
@transform(star, regex(r'3_mapping/(.*)/(.*).bam'), r'4_mapping_qc/\1/\2.idxstat')
def samtools_idxstat(input_file, output_file):
    statement = '''samtools idxstats %(input_file)s > %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '100M')

@follows(bam_index, samtools_idxstat)
@transform(star, regex(r'3_mapping/(.*)/(.*).bam'), r'4_mapping_qc/\1/\2.flagstat')
def samtools_flagstat(input_file, output_file):
    statement = '''samtools flagstat %(input_file)s > %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '100M')


##### Featurecounts #####
    
@follows(bam_index)
@follows(mkdir('5_featurecounts'))
@merge(star, '5_featurecounts/counts.txt')
def count_reads(input_files, output_file):
    input_files_string = ' '.join(input_files)
    statement = '''featureCounts -T 12 -t exon -g gene_id
    -a %(featurecounts_gtf)s -o %(output_file)s %(input_files_string)s
    %(featurecounts_options)s
    2> 5_featurecounts/featurecounts_log.txt'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')

## Multiqc report ##

@follows(count_reads, samtools_idxstat)
#count_reads only there to specify it should be done at end
@merge(samtools_flagstat, '6_multiqc/multiqc_report.html')
def multiqc(input_file, output_file):
    statement = '''multiqc . -o 6_multiqc'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G')


@follows(fastqc, trimming_SE, trimming_PE)
def Trimming():
    pass
    
@follows(star, samtools_flagstat, samtools_idxstat)
def Mapping_qc():
    pass
    
@follows(Mapping_qc, count_reads, multiqc)
def full():
    pass

if __name__=="__main__":
    sys.exit(P.main(sys.argv))
    
    
    
