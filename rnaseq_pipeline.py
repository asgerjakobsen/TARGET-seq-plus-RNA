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
@transform('fastq/*.fastq.gz', regex(r'fastq/(.*).fastq.gz'), r'1_fastqc/\1_fastqc.zip')
def fastqc(input_file, output_file):
    statement = 'fastqc -t 8 %(input_file)s -o 1_fastqc'
    P.run(statement, job_queue=PARAMS['q'], job_threads=8)


##### Trimming #####
    
@follows(mkdir('2_trimmed'))
@follows(mkdir('3_trimmed_fastqc'))
@transform('fastq/*.fastq.gz', regex(r'fastq/(.*)_(.*)_R1_001.fastq.gz'), r'2_trimmed/\1_trimmed.fq.gz')
def trimming(input_file, output_file):
    basename = P.snip(os.path.basename(input_file),"_R1_001.fastq.gz").split("_")[0]
    statement = '''trim_galore --cores 4 -q %(trimgalore_q)s -a "A{100}"
    %(trimgalore_options)s
    %(input_file)s  --basename %(basename)s -o 2_trimmed
    --fastqc_args "-t 4 -o 3_trimmed_fastqc" '''
    P.run(statement, job_queue=PARAMS['q'], job_threads=4)


##### Mapping #####
    
@follows(mkdir('4_mapping'))
@transform(trimming, regex(r"2_trimmed/(.*)_trimmed.fq.gz"), r"4_mapping/\1.bam") 
def star(input_file, output_file):
    outprefix = P.snip(output_file, ".bam")
    statement = '''STAR 
    --runThreadN %(star_threads)s
    --genomeDir %(star_ref)s
    --readFilesIn %(input_file)s 
    --readFilesCommand zcat
    --outStd SAM
    --outSAMunmapped Within
    --outFileNamePrefix %(outprefix)s_
    | samtools view -bu | samtools sort -@ %(star_threads)s -o %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_memory = '8G')
    if P.get_params()["zap_files"]==1:
        IOTools.zap_file(input_file)
       
@transform(star, suffix('.bam'), '.bam.bai')     
def bam_index(input_file, output_file):
    statement = 'samtools index %(input_file)s -@ 12 2> %(output_file)s.log'
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')  

##### Mapping QC #####
    
@follows(bam_index)
@follows(mkdir('5_mapping_qc'))  
@transform(star, regex(r'4_mapping/(.*).bam'), r'5_mapping_qc/\1.idxstat')     
def samtools_idxstat(input_file, output_file):
    statement = '''samtools idxstats %(input_file)s > %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G') 

@follows(bam_index, samtools_idxstat)    
@transform(star, regex(r'4_mapping/(.*).bam'), r'5_mapping_qc/\1.alignment_metrics.txt')     
def alignment_summary_metrics(input_file, output_file):
    statement = '''picard CollectAlignmentSummaryMetrics R=%(picard_ref)s I=%(input_file)s O=%(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '16G') 

@follows(bam_index, samtools_idxstat)
@transform(star, regex(r'4_mapping/(.*).bam'), r'5_mapping_qc/\1.flagstat')     
def samtools_flagstat(input_file, output_file):
    statement = '''samtools flagstat %(input_file)s > %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G')


##### Featurecounts #####
    
@follows(bam_index)
@follows(mkdir('6_featurecounts'))
@merge(star, '6_featurecounts/featurecounts.txt')     
def count_reads(input_files, output_file):
    input_files_string = ' '.join(input_files)
    statement = '''featureCounts -T 12 -t exon -g gene_id
    -a %(featurecounts_gtf)s -o %(output_file)s %(input_files_string)s
    %(featurecounts_options)s'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')

  
#count_reads only there to specify it should be done at end
@follows(count_reads, samtools_idxstat, alignment_summary_metrics, fastqc)
@merge(samtools_flagstat, '7_multiqc/multiqc_report.html')    
def multiqc(input_file, output_file):
    statement = '''multiqc . -o 7_multiqc'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')



@follows(samtools_flagstat, samtools_idxstat, alignment_summary_metrics, fastqc)
def Mapping_qc():
    pass

if __name__=="__main__":
    sys.exit(P.main(sys.argv))
    
    
    
