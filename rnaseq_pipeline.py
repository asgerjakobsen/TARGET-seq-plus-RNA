#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 17:09:47 2020

@author: jakobsen
"""

from ruffus import *
from cgatcore import pipeline as P
import sys

PARAMS = P.get_parameters("pipeline.yml")


@follows(mkdir('fastqc'))
@transform('fastq/*.fastq.gz', regex(r'fastq/(.*).fastq.gz'), r'fastqc/\1_fastqc.zip')
def fastqc(input_file, output_file):
    statement = 'fastqc -t 8 %(input_file)s -o fastqc'
    P.run(statement, job_queue=PARAMS['q'], job_threads=8)

##### Mapping #####
    
@follows(mkdir('hisat2'))
@transform("fastq/*_R1_001.fastq.gz", regex(r"fastq/(.*)_R1_001.fastq.gz"), r"hisat2/\1.bam") 
def hisat2(input_file, output_file):
    outprefix = P.snip(output_file, ".bam")
    statement = '''hisat2 -x %(hisat2_ref)s 
    -U %(input_file)s -p %(hisat2_t)s 2> %(outprefix)s.log 
    | samtools sort -@ %(hisat2_t)s -o %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['hisat2_t'], job_memory = '8G')
       
@transform(hisat2, suffix('.bam'), '.bam.bai')     
def bam_index(input_file, output_file):
    statement = 'samtools index %(input_file)s -@ 12 2> %(output_file)s.log'
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')  

##### Mapping QC #####
    
@follows(bam_index)
@follows(mkdir('mapping_qc'))  
@transform(hisat2, regex(r'hisat2/(.*).bam'), r'mapping_qc/\1.idxstat')     
def samtools_idxstat(input_file, output_file):
    statement = '''samtools idxstats %(input_file)s > %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G') 

@follows(bam_index, samtools_idxstat)    
@transform(hisat2, regex(r'hisat2/(.*).bam'), r'mapping_qc/\1.alignment_metrics.txt')     
def alignment_summary_metrics(input_file, output_file):
    statement = '''picard CollectAlignmentSummaryMetrics R=%(picard_ref)s I=%(input_file)s O=%(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '16G') 

@follows(bam_index, samtools_idxstat)
@transform(hisat2, regex(r'hisat2/(.*).bam'), r'mapping_qc/\1.flagstat')     
def samtools_flagstat(input_file, output_file):
    statement = '''samtools flagstat %(input_file)s > %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G')


##### Featurecounts #####
    
@follows(bam_index)
@merge(hisat2, 'counts.txt')     
def count_reads(input_files, output_file):
    input_files_string = ' '.join(input_files)
    statement = '''featureCounts -T 12 -t exon -g gene_id --primary
    -a %(featurecounts_gtf)s -o %(output_file)s %(input_files_string)s'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')

  
#count_reads only there to specify it should be done at end
@follows(samtools_flagstat, samtools_idxstat, alignment_summary_metrics, fastqc)
@transform(count_reads, suffix('.txt'), '.html')    
def multiqc(input_file, output_file):
    statement = '''multiqc -o %(output_file)s .'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')



@follows(samtools_flagstat, samtools_idxstat, alignment_summary_metrics, fastqc)
def Mapping_qc():
    pass

if __name__=="__main__":
    sys.exit(P.main(sys.argv))
    
    
    