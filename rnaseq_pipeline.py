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

PARAMS = P.get_parameters("pipeline.yml")


@follows(mkdir('fastqc'))
@transform('fastq/*.fastq.gz', regex(r'fastq/(.*).fastq.gz'), r'fastqc/\1_fastqc.zip')
def fastqc(input_file, output_file):
    statement = 'fastqc -t 8 %(input_file)s -o fastqc'
    P.run(statement, job_queue=PARAMS['q'], job_threads=8)


##### Trimming #####
    
@follows(mkdir('trimmed'))
@follows(mkdir('trimmed_fastqc'))
@transform('fastq/*.fastq.gz', regex(r'fastq/(.*)_(.*)_R1_001.fastq.gz'), r'trimmed/\1_trimmed.fq.gz')
def trimming(input_file, output_file):
    basename = P.snip(os.path.basename(input_file),"_R1_001.fastq.gz").split("_")[0]
    statement = '''trim_galore --cores 4 -q %(trimgalore_q)s -a "A{100}" 
    %(input_file)s  --basename %(basename)s -o trimmed
    --fastqc_args "-t 4 -o trimmed_fastqc" '''
    P.run(statement, job_queue=PARAMS['q'], job_threads=4)


##### Mapping #####
    
@follows(mkdir('STAR'))
@transform(trimming, regex(r"trimmed/(.*)_trimmed.fq.gz"), r"STAR/\1.bam") 
def star(input_file, output_file):
    outprefix = P.snip(output_file, ".bam")
    statement = '''STAR 
    --runThreadN %(star_threads)s
    --genomeDir %(star_ref)s
    --readFilesIn %(input_file)s 
    --readFilesCommand zcat
    --outSAMunmapped Within
    --outFileNamePrefix %(outprefix)s_
    | samtools view -bu | samtools sort -@ %(star_threads)s -o %(output_file)s'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['hisat2_t'], job_memory = '8G')
       
@transform(star, suffix('.bam'), '.bam.bai')     
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
@follows(count_reads, samtools_idxstat, alignment_summary_metrics, fastqc)
@merge(samtools_flagstat, 'multiqc/multiqc_report.html')    
def multiqc(input_file, output_file):
    statement = '''multiqc . -o multiqc'''
    P.run(statement, job_queue=PARAMS['q'], job_threads=12, job_memory = '8G')



@follows(samtools_flagstat, samtools_idxstat, alignment_summary_metrics, fastqc)
def Mapping_qc():
    pass

if __name__=="__main__":
    sys.exit(P.main(sys.argv))
    
    
    