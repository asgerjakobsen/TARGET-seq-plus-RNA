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
@active_if(P.get_params()['lanes'] == "split")
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
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=4)


## PE_barcoded merged lanes: paired-end reads with file extension in the format: _R1_001.fastq.gz
@follows(mkdir('2_trimmed'))
@active_if(P.get_params()['input'] == "PE_barcoded" and P.get_params()['lanes'] == "merged")
@subdivide('fastq/*_R1_001.fastq.gz', regex(r"fastq/(.*)_R1_001.fastq.gz"), output = r"2_trimmed/\1/\1_HT*.trimmed_R1.fq.gz")
def trimming_PE_merged(input_file, output_files):
    input_file2 = input_file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
    basename = P.snip(os.path.basename(input_file),"_R1_001.fastq.gz").split("_")[0]
    statement = '''mkdir 2_trimmed/%(basename)s &&
    mkdir -p 3_mapping/%(basename)s &&
    mkdir -p 4_mapping_qc/%(basename)s &&
    cutadapt --cores=4
    -g file:%(cutadapt_barcodes)s
    --nextseq-trim=%(cutadapt_q)s
    -A "A{20}N{80};min_overlap=3" -A AGCAACTCTGCGTTGATACCACTGCTT
    %(cutadapt_options)s
    -o 2_trimmed/%(basename)s/%(basename)s_{name}.trimmed_R1.fq.gz
    -p 2_trimmed/%(basename)s/%(basename)s_{name}.trimmed_R2.fq.gz
    %(input_file)s
    %(input_file2)s
    > 2_trimmed/%(basename)s_trimming_report.txt'''
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=4)


## PE_barcoded split lanes: paired-end reads with file extension in the format: _L001_R1_001.fastq.gz
@follows(mkdir('2_trimmed'))
@follows(merge_lanes)
@active_if(P.get_params()['input'] == "PE_barcoded" and P.get_params()['lanes'] == "split")
@subdivide('2_trimmed/*merged_R1.fq.gz', regex(r"2_trimmed/(.*).merged_R1.fq.gz"), output = r"2_trimmed/\1/\1_HT*.trimmed_R1.fq.gz")
def trimming_PE_split(input_file, output_files):
    input_file2 = input_file.replace("_R1.fq.gz", "_R2.fq.gz")
    outfile_prefix = P.snip(input_file,".merged_R1.fq.gz").split(".")[0]
    basename = P.snip(os.path.basename(input_file),".merged_R1.fq.gz").split("_")[0]
    #lane_suffix = input_file,r"_(.*)_L00\d_R1.fq.gz").split("_")[0]
    statement = '''mkdir 2_trimmed/%(basename)s &&
    mkdir -p 3_mapping/%(basename)s &&
    mkdir -p 4_mapping_qc/%(basename)s &&
    cutadapt --cores=4
    -g file:%(cutadapt_barcodes)s
    --nextseq-trim=%(cutadapt_q)s
    -A "A{20}N{80};min_overlap=3" -A AGCAACTCTGCGTTGATACCACTGCTT
    %(cutadapt_options)s
    -o 2_trimmed/%(basename)s/%(basename)s_{name}.trimmed_R1.fq.gz
    -p 2_trimmed/%(basename)s/%(basename)s_{name}.trimmed_R2.fq.gz
    %(input_file)s
    %(input_file2)s
    > %(outfile_prefix)s_trimming_report.txt'''
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=4)
    if P.get_params()["zap_files"]==1:
        IOTools.zap_file(input_file)
        IOTools.zap_file(input_file2)



##### Mapping #####

@follows(mkdir('3_mapping'))
@follows(trimming_SE)
@follows(trimming_PE_merged)
@follows(trimming_PE_split)
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
        job_options = " -t 24:00:00"
        P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_total_memory = '30G')
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
        job_options = " -t 24:00:00"
        P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_total_memory = '30G')
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
        job_options = " -t 24:00:00"
        P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_total_memory = '30G')
        if P.get_params()["zap_files"]==1:
            IOTools.zap_file(input_file)
            IOTools.zap_file(read1)
            
    else:
        print("Incorrect FASTQ input parameter")


@transform(star, suffix('.bam'), '.bam.bai')     
def bam_index(input_file, output_file):
    statement = 'samtools index %(input_file)s -@ 1 2> %(output_file)s.log'
    job_options = " -t 12:00:00"
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
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=4, job_memory =PARAMS['featurecounts_memory'])


##### STARsolo #####
    
@follows(mkdir('7_starsolo'))
@follows(trimming_SE)
@follows(trimming_PE_merged)
@follows(trimming_PE_split)
@follows(merge_lanes)
@merge('2_trimmed/*/*R2.fq.gz', '7_starsolo/Aligned.out.bam')
def STARsolo(input_files, output_file):
    input_files_string = ','.join(input_files)
    with IOTools.open_file('7_starsolo/manifest.tsv', "w") as outf:
        for x in range (0,len(input_files)):
            filename=str(input_files[x].split("/")[2])
            sample = P.snip(filename, ".trimmed_R2.fq.gz")
            line = input_files[x]  + "\t-\t" + sample + "\n" 
            outf.write(line)
    
    outprefix =  "7_starsolo/" #P.snip(output_file, ".bam")
    statement = '''STAR
    --runThreadN %(star_threads)s
    --genomeDir %(star_ref)s
    --readFilesManifest 7_starsolo/manifest.tsv
    --readFilesCommand zcat
    --soloType SmartSeq   
    --soloFeatures Gene GeneFull_Ex50pAS SJ
    --soloUMIdedup NoDedup
    --soloStrand Forward
    --soloMultiMappers Unique
    --limitOutSJcollapsed 10000000
    --outSAMtype BAM Unsorted
    --outSAMattributes NH HI AS nM RG GX GN
    --outFileNamePrefix %(outprefix)s
    && gzip 7_starsolo/Solo.out/Gene*/*/*
    && gzip 7_starsolo/Solo.out/SJ/*/[bm]*
    '''
    job_options = " -t 72:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_total_memory = '50G')
    if P.get_params()["zap_files"]==1:
        for x in range (0,len(input_files)):
            IOTools.zap_file(input_files[x])
        
##     --readFilesIn %(input_files_string)s
##     2> 7_starsolo/starsolo_log.txt'''
##    --soloCellFilter None

@transform(STARsolo, suffix('.out.bam'), '.sortedByCoord.out.bam')     
def bam_sort_starsolo(input_file, output_file):
    statement = 'samtools sort -@ %(star_threads)s %(input_file)s  -o %(output_file)s'
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_memory = '10G')
    if P.get_params()["zap_files"]==1:
        IOTools.zap_file(input_file)

@transform(bam_sort_starsolo, suffix('.bam'), '.bam.bai')     
def bam_index_starsolo(input_file, output_file):
    statement = 'samtools index %(input_file)s -@ 1'
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=1, job_memory = '100M')


##### Velocyto #####

@follows(mkdir('7_starsolo/bam_split_by_cell'))
@split(bam_sort_starsolo, '7_starsolo/bam_split_by_cell/*.bam')     
def bam_split(input_file, output_files):
    statement = '''samtools split -@ %(star_threads)s %(input_file)s  
    -f "7_starsolo/bam_split_by_cell/%%!.bam"
    '''
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=PARAMS['star_threads'], job_memory = '10G')

@transform(bam_split, suffix('.bam'), '.bam.bai')     
def bam_split_index(input_file, output_file):
    statement = 'samtools index %(input_file)s -@ 1'
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_threads=1, job_memory = '100M')




## Multiqc report ##

@follows(count_reads, samtools_idxstat)
#count_reads only there to specify it should be done at end
@merge(samtools_flagstat, '6_multiqc/multiqc_report.html')
def multiqc(input_file, output_file):
    statement = '''multiqc . -o 6_multiqc'''
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G')


@follows(STARsolo)
#count_reads only there to specify it should be done at end
@merge(STARsolo, '8_multiqc_starsolo/multiqc_report.html')
def multiqc2(input_file, output_file):
    statement = '''multiqc . -o 8_multiqc_starsolo'''
    job_options = " -t 24:00:00"
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G')
    
@follows(mkdir('9_RSeQC'))
@follows(bam_index_starsolo)
@merge(STARsolo, '9_RSeQC/RSeQC.geneBodyCoverage.curves.pdf')
def RSeQC(input_file, output_file):
    statement = '''geneBody_coverage.py 
    -r  %(rseqc_refgene)s
    -i %(input_file)s
    -o 9_RSeQC/RSeQC  
    '''
    P.run(statement, job_queue=PARAMS['q'], job_memory = '8G')    

#     &&
# junction_saturation.py 
# -r %(featurecounts_gtf)s
# -i %(input_file)s      
# -o 9_RSeQC/RSeQC   &&
# read_distribution.py 
# -r %(featurecounts_gtf)s
# -i %(input_file)s      
# > 9_RSeQC/RSeQC.read_distribution.txt


@follows(fastqc, trimming_SE, trimming_PE_merged, trimming_PE_split)
def Trimming():
    pass
    
@follows(star, samtools_flagstat, samtools_idxstat)
def Mapping_qc():
    pass
    
@follows(Mapping_qc, count_reads, multiqc)
def full_featurecounts():
    pass
    
@follows(Trimming, STARsolo)
def starsolo():
    pass

@follows(Trimming, STARsolo, multiqc2)
def full_starsolo():
    pass

if __name__=="__main__":
    sys.exit(P.main(sys.argv))
    
    
    
