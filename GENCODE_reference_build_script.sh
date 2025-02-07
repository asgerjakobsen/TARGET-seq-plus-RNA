#!/bin/bash
#SBATCH --ntasks=10
#SBATCH --mem=40G

# Build notes for GENCODE reference
# June 2021

# Genome metadata
genome="GRCh38"

# Folders
mkdir WholeGenomeFasta
mkdir Annotation
sequence="./WholeGenomeFasta"
annotation="./Annotation"

# ERCC
ercc="ERCC92/ERCC92.fa"
ercc_gtf="ERCC92/ERCC92.gtf"



# Download GENCODE fasta
fasta_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz"
fasta_zip="${sequence}/GRCh38.primary_assembly.genome.fa.gz"
fasta_in="${sequence}/GRCh38.primary_assembly.genome.fa"
fasta_ercc="${sequence}/GRCh38.primary_assembly_with_ERCC92.genome.fa"

if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" > "$fasta_zip"
    zcat "$fasta_zip" > "$fasta_in"
fi


# Append the ERCC annotations to the fasta
cat "$fasta_in" > "$fasta_ercc"
cat "$ercc" >> "$fasta_ercc"


# Download GENCODE GTF
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz"
gtf_zip="${annotation}/gencode.v38.primary_assembly.annotation.gtf.gz"
gtf_in="${annotation}/gencode.v38.primary_assembly.annotation.gtf"


if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" > "$gtf_zip"
    zcat "$gtf_zip" > "$gtf_in"
fi


# Filter the GTF file as for 10X

# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.

BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"

GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag

# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_in" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${annotation}/gene_allowlist"

# Filter the GTF file based on the gene allowlist
gtf_filtered="${annotation}/gencode.v38.primary_assembly.annotation.filtered.gtf"

# Copy header lines beginning with "#"
grep -E "^#" "$gtf_in" > "$gtf_filtered"

# Filter to the gene allowlist
grep -Ff "${annotation}/gene_allowlist" "$gtf_in" \
    >> "$gtf_filtered"

# Append the ERCC annotations

# Append the ERCC annotations to the gtf
gtf_filtered_ercc="${annotation}/gencode.v38.primary_assembly.annotation.filtered.with_ERCC92.gtf"

cat "$gtf_filtered" > "$gtf_filtered_ercc"
cat "$ercc_gtf" >> "$gtf_filtered_ercc"

# Build genome
STAR   --runMode genomeGenerate   --runThreadN 10 --genomeDir STAR_2.7.9a_GENCODE_with_ERCC_filtered_GTF/ \
--genomeFastaFiles WholeGenomeFasta/GRCh38.primary_assembly_with_ERCC92.genome.fa     --sjdbGTFfile Annotation/gencode.v38.primary_assembly.annotation.filtered.with_ERCC92.gtf   --sjdbOverhang 149

