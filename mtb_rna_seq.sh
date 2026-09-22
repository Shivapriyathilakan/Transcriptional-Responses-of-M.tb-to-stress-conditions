#!/bin/bash
set -e
set -o pipefail 

SECONDS=0

cd raw_fastq || exit 1

## Run fastqc on all raw fastq files
fastqc -o "fastqc" *.fastq 

## Run MultiQC on the FastQC output directory
multiqc "fastqc" -o "fastqc"

## trim the raw fastq files
## Define paths and variables
RAW_DIR="raw_fastq"
TRIMMED_DIR="trimmed_fastq"
ADAPTERS="TruSeq3-PE.fa"
THREADS=4 # Defined thread count

## Create output directory
mkdir -p "$TRIMMED_DIR"

GROUP_A_SAMPLES=(
    "SRR13687749" "SRR13687750" "SRR13687751" "SRR13687767" "SRR13687768" "SRR13687769"
    "SRR3725585"  "SRR3725586"  "SRR3725587"  "SRR3725597"  "SRR3725598"  "SRR3725599")

echo "Processing Group A (Longer Reads)..."

for SAMPLE in "${GROUP_A_SAMPLES[@]}"; do
    echo "Trimming $SAMPLE..."
    
    trimmomatic PE -threads $THREADS \
      "${RAW_DIR}/${SAMPLE}_1.fastq" "${RAW_DIR}/${SAMPLE}_2.fastq" \
      "${TRIMMED_DIR}/${SAMPLE}_1_paired.fq.gz" "${TRIMMED_DIR}/${SAMPLE}_1_unpaired.fq.gz" \
      "${TRIMMED_DIR}/${SAMPLE}_2_paired.fq.gz" "${TRIMMED_DIR}/${SAMPLE}_2_unpaired.fq.gz" \
      ILLUMINACLIP:${ADAPTERS}:2:30:10 \
      LEADING:3 TRAILING:3 \
      SLIDINGWINDOW:4:15 \
      MINLEN:36
 done

GROUP_B_SAMPLES=(
    "SRR1917709" "SRR1917710" "SRR1917711" "SRR1917703" "SRR1917704" "SRR1917705" )

echo "Processing Group B (Short Reads - Paired)..."

for SAMPLE in "${GROUP_B_SAMPLES[@]}"; do
    echo "Trimming $SAMPLE..."
    
    trimmomatic PE -threads $THREADS \
      "${RAW_DIR}/${SAMPLE}_1.fastq" "${RAW_DIR}/${SAMPLE}_2.fastq" \
      "${TRIMMED_DIR}/${SAMPLE}_1_paired.fq.gz" "${TRIMMED_DIR}/${SAMPLE}_1_unpaired.fq.gz" \
      "${TRIMMED_DIR}/${SAMPLE}_2_paired.fq.gz" "${TRIMMED_DIR}/${SAMPLE}_2_unpaired.fq.gz" \
      ILLUMINACLIP:${ADAPTERS}:2:30:10 \
      LEADING:3 TRAILING:3 \
      SLIDINGWINDOW:4:12 \
      MINLEN:18
 done

# ---------- QC after trimming ----------
fastqc -o trimmed_fastqc trimmed_fastq/*_paired.fq.gz
multiqc trimmed_fastqc -o trimmed_fastqc

## mapping with Bowtie2
 mkdir -p results

## Allocate CPU cores (20 total)
BOWTIE_THREADS=16
SAMTOOLS_THREADS=4

for r1 in trimmed_fastq/*_1_paired.fq.gz; do
    Skip loop if no matching fastq files exist
    [ -e "$r1" ] || continue

    sample=$(basename "$r1" _1_paired.fq.gz)

     CHECK: Skip if the sorted BAM and its index already exist
     if [ -f "results/${sample}_aligned.sorted.bam" ] && [ -f "results/${sample}_aligned.sorted.bam.bai" ]; then
        echo "Skipping ${sample} (BAM and index already exist)."
        continue
    

    echo "Aligning sample: ${sample}..."

    bowtie2 -p "$BOWTIE_THREADS" -x reference/Mtb_index \
        -1 "$r1" \
        -2 "trimmed_fastq/${sample}_2_paired.fq.gz" \
        | samtools sort -@ "$SAMTOOLS_THREADS" -o "results/${sample}_aligned.sorted.bam"

    # Index using multithreading
    samtools index -@ "$SAMTOOLS_THREADS" "results/${sample}_aligned.sorted.bam"

    echo "Finished ${sample}."
done

echo "All alignments complete!"

 # Directories
BAM_DIR="results"
STATS_DIR="results/qc_stats"
THREADS=4  # Adjust threads for samtools execution

mkdir -p "$STATS_DIR"

echo "=========================================="
echo " Running Samtools QC (flagstat & idxstats)"
echo "=========================================="

for bam in "${BAM_DIR}"/*_aligned.sorted.bam; do
     Skip if no BAM files exist
    [ -e "$bam" ] || continue

    sample=$(basename "$bam" _aligned.sorted.bam)
    echo "Processing QC for: ${sample}..."

     1. Flagstat: High-level alignment metrics (mapping %, duplicates, paired status)
    samtools flagstat -@ "$THREADS" "$bam" > "${STATS_DIR}/${sample}_flagstat.txt"

     2. Idxstats: Chromosome/Contig-level mapped & unmapped read counts
    samtools idxstats -@ "$THREADS" "$bam" > "${STATS_DIR}/${sample}_idxstats.txt"
done

echo "QC complete! Stats files saved to: ${STATS_DIR}/"

# Run MultiQC 
multiqc . -o results/multiqc_report

infer_experiment.py \
  -i results/SRR1917703_aligned.sorted.bam \
  -r reference/NC_000962.3.bed
 
 mkdir -p results/counts

featureCounts -p --countReadPairs \
  -s 2 \
  -T 8 \
  -a reference/NC_000962.3.gff \
  -t gene \
  -g locus_tag \
  -o results/counts/mtb_gene_counts.txt \
  results/*.sorted.bam

duration=$SECONDS
echo "$(($duration / 60)) minutes $((duration % 60)) seconds elapsed."
