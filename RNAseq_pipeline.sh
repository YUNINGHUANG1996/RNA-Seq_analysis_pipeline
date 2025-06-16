#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00

# Usage: ./RNAseq_pipeline.sh pair fastq_accession.txt /path/to/fastq_files mice /path/to/output

# Parse input arguments
TYPE=$1            # pair or single
ACCESSION_FILE=$2  # txt file containing list of sample names
FASTQ_DIR=$3       # path to raw FASTQ files
SPECIES=$4         # species (currently only 'mice' supported)
OUTPUT_DIR=$5      # output path

# Paths for mice
if [ "$SPECIES" == "mice" ]; then
    REF_DIR="/scratch1/yuninghu/PreventE4_transcriptomics/ref/Mice"
    STAR_INDEX="${REF_DIR}/STAR_index"
    GTF="${REF_DIR}/Mus_musculus.GRCm39.110.gtf"
    GENOME_FA="${REF_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
else
    echo "Species $SPECIES not supported yet."
    exit 1
fi

# Save command to README (for reproducibility)
mkdir -p ${OUTPUT_DIR}
echo "Command executed on $(date):" >> ${OUTPUT_DIR}/README.txt
echo "$0 $@" >> ${OUTPUT_DIR}/README.txt
echo " " >> ${OUTPUT_DIR}/README.txt

# Load modules
module load gcc/13.3.0 fastqc/0.12.1 samtools/1.19.2 picard/3.1.1

# Process each sample
while read SAMPLE; do
    echo "Processing sample: $SAMPLE"
    mkdir -p ${OUTPUT_DIR}/${SAMPLE}
    cd ${OUTPUT_DIR}/${SAMPLE}

    if [ "$TYPE" == "pair" ]; then
        FQ1="${FASTQ_DIR}/${SAMPLE}_1.fastq"
        FQ2="${FASTQ_DIR}/${SAMPLE}_2.fastq"
    else
        FQ1="${FASTQ_DIR}/${SAMPLE}.fastq"
    fi

    #### 1. Run FastQC
    if [ "$TYPE" == "pair" ]; then
        fastqc $FQ1 $FQ2 -o .
    else
        fastqc $FQ1 -o .
    fi

    #### 2. Run fastp (Trimming)
    /scratch1/yuninghu/fastp/fastp \
        -i $FQ1 \
        ${FQ2:+-I $FQ2} \
        -o ${SAMPLE}_R1.trimmed.fastq.gz \
        ${FQ2:+-O ${SAMPLE}_R2.trimmed.fastq.gz} \
        --detect_adapter_for_pe --thread 8 --html fastp_report.html

    module load gcc/13.3.0 star/2.7.11b
    #### 3. Run STAR alignment
    if [ "$TYPE" == "pair" ]; then
        STAR --runThreadN 8 \
            --genomeDir $STAR_INDEX \
            --readFilesIn ${SAMPLE}_R1.trimmed.fastq.gz ${SAMPLE}_R2.trimmed.fastq.gz \
            --readFilesCommand zcat \
            --outFileNamePrefix ${SAMPLE}_ \
            --outSAMtype BAM SortedByCoordinate
    else
        STAR --runThreadN 8 \
            --genomeDir $STAR_INDEX \
            --readFilesIn ${SAMPLE}_R1.trimmed.fastq.gz \
            --readFilesCommand zcat \
            --outFileNamePrefix ${SAMPLE}_ \
            --outSAMtype BAM SortedByCoordinate
    fi

    #### 4. Post-alignment QC
    module load gcc/13.3.0 fastqc/0.12.1 samtools/1.19.2 picard/3.1.1
    samtools flagstat ${SAMPLE}_Aligned.sortedByCoord.out.bam > alignment_flagstat.txt
    picard CollectAlignmentSummaryMetrics R=$GENOME_FA I=${SAMPLE}_Aligned.sortedByCoord.out.bam O=alignment_metrics.txt

    #### 5. Run HTSeq
    module load legacy/CentOS7 gcc/11.3.0 openblas/0.3.20 py-htseq/0.11.2

    htseq-count \
      -f bam \
      -r pos \
      -s reverse \
      -i gene_id \
      ${SAMPLE}_Aligned.sortedByCoord.out.bam \
      $GTF > ${SAMPLE}.counts.txt

done < $ACCESSION_FILE

echo "All samples processed successfully."

