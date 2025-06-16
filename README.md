# RNA-Seq Analysis Pipeline

This repository contains a shell-based RNA-Seq processing pipeline for bulk gene expression analysis starting from raw FASTQ files. The pipeline includes quality control, adapter trimming, alignment, post-alignment QC, and read counting.

Currently, the pipeline supports **mouse (Mus musculus)** samples aligned to GRCm39.

---

## Features

- Quality control using **FastQC**
- Adapter trimming using **fastp**
- Alignment using **STAR**
- Post-alignment QC using **samtools** and **Picard**
- Gene-level read quantification using **HTSeq-count**
- SLURM-based job configuration for HPC environments

---

## Pipeline Usage

```bash
./RNAseq_pipeline.sh [pair|single] fastq_accession.txt /path/to/fastq_files mice /path/to/output
```
---

## Arguments

- pair or single
  Specify whether the input data are paired-end or single-end FASTQ files.
- fastq_accession.txt
  A text file listing sample identifiers (one per line). For paired-end data, files should follow the naming convention:
  {SAMPLE}_1.fastq
  {SAMPLE}_2.fastq
- /path/to/fastq_files
  Directory containing raw FASTQ files.
-  mice
  Species indicator. Currently, only mice is supported.
-  /path/to/output
  Directory where all output files will be stored.

---

## Example Command

```bash
./RNAseq_pipeline.sh pair fastq_accession.txt /scratch/data/fastq_files mice /scratch/output
```

---

## Reference Genome (Mouse)

- Genome: Mus musculus GRCm39
- Reference files (STAR index, GTF, genome fasta) must be pre-built and configured in the script.

---

## Output

For each sample, the pipeline generates:
- Quality control reports (FastQC, fastp)
- Aligned BAM files (sorted)
- Alignment metrics (flagstat, picard)
- Gene count tables (HTSeq)




