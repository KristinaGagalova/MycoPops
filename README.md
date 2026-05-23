# Population Admixture analysis of fungal WGS data

A Nextflow pipeline for population genetics analysis of whole-genome resequencing data.

## Pipeline Steps

| Step | Tool | Description |
|------|------|-------------|
|  |  |  |

## Quick Start

```bash
nextflow run MycoPop \
    --input samplesheet.csv \
    --fasta /path/to/reference.fasta \
    --genome_size 139956545 \
    --outdir ./results \
    -profile singularity,slurm
```

## Input: Sample Sheet

A CSV file with a header row and three columns:

| Column | Description |
|--------|-------------|
| `sample` | Unique sample identifier |
| `fastq_1` | Full path to R1 FASTQ (gzipped) |
| `fastq_2` | Full path to R2 FASTQ (gzipped) |

Example:

```csv
sample,fastq_1,fastq_2
395718,/data/fastp/395718_R1.fastp.fq.gz,/data/fastp/395718_R2.fastp.fq.gz
395719,/data/fastp/395719_R1.fastp.fq.gz,/data/fastp/395719_R2.fastp.fq.gz
```

> **TSV support:** Set `--manifest_sep '\t'` if your manifest uses tab separators (matching the original `sample_list_bwa` format).

## Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--input` | ✅ | - | Path to samplesheet CSV/TSV |
| `--fasta` | ✅ | - | Reference genome FASTA (must be pre-indexed with `bwa-mem2 index`) |
| `--genome_size` | ✅ | - | Reference genome size in bp (for coverage calculation) |
| `--outdir` | - | `./results` | Output directory |
| `--manifest_sep` | - | `,` | Manifest separator (`,` for CSV, `\t` for TSV) |

## Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers |
| `singularity` | Run with Singularity containers |
| `conda` | Run with Conda environments |
| `slurm` | Submit jobs to Slurm scheduler (e.g. Pawsey/Setonix) |
| `test` | Small test dataset with minimal resources |

## Output Structure

```
results/
├── bwamem2/
│   └── <sample>.bwa_alnPE.sorted.bam
├── picard/
│   ├── <sample>.markdup.bam
│   └── <sample>.marked_dup_metrics.txt
├── coverage/
│   ├── <sample>.coverage.tsv
│   └── coverage_summary.tsv
└── pipeline_info/
    ├── timeline_*.html
    ├── report_*.html
    ├── trace_*.txt
    └── dag_*.html
```

## Software Versions

- bwa-mem2 2.3
- samtools 1.23.1
- picard 3.4.0

## Credits

Pipeline adapted from bash scripts by Grace Fang (Curtin University / Pawsey).
Original scripts: https://github.com/fc87290118/WPM_population_analysis

## License

MIT
