#!/usr/bin/bash --login

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB #requirements may be high in case user downloads large images
#SBATCH --partition=work
#SBATCH --time=1-00:00:00
#SBATCH --account=pawsey1142
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

module load nextflow/25.04.6
module load singularity/4.1.0-slurm

REF="/scratch/pawsey1142/gracefang/WPM_project/WPM_genome/Bgt_USA_2_genome_v1.fasta"

nextflow run ./main.nf \
    --input "/scratch/pawsey1142/gracefang/WPM_project/MycoPops/WPM_38_Australian_fastq_manifest.csv" \
    --fasta "$REF" \
    --genome_size 140000000\
    --outdir ./results \
    --min_genotype_rate 0.8 \
    -profile pawsey_setonix,singularity \
    -resume 
