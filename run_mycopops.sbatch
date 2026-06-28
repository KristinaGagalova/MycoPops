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

REF=/scratch/pawsey1142/kgagalova/pop_analysis/reference/GCA_900519115.1_Bgt_genome_v3.16_genomic.fna.gz

nextflow run ./main.nf \
    --input samplesheet.csv \
    --fasta $REF \
    --genome_size 140000000\
    --outdir ./results \
    -profile pawsey_setonix,singularity \
    -resume
