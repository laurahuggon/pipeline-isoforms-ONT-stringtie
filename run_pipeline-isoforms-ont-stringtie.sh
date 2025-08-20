#!/bin/bash -l
#SBATCH --job-name=ont_barcode01_no-annotation
#SBATCH --output=/scratch/users/%u/%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=128G
#SBATCH --time=48:00:00

# Set sample barcode
BARCODE="barcode01"

# Copy input (FASTQ) files for this sample to scratch
cp -n /rds/prj/bcn_synaptopathy_in_ftd_als/Emma_Clayton_data/pool1-16/*/fastq_pass/${BARCODE}/*.fastq.gz /scratch/users/$USER/input_fastq/
# Print number of FASTQ files in input_fastq directory
ls /scratch/users/$USER/input_fastq/*.fastq.gz | wc -l

# Load anaconda
module load anaconda3/2022.10-gcc-13.2.0
# Make job environment behave like login shell
source /users/${USER}/.bashrc
# Activate conda environment
source activate ont_stringtie

# Go to repo directory
cd /users/$USER/pipeline-isoforms-ONT-stringtie
# Test run (prints the commands it would run but does not execute anything)
snakemake --use-conda -n all
# Run snakemake
snakemake --use-conda -j 30 all

# Copy output files to RDS
cp -r /scratch/users/$USER/${BARCODE} /rds/prj/bcn_synaptopathy_in_ftd_als/synaptosome_long-read/data/processed/
cp /scratch/users/%u/%j.out /rds/prj/bcn_synaptopathy_in_ftd_als/synaptosome_long-read/data/processed/${BARCODE}/

# Remove input (FASTQ) files from scratch

# Remove output files from scratch