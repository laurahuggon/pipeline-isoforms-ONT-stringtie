#!/bin/bash -l
#SBATCH --job-name=ont_barcode04_annotation
#SBATCH --output=/scratch/prj/bcn_synaptopathy/%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=k21224575@kcl.ac.uk

# Check number of files
ls /scratch/prj/bcn_synaptopathy/input_fastq/*.fastq.gz | wc -l
ls /scratch/prj/bcn_synaptopathy/input_fastq/*barcode04*.fastq.gz | wc -l

# Load anaconda
module load anaconda3/2022.10-gcc-13.2.0
# Make job environment behave like login shell
source /users/${USER}/.bashrc
# Activate conda environment
source activate ont_stringtie

# Go to repo directory
cd /scratch/prj/bcn_synaptopathy/pipeline-isoforms-ONT-stringtie
# Run snakemake
snakemake --use-conda -j 30 all