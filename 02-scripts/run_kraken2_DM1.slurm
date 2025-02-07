#!/bin/bash
#SBATCH --job-name=kraken2_analysis         # Job name
#SBATCH --output=kraken2_%j.out            # Standard output log (%j = job ID)
#SBATCH --error=kraken2_%j.err             # Standard error log
#SBATCH --ntasks=1                         # Run on a single task
#SBATCH --cpus-per-task=2                  # Number of CPU cores
#SBATCH --mem=16G                          # Memory allocation
#SBATCH --partition=standard               # Partition to submit to (adjust as needed)

# Load any necessary modules (if kraken2 is available as a module, load it)
#module load kraken2                         # Adjust based on your cluster's configuration
eval "$(conda shell.bash hook)"
conda activate kraken2

# Run Kraken2
kraken2 --db /data/databases/kraken2/16S_Silva138_20200326/16S_SILVA138_k2db \
        --threads 2 \
        --output /data/projects/2024/Effenberger-Diabetes/out/kraken2/kraken2_output_DM1.txt \
        /data/projects/2024/Effenberger-Diabetes/data/20011/Fastq/DM/DM1/*.fastq.gz
