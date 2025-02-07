#!/bin/bash
#SBATCH --job-name=krrona_report         # Job name
#SBATCH --output=krona_%j.out            # Standard output log (%j = job ID)
#SBATCH --error=krona_%j.err             # Standard error log
#SBATCH --ntasks=1                         # Run on a single task
#SBATCH --cpus-per-task=2                  # Number of CPU cores
#SBATCH --mem=16G                          # Memory allocation
#SBATCH --partition=standard               # Partition to submit to (adjust as needed)

# Load any necessary modules (if kraken2 is available as a module, load it)
#module load kraken2                         # Adjust based on your cluster's configuration
eval "$(conda shell.bash hook)"
conda activate krona_2025
# Run Kraken2
ktImportTaxonomy -t 5 -m 3 kraken2_output_PDM.txt -o krona_report_PDM.html
ktImportTaxonomy -t 5 -m 3 kraken2_output_DM.txt -o krona_report_DM.html
