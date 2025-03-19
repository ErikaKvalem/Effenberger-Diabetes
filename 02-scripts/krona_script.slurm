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
#ktImportTaxonomy -t 5 -m 3 /data/projects/2024/Effenberger-Diabetes/out/kraken2/no_PDM3/PDM/kraken2_output_PDM.txt -o krona_report_PDM.html
#ktImportTaxonomy -t 5 -m 3 /data/projects/2024/Effenberger-Diabetes/out/kraken2/no_PDM3/DM/kraken2_output_DM.txt -o krona_report_DM.html
#ktImportTaxonomy -t 5 -m 3 /data/projects/2024/Effenberger-Diabetes/out/kraken2/no_PDM3/nikotin/kraken2_output_nikotin.txt -o krona_report_nikotin.html
ktImportTaxonomy -t 5 -m 3 /data/projects/2024/Effenberger-Diabetes/out/kraken2/no_PDM3/no_nikotin/kraken2_output_no_nikotin.txt -o krona_report_no_nikotin.html
#ktImportTaxonomy -t 5 -m 3 /data/projects/2024/Effenberger-Diabetes/out/kraken2/no_PDM3/male/kraken2_output_male.txt -o krona_report_male.html
#ktImportTaxonomy -t 5 -m 3 /data/projects/2024/Effenberger-Diabetes/out/kraken2/no_PDM3/female/kraken2_output_female.txt -o krona_report_female.html
ktImportTaxonomy -t 5 -m 3 /data/projects/2024/Effenberger-Diabetes/out/kraken2/no_PDM3/control/kraken2_output_control.txt -o krona_report_control.html
