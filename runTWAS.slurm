#!/bin/sh
#SBATCH --array=1-4
#SBATCH --job-name=TWAS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=0-05:00:00
#SBATCH --mem-per-cpu=16GB
#SBATCH --partition=jclarke,batch
#SBATCH --output=output.log.TWAS/%x_%a.out
#SBATCH --mail-type=ALL

# First we ensure a clean running environment:
module purge

# Load R
ml anaconda/4.12 R/4.3 cmake/3.20

conda activate GWASrMVP

### Initialization
# Get Array ID
i=${SLURM_ARRAY_TASK_ID}

### Parameter file to read 
pheno="pheno.txt"
echo "Parameter file: ${pheno}"
echo

# Read line #i from the parameter file
PHENO=$(sed "${i}q;d" ${pheno})
echo "Parameters are: ${PHENO}"
echo

# Output file
outFile="output_${PHENO}.txt"


# Pass line #i to a R script 
Rscript --verbose TWAS.R ${PHENO} ${outFile}
