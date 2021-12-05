#!/bin/bash
#SBATCH --job-name ChIA-PIPE_nxf
#SBATCH --mail-type=END
#SBATCH --mail-user=${USER}@jax.org   # Your email address
#SBATCH -p compute -q batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 10G
#SBATCH -t 10:00:00


module load singularity

# Specify a location for the installed pipeline
export projectDir=/path/to/Installed_ChIA-PIPE_nextflow_pipeline

nextflow run ${projectDir}/main.nf -profile sumner -c params.config -ansi-log true -resume


exit 0
