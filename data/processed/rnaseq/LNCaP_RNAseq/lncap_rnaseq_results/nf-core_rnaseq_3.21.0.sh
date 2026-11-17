#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=80G
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@address.com

# Load necessary modules
module purge
module add miniconda3/v4
module add singularity/3.10.4

# Initialize Conda and activate the nf-core environment
source /home/software/conda/miniconda3/bin/condainit
conda activate nf-core

# Force a specific Nextflow version for reproducibility
export NXF_VER=25.10.2

# Execute Nextflow nf-core/rnaseq pipeline
# Custom configurations are applied via 'nextflow.config'
nextflow run nf-core/rnaseq -r 3.21.0 \
-c nextflow.config \
-profile slurm,singularity \
--input RNA_samplesheet.csv \
--fasta /path/to/your/reference/GRCh38.primary_assembly.genome.fa \
--gtf /path/to/your/reference/gencode.v49.primary_assembly.basic.annotation.gtf \
--with_umi \
--skip_umi_extract \
--umitools_umi_separator "_" \
--extra_star_align_args "--alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --outFilterMismatchNmax 999 --outFilterMultimapNmax 20 --outFilterType BySJout --outFilterMismatchNoverLmax 0.1 --clip3pAdapterSeq AAAAAAAA" \
--multiqc_title 'multiqc_title' \
--aligner 'star_salmon' \
--save_align_intermeds \
--pseudo_aligner 'salmon' \
--skip_stringtie \
--outdir Update 'outdir' \
--save_reference \
-resume
