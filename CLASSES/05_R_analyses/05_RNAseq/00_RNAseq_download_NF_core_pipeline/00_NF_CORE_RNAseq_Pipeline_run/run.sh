#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=2023_HEPG2_rna_seq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=luva6824@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=20:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err

pwd; hostname; date
echo "Here we go You've requested $SLURM_CPUS_ON_NODE core."

module load singularity/3.1.1

nextflow run nf-core/rnaseq -r 1.4.2 \
-resume \
-profile singularity \
--reads '/scratch/Shares/rinnclass/CLASS_2023/lucyvk/CLASS_2023/CLASSES/05_R_analyses/05_RNAseq/00_RNAseq_download_NF_core_pipeline/data/fastq2/*{_read1,_read2}.fastq.gz' \
--fasta /scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/GRCh38.p13.genome.fa \
--gtf /scratch/Shares/rinnclass/CLASS_2023/data/data/genomes/gencode.v32.annotation.gtf \
--pseudo_aligner salmon \
--gencode \
--email luva6824@colorado.edu \
-c nextflow.config

date
