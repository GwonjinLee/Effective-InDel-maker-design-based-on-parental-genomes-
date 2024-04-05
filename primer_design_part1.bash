#!/bin/sh -l

#SBATCH --job-name primer_design 
#SBATCH -A meixiazhao
#SBATCH --mem=36G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --time=3-0:00
#SBATCH -o /blue/meixiazhao/lee.gwonjin/Maize_GBS/scripts/outtext/primer_design_B97.out
#SBATCH -e /blue/meixiazhao/lee.gwonjin/Maize_GBS/scripts/outtext/primer_design_B97.err
#SBATCH --mail-type=all          
#SBATCH --mail-user=lee.gwonjin@ufl.edu

# Load modules
module load bedtools
module load picard
module load primer3
module load R/4.2
module load ncbi_blast
module load seqkit
module load muscle/3.8.31


#working space
cd /blue/meixiazhao/lee.gwonjin/Maize_GBS/Indel/B97

#gzip -d Zm-B97-REFERENCE-NAM-1.0.fa.gz

# Generate variation file from two fasta files
#~/GSAlign/bin/GSAlign -ind 100 \
#-r ../B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa \
#-q Zm-B97-REFERENCE-NAM-1.0.fa \
#-o B97

# Extract Indel region
#python InDel-1.py B97.vcf Indel.txt

# Select Indel region based on the pre-made list using 10_chrs_markers (+- 1MB)
Rscript Extract_indel_range_chr5910.R

# Extract sequences from selected Indel region
bedtools getfasta -fi ../B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa -bed marker_pos.bed -name -fo All_marker.fasta




###########################################################################################################################
# Design primers 

#Use the primer3plus website for a manual primer design istead (https://www.bioinformatics.nl/cgi-bin/primer3plus/primer3plus.cgi)

###########################################################################################################################