#!/bin/sh -l

#SBATCH --job-name primer_design 
#SBATCH -A meixiazhao
#SBATCH --mem=36G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --time=3-0:00
#SBATCH -o /blue/meixiazhao/lee.gwonjin/Maize_GBS/scripts/outtext/primer_design.out
#SBATCH -e /blue/meixiazhao/lee.gwonjin/Maize_GBS/scripts/outtext/primer_design.err
#SBATCH --mail-type=end          
#SBATCH --mail-user=lee.gwonjin@ufl.edu

# Load modules
module load bedtools
module load picard
module load gatk
#module load samtools
module load primer3
module load R/4.2
module load ncbi_blast
module load seqkit
module load muscle/3.8.31



#working space
cd /blue/meixiazhao/lee.gwonjin/Maize_GBS/Indel/A344

#gzip -d A344/Zm-A344-REFERENCE-NRGENE-2.0.fa.gz

#################################################################################################
# In case there is no alternative genome but a vcf file create an alternative genome fa file

#Make an alternative reference fasta file
## Create dictionary file for fa
#picard CreateSequenceDictionary \ 
#R=/blue/meixiazhao/lee.gwonjin/Maize_GBS/Indel/B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa \
#O=/blue/meixiazhao/lee.gwonjin/Maize_GBS/Indel/B73_ref/Zm-B73-REFERENCE-NAM-5.0.dict

## Index vcf
#gatk IndexFeatureFile -I A344_parent.vcf

## Make alt fa
#gatk FastaAlternateReferenceMaker \
#-R /blue/meixiazhao/lee.gwonjin/Maize_GBS/Indel/B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa \
#-V A344_parent.vcf -O ref/Zm-A344-REFERENCE-NRGENE-2.0.fa

## Rename headers of the fa file for the chrom

#################################################################################################

# Generate variation file from two fasta files
#~/GSAlign/bin/GSAlign -ind 100 \
#-r ../B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa \
#-q A344/Zm-A344-REFERENCE-NRGENE-2.0.fa \
#-o A344.vcf

# Extract Indel region
#python InDel-1.py A344_indel.recode.vcf Indel.txt

# Select Indel region based on the pre-made list using 10_chrs_markers (+- 1MB)
#Rscript Extract_indel_range.R

# Extract sequences from selected Indel region
#bedtools getfasta -fi ../B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa -bed marker_pos.bed -name -fo All_marker.fasta




###########################################################################################################################
# Design primers (not working)
#python fasta2primer3.py All_marker.fasta

#Use the primer3plus website for a manual primer design istead (https://www.bioinformatics.nl/cgi-bin/primer3plus/primer3plus.cgi)
#Find the Primer3plus_Settings file and use for primer design setting

###########################################################################################################################


#Remove the header for blast
#sed 's/>[^~]*:://g' Primers_A344_name.fasta > Primers_A344.fasta
#sed -i 's/chr/>chr/g' Primers_A344.fasta
#sed -i 's/     |X|11.01.2023//g' Primers_A344.fasta  #change the date for today

#If you want to add some more primers, then add them to Primers_A344.fasta

# Check specificity of primers
##For the reference genome
#makeblastdb -in ../B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa -parse_seqids -dbtype nucl 
#blastn -word_size 7 -evalue 1 -query Primers_A344.fasta -db ../B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa -dust no -num_threads 8 -outfmt 6 -out blast_A344_ref_primers.txt

## Extract only uniquely blasted primers (1 hit)
#awk 'NR == FNR {count[$1]++; next} count[$1]==1' blast_A344_ref_primers.txt blast_A344_ref_primers.txt | awk '$3 == "100.000" {print $1}' > blast_A344_primers_ref_filtered.txt

#makeblastdb -in ref/Zm-A344-REFERENCE-NRGENE-2.0.fa -parse_seqids -dbtype nucl 
#blastn -word_size 7 -evalue 1 -query Primers_A344.fasta -db ref/Zm-A344-REFERENCE-NRGENE-2.0.fa -dust no -num_threads 8 -outfmt 6 -out blast_A344_alt_primers.txt
#awk 'NR == FNR {count[$1]++; next} count[$1]==1' blast_A344_alt_primers.txt blast_A344_alt_primers.txt | awk '$3 == "100.000" {print $1}'> blast_A344_primers_alt_filtered.txt

# Retain uniquely blasted primers from both blast results using ref and alt genomes
#grep -xFf blast_A344_primers_ref_filtered.txt blast_A344_primers_alt_filtered.txt > blast_A344_primers_merged.txt #add any addtional primers that either of F or R is ok here!

# Extract only paired primers and make the bed file
#Rscript Primer_processing.R # Check the sample name in the script

# Re-add sequences to primer names
#seqkit grep -n -f blast_A344_primers_paired.txt Primers_A344.fasta > blast_A344_primers_paired.fasta


# Extract sequences from alternative genome (normalize fa file using picard if necessary)
#picard NormalizeFasta I=ref/Zm-A344-REFERENCE-NRGENE-2.0.fa O=ref/norm_Zm-A344-REFERENCE-NRGENE-2.0.fa
## replace the header for the selected chromosome (the alt fa file has different headers)
#sed -i 's/>4 chr4:1-250330460/>chr4/g' ref/norm_Zm-A344-REFERENCE-NRGENE-2.0.fa #check also the fai file

#bedtools getfasta -fi ref/norm_Zm-A344-REFERENCE-NRGENE-2.0.fa -bed A344_primers_alt.bed -s -name -fo A344_primers_altseq.fasta
#sed -i 's/:://g' A344_primers_altseq.fasta # remove ::
#sed -i 's/_Fch/_F ch/g' A344_primers_altseq.fasta # make space
#sed -i 's/_Rch/_R ch/g' A344_primers_altseq.fasta # make space

# Merge primers and alt sequences (blast_A344_primers_paired.fasta + A344_primers_altseq.fasta)
## Sort the input fasta files alphabetically
#module load R/3.2.2 #need to load the old version R for bioawk
#module load gcc/5.2.0
#module load bioawk
#bioawk -c fastx '{print}' blast_A344_primers_paired.fasta | sort -k1,1V | awk '{print ">"$1;print $2}' > blast_A344_primers_paired_alpha.fasta
#bioawk -c fastx '{print}' A344_primers_altseq.fasta | sort -k1,1V | awk '{print ">"$1;print $2}' > A344_primers_altseq_alpha.fasta

#combine two fasta file for input
#paste blast_A344_primers_paired_alpha.fasta A344_primers_altseq_alpha.fasta > primers_cand.fa 

#module load R/4.2 #Re-load the new version
#Rscript Primer_processing2.R

# Get degenerate primers
#for infile in *_primer_input.fa
#	do base=$(basename ${infile} _primer_input.fa) 
#	muscle -in ${infile} -out ${base}.degen_primer.txt -clw 
#	done

#merge all muscle output files
#ls -v *degen_primer.txt | xargs cat > degen_primers.txt

#Remove intermediate files
rm *_primer_input.fa
rm *degen_primer.txt
rm primers_cand.fa
rm blast_A344_primers_paired.fasta
rm A344_primers_altseq_alpha.fasta
rm blast_A344_primers_paired.txt
rm A344_primers_altseq.fasta
rm A344_primers_alt.bed

# Select final primers manually

# Generate the primer list containing the primer name, position, degenrate sequence, and etc.
## Use "degen_primers.txt" for primer sequences and "All_marker.fasta" for PCR products info
## check the PCR product size in Primer3plus