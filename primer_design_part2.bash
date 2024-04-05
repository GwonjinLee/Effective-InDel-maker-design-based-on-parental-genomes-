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


###########################################################################################################################
# Design primers (not working)
#python fasta2primer3.py All_marker.fasta

#Use the primer3plus website for a manual primer design istead (https://www.bioinformatics.nl/cgi-bin/primer3plus/primer3plus.cgi)

###########################################################################################################################


#Remove the header for blast
sed 's/>[^~]*:://g' Primers_B97_name.fasta > Primers_B97.fasta
sed -i 's/chr/>chr/g' Primers_B97.fasta
sed -i 's/     |X|01.04.2024//g' Primers_B97.fasta  #change the date for today!!

#If you want to add some more primers, then add them to Primers_B97.fasta

# Check specificity of primers
##For the B73 reference genome 
###makeblastdb -in ../B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa -parse_seqids -dbtype nucl #(skip if you already have)
blastn -word_size 7 -evalue 1 -query Primers_B97.fasta -db ../B73_ref/Zm-B73-REFERENCE-NAM-5.0.fa -dust no -num_threads 8 -outfmt 6 -out blast_B97_REF_primers.txt

## Extract only uniquely blasted primers (1 hit)
awk 'NR == FNR {count[$1]++; next} count[$1]==1' blast_B97_REF_primers.txt blast_B97_REF_primers.txt | awk '$3 == "100.000" {print $1}' > blast_B97_primers_REF_filtered.txt

##For the alternative genome
#makeblastdb -in Zm-B97-REFERENCE-NAM-1.0.fa -parse_seqids -dbtype nucl 
blastn -word_size 7 -evalue 1 -query Primers_B97.fasta -db Zm-B97-REFERENCE-NAM-1.0.fa -dust no -num_threads 8 -outfmt 6 -out blast_B97_alt_primers.txt
awk 'NR == FNR {count[$1]++; next} count[$1]==1' blast_B97_alt_primers.txt blast_B97_alt_primers.txt | awk '$3 == "100.000" {print $1}'> blast_B97_primers_alt_filtered.txt

# Retain uniquely blasted primers from both blast results using REF and alt genomes
grep -xFf blast_B97_primers_REF_filtered.txt blast_B97_primers_alt_filtered.txt > blast_B97_primers_merged.txt #add any addtional primers that either of F or R is ok here!

# Extract only paired primers and make the bed file
Rscript Primer_processing.R

# Re-add sequences to primer names
seqkit grep -n -f blast_B97_primers_paired.txt Primers_B97.fasta > blast_B97_primers_paired.fasta

# Extract sequences from alternative genome (normalize fa file using picard if necessary)
#picard NormalizeFasta I=Zm-B97-REFERENCE-NAM-1.0.fa O=norm_Zm-B97-REFERENCE-NAM-1.0.fa

########## Only if necessary
########## replace the header for the selected chromosome (the alt fa file has different headers)
########## sed -i 's/>4 chr4:1-250330460/>chr4/g' ref/norm_Zm-A344-REFERENCE-NRGENE-2.0.fa #check also the fai file

# Get fasta file for alt primers
bedtools getfasta -fi norm_Zm-B97-REFERENCE-NAM-1.0.fa -bed B97_primers_alt.bed -s -name -fo B97_primers_altseq.fasta
sed -i 's/:://g' B97_primers_altseq.fasta # remove ::
sed -i 's/_Fch/_F ch/g' B97_primers_altseq.fasta # make space
sed -i 's/_Rch/_R ch/g' B97_primers_altseq.fasta # make space

# Merge primers and alt sequences (blast_B97_primers_paired.fasta + B97_primers_altseq.fasta)
## Sort the input fasta files alphabetically
module load R/3.2.2 #need to load the old version R for bioawk
module load gcc/5.2.0
module load bioawk
bioawk -c fastx '{print}' blast_B97_primers_paired.fasta | sort -k1,1V | awk '{print ">"$1;print $2}' > blast_B97_primers_paired_alpha.fasta
bioawk -c fastx '{print}' B97_primers_altseq.fasta | sort -k1,1V | awk '{print ">"$1;print $2}' > B97_primers_altseq_alpha.fasta

#combine two fasta file for input
paste blast_B97_primers_paired_alpha.fasta B97_primers_altseq_alpha.fasta > primers_cand.fa 

module load R/4.2 #Re-load the new version
Rscript Primer_processing2.R

# Get degenerate primers
for infile in *_primer_input.fa
do	base=$(basename ${infile} _primer_input.fa) 
	muscle -in $infile -clw -out $base.degen_primer.txt
	done

#merge all muscle output files
ls -v *degen_primer.txt | xargs cat > degen_primers.txt

#Remove intermediate files
rm *_primer_input.fa
rm *degen_primer.txt
rm primers_cand.fa
rm blast_B97_primers_paired.fasta
rm B97_primers_altseq_alpha.fasta
rm blast_B97_primers_paired.txt
rm B97_primers_altseq.fasta
rm B97_primers_alt.bed


# Select final primers manually

# Generate the primer list containing the primer name, position, degenrate sequence, and etc.
## Use "degen_primers.txt" for primer sequences and "All_marker.fasta" for PCR products