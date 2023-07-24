import os
fh=open("W22_test.vcf","r")
fh2=open("Extracted_sequence.txt","w")
for line in fh:
	data=line.strip().split("\t")
	start=(int(data[1]))-500
	end=(int(data[1]))+500
	chr=data[0]
	os.system("samtools faidx ref/Zm-B73-REFERENCE-NAM-5.0.fa {chr}:{start}-{end} >> Extracted_sequence.txt")
