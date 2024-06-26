#!/bin/bash
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=5:00:00
 
#SBATCH --output=script.out
#SBATCH --account=416a516a
#SBATCH --mail-type=END
 
### fastp job part
InputDir = "rawReads"
file_list=$(ls $InputDir/*SRR* | sort)

# set the dir to store output files
mkdir -p fastp 
 
# run your executable program with begin and end date and time output
cd fastp
date
 
# only 1 read for SE data
for Read1 in $file_list
do
echo 'currently processing '$Read
modules/fastp -i $InputDir/$Read1 -o $Read1 --detect_adapter_for_pe -5 -r
done
 
date 
##########################
### fastQC job part
# set the dir to store output files
mkdir -p fastQC

cd fastQC
module load fastqc  	### Load the HPC-installed module/program
 
# run your executable program with begin and end date and time output
date
 
for file in $(ls $InputPath*.gz)
do
echo 'currently processing' $file
fastqc $file -o .
done
 
date
 
##########################
### multiQC job part 
cd fastp
 
module load contrib
module load bjoyce3/sarawillis/multiqc/1.20
 
multiqc .
##########################
### STAR mapping job part
# set input file dir
InputDir="fastp"
 
# generate the input file list (SRR accession numbers)
cd $InputDir
fileList=$(ls *SRR* | awk -F. '{print $1}' | sort | uniq)
 
# set output file dir
mkdir -p STAR
Fo="STAR"
module load star
 
# the for loop below assigns each item in the list ($fileList) to the variable (RSSnumber) one by one and executes the commands between do and done for eachi item.
for RSSnumber in $fileList; do 
	date
	STAR \
	--runThreadN 4 \
	--genomeDir human_annotation/gencode.v45.primary_assembly.basic.annotation.gtf \
	--readFilesIn $RSSnumber".fastq.gz" \
	--readFilesCommand zcat \
	--outFileNamePrefix $Fo/$RSSnumber"_" \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMstrandField intronMotif
	date
done

# not sure if primary assembly or basic notation must check to see which one for the genomeDir flag
#line 51: add original RSSnumber_ before the default output file name (Aligned.sortedByCoord.out.bam)
 
 
##########################
### Featurecount job part
courseData="/groups/guangyao/416a516a"
MP="$courseData/modules/subread-2.0.3-Linux-x86_64/bin"
annotation="$courseData/hs.genome.annotation/gencode.v26.primary_assembly.annotation.gtf"
IP=$(ls /xdisk/guangyao/416a516a/516_projects/Group_2/Amy/STAR/*.bam)
 
### set output file dir
mkdir -p /xdisk/guangyao/416a516a/516_projects/Group_2/Amy/featureCounts
OP="/xdisk/guangyao/416a516a/516_projects/Group_2/Amy/featureCounts"
 
date
$MP/featureCounts \
-T 2 \
-a $annotation \
-o $OP/featureCounts_PE.txt \
$IP
date
 
echo 'finished featureCounts using course data'
 
 
echo 'started on online retrieved data'
# website: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/
 
annotation2='/xdisk/guangyao/416a516a/516_projects/Group_2/Amy/human_annotation/gencode.v45.primary_assembly.basic.annotation.gtf'
date
$MP/featureCounts \
-T 2 \
-a $annotation2 \
-o $OP/featureCounts_PE_2.txt \
$IP
date
echo 'finished featureCounts using online data'
 
 

