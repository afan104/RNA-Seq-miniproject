# RNA-Seq-miniproject
Mini RNA-Seq project using differential analysis (DESeq2) and pathway analysis (GSEA, g:Profiler).

# Description of Data
11 samples

# Accessing Data
Accessed data from the GEO database using accession number GSE38006. From the GEO Database, I accessed the SRP number from the "SRA Run Selector" and input SRP013224 into the ENA database. Downloaded the tsv file containing the ftp directories. 

Created acc_list.txt file containing all the values under "fastq_ftp" column and pasted directly into the HPC using vim. Then ran in vim `:%s/ftp/ftp\:\/\/ftp/g` to globally replace ftp with the directory address ftp://ftp (using escape characters). 

Ran wget -i acc_list.txt to download raw reads into rawReads directory.

# Data Pipeline I (Processing Raw Reads)
FastQC (+multiQC) -> fastp -> STAR (+multiQC) -> featureCounts (+multiQC)

STAR requires the human index

featureCounts uses human annotation (https://www.gencodegenes.org/human/release_38.html)

# Data Pipeline II (DE Analysis)


# Data Pipeline III (Pathway Analysis - GSEA, g:Profiler)
