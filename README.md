# RNA-Seq-miniproject
Mini RNA-Seq project using differential analysis (DESeq2) and pathway analysis (GSEA, g:Profiler). Comparing gene expression between T cells infected with HIV 12 hours and 24 hours post infection. Interested in how the HIV virus influences gene expression within T cells over time (12h v 24h). Data based on this experiment: https://pubmed.ncbi.nlm.nih.gov/21933919/.

# Description of Data
11 samples of single end mRNA sequences.

# Accessing Data
Accessed data from the GEO database using accession number GSE38006. From the GEO Database, I accessed the SRP number from the "SRA Run Selector" and input SRP013224 into the ENA database. Downloaded the tsv file containing the ftp directories.

Created acc_list.txt file containing all the values under "fastq_ftp" column and pasted directly into the HPC using vim. Then ran in vim `:%s/ftp/ftp\:\/\/ftp/g` to globally replace ftp with the directory address ftp://ftp (using escape characters). 

Ran `wget -P parentdir/rawReads -i acc_list.txt` to download raw reads into rawReads directory.

# Data Pipeline I (Processing Raw Reads)
FastQC (+multiQC) -> fastp -> STAR (+multiQC) -> featureCounts (+multiQC)

STAR requires the human index

featureCounts uses human annotation (https://www.gencodegenes.org/human/release_38.html)

# Data Pipeline II (DE Analysis)
- Performed quality control using correlation plot and PCA plot
- Rlog for normalization between samples
- Low expression managed using shrinkage
- MA and volcano plots to visualize the differential expression
- Heatmap to compare expression between different experimental groups (12h v 24h)

# Data Pipeline III (Pathway Analysis - GSEA, g:Profiler)
GSEA is a rank-based system, so the entire list of genes are provided to determine differentially regulated gene sets.
g:Profiler is a threshold based system, so I used the DE genes list generated from Pipeline II to determine the differentially regulated gene sets.

Cytoscape can be used to visualize the g:Profiler results

**Note: since data files can be very large, some folders are skeleton folders that can be filled up by running the script.slurm file**
