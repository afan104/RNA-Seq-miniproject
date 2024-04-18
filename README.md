# RNA-Seq-miniproject
Mini RNA-Seq project using differential analysis (DESeq2) and pathway analysis (GSEA, g:Profiler).

# Description of Data

# Accessing Data
Accessed data from the GEO database using accession number GSE38006. From the GEO Database, I accessed the SRP number from the "SRA Run Selector" and input SRP013224 into the ENA database. Downloaded the tsv file containing the ftp directories.

# Data Pipeline I (Processing Raw Reads)
FastQC (+multiQC) -> fastp -> STAR (+multiQC) -> featureCounts (+multiQC)

# Data Pipeline II (DE Analysis)


# Data Pipeline III (Pathway Analysis - GSEA, g:Profiler)
