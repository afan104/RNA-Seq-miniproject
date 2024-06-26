# DESCRIPTION: This file processes the raw counts.
# Outputs: rlog counts, DE genes, heatmap
# Note: DE genes for 
#       1. R1 v C1 (infected v control at 12 h)
#       2. R2 v C2 (infected v control at 24 h)

# steps:
# 1. obtain count.table (replace column names with condition e.g. SRR** > R1.1, remove version names in rows)
#     a. removed duplicate rows (all had 0 counts anyways)
# 2. obtain colData objects: data.frame(row.names, condition, state)
#     a. row.names as gene names (cols of count.table)
#     b. condition as experimental condition (R1, R2, C1, C2)
#     c. state as time point (infected or mock)
# 3. Define 2 DESeq Objects (1 with design as ~condition, 1 with design as ~state)
# 4. rlog the DESeq object (set blind = TRUE)
#     a. save rlog as csv for GSEA downstream
# 5. DESeq analysis on rawcounts DESeq object
# 6. Get the results and DEgenes using each of 3 contrasts:
#     a. DEgenes_1: (R1+R2) v (C1+C2)
#     b. DEgenes_2: R1 v C1 (important)
#     c. DEgenes_3: R2 v C2 (important)
#     note: put the DE genes into a csv for hypergeometric test downstream (g:Profiler)
# 7. Plot the paper-focused heatmap using the combined DEgenes_2 and DEgenes_3 lists

library(DESeq2)

setwd("~/516a/516a/Project/")

##########################################
######step 1: countData
######obtain count.table (remove version names in rows, replace column names with condition e.g. SRR** > R1.1)

# download data
count.table <- read.csv("./featureCounts/featureCounts_PE_Giang.txt", header = T, sep = "\t", skip = 1)
View(count.table)

# there are no duplicates
for (i in 1:nrow(count.table)){
  count.table$Geneid[i] <- gsub("\\..*", "", count.table$Geneid[i])
}  
TRUE %in% duplicated(count.table$Geneid)
length(count.table$Geneid[duplicated(count.table$Geneid)]) # 45 duplicated genes in the 58278 rows

# keep the copy with the largest total count number
duplicated_genes <- unique(count.table$Geneid[duplicated(count.table$Geneid)])

# check if any non NA non 0 values exist: all are 0
View(count.table[count.table$Geneid%in%duplicated_genes, 7:17])
nrow(count.table[count.table$Geneid%in%duplicated_genes, 7:17]) # 90 duplicates total

# initialize list of all duplicate indices to remove
remove_indices <- c()

# find indices to keep
for (gene in duplicated_genes) {
  counts <- count.table[count.table$Geneid==gene, 7:17]
  # create list of duplicate indices for this gene
  duplicate_indices <- which(count.table$Geneid==gene)
  # initialize list of summed counts
  num_duplicates <- nrow(counts)
  count_sums = rep(0, num_duplicates)
  largest <- 0
  largest_i <- 1
  # total counts for each duplicate
  for (i in num_duplicates) {
    curr_sum <- sum(counts[i,], na.rm=TRUE) # count NAs as 0's
    count_sums[i] <- curr_sum
    print(paste(gene, " sum: ", curr_sum))
    if (curr_sum > largest) {
      largest <- curr_sum
      largest_i <- i
      print(paste("largest sum: ", largest))
    }
  }
  # remove the largest duplicate index from duplicate indices for this gene
  duplicate_indices_remove <- duplicate_indices[-largest_i]
  # add all remaining indices to remove_indices list
  remove_indices <- c(remove_indices, duplicate_indices_remove)
}

length(remove_indices) # how many did we remove? 45 out of 58278
remove_indices <- sort(remove_indices, decreasing = T)
count.table <- count.table[-remove_indices,] # remove 45

# check if there are duplicates left
TRUE %in% duplicated(count.table$Geneid)


# set rownames to geneID
rownames(count.table) <- count.table$Geneid

# remove unneccessary cols
colnames(count.table)
count.table <- count.table[,-c(1:6)]
ncol(count.table)
colnames(count.table)

# add column names; refactor to hpi
colnames(count.table) <- c("R1.1", "R1.2", "R1.3", "R2.1", "R2.2", "R2.3",
                           "C1.1", "C1.2", "C1.3", "C2.1", "C2.3")
colnames(count.table)

# write all genes list for GSEA downstream
write.csv(rownames(count.table), "./gene_lists/full_gene_list.csv")
##########################################
######step 2: colData
colData <- data.frame(row.names=colnames(count.table), 
                      condition=c(rep("R1", 3), rep("R2", 3), rep("C1", 3), rep("C2", 2)),
                      state = c( rep("infection", 6), rep("ctrl", 5) )
)
View(colData)

##########################################
######step 3: create DESeq object and save
DESeq.ds1 <- DESeqDataSetFromMatrix (countData = count.table, colData = colData, design = ~condition) 
saveRDS(DESeq.ds1, file = "DESeq.ds1")
DESeq.ds2 <- DESeqDataSetFromMatrix (countData = count.table, colData = colData, design = ~state) 
saveRDS(DESeq.ds2, file = "DESeq.ds2")

##########################################
######step 4: rlog
DESeq.rlog1 <- rlog(DESeq.ds1, blind = TRUE)
rlog.counts_condition <- assay(DESeq.rlog1)
DESeq.rlog2 <- rlog(DESeq.ds2, blind = TRUE)
rlog.counts_state <- assay(DESeq.rlog2)

# save rlog RDS object
saveRDS(rlog.counts_condition, file = "rlog.counts_condition")
saveRDS(rlog.counts_state, file = "rlog.counts_state")
write.csv(rlog.counts_condition, file = "./GSEA_files/rlog.counts.csv") # for GSEA usage later
# doesn't matter which one is used, weird DESEQ encoding defines design which doesn't show up in csv
##########################################
######step 5:DE analysis (can act on raw counts)
DESeq.ds1 <- DESeq(DESeq.ds1)
DESeq.ds2 <- DESeq(DESeq.ds2)

# see results
# a. infected v ctrl
DE.results1<- results (DESeq.ds2, contrast = c('state', 'infection', 'ctrl'), alpha = 0.05)
# b. R1 v C1
DE.results2<- results (DESeq.ds1, contrast = c('condition','R1','C1'), alpha = 0.05)
# c. R2 v C2
DE.results3<- results (DESeq.ds1, contrast = c('condition','R2','C2'), alpha = 0.05)

##########################################
######step 6: Extract DEgenes with FDR<0.05 and FC>1.5
# Produce table of genes differentially expressed v. not
table(DE.results1$padj < 0.05 & abs(DE.results1$log2FoldChange) > 1.5) # (infected v ctrl)
table(DE.results2$padj < 0.05 & abs(DE.results2$log2FoldChange) > 1.5) # in paper (R1 v C1)
table(DE.results3$padj < 0.05 & abs(DE.results3$log2FoldChange) > 1.5) # in paper (R2 v C2)

# DEgenes on padj and logfc
DEgenes1 <- row.names (subset(DE.results1, padj < 0.05 & abs(log2FoldChange) > 1.5)) #infected v ctrl
DEgenes2 <- row.names (subset(DE.results2, padj < 0.05 & abs(log2FoldChange) > 1.5)) # R1 v C1
DEgenes3 <- row.names (subset(DE.results3, padj < 0.05 & abs(log2FoldChange) > 1.5)) # R2 v C2

# place differentially expressed genes into a list for g:Profiler
write.csv(DEgenes2, file = "./gene_lists/DEgenes_12hpi.txt", quote = F, row.names = F) # R1 v C1
write.csv(DEgenes3, file = "./gene_lists/DEgenes_24hpi.txt", quote = F, row.names = F) # R2 v C2
##########################################
######step 7: create heatmap of DEgenes
# just the up/downregulated in infected
png("./heatmap_final.png")
DE_all <- union(DEgenes2, DEgenes3)
mat_DE_all <- rlog.counts_condition[DE_all,c('R1.1', 'R1.2', 'R1.3', 'R2.1', 'R2.2', 'R2.3')]
aheatmap (mat_DE_all,
          Rowv = TRUE , Colv = TRUE,  
          distfun = "euclidean", hclustfun = "median",
          scale="row")
dev.off()