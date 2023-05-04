library("tximport")
library("readr")
library("DESeq2")

# Adapted from http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#transcript-abundance-files-and-tximport-tximeta

samples <- read.table("samples.txt", header=TRUE)

dir_names = paste(samples$sample, ".fastq.quant", sep="")

files <- file.path("quant", dir_names, "quant.sf")
names(files) <- samples$sample
txi <- tximport(files, type="salmon")

# In our case each entry is a gene, but we still need a map of genes to transcripts
# loop over all files and create a vector of the unique Names
gene_names = c()
for (i in 1:length(files)) {
	q = read.table(files[[i]], header=TRUE, sep="\t")
	gene_names = c(gene_names, q$Name)
}
# remove duplicates
gene_names = unique(gene_names)
# sort the vector
gene_names = sort(gene_names)

# create a data frame with two columns, gene and transcript, where each row is a gene
tx2gene = data.frame(gene=gene_names, transcript=gene_names, stringsAsFactors=FALSE)

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsNanomia <- DESeqDataSetFromTximport(txi,
	colData = samples,
	design = ~ samples$tissue)
