library(tximport)
library(rhdf5)

# Read configuration from Snakemake
tx2gene_file <- snakemake@input[["tx2gene"]]
abundance_files <- snakemake@input[["abundance_files"]]
output_file <- snakemake@output[["expression"]]

# Read the annotation file
tx2gene <- read.csv(tx2gene_file)

# Name the abundance files
names(abundance_files) <- basename(dirname(abundance_files))

# Import Kallisto results
txi.kallisto.g <- tximport(abundance_files, type = "kallisto", tx2gene = tx2gene)

# Write expression data
write.table(txi.kallisto.g$abundance, file = output_file,
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE)
