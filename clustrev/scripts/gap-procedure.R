#!/usr/bin/Rscript

args <- commandArgs(TRUE)
if (length(args) != 2) {
	stop("Usage: Rscript gap-procedure.R [input FASTA] [output CSV]; ", length(args), " args given");
	return;
}

## Test values
#fasta.file <- '~/git/papers/clustrev/data/Both.0_TRUE.fas'
#output.csv <- '~/git/papers/clustrev/results/test.csv'

fasta.file <- args[1]  # path to FASTA-format file containing aligned sequences
output.csv <- args[2]

require(ape, quietly=TRUE)
require(GapProcedure, quietly=TRUE)


aln <- read.dna(fasta.file, format="fasta", comment.char=" ", as.matrix=FALSE)

gp <- GapProcedure(as.matrix(aln), 'aK80')
group.counts <- table(gp$class)
labels <- gsub("[ ]+$", "", names(aln))

res <- data.frame(label=labels, group=gp$class, in.cluster=grepl("_1_", names(aln)), predict=(group.counts[gp$class]>1))

write.csv(res, file=output.csv, quote=FALSE, row.names=FALSE)
