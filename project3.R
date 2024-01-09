library(seqinr)
library(ape)
cleanAlignment <- function(alignment, minpcnongap, minpcid)
{
 # make a copy of the alignment to store the new alignment in:
 newalignment <- alignment
 # find the number of sequences in the alignment
 numseqs <- alignment$nb
 # empty the alignment in "newalignment")
 for (j in 1:numseqs) { newalignment$seq[[j]] <- "" }
 # find the length of the alignment
 alignmentlen <- nchar(alignment$seq[[1]])
 # look at each column of the alignment in turn:
 for (i in 1:alignmentlen)
 {
	# see what percent of the letters in this column are non-gaps:
	nongap <- 0
	for (j in 1:numseqs)
	{
	   seqj <- alignment$seq[[j]]
	   letterij <- substr(seqj,i,i)
	   if (letterij != "-") { nongap <- nongap + 1}
	}
	pcnongap <- (nongap*100)/numseqs
	# Only consider this column if at least minpcnongap % of the letters are not gaps:
	if (pcnongap >= minpcnongap)
	{
	   # see what percent of the pairs of letters in this column are identical:
	   numpairs <- 0; numid <- 0
	   # find the letters in all of the sequences in this column:
	   for (j in 1:(numseqs-1))
	   {
		  seqj <- alignment$seq[[j]]
		  letterij <- substr(seqj,i,i)
		  for (k in (j+1):numseqs)
		  {
			 seqk <- alignment$seq[[k]]
			 letterkj <- substr(seqk,i,i)
			 if (letterij != "-" && letterkj != "-")
			 {
				numpairs <- numpairs + 1
				if (letterij == letterkj) { numid <- numid + 1}
			 }
		  }
	   }
	   pcid <- (numid*100)/(numpairs)
	   # Only consider this column if at least %minpcid of the pairs of letters are identical:
	   if (pcid >= minpcid)
	   {
		   for (j in 1:numseqs)
		   {
			  seqj <- alignment$seq[[j]]
			  letterij <- substr(seqj,i,i)
			  newalignmentj <- newalignment$seq[[j]]
			  newalignmentj <- paste(newalignmentj,letterij,sep="")
			  newalignment$seq[[j]] <- newalignmentj
		   }
	   }
	}
 }
 return(newalignment)
}
alignment <- read.alignment(file="/work2/09269/saip/stampede2/projects/project3/humanSeqAligned.fasta", format="fasta")
clean_alignment <- cleanAlignment(alignment, minpcnongap = 30, minpcid = 30)
DNAbin_alignment <- as.DNAbin(clean_alignment)
dist_matrix <- dist.dna(DNAbin_alignment)
NJ_tree <- nj(dist_matrix)
rooted_tree <- root(NJ_tree, "NC_001643.1", r = TRUE)
pdf("/work2/09269/saip/stampede2/projects/project3/project3tree.pdf")
plot.phylo(rooted_tree, type = "p")
dev.off()
pca_output <- prcomp((dist_matrix), scale=T)
pdf("/work2/09269/saip/stampede2/projects/project3/project3pca.pdf")
plot(pca_output$x[, 1], pca_output$x[, 2], xlab = "PC1", ylab = "PC2", main = "PCA Plot",xlim=c(-35,5))
text(pca_output$x[, 1], pca_output$x[, 2], labels = row.names(pca_output$x), pos = 4)
dev.off()
