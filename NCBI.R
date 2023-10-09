# dependencies
library(reticulate)
library(Biostrings)
library(DECIPHER)
library(MSA2dist)
library(dplyr)
library(ggplot2)

# read mRNA seqs and trim non-coding region
read.seqs <- function(seq.path) {
  seq <- readDNAStringSet(seq.path) |> TrimDNA(leftPatterns = "ATGGCGGGCC", rightPatterns = "TAACTGCAT", type = "sequences")
  # seq.aln <- AlignTranslation(seq, readingFrame = NA)
  seq.aln <- AlignSeqs(seq)
}

# calculating and plotting conservation matrix
csvn.analysis <- function(seq){
  aln.score <- stringDist(lhcgr.e, method = "levenshtein")
  csvn.matrix <- msaConservationScore(aln.score, gapvsGap = 0)
}

# calculating and plotting nonsynonymous v synonymous ratio
syn.nonsyn.analysis <- function(seq){
  seq.xy <- seq |> dnastring2codonmat() |> codonmat2xy()
  seq.xy.m <- mutate(seq.xy, NonSynRatio = NonSynMean/SynMean) |> filter(is.finite(NonSynRatio)) |> select(Codon, NonSynRatio)
  ggplot(seq.xy.m, mapping = aes(x = Codon, y = NonSynRatio)) +
    geom_line() +
    ylab("Nonsynonymous/Synonymous Substitution Ratio")
}

# lhcgr.e.path <- "Eutherian_LHCGR_refseq.fasta"
# lhcgr.e <- read.seqs(lhcgr.e.path)
# syn.nosyn.analysis(lhcgr.e)
# fshr.e.path <- "Eutherian_FSHR_refseq.fasta"
# fshr.e <- read.seqs(fshr.e.path)
# tshr.e.path <- "Eutherian_TSHR_refseq.fasta"
# tshr.e <- read.seqs(tshr.e.path)

lhcgr.p.path <- "Primates_LHCGR_refseq.fasta"
lhcgr.p <- read.seqs(lhcgr.p.path)
# syn.nonsyn.analysis(lhcgr.p)
# fshr.p.path <- "Primates_FSHR_refseq.fasta"
# fshr.p <- read.seqs(fshr.p.path)
# syn.nonsyn.analysis(fshr.p)
# tshr.p.path <- "Primates_TSHR_refseq.fasta"
# tshr.p <- read.seqs(fshr.p.path)
# syn.nonsyn.analysis(tshr.p)