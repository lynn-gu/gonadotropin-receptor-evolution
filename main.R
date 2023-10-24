# dependencies
library(sys)
library(Biostrings)
library(DECIPHER)
library(MSA2dist)
library(dplyr)
library(ggplot2)

# read mRNA seqs and trim non-coding region with clipkit
read.seqs <- function(seq.path) {
  # exec_wait(
    # "~/anaconda3/bin/clipkit",
    # args = seq.path,
    # std_out = TRUE
  # )
  # seq.path <- base::paste(seq.path, ".clipkit", sep = "")
  seq <- readDNAStringSet(seq.path, format = "FASTA")
  return(seq)
}

# calculating and plotting conservation matrix
csvn.analysis <- function(seq){
  aln.score <- stringDist(lhcgr.e, method = "levenshtein")
  csvn.matrix <- msaConservationScore(aln.score, gapvsGap = 0)
}

# calculating nonsynonymous v synonymous ratio
syn.nonsyn.analysis <- function(seq){
  # seq.xy <- seq |> dnastring2codonmat() |> codonmat2xy()
  seq.xy <- seq |> 
    # RemoveGaps(removeGaps = "all") |> AlignTranslation() |> 
    dnastring2codonmat(shorten=T) 
  seq.xy <- seq.xy[apply(seq.xy, MARGIN = 1, function(x) sum(x != "---") > 0.9*length(x)),] |> codonmat2xy()
  seq.xy.m <- mutate(seq.xy, NonSynRatio = NonSynMean/SynMean) |> filter(is.finite(NonSynRatio)) |> select(Codon, NonSynRatio)
  return(seq.xy.m)
}

lhcgr.e.path <- "Eutherian_LHCGR_orthologues.fa"
fshr.e.path <- "Eutherian_FSHR_orthologues.fa"
tshr.e.path <- "Eutherian_TSHR_orthologues.fa"
lhcgr.e <- read.seqs(lhcgr.e.path)
fshr.e <- read.seqs(fshr.e.path)
tshr.e <- read.seqs(tshr.e.path)
lhcgr.e.xy <- syn.nonsyn.analysis(lhcgr.e)
fshr.e.xy <- syn.nonsyn.analysis(fshr.e)
tshr.e.xy <- syn.nonsyn.analysis(tshr.e)

lhcgr.p.path <- "Primates_LHCGR_orthologues.fa"
fshr.p.path <- "Primates_FSHR_orthologues.fa"
tshr.p.path <- "Primates_TSHR_orthologues.fa"
lhcgr.p <- read.seqs(lhcgr.p.path)
fshr.p <- read.seqs(fshr.p.path)
tshr.p <- read.seqs(tshr.p.path)
lhcgr.p.xy <- syn.nonsyn.analysis(lhcgr.p)
fshr.p.xy <- syn.nonsyn.analysis(fshr.p)
tshr.p.xy <- syn.nonsyn.analysis(tshr.p)

# plotting eutherian data
ggplot() +
  geom_segment(data = lhcgr.e.xy, mapping = aes(x = Codon, xend = Codon, y = 0, yend = NonSynRatio, color = "LHCGR")) +
  geom_segment(data = fshr.e.xy, mapping = aes(x = Codon, xend = Codon, y = 0, yend = NonSynRatio, color = "FSHR")) +
  geom_segment(data = tshr.e.xy, mapping = aes(x = Codon, xend = Codon, y = 0, yend = NonSynRatio, color = "TSHR")) +
  ylab("Nonsynonymous/Synonymous Substitution Ratio") +
  labs(color = "receptors")

# plotting primates data
ggplot() +
  geom_segment(data = lhcgr.p.xy, mapping = aes(x = Codon, xend = Codon, y = 0, yend = NonSynRatio, color = "LHCGR")) +
  geom_segment(data = fshr.p.xy, mapping = aes(x = Codon, xend = Codon, y = 0, yend = NonSynRatio, color = "FSHR")) +
  geom_segment(data = tshr.p.xy, mapping = aes(x = Codon, xend = Codon, y = 0, yend = NonSynRatio, color = "TSHR")) +
  ylab("Nonsynonymous/Synonymous Substitution Ratio") +
  labs(color = "receptors")
