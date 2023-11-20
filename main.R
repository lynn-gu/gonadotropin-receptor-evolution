# dependencies
library(Biostrings)
library(DECIPHER)
library(MSA2dist)
library(dplyr)
library(ggplot2)
library(ggeasy)

# read DNA seqs
read.seqs <- function(seq.path) {
  seq <- readDNAStringSet(seq.path, format = "FASTA")
  return(seq)
}

# calculating and plotting conservation matrix
csvn.analysis <- function(seq){
  aln.score <- stringDist(lhcgr.e, method = "levenshtein")
  csvn.matrix <- msaConservationScore(aln.score, gapvsGap = 0)
}

codon.position <- function(start.c, mat){
  start.p <- mat |> as.data.frame() |> select ("Human") |> mutate(row_number = row_number()) |> 
    filter(Human != "---") |> slice(start.c) |> 
    pull(row_number)
  return(start.p)
}

# calculating and plotting nonsynonymous v synonymous ratio
syn.nonsyn.analysis <- function(seq, domain.p, recep.name){
  seq.xy <- seq |> dnastring2codonmat(shorten=T) 
  seq.xy <- seq.xy[apply(seq.xy, MARGIN = 1, function(x) sum(x != "---") > 0.9*length(x)),]
  domains.plotting <- as.data.frame(lapply(domain.p$codon, codon.position, seq.xy))
  seq.xy <- codonmat2xy(seq.xy)
  seq.xy.m <- mutate(seq.xy, NonSynRatio = NonSynMean/SynMean) |> filter(is.finite(NonSynRatio)) |> select(Codon, NonSynRatio)
  ggplot(data = seq.xy.m) +
    geom_segment(mapping = aes(x = Codon, xend = Codon, y = 0, yend = NonSynRatio)) +
    geom_segment(mapping = aes(x = 1, xend = domains.plotting[, 2], y = max(NonSynRatio), yend = max(NonSynRatio), size = 2, color = "Extracellular")) +
    geom_segment(mapping = aes(x = domains.plotting[, 2], xend = domains.plotting[, 3], y = max(NonSynRatio), yend = max(NonSynRatio), size = 2, color = "Serpentine")) +
    geom_segment(mapping = aes(x = domains.plotting[, 3], xend = max(Codon), y = max(NonSynRatio), yend = max(NonSynRatio), size = 2, color = "Intracellular")) +
    ylab("Nonsynonymous/Synonymous Substitution Ratio") +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    labs(color = "domain", title = recep.name) + 
    guides(size = none) + theme_bw() +
    easy_center_title() + easy_remove_gridlines()
}

lhcgr.domains <- data.frame(codon = c(25, 362, 627))
fshr.domains <- data.frame(codon = c(25, 362, 638))
tshr.domains <- data.frame(codon = c(34, 361, 632))

lhcgr.e.path <- "Eutherian_LHCGR_orthologues.fa"
fshr.e.path <- "Eutherian_FSHR_orthologues.fa"
tshr.e.path <- "Eutherian_TSHR_orthologues.fa"
lhcgr.e <- read.seqs(lhcgr.e.path)
fshr.e <- read.seqs(fshr.e.path)
tshr.e <- read.seqs(tshr.e.path)
print(syn.nonsyn.analysis(lhcgr.e, lhcgr.domains, "Eutherian LHCGR"))
print(syn.nonsyn.analysis(fshr.e, fshr.domains, "Eutherian FSHR"))
print(syn.nonsyn.analysis(tshr.e, tshr.domains, "Eutherian TSHR"))

lhcgr.p.path <- "Primates_LHCGR_orthologues.fa"
fshr.p.path <- "Primates_FSHR_orthologues.fa"
tshr.p.path <- "Primates_TSHR_orthologues.fa"
lhcgr.p <- read.seqs(lhcgr.p.path)
fshr.p <- read.seqs(fshr.p.path)
tshr.p <- read.seqs(tshr.p.path)
print(syn.nonsyn.analysis(lhcgr.p, lhcgr.domains, "Primate LHCGR"))
print(syn.nonsyn.analysis(fshr.p, fshr.domains, "Primate FSHR"))
print(syn.nonsyn.analysis(tshr.p, tshr.domains, "Primate TSHR"))
