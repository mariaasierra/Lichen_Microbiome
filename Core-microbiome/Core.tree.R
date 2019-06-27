#Core microbiome tree
library(phylotools)
library(ggtree)
library(Biostrings)
library(ape)
#OTUs fasta file
lichen.fasta<-readDNAStringSet("../Data/Lichen_fasta.fasta", format = "fasta")
names(lichen.fasta) <- names(lichen.fasta) %>% # Rename  OTUs
  strsplit(., split="|",fixed=TRUE) %>%
  sapply(.,'[', 2) %>% gsub("\t", "_",.)
fastacore<-lichen.fasta[c(which(names(lichen.fasta) %in% b1$OTU))] ## extract OTUs from main file 
writeXStringSet(fastacore, filepath = "fasta.core",format="fasta")
lichen.core.fasta<-rename.fasta(infile = "fasta.core", ref_table = us, outfile = "core.rename.fasta") ##rename OTUs in fasta file with the OTUs taxonomy 
#Align with MAFFT the lichen core fasta file
#MAFFT-FFT-NS-i 
#mafft --reorder --anysymbol --maxiterate 1000 --6merpair input 
#Then create tree with MAFFT alignment in MAFFT v.7 

tree.core<-read.nexus("core.prev.tree.txt") ##Nexus tree
ggtree(tree.core,layout = "daylight",branch.length = "none")+
  geom_tiplab(size=2,  color="black", hjust = 0)  + 
  geom_hilight_encircle(node = 19, fill = "steelblue", alpha = 0.3) +
  geom_hilight_encircle(node = 24, fill = "khaki3", alpha = 0.3) +
  geom_hilight_encircle(node = 28, fill = "orangered2", alpha = 0.3)
