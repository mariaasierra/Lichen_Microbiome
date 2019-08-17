library(phyloseq)
library(tidyr)
library(plyr)
library(reshape2)
library(microbiome)
library(microbiomeutilities)
library(Biostrings)
library(phylotools)
library(data.table)
library(ape)
library(ggtree)
#Prepare Phyloseq object 
OTU = read.table("../Data/OTU_table.txt", header=TRUE, sep="\t")
tax = read.table("../Data/taxonomy_table.txt", header=TRUE, sep="\t")
row.names(OTU) = OTU$Group
OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group"))]
row.names(tax) = tax$OTU
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]
tax.clean = separate(tax.clean, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep=";")
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("Size", "Strain", "OTU"))]
OTU.UF = otu_table(as.matrix(OTU.clean), taxa_are_rows=F)
tax.UF = tax_table(as.matrix(tax.clean))
meta = read.table("../Data/Metadata.txt", header=TRUE, row.names=1, sep="\t", dec = ".") #Sample names MUST be the same in both files, otherwise is going not going to pair data 
meta<-sample_data(meta)
physeq = phyloseq(OTU.UF, tax.UF, meta)
physeq.lichen = merge_phyloseq(physeq)

#Subset OTUs present in each lichen 
OTU.Usnea = colnames(OTU.UF[meta$Lichen == "Usnea", apply(OTU.UF[meta$Lichen == "Usnea",], MARGIN=2, function(x) any(x >0))])
OTU.Stereo = colnames(OTU.UF[meta$Lichen == "Stereocaulon", apply(OTU.UF[meta$Lichen == "Stereocaulon",], MARGIN=2, function(x) any(x >0))])
OTU.Hyp = colnames(OTU.UF[meta$Lichen == "Hypotrachyna", apply(OTU.UF[meta$Lichen == "Hypotrachyna",], MARGIN=2, function(x) any(x >0))])
OTU.Cora = colnames(OTU.UF[meta$Lichen == "Cora", apply(OTU.UF[meta$Lichen == "Cora",], MARGIN=2, function(x) any(x >0))])
OTU.Sticta = colnames(OTU.UF[meta$Lichen == "Sticta", apply(OTU.UF[meta$Lichen == "Sticta",], MARGIN=2, function(x) any(x >0))])
OTU.Cladonia = colnames(OTU.UF[meta$Lichen == "Cladonia", apply(OTU.UF[meta$Lichen == "Cladonia",], MARGIN=2, function(x) any(x >0))])
OTU.Pelt = colnames(OTU.UF[meta$Lichen == "Peltigera", apply(OTU.UF[meta$Lichen == "Peltigera",], MARGIN=2, function(x) any(x >0))])

#Extract taxonomy of OTUs 

#Usnea
usnea <- subset(tax, subset=OTU %in% OTU.Usnea, select = c("Taxonomy"))
usnea<-as.data.frame(usnea)
setDT(usnea, keep.rownames = "OTU")
#Stereocaulon
stereo <- subset(tax, subset=OTU %in% OTU.Stereo, select = c("Taxonomy"))
stereo<-as.data.frame(stereo)
setDT(stereo, keep.rownames = "OTU")
#Hypotrachina
hypo <- subset(tax, subset=OTU %in% OTU.Hyp, select = c("Taxonomy"))
hypo<-as.data.frame(hypo)
setDT(hypo, keep.rownames = "OTU")
#Cora
cora <- subset(tax, subset=OTU %in% OTU.Cora, select = c("Taxonomy"))
cora<-as.data.frame(cora)
setDT(cora, keep.rownames = "OTU")
#Sticta
sticta <- subset(tax, subset=OTU %in% OTU.Sticta, select = c("Taxonomy"))
sticta<-as.data.frame(sticta)
setDT(sticta, keep.rownames = "OTU")
#Peltigera
peltigera <- subset(tax, subset=OTU %in% OTU.Pelt, select = c("Taxonomy"))
peltigera<-as.data.frame(peltigera)
setDT(peltigera, keep.rownames = "OTU")
#Cladonia
cladonia <- subset(tax, subset=OTU %in% OTU.Cladonia, select = c("Taxonomy"))
cladonia<-as.data.frame(cladonia)
setDT(cladonia, keep.rownames = "OTU")

#We have the OTUs present in each lichen, we now need to extract the sequences of each OTU for each lichen file
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings", version = "3.8")
require(Biostrings)
library(Biostrings)
lichen.fasta<-readDNAStringSet("lichen_fasta.fasta", format = "fasta")
names(lichen.fasta) <- names(lichen.fasta) %>% #Rename  OTUs
  strsplit(., split="|",fixed=TRUE) %>%
  sapply(.,'[', 2) %>% gsub("\t", "_",.)

#Extract sequences of OTUs present in each lichen

#Usnea
fastausnea<-lichen.fasta[c(which(names(lichen.fasta) %in% usnea$OTU))] ## extract OTUs from main file 
writeXStringSet(fastausnea, filepath = "usnea.fasta",format="fasta")
rename.fasta(infile = "usnea.fasta", ref_table = us, outfile = "usnea.rename.fasta") ##rename OTUs in fasta file with the taxonomy 
#Stereocaulon
fastastereo<-lichen.fasta[c(which(names(lichen.fasta) %in% stereo$OTU))]
writeXStringSet(fastastereo, filepath = "stereo.fasta",format="fasta")
rename.fasta(infile = "stereo.fasta", ref_table = ste, outfile = "stereo.rename.fasta") 
#Hypotrachyna
fastahypo<-lichen.fasta[c(which(names(lichen.fasta) %in% hypo$OTU))] 
writeXStringSet(fastahypo, filepath = "hypo.fasta",format="fasta")
rename.fasta(infile = "hypo.fasta", ref_table = hypo, outfile = "hypo.rename.fasta") 
#Cora
fastacora<-lichen.fasta[c(which(names(lichen.fasta) %in% cora$OTU))]
writeXStringSet(fastacora, filepath = "cora.fasta",format="fasta")
rename.fasta(infile = "cora.fasta", ref_table = cora, outfile = "cora.rename.fasta") 
#Sticta
fastasticta<-lichen.fasta[c(which(names(lichen.fasta) %in% sticta$OTU))] 
writeXStringSet(fastasticta, filepath = "sticta.fasta",format="fasta")
rename.fasta(infile = "sticta.fasta", ref_table = sticta, outfile = "sticta.rename.fasta") 
#Sticta
fastasticta<-lichen.fasta[c(which(names(lichen.fasta) %in% sticta$OTU))] 
writeXStringSet(fastasticta, filepath = "sticta.fasta",format="fasta")
rename.fasta(infile = "sticta.fasta", ref_table = sticta, outfile = "sticta.rename.fasta") 
#Peltigera
fastapeltigera<-lichen.fasta[c(which(names(lichen.fasta) %in% peltigera$OTU))] 
writeXStringSet(fastapeltigera, filepath = "Peltigera.fasta",format="fasta")
rename.fasta(infile = "Peltigera.fasta", ref_table = Peltigera, outfile = "Peltigera.rename.fasta")  
#Cladonia
fastacladonia<-lichen.fasta[c(which(names(lichen.fasta) %in% cladonia$OTU))] 
writeXStringSet(fastacladonia, filepath = "Cladonia.fasta",format="fasta")
rename.fasta(infile = "Cladonia.fasta", ref_table = Cladonia, outfile = "Cladonia.rename.fasta") 

#Built tree after aligning with MAFFT each lichen file

#Cladonia
tax<-read.table("../Data/Trees/cladonia.tree.tax.txt") ## Taxonomy of the tree
phylum <- read.table("../Data/Trees/cladonia.phylum.txt") #File with each node phylum
tree<-read.nexus("../Data/Trees/cladonia.tree.txt") ##Newick trees was converted to Nexus tree
dd <- data.frame(taxa = tax,
                 Phylum = phylum)
row.names(dd) <- NULL
ggtree(tree, layout ="fan", aes(color=Phylum), alpha=0.8) %<+% dd + xlim(NA,1) +
  theme_tree(legend.position="bottom")

#Cora
tax<-read.table("../Data/Trees/cora.tree.tax.txt") ## Taxonomy of the tree
phylum <- read.table("../Data/Trees/cora.phylum.txt") #File with each node phylum
tree<-read.nexus("../Data/Trees/cora.tree.txt") ##Newick trees was converted to Nexus tree
dd <- data.frame(taxa = tax,
                 Phylum = phylum)
row.names(dd) <- NULL
ggtree(tree, layout ="fan", aes(color=Phylum), alpha=0.8) %<+% dd + xlim(NA,1) +
  theme_tree(legend.position="bottom")

#Sticta
tax<-read.table("../Data/Trees/sticta.tree.tax.txt") ## Taxonomy of the tree
phylum <- read.table("../Data/Trees/sticta.phylum.txt") #File with each node phylum
tree<-read.nexus("../Data/Trees/sticta.tree.txt") ##Newick trees was converted to Nexus tree
dd <- data.frame(taxa = tax,
                 Phylum = phylum)
row.names(dd) <- NULL
gt<-ggtree(tree, layout ="fan", aes(color=Phylum), alpha=0.8) %<+% dd + xlim(NA,1) +
  theme_tree(legend.position="bottom")

#Peltigera
tax<-read.table("../Data/Trees/peltigera.tree.tax.txt") ## Taxonomy of the tree
phylum <- read.table("../Data/Trees/peltigera.phylum.txt") #File with each node phylum
tree<-read.nexus("../Data/Trees/peltigera.tree.txt") ##Newick trees was converted to Nexus tree
dd <- data.frame(taxa = tax,
                 Phylum = phylum)
row.names(dd) <- NULL
ggtree(tree, layout ="fan", aes(color=Phylum), alpha=0.8) %<+% dd + xlim(NA,1) +
  theme_tree(legend.position="bottom")

#Hypotrachyna
tax<-read.table("../Data/Trees/hypo.tree.tax.txt") ## Taxonomy of the tree
phylum <- read.table("../Data/Trees/hypo.phylum.txt") #File with each node phylum
tree<-read.nexus("../Data/Trees/hypo.tree.txt") ##Newick trees was converted to Nexus tree
dd <- data.frame(taxa = tax,
                 Phylum = phylum)
row.names(dd) <- NULL
gt<-ggtree(tree, layout ="fan", aes(color=Phylum), alpha=0.8) %<+% dd + xlim(NA,1) +
  theme_tree(legend.position="bottom")

#Stereocaulon
tax<-read.table("../Data/Trees/stereo.tree.tax.txt") ## Taxonomy of the tree
phylum <- read.table("../Data/Trees/stereo.phylum.txt") #File with each node phylum
tree<-read.nexus("../Data/Trees/stereo.tree.txt") ##Newick trees was converted to Nexus tree
dd <- data.frame(taxa = tax,
                 Phylum = phylum)
row.names(dd) <- NULL
gt<-ggtree(tree, layout ="fan", aes(color=Phylum), alpha=0.8) %<+% dd + xlim(NA,1) +
  theme_tree(legend.position="bottom")

#Usnea
tax<-read.table("../Data/Trees/usnea.tree.tax.txt") ## Taxonomy of the tree
phylum <- read.table("../Data/Trees/usnea.phylum.txt") #File with each node phylum
tree<-read.nexus("../Data/Trees/usnea.tree.txt") ##Newick trees was converted to Nexus tree
dd <- data.frame(taxa = tax,
                 Phylum = phylum)
row.names(dd) <- NULL
gt<-ggtree(tree, layout ="fan", aes(color=Phylum), alpha=0.8) %<+% dd + xlim(NA,1) +
  theme_tree(legend.position="bottom")






