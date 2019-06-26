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


pr<-plot_richness(physeq.lichen, x="Lichen",  color = "Lichen")+
  scale_color_manual(values =c( "chocolate3","maroon","olivedrab",
                                "mediumaquamarine","steelblue","tan1","pink2"))+geom_point(size=4, alpha=0.7)

pr + geom_boxplot(data = pr$data,color="gray40",size=0.3,fill="gray90", alpha = 0.5,aes(x = Lichen, y = value))+
  theme_test()+ theme(strip.text = element_text(size=20),
                      axis.text.x = element_text(face="bold", size=10, angle = 45, hjust = 1), 
                      axis.text.y = element_text(colour = "grey30", size = 10, face = "italic"),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(), legend.position = "none")
