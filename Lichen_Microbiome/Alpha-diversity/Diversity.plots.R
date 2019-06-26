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

#Diversity plots
plot_bar(physeq.lichen,  x="Sample", y="Abundance", fill = "Phylum") + geom_bar(stat="identity")  + 
  labs(x = "", y = "Relative Abundance") +  
  scale_fill_manual(values=c( "cadetblue3","gold2" ,"salmon4","gray45",
                              "darkorange","goldenrod","navajowhite4",
                              "khaki3","pink4","rosybrown3" ,"darkolivegreen4",
                              "plum3","palevioletred", "mistyrose", "darkorange4",
                              "ivory2","olivedrab3", "darkcyan", "goldenrod3",
                              "indianred2")) + 
  theme(plot.title = element_text(size=20,lineheight=.8, vjust=1, hjust = 0.5,
                                  family="Arial", face = "bold",margin = margin(10, 0, 10, 0)),
        axis.text.y = element_text(color="black", size=12), 
        axis.text.x = element_text(face = "bold"),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 15))

### Most abundant Bacteria
taxic <- as.data.frame(physeq.lichen@tax_table) 
taxic$OTU <- rownames(taxic)  
colnames(taxic)
taxmat <- as.matrix(taxic)
new.tax <- tax_table(taxmat) 
tax_table(physeq.lichen) <- new.tax 
#Phylum level
ps1.com.phy <- microbiome::aggregate_taxa(physeq.lichen, "Phylum", top = 5)
lich.comp <- microbiomeutilities::phy_to_ldf(ps1.com.phy, transform.counts = "compositional")
ggstripchart(lich.comp, "Lichen", "Abundance", 
             facet.by = "Phylum", color = "Lichen",
             size = 2, add = c("boxplot", "mean_sd"),
             add.params = list(color = "gray40", size=0.3,fill="gray90") )+
  theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
                      legend.position = "none") + ggtitle("Phylum")

# Class level
ps1.com.class <- microbiome::aggregate_taxa(physeq.lichen, "Class", top = 5)
lich.comp <- microbiomeutilities::phy_to_ldf(ps1.com.class, transform.counts = "compositional")
ggstripchart(lich.comp, "Lichen", "Abundance", 
             facet.by = "Class", color = "Lichen",
             size = 2, add = c("boxplot", "mean_sd"),
             add.params = list(color = "gray40", size=0.3,fill="gray90") )+
  theme_test()+ theme(axis.text.x = element_blank(),
                      legend.position = "right") + ggtitle("Class")

