library(phyloseq)
library(ALDEx2)
library(ggstance)
#Lichen microbiomes comparison with ALDEx2 
OTU = read.table("../Data/OTU_table.txt", header=TRUE, sep="\t")
tax = read.table("../Data/taxonomy_table.txt", header=TRUE, sep="\t")
row.names(OTU) = OTU$Group
OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group"))]
OTU.UF = otu_table(as.matrix(OTU.clean), taxa_are_rows=F)
t.otus <- t(OTU.UF) #Transpose OTU
conds.lichen <- c(rep("Cladonia", 9), rep("Cora", 7), rep("Hypotrachyna", 7),
                  rep("Peltigera", 3),rep("Stereocaulon", 7),rep("Sticta", 7),rep("Usnea", 7)) #Lichen genera are the conditions to compare
x.clr <- aldex.clr(t.otus, conds.lichen, denom = "all") #clr: centred log-ratio transformation
x.al <- aldex.glm(x.clr, conds.lichen) #test glm: general linear model and Kruskal Wallace
OTU.aldex = x.al[which(x.al$kw.ep < 0.05 | x.al$glm.ep < 0.05),] #OTUs with p<0.05  #glm ANOVA

#Generate new phyloseq file with OTUs differentially abundant between lichen genera p<0.05 from ALDEx2
OTU.aldex <- cbind(rownames(OTU.aldex), data.frame(OTU.aldex, row.names=NULL))
colnames(OTU.aldex)[1] <- "OTU"
taxa.pvalues.aldex <- subset(tax, subset=OTU %in% OTU.aldex$OTU, select = c("Taxonomy")) ##Extract Taxonomy for ALDEx2 OTUs
row.names(taxa.pvalues.aldex) = OTU.aldex$OTU
tax.clean = separate(taxa.pvalues.aldex, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
tax.UF = tax_table(as.matrix(tax.clean))
OTU.clean<-OTU.UF[,OTU.aldex$OTU]
OTU.aldex.uf = otu_table(as.matrix(OTU.clean), taxa_are_rows=F)
meta = read.table("lichen_metadata.txt", header=TRUE, row.names=1, sep="\t", dec = ".") #Sample names MUST be the same in both files, otherwise is going not going to pair data 
meta<-sample_data(meta)
physeq = phyloseq(OTU.aldex.uf, tax.UF, meta)
physeq.lichen.aldex = merge_phyloseq(physeq) ## New phyloseq file with ALDEx2 OTUs

#Lichen PCoA by genera
ord = ordinate(physeq.lichen.aldex, "PCoA", "bray")
(ordplot <- plot_ordination(physeq.lichen.aldex, ord, "Samples", color="Lichen", axes = 1:2))
ordplot +  stat_ellipse(type = "t",linetype = 2,alpha=0.5) + geom_point(size=4, alpha=0.5) +
  geom_point(size=4, alpha=0.5)+ ggtitle("OTUs Aldex") + theme_test() + 
  scale_color_manual(values =c( "chocolate3","maroon","olivedrab",
                                "mediumaquamarine","steelblue","tan1","pink2"))

#Heatmap presence/absence all samples
OTU.aldex.bm[OTU.aldex.uf>0] <-1 #Convert to binary presence/absence
physeq = phyloseq(OTU.aldex.bm, tax.UF, meta) 
physeq.heatmap.aldex = merge_phyloseq(physeq)
OTU2 = as(otu_table(physeq.heatmap.aldex), "matrix")
OTUdf = as.data.frame(OTU2)
OTUfd= t(OTUdf) #Transpose matrix
data.dist <- vegdist(OTUfd, method = "bray")
row.clus <- hclust(data.dist, "mcquitty") #hclust with WPGMA
ann<- read.table("../Data/Lichen.annotations.txt", header = T)
col.pal <- RColorBrewer::brewer.pal(2, "PuBuGn")
ann_colors = list(Lichen = c(Cladonia="chocolate3",
                             Cora="maroon",Hypotrachyna="olivedrab",
                             Peltigera="mediumaquamarine",Stereocaulon="steelblue",
                             Sticta="tan1",Usnea="pink2"))
pheatmap(as.matrix(OTUfd), border_color = "black",fontsize_row=6, 
         fontsize_col = 10, color = col.pal, scale = "none",
         cluster_cols = F, cluster_rows = row.clus,
         annotation_col  = ann, labels_col = NULL,annotation_colors  = ann_colors)

#Tree and relative abundance of ALDEx2 OTUs
tax<-read.table("../Data/aldex.tax.txt") #taxonomy of the tree
phylum <- read.table("../Data/aldex.phylum.txt") #File with each node phylum
tree.aldex<-read.nexus("../Data/aldex.tree.txt") ##Nexus tree
otus<-as.data.frame(read.csv("../Data/relabun.otus.csv"))
df <- data.frame(taxa = tax,
                 Phylum= phylum)
row.names(df) <- NULL
gt<-ggtree(tree.aldex,layout = "rectangular", open.angle = T) %<+% df + geom_tippoint(aes(color=Phylum))+
  geom_tiplab(size=0, align = T, linesize = 0.1, color="black") 
facet_plot(gt, panel = "Rel. Abundance", 
           data= otus,
           geom_barh, mapping = aes(x = value, fill = as.factor(lichen)), stat="identity") +
  scale_fill_manual(values =c( "chocolate3","maroon","olivedrab",
                               "mediumaquamarine","steelblue","tan1","pink2")) + 
  theme(strip.text = element_text(size=15), panel.spacing =unit(0.2, "lines"), legend.position = "left") +
  labs(fill="Lichen") + guides(color = guide_legend(override.aes = list(size=5)))






