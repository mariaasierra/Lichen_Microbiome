library(ape)
library(ggtree)
library(utils)
library(ggplot2)

#mafft --reorder --anysymbol --maxiterate 1000 --6merpair input 

#Read tree from MAFFT
tree<-read.newick("../Data/tree.ref.core.txt")
ring<-read.table("../Data/lichen.ring.txt")
tax<-read.table("../Data/taxa.txt")
phylum<-read.table("../Data/phylum.ring.txt")

data.lich<-data.frame(Taxa = tax,
               Ring= ring,
               Phylum=phylum)
row.names(data.lich) <- NULL


lich.tree<-ggtree(tree,layout = "rectangular", aes(colour=Phylum), size=0.5) %<+% data.lich + 
  geom_aline(size=1)+ geom_tippoint(aes(colour=Phylum))+
  geom_tiplab( offset = 0.03, size=4,hjust = 0, align = T,linetype="dotted",aes(angle=0))+
  scale_colour_brewer(palette =  "Set2",
                      na.value="gray30") 
#Heatmap ring
gheatmap(lich.tree, ring, offset = 0.07,width=0.02,font.size=0,
         colnames_offset_y = 2, 
         colnames_angle=90,
         colnames_position = "top", hjust = 1)+
    scale_fill_manual(values = c("gray90","#055575"),
                      na.value="white",breaks=c("Yes","No")) +
  theme_tree(legend.position="right", 
             legend.text = element_text(colour="gray30", size = 20), 
             legend.title=element_text(size = 25),
             legend.key.size=unit(1.5,"line")) + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(fill = "Lichen host")


