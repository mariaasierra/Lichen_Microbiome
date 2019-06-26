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

#Rarefaction
calculate_rarefaction_curves <- function(physeq.lichen, measures, depths) {
  estimate_rarified_richness <- function(physeq.tree, measures, depth) {
    if(max(sample_sums(physeq.lichen)) < depth) return()
    physeq.tree <- prune_samples(sample_sums(physeq.lichen) >= depth, physeq.tree)
    
    rarified_psdata <- rarefy_even_depth(physeq.lichen, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }  
  
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, physeq.tree = physeq.tree, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(physeq.tree, c('Observed'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 100))
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(physeq.tree)), by.x = 'Sample', by.y = 'row.names')


ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Lichen,
    group = Sample)) + geom_line() +
  xlab("No. of sequences") + ylab("No. of OTUs")+ scale_colour_manual(values =c( "chocolate3","maroon","olivedrab",
                                                                                 "mediumaquamarine","steelblue","tan1","pink2"))
dev.off()
