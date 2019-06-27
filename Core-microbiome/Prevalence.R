### Prevalence Total Lichen Microbiome ###
library(phyloseq)
library(data.table)
#Prepare phyloseq object 
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
meta = read.table("Lichen_metadata.txt", header=TRUE, row.names=1, sep="\t", dec = ".") #Sample names MUST be the same in both files, otherwise is going not going to pair data 
meta<-sample_data(meta)
physeq = phyloseq(OTU.UF, tax.UF, meta)
physeq.lichen = merge_phyloseq(physeq)

#Prevalence 
lichen_melt = function(physeq,
                     includeSampleVars = character(),
                     omitZero = FALSE){
  # supports otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "TaxaID")
  # Enforce character TaxaID key
  otudt[, TaxaIDchar := as.character(TaxaID)]
  otudt[, TaxaID := NULL]
  setnames(otudt, "TaxaIDchar", "TaxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "TaxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID] #Sum of counts
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. 
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "TaxaID")
    # Enforce character TaxaID key
    taxdt[, TaxaIDchar := as.character(TaxaID)]
    taxdt[, TaxaID := NULL]
    setnames(taxdt, "TaxaIDchar", "TaxaID")
    # Join with tax table
    setkey(taxdt, "TaxaID")
    setkey(mdt, "TaxaID")
    mdt <- taxdt[mdt]
  }
  setkey(mdt, "TaxaID")
  return(mdt)
}

mdt = lichen_melt(physeq.lichen)
lichen.prev = mdt[, list(Prevalence = mean(count > 0), #Mean of samples with counts greater than cero
                    TotalCounts = mean(count)),
             by = TaxaID]
addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)])) #Color taxa by phylum
# Join by TaxaID
setkey(lichen.prev, TaxaID)
setkey(addPhylum, TaxaID)
lichen.prev <- addPhylum[lichen.prev]
showPhyla = lichen.prev[, sum(TotalCounts), by = Phylum][order(-V1)][1:10]$Phylum
setkey(lichen.prev, Phylum)
#Prevalence plot of total OTUs
ggplot(lichen.prev[showPhyla], 
       mapping = aes(Prevalence, TotalCounts, color = Phylum)) + 
  geom_jitter(size = 5, alpha = 0.7) + 
  annotate(geom='label', x=0.3, y=9999, label="Pan", fontface=2) +
  annotate(geom='label', x=0.05, y=9999, label="Peripheral", fontface=2) +
  annotate(geom='label', x=0.95, y=9999, label="Core", fontface=2) +
  geom_vline(xintercept=0.25, color='gray20') +
  geom_vline(xintercept=0.9, color='gray20') +
  scale_y_log10() +  theme_minimal() + ylab("Total Counts") +
  theme(strip.text = element_text(size=20),
        axis.text.x = element_text(colour = "black", size=11 ), 
        axis.text.y = element_text(colour = "black", size = 11),
        axis.title.x = element_text(size=13), 
        axis.title.y = element_text(size=13),
        legend.text  = element_text(size=12, colour = "gray40"),
        legend.title = element_text(size = 14))




