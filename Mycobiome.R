# Pre-harvest maize fungal microbiome and mycotoxin contamination: Case of Zambia's Different rainfall patterns
######### Figures and Stats ##################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require (dada2); packageVersion("dada2")
require(phyloseq); packageVersion("phyloseq")
require (Biostrings); packageVersion("Biostrings")

# Sample dataframe:
samdf <- readRDS("https://github.com/bkatati/mmycobiome/samdf.rds")
head(samdf)
# NB: the dataframe includes a fungal mock community represented as 'EQ' under column district.



######################
# [A] MYCOBIOME ####
######################

# NB: The saved 'rds' files for the mycobiome were generated raw through the Divisive Amplicon Denoising Algorithym (DADA) V2.
# Run code for DADA 2 may be found at external link https://benjjneb.github.io/dada2/tutorial_1_8.html

# Reopen our saved rds files:
# Table of non-chimeric sequences:
seqtab.nochim <- readRDS("https://github.com/bkatati/mmycobiome/seqtab.rds")
head(seqtab.nochim)
samples.out <- rownames(seqtab.nochim)
head(samples.out)
# Next, you may use "View()" or "print()" to view the rds contents

# rds of assigned taxa to the seqtab.nochim ASVs, assignment done in DADA2 using the Unite Fungal Database
taxa <- readRDS("https://github.com/bkatati/mmycobiome/taxa_Unite2019.rds")

# Generation of phyloseq object of samples ID (rows) x ASVs (Columns)
# ( NB: an already generated phyloseq object is available here:)
# ps <- readRDS("C:/Users/katat001/OneDrive - Wageningen University & Research/ACADEMIC/Thesis/Publications/Github/mmycobiome/ps.rds")

rownames(samdf) <- samples.out
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# remove mock community:
remove <- c("S3EQA05", "S3EQA15", "S3EQA50", "S3EQF50", "S3EQF85", "S3EQF95", "S3EQX32", "S3EQY32")
pr <- prune_samples(!(sample_names(ps) %in% remove), ps)
pr

# Generate genera abundances based on the new phyloseq object 'pr':
library(plyr)
psg0 <- tax_glom(pr, "Genus")
psg1 <- transform_sample_counts(psg0, function(OTU)100* OTU / sum(OTU))
psg2 <- merge_samples(psg1, "Field")
psg3 <- transform_sample_counts(psg2, function(OTU)100* OTU / sum(OTU))
psg4 <- psg3
otu_table(psg4) <- t(otu_table(psg4))
OTUg <- otu_table(psg4)
TAXg <- tax_table(psg4)[, "Genus"]
mycobiome <- merge(TAXg, OTUg, by=0, all = TRUE)
mycobiome$Row.names = NULL
mycobiome$Mean=rowMeans(mycobiome[,-c(1)], na.rm=TRUE)
mycobiome <- mycobiome[order(desc(mycobiome$Mean)),]
head(mycobiome)



######################################################################
