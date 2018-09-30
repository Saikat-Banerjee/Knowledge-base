#!/usr/bin/env/Rscript
#Steps for microbiome data analysis from raw reads to community analysis
#description of the steps: https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#different_ordination_projections
#check for dependencies
library("knitr")
library("BiocStyle")
.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}
# check packages
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
set.seed(100)
#complementary analysis

.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                     "reshape2", "PMA", "structSSI", "ade4",
                     "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
  devtools::install_github(.github_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}

#CHANGE to the directory containing the fastq files.
#Note: unzip the fastq file before running these steps.
miseq_path <- choose.dir()
print(file.list<- list.files(miseq_path))
# Sort forward/reverse reads in same order
sorted.Fs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
sorted.Rs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names,format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(sorted.Fs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
sorted.Fs.path <- file.path(miseq_path, sorted.Fs)
sorted.Rs.path <- file.path(miseq_path, sorted.Rs)
#inspect Sorted files
print(inspect.FS <- head(sorted.Fs.path))
print(inspect.RS <- head(sorted.Rs.path))
#QC plot
print(QC.FS<- plotQualityProfile(sorted.Fs.path[1:2]))
print(QC.RS<- plotQualityProfile(sorted.Rs.path[1:2]))
# Place filtered files in filtered/ subdirectory
filt_path <- file.path(miseq_path, "filtered") 
if(!file_test("-d", filt_path)) dir.create(filt_path)
#define the filenames for the filtered fastq.gz files
filt.Fs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filt.Rs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
#Filter the forward and reverse reads:
filt.out <- filterAndTrim(sorted.Fs.path, filt.Fs, sorted.Rs.path, filt.Rs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=F)
print(inspect.filt.out <- head(filt.out))
#Dereplication to infer sequence variants (ASVs)
derep.Fs <- derepFastq(filt.Fs, verbose=TRUE)
derep.Rs <- derepFastq(filt.Rs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derep.Fs) <- sampleNames
names(derep.Rs) <- sampleNames
#learn errors
err.F <- learnErrors(filt.Fs, multithread=F)
err.R <- learnErrors(filt.Rs, multithread=F)
#plot errors
plotErrors(err.F,nominalQ = T)
plotErrors(err.R, nominalQ = T)
#unique sequences
dada.Fs <- dada(derep.Fs, err=err.F, multithread=F)
dada.Rs <- dada(derep.Rs, err = err.R, multithread = F)
#Inspecting the dada-class object
dada.Fs[[1]]
#Construct sequence table and remove chimeras
merge.pairs <- mergePairs(dada.Fs, derep.Fs, dada.Rs, derep.Rs)
seq.tab.All <- makeSequenceTable(merge.pairs[!grepl("Mock", names(merge.pairs))])
table(nchar(getSequences(seq.tab.All)))
seqtabNoC <- removeBimeraDenovo(seq.tab.All)

#Assign taxonomy
fastaRef.train <- choose.files()
taxa.Tab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef.train, multithread=F)
unname(head(taxa.Tab))
#Construct phylogenetic tree
seqs <- getSequences(seqtabNoC)
# propagate to the tip labels of the tree
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
#Combine data into a phyloseq object
#import meta file
samdf <- read.csv(choose.files(),header = T)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
# Remove dupicate entries for reverse reads
samdf <- samdf[!duplicated(samdf$SampleID),] 
# Fix discrepancy
rownames(seq.tab.All) <- gsub("124", "125", rownames(seq.tab.All))
all(rownames(seq.tab.All) %in% samdf$SampleID)
rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
               "host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
               "diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seq.tab.All), keep.cols]
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE),sample_data(samdf), tax_table(taxa.Tab),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
# Show available ranks
print(rank_names <- rank_names(ps))
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
#filter uncharacterized taxa
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
#Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Define phyla to filter
filterPhyla = c("Verrucomicrobia", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
# Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 1.5) +  geom_point(size = 1, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
# Define prevalence threshold as 5% of total samples
print (prevalenceThreshold <- (0.05 * nsamples(ps)))
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
h1= 0.4 
ps4 = tip_glom(ps2, h = h1)
#plot trees
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
#Abundance value transformation
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)
#Subset by taxonomy
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)

