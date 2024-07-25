#Script for June analysis
library("ShortRead")
library(dada2)
library(phyloseq)
library(vegan)
library("eulerr")
library("ggplot2")
library("ShortRead")
library(readr)
library(dplyr)
library(ggpubr)
library(scales)
library(grid)
library(reshape2)
library(DECIPHER)
library(phangorn)
library(seqinr)
library(MicEco) #used and mentioned
library(microbiome) #used
library(microbiomeutilities) #used
library(BioStrings) #used


#STARTS HERE

#Sample metadata update
metadata_final <- read.table("~/metadata_final.txt" , header = TRUE, fill= TRUE) 
rownames(metadata_final) <- metadata_final$SampleID
head(metadata_final)
phyloseq::sample_data(metadata_final)
dim(metadata_final)
sample_data(my_physeq) <- metadata_final
save(my_physeq, file="~/my_physeq_final.Rdata")


#Subsetting samples to keep samples from Tidepool 1.
my_physeq.int <- subset_samples(my_physeq,Tidepool%in%("One"))
my_physeq.int
plot_bar(my_physeq.int , fill= "Genus")

#Preliminary data representations before filtering
# Calculate the total number of sequences
total_seqs <- sum (taxa_sums(phyASV.f.c.int_tsn))
print(total_seqs)

#612244 before filtering, 2700 seqs after filtering by Genus

rarecurve(my_physeq, step=50, cex=0.5)
#and
tab <- otu_table(my_physeq)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- t(tab) # transpose observations to rows
library(vegan)
rare <- rarecurve(tab, step=10000, lwd=2, ylab="OTU",  label=F)
rarecurve(tab, step=50, cex=0.5) 

plot_heatmap(my_physeq, taxa.label = "Genus")
plot_heatmap(my_physeq, taxa.label = "Genus")

#Creating filtered file containing only Genera level assignments
gen_keep=c("g__cladeH", "g__Cladocopium","g__Durusdinium", "g__Fugacium", "g__Gerakladium", "g__cladeI" )
my_physeq.int.filt <- subset_taxa(my_physeq.int, Genus%in%gen_keep)

#Adding tree based on current ITS2 alignments and nr28s-rDNA distances: refer to "treee_nr28s.r" script
phy_tree(my_physeq.int.filt) <- phy_tree(uber.tree)

#Visualizing diversity

# We rename this phyloseq object to shorten it for use
z.int <- my_physeq.int.filt 


#Venn diagrams
 
sample_data(my_physeq.int)<- metadata_final
ps_venn(my_physeq.int,
  "Interaction",
  fraction = 0,
  weight = FALSE,
  relative = FALSE,
  plot = TRUE,
 )


#Alpha Diversity: Check script "alpha_intseas.r"
alpha_diversity <- estimate_richness(my_physeq.int, measures = c( "Observed","Chao1","Shannon", "InvSimpson"))
print(alpha_diversity)


plot_richness(my_physeq.int, x="Season", measures=c("Observed","Chao1" ,"Shannon", "InvSimpson")) + geom_boxplot()


sample_data_df <- as.data.frame(sample_data(my_physeq.int))
# Combine alpha diversity metrics with sample data
alpha_diversity_df <- cbind(sample_data_df, alpha_diversity)
group_variable1 <- alpha_diversity_df$Interaction
group_variable2 <- alpha_diversity_df$Season

# Perform pairwise Wilcoxon test for Shannon diversity
#Interaction
pairwise_wilcox_obs <- pairwise.wilcox.test(alpha_diversity_df$Observed, group_variable1, p.adjust.method = "BH")
print(pairwise_wilcox_obs)
pairwise_wilcox_chao1 <- pairwise.wilcox.test(alpha_diversity_df$Chao1, group_variable1, p.adjust.method = "BH")
print(pairwise_wilcox_chao1)
pairwise_wilcox_shannon <- pairwise.wilcox.test(alpha_diversity_df$Shannon, group_variable1, p.adjust.method = "BH")
print(pairwise_wilcox_shannon)
pairwise_wilcox_invsimp <- pairwise.wilcox.test(alpha_diversity_df$InvSimpson, group_variable1, p.adjust.method = "BH")
print(pairwise_wilcox_invsimp)
#Season
pairwise_wilcox2_obs <- pairwise.wilcox.test(alpha_diversity_df$Observed, group_variable2, p.adjust.method = "BH")
print(pairwise_wilcox2_obs)
pairwise_wilcox2_chao1 <- pairwise.wilcox.test(alpha_diversity_df$Chao1, group_variable2, p.adjust.method = "BH")
print(pairwise_wilcox2_chao1)
pairwise_wilcox2_shannon <- pairwise.wilcox.test(alpha_diversity_df$Shannon, group_variable2, p.adjust.method = "BH")
print(pairwise_wilcox2_shannon)
pairwise_wilcox2_invsimp <- pairwise.wilcox.test(alpha_diversity_df$InvSimpson, group_variable2, p.adjust.method = "BH")
print(pairwise_wilcox2_invsimp)


#Beta Diversity
#transformation of counts 
save(z.int, file="~/z.int.Rdata")
z.int_tsn <- transform_sample_counts(z.int, function(ASV) ASV/sum(ASV) *100 )


#Barplot 
plot_bar(z.int_tsn, x = "SampleID", fill = "Genus") +
  facet_wrap(~Interaction, scales = "free_x", nrow = 1) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 15) 
  )

#betadiversity script according to from https://scienceparkstudygroup.github.io/microbiome-lesson/06-beta-diversity/index.html
#https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
#Claar et. al (2020)
library(gridExtra)

theme_classic()
pcoa_zint = ordinate(z.int_tsn, "PCoA", "wunifrac")
pcoa_zintbray = ordinate(z.int_tsn, "PCoA", "bray")

pcoa_int <- plot_ordination(z.int_tsn, pcoa_zint , type="samples", color = "Interaction") 
themed_pcoa1 <- pcoa_int + geom_point(size = 2) + geom_polygon(aes(fill=Interaction ), alpha=0.2)  
themed_pcoa1 + theme_classic() + theme(    legend.position = c(0.90, 0.90), # Move legend to the top right inside the plot
                                           legend.justification = c("right", "top"),
  legend.text = element_text(size = rel(2)), # Increase legend text size by 1.5
  legend.title = element_text(size = rel(2)) # Increase legend title size by 1.5
)
pcoa_seas <- plot_ordination(z.int_tsn, pcoa_zint , type="samples", color = "Season") 
themed_pcoa2 <- pcoa_seas + geom_point(size = 2) + geom_polygon(aes(fill=Season ), alpha=0.2) 
themed_pcoa2 + theme_classic() + theme( legend.position = c(0.90, 0.90), # Move legend to the top right inside the plot
                                        legend.justification = c("right", "top"),
  legend.text = element_text(size = rel(2)), # Increase legend text size by 1.5
  legend.title = element_text(size = rel(2)) # Increase legend title size by 1.5
)

NMDS_zint = ordinate(z.int_tsn, "NMDS", "wunifrac")
NMDS_int <- plot_ordination(z.int_tsn, NMDS_zint , type="samples", color = "Interaction") 
themed_NMDS1 <- NMDS_int + geom_point(size = 2) + geom_polygon(aes(fill=Interaction ), alpha=0.2)  
themed_NMDS1 + theme_classic() + theme(    legend.position = c(0.90, 0.90), # Move legend to the top right inside the plot
                                           legend.justification = c("right", "top"),
                                           legend.text = element_text(size = rel(2)), # Increase legend text size by 1.5
                                           legend.title = element_text(size = rel(2)) # Increase legend title size by 1.5
)
NMDS_seas <- plot_ordination(z.int_tsn, NMDS_zint , type="samples", color = "Season") 
themed_NMDS2 <- NMDS_seas + geom_point(size = 2) + geom_polygon(aes(fill=Season ), alpha=0.2) 
themed_NMDS2 + theme_classic() + theme( legend.position = c(0.90, 0.90), # Move legend to the top right inside the plot
                                        legend.justification = c("right", "top"),
                                        legend.text = element_text(size = rel(2)), # Increase legend text size by 1.5
                                        legend.title = element_text(size = rel(2)) # Increase legend title size by 1.5
)






#tests permanova from https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
set.seed(1)
# Calculate distance matrix
perm.int_bc <- phyloseq::distance(z.int, method = "bray")
perm.int_wu <- phyloseq::distance(z.int, method = "wunifrac")
# make a data frame from the sample_data
sampledzint <- data.frame(sample_data(z.int))

# Adonis test
adonis2(perm.int_wu ~ Interaction, data = sampledzint)
adonis2(perm.int_wu ~ Season, data = sampledzint)

#betadisper PERMADISP for homogeneity of dispersion
beta.int <- betadisper(perm.int_wu, sampledzint$Interaction)
permutest(beta.int, permutations = 999, pairwise= TRUE)
beta.seas <- betadisper(perm.int_wu, sampledzint$Season)
permutest(beta.seas, permutations = 999, pairwise= TRUE)

#pairwise PERMANOVA

set.seed(99)
metadata <- data.frame(sample_data(z.int))
cbn <- combn(x = unique(metadata$Season), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(z.int, Season %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = "wunifrac") ~ Season, 
                                data = metadata_sub)
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}

p.adj1 <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p = p, p.adj1 = p.adj1)
p.table

#----------
  
#ANOSIM
anosimdist = phyloseq::distance(z.int, method="bray")
anosimdist
metaanosim <- data.frame(sample_data(z.int))
anosim(anosimdist, metaanosim$Interaction)
anosim(anosimdist, metaanosim$Season)


#Simper
## Simper analysis, extract abundance matrix from a phyloseq object
zz_OTUs = as(otu_table(z.int), "matrix")
# transpose so we have the OTUs as columns
if(taxa_are_rows(z.int)){zz_OTUs <- t(zz_OTUs)}
# Coerce the object to a data.frame
zz_OTUs_scaled = as.data.frame(zz_OTUs)
# running the simper analysis on the dataframe and the variable of interest "time"
zz_simper <- simper(zz_OTUs_scaled, metadata_27_samples$Interaction, permutations = 999)
# printing the top OTUs
print(zz_simper)
summary(zz_simper)



