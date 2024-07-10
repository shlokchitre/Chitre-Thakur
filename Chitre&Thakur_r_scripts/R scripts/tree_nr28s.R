#This script ws used for phylogenetic tree construction following the R cript included in supplementary data provided with Claar et. al (2020)

library(readr)
library("ggplot2")
library("ShortRead")
library(stringr)
library(reshape2)
library(phyloseq)
library(seqinr)
library(phangorn) 
library(caroline)
library(DECIPHER)
library("ape")


#File in use : Filtered file by genus, 27 samples
my_physeq.int.filt
phyASV.f.c <- my_physeq.int.filt #for purpose of this script

# Write tax table for phylogenetically-informed diversity
write.table(data.frame(tax_table(phyASV.f.c)), "~/tree/taxtable_fortree.txt", row.names=T, quote=F)
#Write otu table for phylogenetically-informed diversity
write.delim(data.frame(otu_table(phyASV.f.c)), "~/tree/otu_table.tsv", quote = FALSE, row.names = T, sep = "\t")

#Lets segregate

# ASV.As <- subset_taxa(phyASV.f.c,Genus=="g__Symbiodinium") no seqs in data
ASV.Bs <- subset_taxa(phyASV.f.c,Genus=="g__Breviolum")
ASV.Cs <- subset_taxa(phyASV.f.c,Genus=="g__Cladocopium")
ASV.Ds <- subset_taxa(phyASV.f.c,Genus=="g__Durusdinium")
# ASV.Es <- subset_taxa(phyASV.f.c,Genus=="g__Effrenium") no seqs in data
ASV.Fs <- subset_taxa(phyASV.f.c,Genus=="g__Fugacium")
ASV.Gs <- subset_taxa(phyASV.f.c,Genus=="g__Gerakladium")
ASV.Hs <- subset_taxa(phyASV.f.c,Genus=="g__cladeH")
ASV.Is <- subset_taxa(phyASV.f.c,Genus=="g__cladeI")

ASV.seqs <- refseq(phyASV.f.c)
min(width(ASV.seqs))
max(width(ASV.seqs))


ASV.B.seqs <- refseq(ASV.Bs)
writeXStringSet(ASV.B.seqs, # write to a new FASTA file
                file="~/ASV_B_tree_seqs.fasta")
ASV.B.seqs.aligned <- refseq(ASV.Bs) # There is only one sequence here! Can't align
writeXStringSet(ASV.B.seqs.aligned, # write the alignment to a new FASTA file
                file="~/tree/ASV_B_tree_seqs_aligned.fasta")

ASV.C.seqs <- refseq(ASV.Cs)
writeXStringSet(ASV.C.seqs, # write to a new FASTA file
                file="~/tree/ASV_C_tree_seqs.fasta")
ASV.C.seqs.aligned <- AlignSeqs(ASV.C.seqs)
writeXStringSet(ASV.C.seqs.aligned, # write the alignment to a new FASTA file
                file="~/tree/ASV_C_tree_seqs_aligned.fasta")

ASV.D.seqs <- refseq(ASV.Ds)
writeXStringSet(ASV.D.seqs, # write to a new FASTA file
                file="~/tree/ASV_D_tree_seqs.fasta")
ASV.D.seqs.aligned <- AlignSeqs(ASV.D.seqs)
writeXStringSet(ASV.D.seqs.aligned, # write the alignment to a new FASTA file
                file="~/tree/ASV_D_tree_seqs_aligned.fasta")

ASV.F.seqs <- refseq(ASV.Fs)
writeXStringSet(ASV.F.seqs, # write to a new FASTA file
                file="~/tree/ASV_F_tree_seqs.fasta")
ASV.F.seqs.aligned <- AlignSeqs(ASV.F.seqs)
writeXStringSet(ASV.F.seqs.aligned, # write the alignment to a new FASTA file
                file="~/tree/ASV_F_tree_seqs_aligned.fasta")

ASV.G.seqs <- refseq(ASV.Gs)
writeXStringSet(ASV.G.seqs, # write to a new FASTA file
                file="~/tree/ASV_G_tree_seqs.fasta")
ASV.G.seqs.aligned <- AlignSeqs(ASV.G.seqs)
writeXStringSet(ASV.G.seqs.aligned, # write the alignment to a new FASTA file
                file="~/tree/ASV_G_tree_seqs_aligned.fasta")

#for H
ASV.H.seqs <- refseq(ASV.Hs)
writeXStringSet(ASV.H.seqs, # write to a new FASTA file
                file="~/tree/ASV_H_tree_seqs.fasta")
ASV.H.seqs.aligned <- AlignSeqs(ASV.H.seqs)
writeXStringSet(ASV.H.seqs.aligned, # write the alignment to a new FASTA file
                file="~/tree/ASV_H_tree_seqs_aligned.fasta")

ASV.I.seqs <- refseq(ASV.Is)
writeXStringSet(ASV.I.seqs, # write to a new FASTA file
                file="~/tree/ASV_I_tree_seqs.fasta")
#ASV.I.seqs.aligned <- AlignSeqs(ASV.I.seqs) #only one seq
ASV.I.seqs.aligned <- refseq(ASV.Is)
writeXStringSet(ASV.I.seqs.aligned, # write the alignment to a new FASTA file
                file="~/tree/ASV_I_tree_seqs_aligned.fasta")

#Save it all
writeXStringSet(c(ASV.B.seqs.aligned,ASV.C.seqs.aligned,ASV.D.seqs.aligned,ASV.F.seqs.aligned,ASV.G.seqs.aligned,ASV.H.seqs.aligned,ASV.I.seqs.aligned), 
                file="~/tree/ASV_ALL_tree_seqs_aligned.fasta")
#all alignment data gathered, moving towards distance matrix based data

#https://rdrr.io/rforge/seqinr/man/dist.alignment.html
#returns sqrt of pairwise genetic distance, then squared the matrices

B.seqs <- read.alignment(file = "~/tree/ASV_B_tree_seqs_aligned.fasta", format= "fasta")
B.dis <- (as.matrix(dist.alignment(B.seqs, matrix = "identity" )))^2
write.csv(B.dis, file="~/tree/ASV_B_dis_matx.csv")

C.seqs <- read.alignment(file = "D:/Metagenomics work/CSS norms/tree/ASV_C_tree_seqs_aligned.fasta", format= "fasta")
C.dis <- (as.matrix(dist.alignment(C.seqs, matrix = "identity" )))^2
write.csv(C.dis, file="D:/~/tree/ASV_C_dis_matx.csv")

D.seqs <- read.alignment(file = "D:/~/tree/ASV_D_tree_seqs_aligned.fasta", format= "fasta")
D.dis <- (as.matrix(dist.alignment(D.seqs, matrix = "identity" )))^2
write.csv(D.dis, file="D:/~/tree/ASV_D_dis_matx.csv")

F.seqs <- read.alignment(file = "D:/Metagenomics work/CSS norms/tree/ASV_F_tree_seqs_aligned.fasta", format= "fasta")
F.dis <- (as.matrix(dist.alignment(F.seqs, matrix = "identity" )))^2
write.csv(F.dis, file="~/tree/ASV_F_dis_matx.csv")

G.seqs <- read.alignment(file = "D:/Metagenomics work/CSS norms/tree/ASV_G_tree_seqs_aligned.fasta", format= "fasta")
G.dis <- (as.matrix(dist.alignment(G.seqs, matrix = "identity" )))^2
write.csv(G.dis, file="D:/~/tree/ASV_G_dis_matx.csv")

H.seqs <- read.alignment(file = "D:/Metagenomics work/CSS norms/tree/ASV_H_tree_seqs_aligned.fasta", format= "fasta")
H.dis <- (as.matrix(dist.alignment(H.seqs, matrix = "identity" )))^2
write.csv(H.dis, file="D:/~/tree/ASV_H_dis_matx.csv")

I.seqs <- read.alignment(file = "D:/Metagenomics work/CSS norms/tree/ASV_I_tree_seqs_aligned.fasta", format= "fasta")
I.dis <- (as.matrix(dist.alignment(I.seqs, matrix = "identity" )))^2
write.csv(I.dis, file="D:/~/tree/ASV_I_dis_matx.csv")          

#give clade distances using average 28s distance from Pochon and Gates 2010

B_C <- matrix(0.114, ncol=ncol(B.dis), nrow=nrow(C.dis), dimnames=list(rownames(C.dis), colnames(B.dis)))
B_D <- matrix(0.1705, ncol=ncol(B.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(B.dis)))
B_F <- matrix(0.1355, ncol=ncol(B.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(B.dis)))
B_G <- matrix(0.21, ncol=ncol(B.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(B.dis)))
B_H <- matrix(0.112, ncol=ncol(B.dis), nrow=nrow(H.dis), dimnames=list(rownames(H.dis), colnames(B.dis)))
B_I <- matrix(0.151, ncol=ncol(B.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(B.dis)))
C_D <- matrix(0.1520, ncol=ncol(C.dis), nrow=nrow(D.dis), dimnames=list(rownames(D.dis), colnames(C.dis)))
C_F <- matrix(0.0945, ncol=ncol(C.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(C.dis)))
C_G <- matrix(0.187, ncol=ncol(C.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(C.dis)))
C_H <- matrix(0.041, ncol=ncol(C.dis), nrow=nrow(H.dis), dimnames=list(rownames(H.dis), colnames(C.dis)))
C_I <- matrix(0.137, ncol=ncol(C.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(C.dis)))
D_F <- matrix(0.1691, ncol=ncol(D.dis), nrow=nrow(F.dis), dimnames=list(rownames(F.dis), colnames(D.dis)))
D_G <- matrix(0.1795, ncol=ncol(D.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(D.dis)))
D_H <- matrix(0.144, ncol=ncol(D.dis), nrow=nrow(H.dis), dimnames=list(rownames(H.dis), colnames(D.dis)))
D_I <- matrix(0.169125, ncol=ncol(D.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(D.dis)))
#for F, first average horizontal, add vertical, divide by 4 , means 0.1785+0.0.1685+0.160+0.1695=0.1691
F_G <- matrix(0.2072, ncol=ncol(F.dis), nrow=nrow(G.dis), dimnames=list(rownames(G.dis), colnames(F.dis)))
F_H <- matrix(0.09125, ncol=ncol(F.dis), nrow=nrow(H.dis), dimnames=list(rownames(H.dis), colnames(F.dis)))
F_I <- matrix(0.14925, ncol=ncol(F.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(F.dis)))
G_H <- matrix(0.176, ncol=ncol(G.dis), nrow=nrow(H.dis), dimnames=list(rownames(H.dis), colnames(G.dis)))
G_I <- matrix(0.194, ncol=ncol(G.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(G.dis)))
H_I <- matrix(0.130, ncol=ncol(H.dis), nrow=nrow(I.dis), dimnames=list(rownames(I.dis), colnames(H.dis)))

#build ACDG matrix
col1 <- rbind(B.dis, B_C, B_D, B_F, B_G, B_H, B_I)
col2 <- rbind(matrix(NA, nrow=nrow(B.dis), ncol=ncol(C.dis), dimnames=list(rownames(B.dis), colnames(C.dis))), C.dis, C_D, C_F, C_G, C_H, C_I)
col3 <- rbind(matrix(NA, nrow=nrow(B.dis)+nrow(C.dis), ncol=ncol(D.dis), dimnames=list(c(rownames(B.dis), rownames(C.dis)), colnames(D.dis))), D.dis, D_F, D_G, D_H, D_I)
col4 <- rbind(matrix(NA, nrow=nrow(B.dis)+nrow(C.dis)+nrow(D.dis), ncol=ncol(F.dis), dimnames=list(c(rownames(B.dis), rownames(C.dis), rownames(D.dis)), colnames(F.dis))), F.dis, F_G, F_H, F_I)
col5 <- rbind(matrix(NA, nrow=nrow(B.dis)+nrow(C.dis)+nrow(D.dis)+nrow(F.dis), ncol=ncol(G.dis), dimnames=list(c(rownames(B.dis), rownames(C.dis), rownames(D.dis), rownames(F.dis)), colnames(G.dis))), G.dis, G_H, G_I)
col6 <- rbind(matrix(NA, nrow=nrow(B.dis)+nrow(C.dis)+nrow(D.dis)+nrow(F.dis)+nrow(G.dis), ncol=ncol(H.dis), dimnames=list(c(rownames(B.dis), rownames(C.dis), rownames(D.dis), rownames(F.dis),  rownames(G.dis)), colnames(H.dis))), H.dis, H_I)
col7 <- rbind(matrix(NA, nrow=nrow(B.dis)+nrow(C.dis)+nrow(D.dis)+nrow(F.dis)+nrow(G.dis)+nrow(H.dis), ncol=ncol(I.dis), dimnames=list(c(rownames(B.dis), rownames(C.dis), rownames(D.dis), rownames(F.dis),  rownames(G.dis), rownames(H.dis)), colnames(I.dis))), I.dis)

#HUSH
ubermatrix <- cbind(col1, col2, col3, col4, col5, col6, col7)
dim(ubermatrix)

#build tree
uber.tree <- phangorn::upgma(ubermatrix)
plot(uber.tree, main="UPGMA")

#write tree to file
write.tree(uber.tree, file="~/tree/uber.tre")

class(uber.tree)
# Join to phyloseq object
phy_tree(my_physeq.int.filt) <- phy_tree(uber.tree)
plot_tree(my_physeq.int.filt, label.tips = "Genus", color = "Genus")

