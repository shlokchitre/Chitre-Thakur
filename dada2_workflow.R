#DADA2 pipeline , modified as per script from github repository: Claar et. al 2020


library("ShortRead")
library(dada2)
library(phyloseq)
library("BioStrings")
library("ggplot2")
library(vegan)
library(readr)

#originated from claar dada_2.r
# https://benjjneb.github.io/dada2/ITS_workflow.html
#https://cutadapt.readthedocs.io/en/stable/installation.html for installation information
cutadapt <- "/usr/bin/cutadapt" # Where is cutadapt located on your machine
system2(cutadapt, args = "--version") # Check that cutadapt works



# Create file lists
path <- "~/sequencing data" 
fnFs <- sort(list.files(path, pattern="_R1.fastq"))	
fnRs <- sort(list.files(path, pattern="_R2.fastq"))	

fnFs <- file.path(path, fnFs) # Append the filepath to the filenames - forward reads
fnRs <- file.path(path, fnRs) # Append the filepath to the filenames - reverse reads
fnFs
fnRs

# Identify primers used for this project # These are  the ITS2 primers used in current study
FWD <- "GAATTGCAGAACTCCGTGAACC"  
REV <- "CGGGTTCWCTTGTYTGACTTCATGC"  

# Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Pre-filter the sequences just to remove those with Ns, but perform no other filtering for now
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, # dada does not allow Ns, so must be zero
              multithread = TRUE, 
              compress = FALSE) # Turn off compression (default is TRUE), cutadapt doesn't work on compressed files

# Count the number of times the primers appear in the forward and reverse read, while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, ShortRead::sread(ShortRead::readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# We did not find any primers and adapters as sequencing company submitted cleared demultiplexed data.


#Workflow
#Data checked for specific ITS primers used, none found.

#SeTTING PATH
path <- "~/ITSfiles"  # Path to where sequences are stored
fnFs <- sort(list.files(path, pattern="_R1.fastq"))	
fnRs <- sort(list.files(path, pattern="_R2.fastq"))	

fnFs <- file.path(path, fnFs) # Append the filepath to the filenames - forward reads
fnRs <- file.path(path, fnRs) # Append the filepath to the filenames - reverse reads
fnFs
fnRs

#FILTER AND TRIM
# Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     maxN = 0, # Standard filtering parameter
                     maxEE = c(2, 4), # Sets the maximum number of âexpected errorsâ allowed in a read
                     trimRight = c(5,60), #Poor quality reverse reads
                     trimLeft = c(1,2),
                     # The number of nucleotides to remove from the end of each read - chosen based on declining quality scors in plot
                     minLen = 170, # Enforced a minLen here, to get rid of spurious short sequences; was 50 before
                     rm.phix = TRUE, # Removing any phiX sequences
                     compress = TRUE, 
                     multithread = FALSE, # on windows, set multithread = FALSE
                     verbose = TRUE) 
head(out)

#Error model generation # Learn error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
# Visualize the estimated error rates as a sanity check
plotErrors(errF, nominalQ = TRUE)

#Dereplication of samples.
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)


# DADA2 Flagship Sample Inference - apply the core sample inference algorithm to the dereplicated data.
dadaFs <- dada(derepFs, err = errF, pool= "pseudo", multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, pool = "pseudo", multithread = TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                      trimOverhang=TRUE, verbose=TRUE)

# Construct amplicon sequence variant table (ASV) table for FORWARD ONLY
#edited for Forward only
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Remove chimeras FORWARD ONLY
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Inspect distribution of sequence lengths:
table(nchar(getSequences(seqtab.nochim)))

#Disribution of Seq length
table(nchar(getSequences(seqtabM.nochim)))
hist(nchar(getSequences(seqtabM.nochim)), main = "Distribution of SEQ lengths")


# Function to get the sum of unique sequences 
getN <- function(x) sum(getUniques(x))

# Combine the output from various steps into a tracking matrix
track <- cbind(out, 
               filtered = sapply(dadaFs, getN), 
               denoisedF = sapply(dadaRs, getN), 
               merged = sapply(mergers, getN), 
               nonchim = rowSums(seqtab.nochim))

# Assign column names to the tracking matrix
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
print(track)



# Track reads through the pipeline - R wizard pending
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

SYM_Refs
sym.ref1 <- "~/ITS2db_claar.fasta" 
sym.ref2 <- "~/SymportalrefSeqDB.fasta" 
sym.ref3 <- "~/AllITS2dbn1_Dinophy_Symbio.fasta" 

# Assign taxonomy RDP classifier
taxaM <- assignTaxonomy(seqtabM.nochim, sym.ref, multithread = TRUE, tryRC = TRUE, minBoot = 50, verbose = TRUE)
taxaM.print <- taxaM  # Removing sequence rownames for display only
rownames(taxaM.print) <- NULL
head(taxaM.print)
tpM <- data.frame(taxaM.print)
write.table(taxaM,file="~/taxatable.txt",sep="\t", quote=F, col.names = NA)

# Import to phyloseq

samdf <- read.table("/metadatax.txt" , header = TRUE, fill= TRUE) # Read in sample data
rownames(samdf) <- samdf$SampleID
head(samdf)

# We now construct a phyloseq object directly from the dada2 outputs.
my_physeq <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                            sample_data(samdf), 
                            tax_table(taxa))

#Make refseq for phyloseq
my_physeq <- add_refseq(my_physeq, tag = "ASV")

library(Biostrings)
dna <- Biostrings::DNAStringSet(taxa_names(my_physeq))
names(dna) <- taxa_names(my_physeq)
my_physeq <- merge_phyloseq(my_physeq, dna)
taxa_names(my_physeq) <- paste0("ASV", seq(ntaxa(my_physeq)))
my_physeq

save(my_physeq, file="~/Phyloseq_Chitre.RDdata")


refseqslist <- refseq(my_physeq, errorIfNULL=TRUE)
write.table(refsasvs ,file="~/taxanames.txt",sep="\t", quote=F, col.names = NA)
taxa_sums1 <- taxa_sums1(my_physeq)
write.table(taxa_sums1 ,file="~/taxasums.txt",sep="\t", quote=F, col.names = NA)

