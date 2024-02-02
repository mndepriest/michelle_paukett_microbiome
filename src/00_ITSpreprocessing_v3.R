##ITS RUN TEST

# Notes and resources -----------------------------------------------------


# R version 4.2.2 on 4.11.23

# Michelle notes that the seq center removed primers, so only tags (defined below)
# remain. 

# Dada2 ITS pipeline
# browse("https://benjjneb.github.io/dada2/ITS_workflow.html")

# Manage packages ------------------------------------------------------

# This just makes it easier to see that all of your packages were loaded successfully:

# Define all necessary packages
list.of.packages <-
    c(library(devtools),
    library(BiocManager),
    library(dada2),
    library(Biostrings),
    library(ShortRead),
    library(igraph),
    library(phyloseq),
    library(dplyr),
    library(ggplot2),
    library(Rcpp),
    library(here),
    library(patchwork))

# Install new packages if necessary
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]

if (length(new.packages))
  install.packages(new.packages)

# Load all packages
packages_load <-
  lapply(list.of.packages, require, character.only = TRUE)

# Check that all packages were loaded successfully
if (any(as.numeric(packages_load) == 0)) {
  warning(paste("Package/s: ", paste(list.of.packages[packages_load != TRUE], sep = ", "), "not loaded!"))
} else {
  print("All packages were successfully loaded.")
}

# Clear global environment
rm(list.of.packages, new.packages, packages_load)

# Managing workplace, path, etc -------------------------------------------


# Set home directory/path location and build/organize your R project files
# Using 'here()' instead of 'setwd()' makes the code more easily reproducible 
# on other computers. 
here("C:/Users/HumanGarbage/OneDrive - The Pennsylvania State University/Esker Lab/Michelle_ITS")
basepath <- "C:/Users/HumanGarbage/OneDrive - The Pennsylvania State University/Esker Lab/Michelle_ITS/"

# Make a list of all folders that need to exist in "Michelle_ITS"
folders <- c("data", "source", "out", "temp")
data.folders <- c("fastqGZ", "filtered", "derep", "dada", "merged", "seqtabs", "nochim", 
                  "aligned")

for (i in folders){
  name <- as.character(i)
  result <- paste(basepath, as.character(i), sep = "")
  assign(name, result)
}

for (i in data.folders){
  name <- as.character(i)
  result <- paste(as.character(data), "/", as.character(i), sep = "")
  assign(name, result)
}

# Define primers to determine presence later
FWD <- "AGCCTCCGCTTATTGATATGCTTAART"
REV <- "AACTTTYRRCAAYGGATCWCT"

# Looking at fastqs --------------------------------------------------------------


# Create a list of all forward sample paths in 'fnFs'
fnFs <- sort(list.files(fastqGZ, pattern = "L001_R1_001.fastq.gz", full.names = TRUE))
readFastq(fnFs[1])

# Create a list of all reverse sample paths in 'fnRs'
fnRs <- sort(list.files(fastqGZ, pattern = "_L001_R2_001.fastq.gz", full.names = TRUE))
readFastq(fnRs[1])


# Checking, sequencing quality
# plotQualityProfile(fnFs[1])
# plotQualityProfile(fnRs[1])
# fnFs are good until around 225, fnRs are SIGNIFIC,ANTLY worse, dropping at 150

# Make base sample names
sample.names <- sapply(strsplit(basename(fnFs), "_L001"), `[`,1)


# Filter for ambiguity in bp ID ####

# Create locations for samples after primer filtering is performed
filtFs <- file.path(data, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(data, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Filter ambiguous (and poor quality) reads from set
# Aggressive trimming makes it possible to merge later. 
output <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(250, 100),
                        maxN = c(1,1), maxEE = c(2,2),
                        compress = TRUE, multithread = FALSE)
output

# Check filtered quality
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])

# Check for stray primers ####

# List all primer seq variations
allOrients <- function(primer){
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}

# Get all possible FWD and REV orientations
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Function that will show presence of all possible orientations and their complements
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Finding # of primer instances in EACH sample
for (i in 1:11){
print(rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[i]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[i]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[i]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[i]])))
}
# Some instances of primer occurrence WITHOUT truncation, but still low. With truncation
# there are no instances of primer occurrence, indicating that ambig bases in the 
# latter part of the rev reads are the source

# Load in filtered files if already processed: (limits extending session temp memory)
# Create a list of all forward sample paths in 'fnFs'
# filtFs <- sort(list.files(filtered, pattern = "F_filt.fastq.gz", full.names = TRUE))
# Create a list of all reverse sample paths in 'fnRs'
# filtRs <- sort(list.files(filtered, pattern = "R_filt.fastq.gz", full.names = TRUE))
# sample.names <- sapply(strsplit(basename(filtFs), "_F_filt.fastq.gz"), `[`,1)


# Error rates ####


# Learning error rates from 990446 bases in 3946 reads in 11 samples
errF <- learnErrors(filtFs, multithread = FALSE)
saveRDS(errF, paste0(out, "/", "errF.RDS"))
errR <- learnErrors(filtRs, multithread = FALSE)
saveRDS(errR, paste0(out, "/", "errR.RDS"))
# The error "Not all seqs are the same length." did NOT occur and it should have
# bc this is ITS, not 16S

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Dereplicate ####

# Create OTU counts for smaller file size
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# derepFs looks normal, but virtually EVERY read in derepRs is a 'unique' read
# Relaxing and restricting maxN and maxEE have virtually no effect, but truncating
# the poor quality parts helps A LOT

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Estimate error rate on dereplicated unique seqs ####

dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)
# Based on this feedback, we COULD get as many successful unique merges as there are
# unique seqs in dadaFs but we rarely do

# Merge #### 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Seqtab ####

# Create sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove chims
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# With trimming at c(250,100), 56 seqs out of 224 are bimeric

# Save as .rds 
saveRDS(seqtab.nochim, paste0(out, "/seqtab250_100.rds"))

# Look at the distribution of seq length. Can be anywhere from 250 (perfect overlap)
# to 338 (longest seq with min overlap)
table(nchar(getSequences(seqtab.nochim)))
# Median at 324, minimum overlap possible at 12bp

# View how many reads we had at each stage
getN <- function(x) sum(getUniques(x))
track <- cbind(output, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# Track output thru processing
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)
# Pretty poor output

# Using UNITE + euks 2022 ####
unite.ref <- paste0(data, "/uniteITS.fasta")
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)


# Save tax file as RDS
saveRDS(taxa, paste0(out, "/taxa250_100.rds"))


