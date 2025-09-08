## DADA2 pipeline was used for filtering and trimming sequences based on quality profiles, Dereplication by combining all identical sequencing reads into into “unique sequences”, merge paired reads, construct an asv table and remove chimeras. Taxonomy assigned with OBItools using ASV table from DADA2

## This pipeline was adapted from a tutorial by the creators of DADA2 (https://benjjneb.github.io/dada2/tutorial.html), we have used the standard default parameters suggested

## You need to run this twice, once for each set of sequences from MiSeq and NovaSeq

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.20")

library(tidyverse) # version 2.0
library(dada2)

###### Arthropod (NovaSeq)

path <- "arthropod_trimmed_24"           #wherever your results from cutadapt were stored
list.files(path) # make sure your path is loaded correctly

#read in the forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

## Extract sample names, assuming filenames have format
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1:2], collapse = "_"))

## check out sequence quality
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

## create output for the filtering process
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Standard filtering parameters
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2),truncQ=2, rm.phix = F,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE)
rownames(out)<-sapply(strsplit(basename(rownames(out)), "_"), function(x) paste(x[1:2], collapse = "_"))

## some files were empty so needed to remove them
remove_names <- rownames(out[out[, 2] == 0, ])
filtFs <- filtFs[!names(filtFs) %in% remove_names]
filtRs <- filtRs[!names(filtRs) %in% remove_names]

## estimating error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

## dereplicating
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names %>% 
  discard(sample.names %in% remove_names)
names(derepRs) <- sample.names %>% 
  discard(sample.names %in% remove_names)

## the core sample inference algorithm (removes sequencing errors)
dadaFs <- dada(derepFs , err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## merge the forward and reverse reads together to obtain the full denoised sequences. 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, trimOverhang = T)

##  construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Checking the frequency of chimeric sequences 
sum(seqtab.nochim)/sum(seqtab)
table(nchar(getSequences(seqtab.nochim)))

out2 <- as.data.frame(out)
out2 <- out2 %>% filter(!rownames(out2) %in% remove_names)

## Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)

#removing asv with less than 10 total hits
seqtab.nochim_filt<-seqtab.nochim %>% as.data.frame() %>%
  select_if(~sum(.) > 10)

dim(seqtab.nochim_filt)

## Making the output from DADA2 compatible with obitools
asv_seqs <- colnames(seqtab.nochim_filt)
asv_headers <- vector(dim(seqtab.nochim_filt)[2], mode="character")

for (i in 1:dim(seqtab.nochim_filt)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## Sequences for attaching taxonomy with obitools
asv_fasta <- c(rbind(asv_headers, asv_seqs))
writeLines(asv_fasta, "arthropod2024_dadaoutput.fasta")

## per sample sequence counts
asv_tab <- t(seqtab.nochim_filt)
row.names(asv_tab) <- sub(">", "", asv_headers)
head(asv_tab)
write.table(asv_tab, "arthropod_ASVs_counts_24.txt", sep="\t", quote=F)


###### Arthropod (MiSeq)

path <- "arthropod_trimmed_23"           #wherever your results from cutadapt were stored
list.files(path) # make sure your path is loaded correctly

#read in the forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

## Extract sample names, assuming filenames have format
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1:2], collapse = "_"))

## check out sequence quality
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

## create output for the filtering process
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Standard filtering parameters
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2),truncQ=2, rm.phix = F,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE)
rownames(out)<-sapply(strsplit(basename(rownames(out)), "_"), function(x) paste(x[1:2], collapse = "_"))

## some files were empty so needed to remove them
remove_names <- rownames(out[out[, 2] == 0, ])
filtFs <- filtFs[!names(filtFs) %in% remove_names]
filtRs <- filtRs[!names(filtRs) %in% remove_names]

## estimating error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

## dereplicating
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names %>% 
  discard(sample.names %in% remove_names)
names(derepRs) <- sample.names %>% 
  discard(sample.names %in% remove_names)

## the core sample inference algorithm (removes sequencing errors)
dadaFs <- dada(derepFs , err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## merge the forward and reverse reads together to obtain the full denoised sequences. 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, trimOverhang = T)

##  construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Checking the frequency of chimeric sequences 
sum(seqtab.nochim)/sum(seqtab)
table(nchar(getSequences(seqtab.nochim)))

out2 <- as.data.frame(out)
out2 <- out2 %>% filter(!rownames(out2) %in% remove_names)

## Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)

#removing asv with less than 10 total hits
seqtab.nochim_filt<-seqtab.nochim %>% as.data.frame() %>%
  select_if(~sum(.) > 10)

dim(seqtab.nochim_filt)

## Making the output from DADA2 compatible with obitools
asv_seqs <- colnames(seqtab.nochim_filt)
asv_headers <- vector(dim(seqtab.nochim_filt)[2], mode="character")

for (i in 1:dim(seqtab.nochim_filt)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## Sequences for attaching taxonomy with obitools
asv_fasta <- c(rbind(asv_headers, asv_seqs))
writeLines(asv_fasta, "arthropod2023_dadaoutput.fasta")

## per sample sequence counts
asv_tab <- t(seqtab.nochim_filt)
row.names(asv_tab) <- sub(">", "", asv_headers)
head(asv_tab)
write.table(asv_tab, "arthropod_ASVs_counts_23.txt", sep="\t", quote=F)

############################################################################################################################


############# Plants (NovaSeq)
# same code as for arthropods

path <- "plant_trimmed_24"             #wherever your results from cutadapt were stored

list.files(path) ## make sure your path is loaded correctly

#read in the forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

## Extract sample names, assuming filenames have format
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1:2], collapse = "_"))

## check out sequence quality
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

## create output for the filtering process
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Standard filtering parameters
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2),truncQ=2, rm.phix = F,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE)

rownames(out)<-sapply(strsplit(basename(rownames(out)), "_"), function(x) paste(x[1:2], collapse = "_"))

## some files were empty so needed to remove them
remove_names <- rownames(out[out[, 2] == 0, ])
filtFs <- filtFs[!names(filtFs) %in% remove_names]
filtRs <- filtRs[!names(filtRs) %in% remove_names]

## estimating error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

## dereplicating
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names %>% 
  discard(sample.names %in% remove_names)
names(derepRs) <- sample.names %>% 
  discard(sample.names %in% remove_names)

## the core sample inference algorithm (removes sequencing errors)
dadaFs <- dada(derepFs , err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, trimOverhang = T)

head(mergers[[1]])

##  construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Checking the frequency of chimeric sequences 
sum(seqtab.nochim)/sum(seqtab)
table(nchar(getSequences(seqtab.nochim)))

out2 <- as.data.frame(out)
out2 <- out2 %>% filter(!rownames(out2) %in% remove_names)

## Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)

#removing asv with less than 10 total hits
seqtab.nochim_filt<-seqtab.nochim %>% as.data.frame() %>%
  select_if(~sum(.) > 10)

dim(seqtab.nochim_filt)

## Making the output from DADA2 compatible with obitools
asv_seqs <- colnames(seqtab.nochim_filt)
asv_headers <- vector(dim(seqtab.nochim_filt)[2], mode="character")

for (i in 1:dim(seqtab.nochim_filt)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## Sequences for attaching taxonomy with obitools
asv_fasta <- c(rbind(asv_headers, asv_seqs))
writeLines(asv_fasta, "plant2024_dadaoutput.fasta")

## per sample sequence counts
asv_tab <- t(seqtab.nochim_filt)
row.names(asv_tab) <- sub(">", "", asv_headers)
tail(asv_tab)
write.table(asv_tab, "plant_ASVs_counts_24.txt", sep="\t", quote=F)

############# Plants (NovaSeq)
# same code as for arthropods

path <- "plant_trimmed_23"             #wherever your results from cutadapt were stored

list.files(path) ## make sure your path is loaded correctly

#read in the forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

## Extract sample names, assuming filenames have format
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[1:2], collapse = "_"))

## check out sequence quality
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

## create output for the filtering process
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Standard filtering parameters
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2),truncQ=2, rm.phix = F,
                     compress=TRUE, multithread=TRUE, matchIDs=TRUE)

rownames(out)<-sapply(strsplit(basename(rownames(out)), "_"), function(x) paste(x[1:2], collapse = "_"))

## some files were empty so needed to remove them
remove_names <- rownames(out[out[, 2] == 0, ])
filtFs <- filtFs[!names(filtFs) %in% remove_names]
filtRs <- filtRs[!names(filtRs) %in% remove_names]

## estimating error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

## dereplicating
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names %>% 
  discard(sample.names %in% remove_names)
names(derepRs) <- sample.names %>% 
  discard(sample.names %in% remove_names)

## the core sample inference algorithm (removes sequencing errors)
dadaFs <- dada(derepFs , err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

## merge the forward and reverse reads together to obtain the full denoised sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, trimOverhang = T)

head(mergers[[1]])

##  construct an amplicon sequence variant table (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Checking the frequency of chimeric sequences 
sum(seqtab.nochim)/sum(seqtab)
table(nchar(getSequences(seqtab.nochim)))

out2 <- as.data.frame(out)
out2 <- out2 %>% filter(!rownames(out2) %in% remove_names)

## Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
head(track)

#removing asv with less than 10 total hits
seqtab.nochim_filt<-seqtab.nochim %>% as.data.frame() %>%
  select_if(~sum(.) > 10)

dim(seqtab.nochim_filt)

## Making the output from DADA2 compatible with obitools
asv_seqs <- colnames(seqtab.nochim_filt)
asv_headers <- vector(dim(seqtab.nochim_filt)[2], mode="character")

for (i in 1:dim(seqtab.nochim_filt)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

## Sequences for attaching taxonomy with obitools
asv_fasta <- c(rbind(asv_headers, asv_seqs))
writeLines(asv_fasta, "plant2023_dadaoutput.fasta")

## per sample sequence counts
asv_tab <- t(seqtab.nochim_filt)
row.names(asv_tab) <- sub(">", "", asv_headers)
tail(asv_tab)
write.table(asv_tab, "plant_ASVs_counts_23.txt", sep="\t", quote=F)

## Move to Obitools/Python