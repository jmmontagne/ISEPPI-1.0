#Load packages
source('http://www.bioconductor.org/biocLite.R')
biocLite('ShortRead')
biocLite('tidyverse')

library('ShortRead') #for readFastq
library('tidyverse') #for read_tsv
library('caroline')

###Read in fq files###
fq_files <- list.files('.', '.fq$', full=TRUE)
names(fq_files) <- sub(".fq", "", basename(fq_files))
#names(fq_files)

#Read in 1MM primer sequences
MM_primers <- sread(readFasta("1MM_primers.fa"))
seqs_MM_primers <- as.character(MM_primers)
###Import 1MM efficiency vector
MM_allele_df <- readRDS('1MM_unitEfficiencies.rds')

for (x in names(fq_files)) {fq <- readFastq(paste(x,".fq",sep=""))
  #Access sequences
  seqs_fq <- sread(fq)
  remove(fq)

  ##Read in MiXCR data
  mixcr_data <- read_tsv(file=paste(x, "_clones.txt", sep=""), col_types = cols_only(nSeqCDR3 = col_character(), clonalSequence = col_character(), aaSeqCDR3 = col_character(), cloneCount = col_number(), allJHitsWithScore = col_character(), reads = col_character()))

  ##Make list of primers used per CDR3
  mixcr_data$reads <- as.character(mixcr_data$reads)
  mixcr_data$reads <- sapply(mixcr_data$reads, function(x) as.numeric(unlist(strsplit(x, split=','))))
  ##add 1 because reads start at 0
  mixcr_data$reads <- sapply(mixcr_data$reads, function(x) x+1)
  mixcr_data$readSeqs <- sapply(mixcr_data$reads, function(x) data.frame(seqs_fq[x,]))
  mixcr_data$primerSeqs <- sapply(mixcr_data$readSeqs, function(x) as.list(substr(x, 1, 20)))
  rm(seqs_fq)

  ##Make primer usage vector
  #add number of times each 20mer appears
  mixcr_data$primerSeqs <- sapply(mixcr_data$primerSeqs, function(x) table(unlist(x)))
  #get an ordered vector for the number of times each primer is used
  mixcr_data$primerCounts <- lapply(mixcr_data$primerSeqs, function(x) as.vector((x[seqs_MM_primers])))
  #change NAs to 0s
  mixcr_data$primerCounts <- lapply(mixcr_data$primerCounts, sapply, function(x) replace(x,is.na(x),0))

  ##Convert primer usage vector to unit vector
  make_unit_vector <- function(x) {x / sqrt(sum(x^2))} #function to generate unit vector
  mixcr_data$unitPrimerCounts <- lapply(mixcr_data$primerCounts, function(x) make_unit_vector(x))

  #compute euclidean distance for each of the primercount unit vectors with each of the 1MM efficiency unit vectors
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
  indx <- expand.grid(allele_effs=1:nrow(MM_allele_df), primerCounts=1:nrow(mixcr_data))
  dists <- as.vector(mapply(euc.dist, mixcr_data$unitPrimerCounts[indx[,2]],MM_allele_df$unit_efficiency_vector[indx[,1]]))
  dists_by_CDR3 <- split(dists, ceiling(seq_along(dists)/47))
  mixcr_data$euc_dist <- dists_by_CDR3
  rm(dists_by_CDR3, indx, dists)

  #determine index of minimum Euclidean distance and match it to unique FR3
  mixcr_data$min_euc_dist <- as.numeric(lapply(mixcr_data$euc_dist, function(x) which(x==min(x))))
  mixcr_data$UPGMA_FR3 <- as.numeric(lapply(mixcr_data$min_euc_dist, function(x) MM_allele_df$UPGMA_order[x]))
  mixcr_data$predictedFR3 <- MM_allele_df$alleles[as.numeric(unlist(mixcr_data$min_euc_dist))]

  saveRDS(mixcr_data, file=paste(x, "_final_clones.RDS", sep=""))
  mixcr_data$primerCounts <- NULL
  mixcr_data$primerSeqs <- NULL
  mixcr_data$readSeqs <- NULL
  mixcr_data$unitPrimerCounts <- NULL
  mixcr_data$euc_dist <- NULL
  mixcr_data$min_euc_dist <- NULL
  mixcr_data$reads <- NULL
  write.delim(mixcr_data, file=paste(x, "_final_clones.txt", sep=""), quote = FALSE, row.names=FALSE, sep = "\t")
}
