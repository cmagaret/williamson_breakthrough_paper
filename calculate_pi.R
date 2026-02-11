

# script to calculate Pi (sequence diversity) from Env amino acid sequences


# refresh our workspace
rm (list=ls ())
set.seed (123)
library (seqinr)
library (divIndex)

# load our data
path.data <- "/path/to/sequence/files"
filenames <- list.files (path.data, pattern="*.fasta", full.names=TRUE)

# initialize results object
results <- data.frame (pubid=substr (filenames, 81, 89), visit=NA, num.seqs=NA, pi=NA)

# take care of business
for (filename.tmp in filenames) {

  # load sequences
  seq.tmp <- read.fasta (filename.tmp, seqtype="AA")
  pubid.tmp <- substr (filename.tmp, 81, 89)

  # parse out the sequence metadata and determine the minimum visit code
  seq.bits <- as.data.frame (matrix (unlist (strsplit (names (seq.tmp), split="_")), ncol=5, byrow=T))
  seq.bits[, 3] <- as.numeric (seq.bits[, 3])
  earliest.visit <- min (seq.bits[, 3])

  # aggregate our sequence inforamtion
  seqname.keepers <- names (seq.tmp)[seq.bits[, 3] == earliest.visit]
  seq.vector <- data.frame (seqname=seqname.keepers, seq=NA)
  for (seqname.tmp in seqname.keepers) {
    seq.vector[seq.vector$seqname == seqname.tmp, "seq"] <- paste (seq.tmp[[seqname.tmp]], collapse="")
  }

  # make our sequence vector and counts required by the Pi function
  seq.table <- table (seq.vector[, 2])
  seq.pi <- data.frame (seq=names (seq.table), counts=as.vector (seq.table))

  # calculate and record Pi for this participant
  results[results$pubid == pubid.tmp, "visit"] <- earliest.visit
  results[results$pubid == pubid.tmp, "num.seqs"] <- length (seq.tmp)
  results[results$pubid == pubid.tmp, "pi"] <- pi (seqs=seq.pi[, 1], count=seq.pi[, 2], weighted=F)
}

# export results
write.csv (results, file="~/Desktop/amp_pi_v1.csv", row.names=F)















