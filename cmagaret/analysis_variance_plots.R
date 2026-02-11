

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
#  x = VRC01 conc at EDDI-7
#
#  y1a = % sequences resistant (IC80 >= 3)
#  y1b = % sequences resistant (IC80 >= 1)
#  y2 = IC80 variance
#  y3 = Pi
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# V703_1240 (high dose) is a dual infection and should be excluded and 
# V703_2805 from placebo group needs to be added.  

# ---------------------------------------------------------------------------- #
# STEP 0:  get everything ready
# ---------------------------------------------------------------------------- #


# initialize our working environment
rm (list=ls ())

# load data
setwd ("/path/to/data")
data.analysis <- read.csv ("data_for_analysis_v2.csv", header=T)
data.pi <- read.csv ("amp_pi_v1.csv", header=T)
data.ic80 <- read.csv ("ic80.csv", header=T)
data.seq.703 <- read.csv ("V703_REN_rm_flagged_seqs_00~count_uncollapsed_seqs_for_TP-lineage (1).csv", header=T)
data.seq.704 <- read.csv ("V704_REN_flagged_seqs_removed_lineages_summary.csv", header=T)

# reshape the analysis dataset so that it only includes the information we need
data.lineages <- data.frame (pubid=paste0 ("V", data.analysis$Trial, "_", sprintf ("%04d", data.analysis$PTID.CAM)),
                             ic80=data.analysis$ug.ml, Lineage=data.analysis$Lineage)
#data.lineages[data.lineages$Lineage == "MxL1.1", "Lineage"] <- "MxL1"
data.lineages <- data.lineages[data.lineages$Lineage %in% c ("MxL1", "MxL1.1", "MxL2", "MxL3", "MxL4"), ]
data.lineages[data.lineages$ic80 == ">100", "ic80"] <- as.character (100 * sqrt (2))
data.lineages[data.lineages$ic80 == "no growth", "ic80"] <- NA
data.lineages$ic80 <- as.numeric (data.lineages$ic80)

# compose analysis dataset:  Pi
data.plot <- data.frame (pubid=sort (unique (substr (data.analysis$plasmid_name, 1, 9)))[-1:-2])
data.plot <- merge (data.plot, data.pi, by=1)

# compose analysis dataset:  predicted VRC01 concentration at EDDI-7
data.ic80$pubid <- paste0 ("V", data.ic80$Trial, "_", sprintf ("%04d", data.ic80$PTID.CAM))
data.plot <- merge (data.plot, data.ic80[, c ("pubid", "Cc_EDDI.7")], by=1)
data.plot$percent.resistant.1 <- NA
data.plot$percent.resistant.3 <- NA
data.plot$ic80.variance <- NA
data.plot <- merge (data.plot, data.ic80[, c ("pubid", "rx_code")], by=1, all.x=T)

# remove participant V703_1240
data.plot <- data.plot[data.plot$pubid != "V703_1240", ]

# remove participant V704_1654
data.plot <- data.plot[data.plot$pubid != "V704_1654", ]


# ---------------------------------------------------------------------------- #
# STEP 1:  take care of business
# ---------------------------------------------------------------------------- #


# compile our sequence counts per lineage
for (pubid.tmp in data.plot$pubid) {

  if (pubid.tmp %in% data.lineages$pubid) {
    lineages.tmp <- data.lineages[data.lineages$pubid == pubid.tmp, "Lineage"]
    lineages.tmp[lineages.tmp == "MxL1.1"] <- "MxL1"
    ic80.tmp <- data.lineages[data.lineages$pubid == pubid.tmp, "ic80"]
    lineages.tmp <- lineages.tmp[!is.na (ic80.tmp)]
    ic80.tmp <- as.vector (na.omit (ic80.tmp))

    # consolidate identical lineages (recording the mean of the IC80)
    if (sum (table (lineages.tmp) > 1) > 0) {
      lineages.unique <- sort (unique (lineages.tmp))
      ic80.unique <- rep (NA, length (lineages.unique))
      for (lineage.num in 1:length (lineages.unique)) {
        ic80.unique[lineage.num] <- mean (ic80.tmp[lineages.tmp == lineages.unique[lineage.num]])
      }
      lineages.tmp <- lineages.unique
      ic80.tmp <- ic80.unique
    }

    # subset the data for the participant, from the proper trial
    if (substr (pubid.tmp, 2, 4) == "703") {
      data.tmp <- data.seq.703[substr (data.seq.703$Id, 1, 9) == pubid.tmp, ]
    } else if (substr (pubid.tmp, 2, 4) == "704") {
      data.tmp <- data.seq.704[substr (data.seq.704$Id, 1, 9) == pubid.tmp, ]
    }

    # capture the data for the earliest visit
    data.tmp <- data.tmp[as.numeric (substr (data.tmp$Id, 11, 13)) == min (as.numeric (substr (data.tmp$Id, 11, 13))), ]
    if (nrow (data.tmp) > 1) {
      print (paste0 ("BARF!  (", pubid.tmp, ")"))
    }

    # capture our primary numbers
    total.tmp <- data.tmp$Total
    counts.tmp <- data.tmp[, paste0 (unique (lineages.tmp), ".c.")]

    #
    data.plot[data.plot$pubid == pubid.tmp, "percent.resistant.1"] <- sum (counts.tmp[ic80.tmp >= 1]) / total.tmp
    data.plot[data.plot$pubid == pubid.tmp, "percent.resistant.3"] <- sum (counts.tmp[ic80.tmp >= 3]) / total.tmp

    # calculate our variances
    conc.tmp <- NULL
    for (lineage.num in 1:length (lineages.tmp)) {
      conc.tmp <- c (conc.tmp, rep (ic80.tmp[lineage.num], counts.tmp[lineage.num]))
    }
    data.plot[data.plot$pubid == pubid.tmp, "ic80.variance"] <- var (log10 (conc.tmp))
  }
}

# assess our missing data
nrow (na.omit (data.plot[, c ("Cc_EDDI.7", "percent.resistant.1")]))
nrow (na.omit (data.plot[, c ("Cc_EDDI.7", "percent.resistant.3")]))
nrow (na.omit (data.plot[, c ("Cc_EDDI.7", "ic80.variance")]))
nrow (na.omit (data.plot[, c ("Cc_EDDI.7", "pi")]))
sum (data.plot$rx_code != "C3" & is.na (data.plot$percent.resistant.3))

# outliers
data.plot[order (data.plot$ic80.variance, decreasing=T), c ("pubid", "rx_code", "Cc_EDDI.7", "ic80.variance")]
data.analysis[, c (2:3, 5:6, 9, 13)]
data.seq.703[55:56, ]
data.seq.703[89:90, ]

# add plot information:  color code by treatment:
#   African sites:
#     Low Dose:   #CC99FF
#     High Dose:  #7030A0
#   N. American sites:
#     Low Dose:   #A9D08E
#     High Dose:  #009900
#   PBO:  #BFBFBF
data.plot$plot.col <- data.plot$rx_code
data.plot$plot.col[data.plot$rx_code == "T1"] <- "#CC99FF"
data.plot$plot.col[data.plot$rx_code == "T2"] <- "#7030A0"
data.plot$plot.col[data.plot$rx_code == "C3"] <- "#BFBFBF"
data.plot$plot.col[data.plot$rx_code == "T1" & substr (data.plot$pubid, 4, 4) == "4"] <- "#A9D08E"
data.plot$plot.col[data.plot$rx_code == "T2" & substr (data.plot$pubid, 4, 4) == "4"] <- "#009900"

# add plot information:  indicate folks with both highly sensitive/resistant 
# lineages
#   -- V703_0514
#   -- V703_1714
#   -- V703_2769
data.plot$plot.pch <- 1
data.plot$plot.pch[data.plot$pubid %in% c ("V703_0514", "V703_1714", "V703_2769")] <- 5

# add plot information:  indicate outliers
data.plot$plot.pch.outliers <- data.plot$plot.pch
data.plot$plot.pch.outliers[data.plot$pubid == "V703_2769"] <- 23


# ---------------------------------------------------------------------------- #
# STEP 2:  make plot
# ---------------------------------------------------------------------------- #

# go!
pdf (file="plot_dreeves_v5_full.pdf", height=4, width=11)
par (mfrow=c (1, 3), mar = c (6, 4, 4, 2) + 0.1)
plot (data.plot$Cc_EDDI.7, data.plot$ic80.variance, main=expression (bold ("A. Variance of IC"[80])), 
      xlab="Predicted VRC01 Concentration at Acquisition (µg/ml)", ylab=expression ("Variance of Log10 IC"[80]),
      col=data.plot$plot.col, pch=data.plot$plot.pch, bg=data.plot$plot.col,
      sub=paste0 ("(p = ", round (cor.test (data.plot$Cc_EDDI.7, data.plot$ic80.variance, method="spearman")$p.value, 3), ")"))
abline (v=20, col="#FF7E79", lty=2)
rect (-2, 0.035, 16, 1.76, lty=2, col=NA, border="gray25")
legend ("topright", pch=c (15, 15, 15, 15, 5), col=c ("#CC99FF", "#7030A0", "#A9D08E", "#009900", "black"), cex=0.8,
        legend=c ("VRC01 Low Dose (HVTN 703)", "VRC01 High Dose (HVTN 703)", "VRC01 Low Dose (HVTN 704)", "VRC01 High Dose (HVTN 704)", "Highly Discordant Participant"))
plot (data.plot$Cc_EDDI.7, data.plot$percent.resistant.3, main="B. Proportion of Resistant Sequences",
      xlab="Predicted VRC01 Concentration at Acquisition (µg/ml)", ylab="Proportion of Resistant Sequences",
      col=data.plot$plot.col, pch=data.plot$plot.pch,
      sub=paste0 ("(p = ", round (cor.test (data.plot$Cc_EDDI.7, data.plot$percent.resistant.3, method="spearman")$p.value, 3), ")"))
abline (v=20, col="#FF7E79", lty=2)
plot (data.plot$Cc_EDDI.7, data.plot$pi, main="C. Nucleotide Diversity", 
      xlab="Predicted VRC01 Concentration at Acquisition (µg/ml)", ylab="Nucleotide Diversity (Pi)",
      col=data.plot$plot.col, pch=data.plot$plot.pch,
      sub=paste0 ("(p = ", round (cor.test (data.plot$Cc_EDDI.7, data.plot$pi, method="spearman")$p.value, 3), ")"))
abline (v=20, col="#FF7E79", lty=2)
dev.off ()


# ---------------------------------------------------------------------------- #
#                                    - 30 -
# ---------------------------------------------------------------------------- #



