

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# source ("analysis_h1_v3a.R")
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Finding:  Across all my 703 and 704 ML participants I found that there was a 
#           higher proportion of individuals with discordant VRC01 phenotype in 
#           the treatment group compared to placebo.
#
# Stats calculation:  I want to know if these differences in numbers/
#                     proportions is significant.
# 
# We could:  master treatment (ie VRC01 low and high dose pooled) vs placebo, 
#            then, 
#            A three level approach looking for an increasing effect (if any)
# ---------------------------------------------------------------------------- #

# for V4:  "V703_1240 (high dose) is a dual infection and should be excluded 
#           and V703_2805 from placebo group needs to be added."

# ---------------------------------------------------------------------------- #
# STEP 0:  get that motor runnin
# ---------------------------------------------------------------------------- #


# initialize our workspace
rm (list=ls ())
set.seed (1479)

# load required packages
library ("DescTools")

# CONFIG:  number of permutation iterations
perm.iter <- 10000


# ---------------------------------------------------------------------------- #
# STEP 1:  load and set up data
# ---------------------------------------------------------------------------- #


# load data
setwd ("/path/to/data")
data <- read.csv ("ptids_w_mli_v5.csv", header=T)

# QA:  confirm that my manual additions are fine
ptid.tmp <- substr (data$plasmid_name, 6, 9)
ptid.tmp[nchar (ptid.tmp) != 4] <- NA
ptid.tmp <- as.numeric (ptid.tmp)
table (data$PTID == ptid.tmp)

# QA:  confirm that we don't have participants with the same "PTID" in different 
# studies
#intersect (unique (data[data$Trial == 703, "PTID"]), unique (data[data$Trial == 704, "PTID"]))

# we need to reformat our ID80 readouts
data$ug.ml[data$ug.ml == "no data"] <- NA
data$ug.ml[data$ug.ml == "no growth"] <- NA
data$ug.ml[data$ug.ml == ">100"] <- 100
data$ug.ml <- as.numeric (data$ug.ml)

# remove participant V703_1240
data <- data[data$PTID.CAM != 1240, ]

# remove participant V704_1654
data <- data[data$PTID.CAM != 1654, ]


# ---------------------------------------------------------------------------- #
# STEP 2:  set up our analysis data frame
# ---------------------------------------------------------------------------- #


# initialize our analysis dataset
deltas <- list ()
data.analysis <- NULL

# loop through each ptid
for (ptid.tmp in unique (data$PTID.CAM)) {

  # capture the data for this ptid and remove rows with missing IC80s
  #   Q:  lineage to include:  start with "M"?
  data.tmp <- data[data$PTID.CAM == ptid.tmp & substr (data$Lineage, 1, 1) == "M" & !is.na (data$ug.ml), ]

  # our participant requires at least two lineages
  if (nrow (data.tmp) >= 2) {

    # reformat the pubid
    pubid.tmp <- paste0 (data.tmp[1, "Trial"], "_", sprintf ("%04d", data.tmp[1, "PTID.CAM"]))

    # get our pairs of IC80 values and their deltas
    pairs.tmp <- CombPairs (log10 ((data.tmp$ug.ml)))
    deltas.tmp <- abs (pairs.tmp[, 1] - pairs.tmp[, 2])
    deltas[[pubid.tmp]] <- deltas.tmp

    # append our results
    data.analysis <- rbind (data.analysis, data.frame (trial=data.tmp[1, "Trial"], pubid=pubid.tmp, mean.delta=mean (deltas.tmp), treatment=data.tmp[1, "VRC01.Treatment.Group"]))
  }
}

# reshape the data for our plot
data.plot.703 <- data.analysis[data.analysis$trial == 703, ][order (data.analysis[data.analysis$trial == 703, "mean.delta"]), ]
data.plot.704 <- data.analysis[data.analysis$trial == 704, ][order (data.analysis[data.analysis$trial == 704, "mean.delta"]), ]
data.plot.pooled <- rbind (data.plot.703, data.plot.704)
data.plot.pooled <- data.plot.pooled[order (data.plot.pooled$mean.delta), ]

# re-org the plot data by treatment group
data.plot.703 <- rbind (data.plot.703[data.plot.703$treatment == "Control", ], 
                        data.plot.703[data.plot.703$treatment == "10mg/kg", ],
                        data.plot.703[data.plot.703$treatment == "30mg/kg", ])
data.plot.704 <- rbind (data.plot.704[data.plot.704$treatment == "Control", ], 
                        data.plot.704[data.plot.704$treatment == "10mg/kg", ],
                        data.plot.704[data.plot.704$treatment == "30mg/kg", ])
data.plot.pooled <- rbind (data.plot.pooled[data.plot.pooled$treatment == "Control", ], 
                           data.plot.pooled[data.plot.pooled$treatment == "10mg/kg", ],
                           data.plot.pooled[data.plot.pooled$treatment == "30mg/kg", ])

# create an orginal equivalent of our pooled variable for the J-T test
data.analysis$treatment.num <- 1
data.analysis[data.analysis$treatment == "10mg/kg", "treatment.num"] <- 2
data.analysis[data.analysis$treatment == "30mg/kg", "treatment.num"] <- 3


# ---------------------------------------------------------------------------- #
# STEP 3:  do our statistical tests
# ---------------------------------------------------------------------------- #


# HVTN 703
p.wilcoxon.v703 <- wilcox.test (x=data.analysis[data.analysis$trial == 703 & data.analysis$treatment == "Control", "mean.delta"],
                                y=data.analysis[data.analysis$trial == 703 & data.analysis$treatment != "Control", "mean.delta"],
                                alternative="two.sided")$p.value
p.jt.v703 <- JonckheereTerpstraTest (x=data.analysis[data.analysis$trial == 703, "mean.delta"],
                                     g=data.analysis[data.analysis$trial == 703, "treatment.num"],
                                     alternative="increasing")$p.value

# HVTN 704
p.wilcoxon.v704 <- wilcox.test (x=data.analysis[data.analysis$trial == 704 & data.analysis$treatment == "Control", "mean.delta"],
                                y=data.analysis[data.analysis$trial == 704 & data.analysis$treatment != "Control", "mean.delta"],
                                alternative="two.sided")$p.value
p.jt.v704 <- JonckheereTerpstraTest (x=data.analysis[data.analysis$trial == 704, "mean.delta"],
                                     g=data.analysis[data.analysis$trial == 704, "treatment.num"],
                                     alternative="increasing")$p.value

# trials pooled
p.wilcoxon.pooled <- wilcox.test (x=data.analysis[data.analysis$treatment == "Control", "mean.delta"], 
                                  y=data.analysis[data.analysis$treatment != "Control", "mean.delta"],
                                  alternative="two.sided")$p.value
p.jt.pooled <- JonckheereTerpstraTest (x=data.analysis$mean.delta,
                                       g=data.analysis$treatment.num,
                                       alternative="increasing")$p.value


#  HVTN 703:
#    p (Wilcoxon) = 0.492
#    p (J-T) = 0.398
#
#  HVTN 704:
#    p (Wilcoxon) = 0.134
#    p (J-T) = 0.053
#
#  Pooled:
#    p (Wilcoxon) = 0.111
#    p (J-T) = 0.116


# ---------------------------------------------------------------------------- #
# STEP 4:  make plot
# ---------------------------------------------------------------------------- #


# set up plot
pdf (file="~/Desktop/avt_cm_h1_deltas_v6.pdf", height=5, width=20)
layout (matrix (c (1, 2, 3), nrow=1), widths=c (0.3, 0.38, 0.6))

# HVTN 703
#   Low Dose:  #CC99FF
#   High Dose:  #7030A0
#   PBO:  #BFBFBF
plot (1:nrow (data.plot.703), data.plot.703$mean.delta, ylim=c (0, 3), 
      xaxt="n", type="n", main="HVTN 703", xlab="Participants", ylab="Pairwise Difference IC80 (Log10)")
axis (side=1, at=c (7.5, 13.5), labels=F)
abline (v=1:nrow (data.plot.703), lty=3, col="lightgray")

# establish a vector of treatment colors
trt.col <- rep ("#7030A0", nrow (data.plot.703))
trt.col[data.plot.703$treatment == "Control"] <- "#BFBFBF"
trt.col[data.plot.703$treatment == "10mg/kg"] <- "#CC99FF"

# plot
for (x in 1:nrow (data.plot.703)) {
  points (rep (x, length (deltas[[data.plot.703[x, "pubid"]]])), deltas[[data.plot.703[x, "pubid"]]], col="darkgray", cex=2)
}
points (1:nrow (data.plot.703), data.plot.703$mean.delta, pch="-", col=trt.col, cex=3, ylim=c (1, 3))
mtext (c ("Placebo", "Low Dose", "High Dose"), side=1, line=1, col=c ("#BFBFBF", "#CC99FF", "#7030A0"), at=c (4, 10.5, 16.5), cex=0.7)

# HVTN 704
#   Low Dose:  #A9D08E
#   High Dose:  #009900
#   PBO:  #BFBFBF
plot (1:nrow (data.plot.704), data.plot.704$mean.delta, ylim=c (0, 3), 
      xaxt="n", type="n", main="HVTN 704", xlab="Participants", ylab="Pairwise Difference IC80 (Log10)")
axis (side=1, at=c (11.5, 20.5), labels=F)
abline (v=1:nrow (data.plot.704), lty=3, col="lightgray")

# establish a vector of treatment colors
trt.col <- rep ("#009900", nrow (data.plot.704))
trt.col[data.plot.704$treatment == "Control"] <- "#BFBFBF"
trt.col[data.plot.704$treatment == "10mg/kg"] <- "#A9D08E"

# plot
for (x in 1:nrow (data.plot.704)) {
  points (rep (x, length (deltas[[data.plot.704[x, "pubid"]]])), deltas[[data.plot.704[x, "pubid"]]], col="darkgray", cex=2)
}
points (1:nrow (data.plot.704), data.plot.704$mean.delta, pch="-", col=trt.col, cex=3, ylim=c (1, 3))
mtext (c ("Placebo", "Low Dose", "High Dose"), side=1, line=1, col=c ("#BFBFBF", "#A9D08E", "#009900"), at=c (6, 16, 24), cex=0.7)

# trials pooled
#   Low Dose:  #F4B084
#   High Dose:  #C65911
#   PBO:  #BFBFBF
plot (1:nrow (data.plot.pooled), data.plot.pooled$mean.delta, ylim=c (0, 3), 
      xaxt="n", type="n", main="AMP Trials Pooled", xlab="Participants", ylab="Pairwise Difference IC80 (Log10)")
axis (side=1, at=c (18.5, 33.5), labels=F)
abline (v=1:nrow (data.plot.pooled), lty=3, col="lightgray")

# establish a vector of treatment colors
trt.col <- rep ("#C65911", nrow (data.plot.pooled))
trt.col[data.plot.pooled$treatment == "Control"] <- "#BFBFBF"
trt.col[data.plot.pooled$treatment == "10mg/kg"] <- "#F4B084"

# plot
for (x in 1:nrow (data.plot.pooled)) {
  points (rep (x, length (deltas[[data.plot.pooled[x, "pubid"]]])), deltas[[data.plot.pooled[x, "pubid"]]], col="darkgray", cex=2)
}
points (1:nrow (data.plot.pooled), data.plot.pooled$mean.delta, pch="-", col=trt.col, cex=3, ylim=c (1, 3))
mtext (c ("Placebo", "Low Dose", "High Dose"), side=1, line=1, col=c ("#BFBFBF", "#F4B084", "#C65911"), at=c (9.5, 26, 40), cex=0.7)

# done
dev.off ()


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #






# ---------------------------------------------------------------------------- #
#                                    - 30 -
# ---------------------------------------------------------------------------- #

