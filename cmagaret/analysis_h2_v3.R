

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#
# source ("analysis_h2_v3.R")
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Finding:  between 2 timepoints there seemed to be an increase in lineages 
#           with a VRC01 resistant clone and decrease in lineages with a 
#           Sensitive clone.
#
# Stats:  I want to know if this expansion of resistant lineages over time was 
#         significantly higher in treated vs placebo.  For this calculation I 
#         don't think we should just focus on the participants with discordant 
#         phenotype -- rather we should look at all multilineage participants 
#         (sticking with those that met primary endpoint only).
#
# We could :  
#
# -- Compare all placebo to all treated -- among ML participants with a 
#    discordant phenotype per our definition.
# -- Compare all placebo to all treated -- among all ML participants
# ---------------------------------------------------------------------------- #


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

# we need to reformat our ID80 readouts
data$ug.ml[data$ug.ml == "no data"] <- NA
data$ug.ml[data$ug.ml == "no growth"] <- NA
data$ug.ml[data$ug.ml == ">100"] <- 100
data$ug.ml <- as.numeric (data$ug.ml)

# correct faults with essential data
data$Lineage.frequency.tp1[data$Lineage.frequency.tp1 == ""] <- NA
data$Lineage.frequency.tp1[data$Lineage.frequency.tp1 == "add data"] <- NA
data$Lineage.frequency.tp1[data$Lineage.frequency.tp1 == "query "] <- NA
data$Lineage.frequency.tp2[data$Lineage.frequency.tp2 == ""] <- NA
data$Lineage.frequency.tp2[data$Lineage.frequency.tp2 == "add data"] <- NA
data$Lineage.frequency.tp2[data$Lineage.frequency.tp2 == "query"] <- NA
data$Lineage.frequency.tp1 <- as.numeric (data$Lineage.frequency.tp1)
data$Lineage.frequency.tp2 <- as.numeric (data$Lineage.frequency.tp2)

# expand the discordant phenotype variable
data$discordant.phenotype.cam <- NA
for (ptid in unique (data$PTID.CAM)) {
  data.tmp <- data[data$PTID.CAM == ptid, "Discordant.phenotype...3."]
  data.tmp <- data.tmp[data.tmp != ""]
  data[data$PTID.CAM == ptid, "discordant.phenotype.cam"] <- data.tmp
}

# dateify our dates
data[, 28] <- as.Date (data[, 28], tryFormats="%m/%d/%Y")
data[, 29] <- as.Date (data[, 29], tryFormats="%m/%d/%Y")
data[, 30] <- as.Date (data[, 30], tryFormats="%m/%d/%Y")
data[, 34] <- as.Date (data[, 34], tryFormats="%m/%d/%Y")

# QA:  which participants have multiple rows for the same lineage?
#for (ptid.tmp in sort (unique (data$PTID.CAM))) {
#  data.tmp <- data[data$PTID.CAM == ptid.tmp & substr (data$Lineage, 1, 1) == "M" & !is.na (data$ug.ml), ]
#  if (length (unique (data.tmp$Lineage)) != nrow (data.tmp)) {
#    print (c (ptid.tmp, data.tmp$Lineage))
#  }
#}

# let's remove the super-infected participant:  V703_1240
data <- data[data$PTID.CAM != 1240, ]

# let's keep only those with discordant phenotypes and dichotomize the 
# phenotype labels
#data.discordant <- data[data$discordant.phenotype.cam == "Y", ]
data.discordant <- data
data.discordant$lineage.phenotype.cam <- data.discordant$Lineage.VRC01.phenotype
data.discordant[data.discordant$lineage.phenotype.cam == "hypersensitive", "lineage.phenotype.cam"] <- "sensitive"
data.discordant[data.discordant$lineage.phenotype.cam == "intermediate resistance", "lineage.phenotype.cam"] <- "resistant"

# make an ordinal phenotype category
data.discordant$lineage.phenotype.num <- rep (1, nrow (data.discordant))
data.discordant[data.discordant$Lineage.VRC01.phenotype == "sensitive", "lineage.phenotype.num"] <- 2
data.discordant[data.discordant$Lineage.VRC01.phenotype == "intermediate resistance", "lineage.phenotype.num"] <- 3
data.discordant[data.discordant$Lineage.VRC01.phenotype == "resistant", "lineage.phenotype.num"] <- 4


# ---------------------------------------------------------------------------- #
# STEP 2:  set up our analysis data frame
# ---------------------------------------------------------------------------- #


# initialize our analysis dataset
data.analysis <- NULL

# loop through each ptid
for (ptid.tmp in sort (unique (data.discordant$PTID.CAM))) {
  data.tmp <- na.omit (data.discordant[data.discordant$PTID.CAM == ptid.tmp & substr (data.discordant$Lineage, 1, 1) == "M" & !is.na (data.discordant$ug.ml), 
                            c ("Trial", "PTID.CAM", "VRC01.Treatment.Group", "Lineage", 
                               "Lineage.VRC01.phenotype", "ug.ml", "Lineage.frequency.tp1", "Lineage.frequency.tp2",
                               "Date.timepoint.1", "Date.timepoint.2")])
  trial.tmp <- unique (data.tmp$Trial)

  # our participant requires at least one row and different dates
  if (nrow (data.tmp) >= 1 & (data.tmp[1, "Date.timepoint.1"] != data.tmp[1, "Date.timepoint.2"])) {

    # reformat the pubid and note
    pubid.tmp <- paste0 (data.tmp[1, "Trial"], "_", sprintf ("%04d", data.tmp[1, "PTID.CAM"]))

    # loop through each lineage
    for (lineage.tmp in sort (unique (data.tmp$Lineage))) {

      # capture our mean IC80
      data.lineage.tmp <- data.tmp[data.tmp$Lineage == lineage.tmp, ]
      ic80.mean <- mean (log10 (data.tmp$ug.ml))

      # set up an entry to record our analysis data
      results.tmp <- data.frame (trial=trial.tmp,
                                 pubid=pubid.tmp,
                                 treatment=data.lineage.tmp$VRC01.Treatment.Group,
                                 lineage=lineage.tmp,
                                 phenotype=data.lineage.tmp$Lineage.VRC01.phenotype,
                                 prop.t1=mean (data.lineage.tmp$Lineage.frequency.tp1),
                                 prop.t2=mean (data.lineage.tmp$Lineage.frequency.tp2),
                                 prop.delta=mean (data.lineage.tmp$Lineage.frequency.tp2) - mean (data.lineage.tmp$Lineage.frequency.tp1),
                                 num.days=data.lineage.tmp$Date.timepoint.2 - data.lineage.tmp$Date.timepoint.1,
                                 prop.slope=(mean (data.lineage.tmp$Lineage.frequency.tp2) - mean (data.lineage.tmp$Lineage.frequency.tp1)) / as.numeric (data.lineage.tmp$Date.timepoint.2 - data.lineage.tmp$Date.timepoint.1),
                                 ic80.mean=ic80.mean)
      # append our results
      data.analysis <- rbind (data.analysis, results.tmp)
    }
  }
}

# create a new ordinal variable for our treatment assignment
data.analysis$treatment.num <- 1
data.analysis[data.analysis$treatment == "10mg/kg", "treatment.num"] <- 2
data.analysis[data.analysis$treatment == "30mg/kg", "treatment.num"] <- 3


# ---------------------------------------------------------------------------- #
# STEP 3:  make plots
# ---------------------------------------------------------------------------- #


# make plot, one panel per treatment assignment
pdf ("~/h2_v4_all_pooled.pdf", height=10, width=6)
par (mfrow=c (3, 1))

# treatment = placebo
data.plot <- data.analysis[data.analysis$treatment.num == 1, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="Placebo",
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# treatment = low dose VRC01
data.plot <- data.analysis[data.analysis$treatment.num == 2, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="10mg/kg", 
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# treatment = high dose VRC01
data.plot <- data.analysis[data.analysis$treatment.num == 3, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="30mg/kg",
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# done
dev.off ()






# make plot, one panel per treatment assignment
pdf ("~/h2_v4_all_v703.pdf", height=10, width=6)
par (mfrow=c (3, 1))

# treatment = placebo
data.plot <- data.analysis[data.analysis$treatment.num == 1 & data.analysis$trial == 703, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="Placebo",
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# treatment = low dose VRC01
data.plot <- data.analysis[data.analysis$treatment.num == 2 & data.analysis$trial == 703, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="10mg/kg", 
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# treatment = high dose VRC01
data.plot <- data.analysis[data.analysis$treatment.num == 3 & data.analysis$trial == 703, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="30mg/kg",
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# done
dev.off ()







# make plot, one panel per treatment assignment
pdf ("~/Desktop/h2_v4_all_v704.pdf", height=10, width=6)
par (mfrow=c (3, 1))

# treatment = placebo
data.plot <- data.analysis[data.analysis$treatment.num == 1 & data.analysis$trial == 704, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="Placebo",
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# treatment = low dose VRC01
data.plot <- data.analysis[data.analysis$treatment.num == 2 & data.analysis$trial == 704, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="10mg/kg", 
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# treatment = high dose VRC01
data.plot <- data.analysis[data.analysis$treatment.num == 3 & data.analysis$trial == 704, ]
plot.pch <- as.integer (as.factor (data.plot$pubid))
plot.col <- rep ("blue", nrow (data.plot))
plot.col[data.plot$phenotype == "intermediate resistance" | data.plot$phenotype == "resistant"] <- "red"
plot (rep (0, nrow (data.plot)), data.plot$prop.t1, pch=plot.pch, col=plot.col, main="30mg/kg",
      xlim=c (0, 30), ylim=c (0, 1), xlab="Days After T1", ylab="Proportion of Lineage")
points (data.plot$num.days, data.plot$prop.t2, pch=plot.pch, col=plot.col)
segments (rep (0, nrow (data.plot)), data.plot$prop.t1, data.plot$num.days, data.plot$prop.t2, col=plot.col)

# done
dev.off ()


# ---------------------------------------------------------------------------- #
# STEP 4:  do our statistical tests
# ---------------------------------------------------------------------------- #


# "First, for each of the 2 variables (X on the x-axis, Y on the y-axis), scale 
#  them to have mean zero and standard deviation one.  Call these 2 new 
#  variables X* and Y*"
#
rescale <- function (x) {
  x.tmp <- x - mean (x)
  x.tmp <- x.tmp / sd (x.tmp)
  return (x.tmp)
}
x.star <- rescale (data.analysis$prop.delta)
y.star <- rescale (data.analysis$ic80.mean)

# "2. Define a new variable as Z* = X* + Y*.  The scaling in step 1 was needed 
#  to have this sum make sense."
z.star <- x.star + y.star

# "3. Do a Jonckheere-Terpstra test to see if Z* tends to increase going from 
#  placebo to low dose VRC01 to high dose VRC01 (one-sided test).
#
# "The logic here is Z* will be low if both X and Y are low, and X* will be 
#  high if both X and Y are high, so this fits your description of what we're 
#  looking a test to pick up if its there."

# HVTN 703
p.jt.v703 <- JonckheereTerpstraTest (x=z.star[data.analysis$trial == 703],
                                     g=data.analysis[data.analysis$trial == 703, "treatment.num"],
                                     alternative="increasing")$p.value

# HVTN 704
p.jt.v704 <- JonckheereTerpstraTest (x=z.star[data.analysis$trial == 704],
                                     g=data.analysis[data.analysis$trial == 704, "treatment.num"],
                                     alternative="increasing")$p.value

# pooled
p.jt.pooled <- JonckheereTerpstraTest (x=z.star,
                                       g=data.analysis$treatment.num,
                                       alternative="increasing")$p.value




