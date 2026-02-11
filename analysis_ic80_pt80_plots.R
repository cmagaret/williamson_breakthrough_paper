

# script to visualize IC80 and PT80 among participants with discordant lineages

# load data
setwd ("/path/to/data")
data.ic80 <- read.csv ("ic80_20240912.csv", header=T)
data.pt80 <- read.csv ("pt80_20240912.csv", header=T)

# filter
data.pt80 <- data.pt80[data.pt80$Clone.used.for.fold.calculation == "y", ]
data.pt80 <- data.pt80[!duplicated (data.pt80), ]
data.pt80 <- data.pt80[data.pt80$PTID != 132, ]

# refactor into plotting info
data.plot <- data.frame (ptid=sort (unique (data.pt80$PTID)), ic80.sens=NA, ic80.resis=NA)
for (ptid in data.plot$ptid) {
  data.tmp <- data.pt80[data.pt80$PTID == ptid, ]
  data.plot[data.plot$ptid == ptid, "ic80.sens"] <- min (data.tmp$ug.ml)
  data.plot[data.plot$ptid == ptid, "ic80.resis"] <- max (data.tmp$ug.ml)
}

# set up some plot parameters
plot.col <- rep ("red", nrow (data.pt80))
plot.col[data.pt80$VRC01.Treatment.Group == "10mg/kg"] <- "blue"

# plot
pdf ("amp_discordant_ic80_pt80_v1.pdf", height=6, width=6)
plot (log10 (data.pt80$ug.ml), data.pt80$PT80.at.EDDI.7, xlim=c (-1, 2),
      xlab="IC80", ylab="PT80", xaxt="n", col=plot.col)
axis (side=1, at=seq (-1, 2, 0.5), labels=10^seq (-1, 2, 0.5))
abline (v=0, lty=2, col="darkgray")
dev.off ()

# plot something else
pdf ("amp_discordant_ic80_v1.pdf", height=6, width=6)
plot (log10 (data.plot$ic80.sens), log10 (data.plot$ic80.resis), xlim=c (-1, 2), ylim=c (-1, 2),
      xlab="IC80 of Sensitive Virus", ylab="IC80 of Resistant Virus", xaxt="n", yaxt="n", col=plot.col)
axis (side=1, at=seq (-1, 2, 0.5), labels=10^seq (-1, 2, 0.5))
axis (side=2, at=seq (-1, 2, 0.5), labels=10^seq (-1, 2, 0.5))
#abline (a=0, b=1, lty=2, col="darkgray")
dev.off ()




