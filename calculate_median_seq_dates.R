

# initialize our workspace
rm (list=ls ())
setwd ("/path/to/data")
data <- read.csv ("amp_pubids_visits_dates.csv", header=T)

# initialize results dataframe
pubids <- sort (unique (data$pub_id))
results <- data.frame (pubid=pubids, num.visits=NA, visit.1=NA, visit.2=NA, visit.3=NA, visit.4=NA, 
                       date.1=NA, date.2=NA, date.3=NA, date.4=NA, delta.date.1.date.2=NA)

# loop through each pubid
for (pubid.tmp in pubids) {
  data.tmp <- data[data$pub_id == pubid.tmp, ]
  data.tmp$drawdt <- as.Date (data.tmp$drawdt)
  results[results$pubid == pubid.tmp, 2] <- nrow (data.tmp)

  # number of visits == 1
  if (nrow (data.tmp) == 1) {
    results[results$pubid == pubid.tmp, 3] <- data.tmp$visitno_craig
    results[results$pubid == pubid.tmp, 7] <- as.character (data.tmp$drawdt)
    results[results$pubid == pubid.tmp, 11] <- 0

  # number of visits == 2
  } else if (nrow (data.tmp) == 2) {
    results[results$pubid == pubid.tmp, 3:4] <- data.tmp$visitno_craig
    results[results$pubid == pubid.tmp, 7:8] <- as.character (data.tmp$drawdt)
    results[results$pubid == pubid.tmp, 11] <- sort (data.tmp$drawdt)[2] - sort (data.tmp$drawdt)[1]

  # number of visits == 3
  } else if (nrow (data.tmp) == 3) {
    results[results$pubid == pubid.tmp, 3:5] <- data.tmp$visitno_craig
    results[results$pubid == pubid.tmp, 7:9] <- as.character (data.tmp$drawdt)
    results[results$pubid == pubid.tmp, 11] <- sort (data.tmp$drawdt)[2] - sort (data.tmp$drawdt)[1]

  # number of visits == 4
  } else if (nrow (data.tmp) == 4) {
    results[results$pubid == pubid.tmp, 3:6] <- data.tmp$visitno_craig
    results[results$pubid == pubid.tmp, 7:10] <- as.character (data.tmp$drawdt)
    results[results$pubid == pubid.tmp, 11] <- sort (data.tmp$drawdt)[2] - sort (data.tmp$drawdt)[1]
  }
}

# export results
write.csv (results, file="amp_date_delta_seq_visits_v1.csv", row.names=F)





