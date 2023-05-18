# This data was obtained from:
# https://doi.org/10.1371/journal.pone.0114947.s001
FILE <- "data-raw/duong2015.csv"

# --------------------------------------------------------------------
# DATA GENERATING MECHANISM ------------------------------------------
# --------------------------------------------------------------------

process.study <- function(file){
  # Read in data
  df <- fread(file)

  # Data pre-processing
  df <- df[!is.na(days)]

  ids <- df[, .(id1)]
  ids <- unique(ids)
  ids[, id.key := .I]

  df <- merge(df, ids, by=c("id1"))
  setorder(df, id.key, days)
  df <- df[, .(id.key, days)]
  df[, samp := 1:.N, by="id.key"]

  # Calculate gap days between samples
  df[, last.time := shift(days), by="id.key"]
  df[, gap := days - last.time]

  # Calculate total number of samples
  df[, num.samples := .N, by="id.key"]

  # Calculate first day
  df[, first.samp := lapply(.SD, min), .SDcols="days", by="id.key"]

  return(df)
}

# Raw Data Analysis --------------------------------------------
# Save for internal package functions
duong <- process.study(file=FILE)
cephia <- read.csv("data-raw/cephia_public_use_dataset_20210604.csv")
usethis::use_data(duong, cephia, internal=TRUE, overwrite=TRUE)
