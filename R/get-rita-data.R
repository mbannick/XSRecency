
#' Download latest CEPHIA dataset
#'
#' Downloads the latest CEPHIA dataset from concept DOI
#' 10.5281/zenodo.4900633, and returns a local filepath where it can be found.
#' This function can be used to download the latest version, and use it in the
#' [getRitaOptions()] and [createRitaCephia()] functions.
#'
#' @param path A folder to store the CEPHIA dataset
#' @param timeout The amount of time to allow R to download file
#' @export
#' @import zen4R
#' @return A full filepath to the downloaded CEPHIA dataset
downloadCephia <- function(path, timeout=6000){

  # Set timeout, otherwise the dataset does not download
  options(timeout=timeout)

  cat("Getting latest CEPHIA dataset using
      the concept DOI 10.5281/zenodo.4900633...")

  zenodo <- ZenodoManager$new()
  rec <- zenodo$getRecordByConceptDOI("10.5281/zenodo.4900633")
  lastDOI <- rec$getLastDOI()
  rec <- zenodo$getRecordByDOI(lastDOI)

  files <- sapply(rec$files, function(x) x$filename)
  files <- files[grep("^cephia_public_use_dataset", files)]
  csvfile <- files[grep(".csv$", files)]

  if(length(csvfile) > 1) stop("Too many files found.")
  if(length(csvfile) < 1) stop("No csv files found with correct name.")

  rec$downloadFiles(files=csvfile, path=path)

  return(paste0(path, "/", csvfile))
}

get.cephia.use <- function(filepath=NULL){

  if(is.null(filepath)){
    cephia <- XSRecency:::cephia
  } else {
    cephia <- read.csv(filepath)
  }

  return(cephia %>%
           dplyr::filter(
             cephia_panel == "CEPHIA 1 Evaluation Panel",
             assay_result_field == "final_result",
             hiv_subtype %in% c("A1", "B", "C", "D"),
             !is.na(days_since_eddi)))
}

#' See which recency assays are available by HIV subtype
#'
#' @import tidyr
#' @import dplyr
#'
#' @param filepath Optional full filepath to downloaded CEPHIA dataset (possibly obtained via the `downloadCephia` function).
#'                 If NULL, then uses CEPHIA version in this package.
#'
#' @export
#' @examples
#' getRitaOptions()
getRitaOptions <- function(filepath=NULL){
  cephia <- get.cephia.use(filepath=filepath)
  cephia <- cephia %>% dplyr::filter(
    cephia_panel == "CEPHIA 1 Evaluation Panel",
    assay_result_field == "final_result",
    hiv_subtype %in% c("A1", "B", "C", "D"),
    !is.na(days_since_eddi))
  with(cephia, table(assay, hiv_subtype))
}

#' Create a dataset based on an algorithm and CEPHIA recency testing data.
#'
#' Get data frame using CEPHIA data for a given recency
#' algorithm specified by the user. Optional subtype and ART exclusions.
#'
#' @param assays A vector of assays to include
#' @param algorithm A function that defines the recency indicator
#'                  with arguments in the same order as the `assays` vector.
#'                  Arguments do not need to have the same name as `assays`.
#'                  E.g., if you have `assays = c("BED", "viral_load")`,
#'                  you can have `algorithm = function(b, v) ...` where `b`
#'                  indicates BED and `v` indicates viral load.
#' @param subtype HIV subtypes to include (one of "A1", "B", "C", "D").
#'                By default includes everyone.
#' @param ever_art Subset data to only those who have used ARTs or
#'                 have not. By default includes everyone.
#' @param filepath Optional full filepath to downloaded CEPHIA dataset (possibly obtained via the `download_cephia`).
#'                 If NULL, then uses CEPHIA version in this package.
#'
#' @return A data frame with the following columns:
#' \item{id}{ID column (there can be multiple measurements per individual)}
#' \item{ui}{infection duration in **years**
#' \item{ri}{identified as recent based on algorithm}
#'
#' @export
#' @import tidyr
#' @import dplyr
#' @import purrr
#' @importFrom data.table data.table setnames
#'
#' @examples
#' f <- function(b, l, v){
#'   ifelse((b > 1 & l < 3 & !is.na(v)), 1, 0)
#' }
#' createRitaCephia(assays=c("BED", "LAg-Sedia", "viral_load"),
#'              algorithm=f)
#' f <- function(l, v){
#'   v <- ifelse(l > 1.5, 0, v)
#'   return(
#'     ifelse((l <= 1.5) & (v > 1000), 1, 0)
#'   )
#' }
#' test <- createRitaCephia(assays=c("LAg-Sedia", "viral_load"), algorithm=f)
createRitaCephia <- function(assays, algorithm, subtype=NULL, ever_art=NULL, filepath=NULL){
  cephia <- get.cephia.use(filepath=filepath)

  # APPLY FILTERS
  if(!is.null(subtype)){
    cephia <- cephia %>%
      dplyr::filter(hiv_subtype %in% subtype)
  }
  if(!is.null(ever_art)){
    if(ever_art == 1){
      cephia <- cephia %>%
        dplyr::filter(!is.na(days_since_first_art))
    } else {
      cephia <- cephia %>%
        dplyr::filter(is.na(days_since_first_art))
    }
  }

  # FUNCTION TO GET A SINGLE ASSAY DATA FRAME BASED
  # ON LATEST TESTING DATE IF DUPES
  get.assay.df <- function(assay_name){

    if(assay_name == "viral_load"){
      assay.1 <- cephia %>%
        select(participant_identifier, visit_identifier, viral_load_closest_to_visit) %>%
        unique() %>%
        rename(viral_load=viral_load_closest_to_visit)
      num.missing <- nrow(assay.1 %>% filter(is.na(viral_load)))
    } else if(assay_name == "cd4"){
      assay.1 <- cephia %>%
        select(participant_identifier, visit_identifier, cd4_count_at_visit) %>%
        unique() %>%
        rename(cd4=cd4_count_at_visit)
      num.missing <- nrow(assay.1 %>% filter(is.na(cd4)))
    } else {
      assay.1 <- cephia %>%
        filter(assay == assay_name)
      assay.1 <- assay.1 %>%
        group_by(participant_identifier, visit_identifier) %>%
        summarize(assay_result_value=median(assay_result_value)) %>%
        ungroup()

      assay.1 <- assay.1 %>%
        select(participant_identifier, visit_identifier, assay_result_value)

      num.missing <- nrow(assay.1 %>% filter(is.na(assay_name)))
      setnames(assay.1, "assay_result_value", assay_name)
    }

    if(num.missing > 0) cat("There are", num.missing, "missing values for", assay_name, "\n")

    return(assay.1)
  }

  # GET LIST OF ASSAY DATA FRAMES AND MERGE THEM
  assay_dfs <- lapply(assays, get.assay.df)
  all_assays <- assay_dfs %>%
    reduce(full_join, by=c("participant_identifier", "visit_identifier"))

  # GET INFECTION DURATION
  cephia.dur <- cephia %>%
    select(participant_identifier, visit_identifier, days_since_eddi) %>%
    unique()

  # MERGE INFECTION DURATION AND ASSAY VALUES
  df <- full_join(cephia.dur, all_assays,
                  by=c("participant_identifier", "visit_identifier")) %>%
    rename(ui=days_since_eddi)
  df <- data.table::as.data.table(data)

  # APPLY ALGORITHMS
  setnames(df, assays, formalArgs(algorithm))
  arglist <- as.list(df[, formalArgs(algorithm), with=F])
  applyit <- function(...) mapply(algorithm, ...)
  df$ri <- do.call(applyit, arglist)

  # CLEAN FINAL DATA FRAME
  final_df <- df %>%
    select(participant_identifier, ui, ri) %>%
    rename(id=participant_identifier)

  # NUMBER OF MISSING RI
  num.missing <- nrow(final_df %>% filter(is.na(ri)))
  if(num.missing > 0) cat(
    "Removing", num.missing,
    "observations with missing recency indicator after application of the algorithm.")

  # CONVERT TIME UNIT TO YEARS, SINCE FROM THE CEPHIA DATA IT IS IN DAYS
  final_df$ui <- final_df$ui / 365.25
  
  return(final_df)
}
