

#' @import tidyr
#' @import dplyr
rita.subtype <- function(){
  cephia <- XSRecency:::cephia
  cephia <- cephia %>% dplyr::filter(
    cephia_panel == "CEPHIA 1 Evaluation Panel",
    assay_result_field == "final_result",
    hiv_subtype %in% c("A1", "B", "C", "D"),
    !is.na(days_since_eddi))
  with(cephia, table(assay, hiv_subtype))
}

get.cephia.use <- function(){
  return(XSRecency:::cephia %>%
           dplyr::filter(
             cephia_panel == "CEPHIA 1 Evaluation Panel",
             assay_result_field == "final_result",
             hiv_subtype %in% c("A1", "B", "C", "D"),
             !is.na(days_since_eddi)))
}


#' Get external data frame using CEPHIA data for a given recency
#' algorithm specified by the user.
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
#' @export
#' @import tidyr
#' @import dplyr
#' @import purrr
#'
#' @examples
#' f <- function(b, l, v){
#'   ifelse((b > 1 & l < 3 & !is.na(v)), 1, 0)
#' }
#' get.assay.df(assays=c("BED", "LAg-Sedia", "viral_load"),
#'              algorithm=f)
#'
#' f <- function(l, v){
#'   v <- ifelse(l > 1.5, 0, v)
#'   return(
#'     ifelse((l <= 1.5) & (v > 100), 1, 0)
#'   )
#' }
#' test <- get.assay.df(assays=c("LAg-Sedia", "viral_load"), algorithm=f)
get.assay.df <- function(assays, algorithm, subtype=NULL, ever_art=NULL){
  cephia <- get.cephia.use()

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
      assay.1.latest <- assay.1 %>%
        group_by(participant_identifier, visit_identifier) %>%
        summarize(test_date=max(test_date)) %>%
        ungroup()

      assay.1 <- merge(assay.1, assay.1.latest,
                       by=c("participant_identifier", "visit_identifier", "test_date"))
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

  return(final_df)
}
