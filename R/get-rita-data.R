
rita.subtype <- function(){
  cephia <- XSRecency:::cephia
  cephia <- cephia %>% dplyr::filter(
    cephia_panel == "CEPHIA 1 Evaluation Panel",
    assay_result_field == "final_result")
  with(cephia, table(assay, hiv_subtype))
}

#' Get recency assay data for a phi dataset
#'
#' @param assay Recent infection testing algorithm to use (see `rita.subtype()` for options)
#' @param subtype Optional HIV subtypes as a vector (see `rita.subtype()` for options)
#' @returns A dataset with `ri` (recency indicator), `ui` (infection duration [in]) and `id` (participant ID)
#'
#' @import tidyr
#' @import dplyr
get.rita.data <- function(assay_name, threshold, subtype=NULL, ever_art=NULL){

  cephia <- XSRecency:::cephia %>%
    dplyr::filter(assay == assay_name,
           assay_result_field == "final_result",
           cephia_panel == "CEPHIA 1 Evaluation Panel") %>%
    tidyr::pivot_wider(names_from = assay_result_field, values_from = assay_result_value)

  if(!is.null(subtype)){
    cephia <- cephia %>%
      dplyr::filter(hiv_subtype %in% subtype)
  }

  cephia <- cephia %>%
    dplyr::select(
      days_since_eddi,
      final_result,
      viral_load_closest_to_visit,
      participant_identifier) %>%
    dplyr::rename(
      LAg=final_result,
      viral_load=viral_load_closest_to_visit,
      id=participant_identifier,
      ui=days_since_eddi) %>%

    # TODO: Build in process for other algorithms?
    # Add a cutoff value for the final_result column, then add optional viral load logic?
    dplyr::mutate(viral_load=replace(viral_load, is.na(viral_load), 0)) %>%
    dplyr::filter(!is.na(ui)) %>%
    dplyr::mutate(ri=as.numeric((LAg <= 1.5) & (viral_load > 1000)),
           ui=ui/365.25) %>%
    dplyr::select(id, ui, ri)

  return(cephia)
}
