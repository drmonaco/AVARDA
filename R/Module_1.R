#' Perform AVARDA Module 1 Analysis - Accounting for peptide - virus xreactivity
#'
#' @param case_data A data.frame
#' @param blast A data.frame
#' @param total_prob A data.frame
#' @param threshold A number
#' @return Returns a three element. Element 1 is matrix of all peptide enriched above \code{threshold} in \code{case_data} showing all alignments defined in \code{blast}, thresholded by a predetermined value. Element 2 shows a list of viruses from \code{blast} with potential reactivity. Element 3 is similar to element 1 but not thresholded.
#' @examples
#' \dontrun{}
#' @export
#' @import dplyr
#' @import parallel
#' @import tibble

Module_1 = function(case_data,blast,total_prob,threshold = 5){

case_data.R = case_data %>% filter(.[[2]] >= threshold) %>% mutate(V1 = as.character(V1))%>% inner_join(blast,by = "V1") %>% select(-2) # threshold hits matrix and create blast matrix

case.binary = binary(case_data.R,80) # create binary matrix of evidence alignments

csums_index = which(case.binary[-1] %>% colSums()  > 3) + 1 # get just viruses that are enriched
# think about thresholding
case.binary.sub.evidence = case.binary %>% select(c(1,csums_index)) # subset evidence peptides above index - binary matrix of evidence peptides by alignnment

case.binary.sub.2 = (case.binary.sub.evidence [-1]%>% colSums()) %>% as.data.frame() %>% rename(peps = ".") %>%
  rownames_to_column() %>% mutate(all_peptides = dim(case.binary.sub.evidence)[1])  %>% left_join(total_prob, "rowname")

  case.binary.sub.2 = case.binary.sub.2 %>% mutate(pVal = mapply(bt, case.binary.sub.2$peps, case.binary.sub.2$all_peptides,case.binary.sub.2$V1)) %>%
  arrange(pVal)

case.sub.all = case_data.R %>% filter(.[[1]] %in% case.binary.sub.evidence[[1]]) %>% select(c(1,csums_index))

output = list(case.binary.sub.evidence,case.binary.sub.2,case.sub.all)

return(output)
}


