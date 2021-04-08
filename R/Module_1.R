Module_1 = function(case_data,blast,total_prob,dict){

case_data.R = case_data %>% filter(.[[2]] > 5) %>% left_join(blast) %>% select(-2) # threshold hits matrix and create blast matrix

case.binary = binary(case_data.R,80) # create binary matrix of evidence alignments

csums_index = which(case.binary[-1] %>% colSums()  > 5) + 1 # get just viruses that are enriched

case.binary.sub.evidence = case.binary %>% select(c(1,csums_index)) # subset evidence peptides above index - binary matrix of evidence peptides by alignnment

case.binary.sub.2 = (case.binary.sub.evidence [-1]%>% colSums()) %>% as.data.frame() %>% rename(peps = ".") %>% mutate(all_peptides = dim(case.binary.sub.evidence)[1]) %>%
  rownames_to_column() %>% left_join(total_prob)

  case.binary.sub.2 = case.binary.sub.2 %>% mutate(pVal = mapply(bt, case.binary.sub.2$peps, case.binary.sub.2$all_peptides,case.binary.sub.2$V1)) %>%
  arrange(pVal)

case.binary.sub.all = case_data.R %>% filter(.[[1]] %in% case.binary.sub.evidence[[1]]) # unbinaried data for all viruses

output = list(case.binary.sub.evidence,case.binary.sub.2,case.binary.sub.all)

return(output)
}


