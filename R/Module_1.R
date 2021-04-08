Module_1 = function(case_data,blast,total_prob,dict){

case_data.R = case_data %>% filter(.[[2]] > 5) %>% left_join(blast) %>% select(-2)

case.binary = binary(case_data.R,80)

csums_index = which(case.binary[-1] %>% colSums()  > 5) +1

case.binary.sub = case.binary %>% select(c(1,csums_index))

case.binary.sub.2 = (case.binary.sub [-1]%>% colSums()) %>% as.data.frame()

return(case.binary.sub)
}
