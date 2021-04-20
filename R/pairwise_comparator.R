#' @export
#' @import igraph

pairwise_comparator = function(comp,pairwise,dict){
  # comp = x[[3]][,c(1,3,7)]
  comp2 = comp %>% filter(xor(.[[2]] == 0 , .[[3]] == 0 ))
  evidence = binary(comp2,80)
  all = binary(comp2,5)
  vi = colnames(comp)[2]
  vj = colnames(comp)[3]
  probability_i_j = pairwise %>% filter(.[[1]] == vj) %>% select(vi) %>% unlist()
  probability_j_i = pairwise %>% filter(.[[1]] == vi) %>% select(vj) %>% unlist()

  vi_f = evidence %>% filter(.[[2]] >0) %>% select(1) %>% filter_avarda(edge = dict)
  vj_f = evidence %>% filter(.[[3]] >0) %>% select(1) %>% filter_avarda(edge = dict)

  if(vi_f >= 3 & vj_f >= 3){

  all_i_f = comp %>% select(1) %>%
    filter(!(V1 %in% (comp %>% filter(.[[3]] >0) %>% select(1) %>% unlist))) %>%  #remove all v_j peptides xr and evidence
    filter(!(V1 %in% (comp %>% filter(.[[2]] >0 & .[[2]] <80) %>% select(1) %>% unlist))) %>%  #remove v_i xr peptides
    filter_avarda(edge = dict)
  all_j_f = comp %>% select(1) %>%
    filter(!(V1 %in% (comp %>% filter(.[[2]] >0) %>% select(1) %>% unlist))) %>%  #remove all v_i peptides xr and evidence
    filter(!(V1 %in% (comp %>% filter(.[[3]] >0 & .[[3]] <80) %>% select(1) %>% unlist))) %>%  #remove v_j xr peptides
    filter_avarda(edge = dict)

  prob_vi = binom.test(vi_f,all_i_f,probability_i_j)$p.value
  prob_vj = binom.test(vj_f,all_j_f,probability_j_i)$p.value
  output = list(prob_vi,prob_vj)
  return(output)
  }
  if(vi_f >= 3 & vj_f <3){
    all_i_f = comp %>% select(1) %>%
      filter(!(V1 %in% (comp %>% filter(.[[3]] >0) %>% select(1) %>% unlist))) %>%  #remove all v_j peptides xr and evidence
      filter(!(V1 %in% (comp %>% filter(.[[2]] >0 & .[[2]] <80) %>% select(1) %>% unlist))) %>%  #remove v_j xr peptides
      filter_avarda(edge = dict)
    prob_vi = binom.test(vi_f,all_i_f,probability_i_j)$p.value
    output = list(prob_vi,1)
    return(output)
  }
  if(vi_f < 3 & vj_f >=3){
    all_j_f = comp %>% select(1) %>%
      filter(!(V1 %in% (comp %>% filter(.[[2]] >0) %>% select(1) %>% unlist))) %>%  #remove all v_j peptides xr and evidence
      filter(!(V1 %in% (comp %>% filter(.[[3]] >0 & .[[3]] <80) %>% select(1) %>% unlist))) %>%  #remove v_j xr peptides
      filter_avarda(edge = dict)
    prob_vj = binom.test(vj_f,all_j_f,probability_j_i)$p.value
    output = list(1,prob_vj)
    return(output)
  }
  if(vi_f < 3 & vj_f <3){
    output = list(1,1)
    return(output)
  }
}
