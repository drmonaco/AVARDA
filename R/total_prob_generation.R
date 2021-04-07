total_prob_generation = function(avarda,threshold = 80, col_cut = 1, row_cut = 1){
  avarda2 <- avarda[,-1]
  rownames(avarda2) <- avarda[,1] %>% unlist()
  avarda2[avarda2 <threshold] = 0
  avarda2[avarda2 >=threshold] = 1
  avarda2 = avarda2 %>% mutate(rsum = rowSums(.)) %>% filter(rsum>=row_cut) %>% select(-rsum)
  cSum = cbind(colSums(avarda2)) %>% as.data.frame() %>% rownames_to_column() %>% filter(V1 >= col_cut)
  cSum$V1 = cSum$V1/dim(avarda2)[1]
  return(cSum)
}
