#' @export

Module_2 = function(x,dict,total_prob){
  all = filter_avarda(edge = dict,vertex = x %>% select(1))
  df = data.frame(matrix(ncol = 2,nrow = dim(x)[2]))
  colnames(df) = c("virus","filtered")
  x = x %>% binary(threshold = 80)
  for(R in 2:dim(x)[2]){
  if(dim(x %>% filter(.[[R]] > 0))[1]>=3){
  x.R = x %>%  filter(.[[R]] > 0) %>% select(1) %>% unlist()
  df[R,1] = colnames(x)[R]
  df[R,2] = filter_avarda(edge = dict,vertex = x.R)
  }
  }
  df$all_f = all
  df = df %>% left_join(total_prob %>% rename(virus = "rowname"),by = "virus") %>% filter(!is.na(virus))%>% filter(filtered>2)

  df = df %>% filter(!is.null(virus)) %>% mutate(pVal_f = mapply(bt, df$filtered, df$all_f,df$V1)) %>% arrange(pVal_f)
  return(df)
}
