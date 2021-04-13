#' @export

dict= fread("~/Desktop/bin2/Phageome/blast_dictionary_phageome.csv")

Module_2 = function(x,dict,total){
  all = filter_avarda(edge = dict,vertex = x %>% select(1))
  df = data.frame(matrix(ncol = 2,nrow = dim(x)[2]))
  colnames(df) = c("virus","filtered")
  for(R in 2:dim(x)[2]){
  x.R = x %>% binary(threshold = 80) %>%  filter(.[[R]] > 0) %>% select(1)
  df[R,1] = colnames(x)[R]
  df[R,2] = filter_avarda(edge = dict,vertex = x.R)
  }
  df$all_f = all
  df = df %>% left_join(total %>% rename(virus = "rowname")) %>% filter(!is.na(virus))

  df = df %>% mutate(pVal_f = mapply(bt, df$filtered, df$all_f,df$V1))

  df %>% left_join(x[2] %>% as.data.frame()%>% rename(virus = rowname) %>% select(-V1)) %>% arrange(pVal_f)
  return(df)
}


