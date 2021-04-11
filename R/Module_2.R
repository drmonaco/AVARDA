#' @export

dict= fread("~/Desktop/bin2/Phageome/blast_dictionary_phageome.csv")

Module_2 = function(x,dict,total){
  all = filter_avarda(edge = dict,vertex = x[[1]][[1]])
  df = data.frame(matrix(ncol = 2,nrow = dim(x[[1]])[2]))
  colnames(df) = c("virus","filtered")
  for(R in 2:dim(x[[1]])[2]){
  x.R = x[[1]] %>% filter(.[[R]] > 0)
  df[R,1] = colnames(x[[1]])[R]
  df[R,2] = filter_avarda(edge = dict,vertex = x.R[[1]])
  }
  df$all = all
  df = df %>% left_join(total %>% rename(virus = "rowname")) %>% filter(!is.na(virus))

  df = df %>% mutate(pVal = mapply(bt, df$filtered, df$all_filtered,df$V1)) %>%
    arrange(pVal)

  return(df)
}


