#' @export
#' @import progress

total_probs_table = function(blast,dict){
  df = data.frame(matrix(ncol = 2,nrow = dim(blast)[2]))
  colnames(df) = c("virus","filtered")
  pb <- progress_bar$new(total = dim(blast)[2])
  blast.bi = blast %>% binary(threshold = 80)
  blast.bi = blast.bi %>% filter(rowSums(blast.bi[-1])>0)
  for(R in 2:dim(blast)[2]){
    x.R = blast.bi %>% select(1,R) %>%  filter(.[[2]] > 0) %>% select(1)
    print(c(R,dim(x.R)[1]))
    df[R,1] = colnames(blast)[R]
    df[R,2] = filter_avarda(edge = dict,vertex = x.R)
    pb$tick()
  }
  df$all = dim(blast.bi)[1]
  df$prob = df$filtered/df$all
  return(df)
}
