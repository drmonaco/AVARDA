#' @export

binary = function(matrix,threshold){
  g2 <- matrix %>% as.data.frame()%>% mutate_if(is.numeric, (funs(as.integer(./threshold)))) %>% as.data.frame()
  g2 <-  g2%>% mutate_if(is.numeric, (funs(ceiling((.)/(.+1))))) %>% as.data.frame()
  return(g2)
  }


