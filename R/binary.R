#' @export

binary = function(matrix,threshold){
  g2 <- matrix %>% as.data.frame()%>% mutate_if(is.numeric, (funs(as.integer(./threshold)))) %>% as.data.frame()
  return(g2)
  }


