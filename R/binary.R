#' Take a matrix and binary it
#'
#' @param matrix A matrix
#' @param threshold A number.
#' @return The binaried version of \code{matrix} setting things above \code{threshold} to 1 and below to zero.
#' @examples
#' binary(matrix(c(100,50,0,100),nrow = 2), 70)
#' @export

binary = function(matrix,threshold){
  g2 <- matrix %>% as.data.frame()%>% mutate_if(is.numeric, (funs(as.integer(./threshold)))) %>% as.data.frame()
  g2 <-  g2%>% mutate_if(is.numeric, (funs(ceiling((.)/(.+1))))) %>% as.data.frame()
  return(g2)
  }


