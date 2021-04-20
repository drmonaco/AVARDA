#' Perform a binomial test
#'
#' @param case_data A data.frame
#' @param blast A data.frame
#' @param total_prob A data.frame
#' @param threshold A number
#' @return Returns a three element. Element 1 is matrix of all peptide enriched above \code{threshold} in \code{case_data} showing all alignments defined in \code{blast}, thresholded by a predetermined value. Element 2 shows a list of viruses from \code{blast} with potential reactivity. Element 3 is similar to element 1 but not thresholded.
#' @examples
#' \dontrun{}
#' @export
#' @import dplyr
#' @import parallel
#' @import tibble
#' @import igraph
#'

filter_avarda2  = function(edge,vertex){ #independence filter that takes a dictionary (defined above) and a set of nodes and tells the minimal number of unique epitopes
  nodes = unlist(vertex)
  links_filtered = edge %>% filter(V1 %in% nodes)
  links_filtered = links_filtered %>% filter(V2 %in% nodes)
  x_1_sum = length(nodes)
  if(dim(links_filtered)[1]!=0){
    net <- simplify(as.undirected(graph_from_data_frame(d=links_filtered,vertices=nodes, directed=F) ))
    x = decompose.graph(net)
    x_1 = x[sapply(x,vcount)<30]
    x_1_sum  = sum(unlist(lapply(x_1,independence.number)))
    x_2 = x[sapply(x,vcount)>=30]
    temp = c()
    #x_2 = x
    if(length(x_2) >0){
      for(R in 1:length(x_2)){
        x_2_r = x_2[[R]]
        while(max(degree(x_2_r))>10){
          toss = degree(x_2_r)==max(degree(x_2_r))
          drop = V(x_2_r)[which(toss == TRUE)]
          x_2_r = delete_vertices(x_2_r, drop)
        }
        while(max(degree(x_2_r))>0){
          toss = degree(x_2_r)==max(degree(x_2_r))
          drop = V(x_2_r)[which(toss == TRUE)[1]]
          x_2_r = delete_vertices(x_2_r, drop)
        }
        #x_l = decompose.graph(x_2_r)
        temp[R] = vcount(x_2_r)
      }
    }
    #print(sum(x_1_sum)+sum(temp))
    return(sum(x_1_sum)+sum(temp))
  }
  if(dim(links_filtered)[1]==0){
    return(length(nodes))
  }

}
