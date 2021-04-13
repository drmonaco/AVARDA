#' @export
#' @import igraph

filter_avarda  = function(edge,vertex){ #independence filter that takes a dictionary (defined above) and a set of nodes and tells the minimal number of unique epitopes
  nodes = unlist(vertex)
  links_filtered = edge %>% filter(V1 %in% nodes)
  links_filtered = links_filtered %>% filter(V2 %in% nodes)
  x_1_sum = length(nodes)
  if(dim(links_filtered)[1]!=0){
    net <- simplify(as.undirected(graph_from_data_frame(d=links_filtered,vertices=nodes, directed=F) ))
    x = decompose.graph(net)
    x_1_sum = sum(unlist(lapply(x,independence.number)))
  }
  return(x_1_sum)
}
