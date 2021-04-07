filter.R  = function(edge,vertex){ #independence filter that takes a dictionary (defined above) and a set of nodes and tells the minimal number of unique epitopes
  nodes = unlist(vertex)
  links_filtered = subset(edge,unlist(edge[,1]) %in% nodes)
  links_filtered = subset(links_filtered,links_filtered[,2] %in% nodes)
  if(dim(links_filtered)[1]!=0){

    net <- as.undirected(graph_from_data_frame(d=links_filtered,vertices=nodes, directed=F) )
    x = decompose.graph(net)
    x_1 = x[sapply(x,vcount)<30]
    x_1_sum  = sum(unlist(lapply(x_1,independence.number)))
    x_2 = x[sapply(x,vcount)>=30]
    temp = c()
    #x_2 = x
    if(length(x_2) >0){
      for(R in 1:length(x_2)){
        x_2_r = x_2[[R]]
        while(max(degree(x_2_r)>5)){

          toss = degree(x_2_r)==max(degree(x_2_r))
          x_2_r = delete_vertices(x_2_r, V(x_2_r)[toss])
        }
        x_l = decompose.graph(x_2_r)
        temp[R] = sum(unlist(lapply(x_l,independence.number)))
      }
    }
    return(sum(x_1_sum)+sum(temp))
  }
  if(dim(links_filtered)[1]==0){
    return(length(nodes))
  }
}
