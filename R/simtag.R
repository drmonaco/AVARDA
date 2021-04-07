simtag  = function(last){
  results = as.data.frame(last[1]) #take the virus of two
  results = cbind(results,0) #
  table = as.data.frame(last[2]) #
  names = subset(results[,1],results[,2]<.05)
  table = subset(table,table[,1] %in% names)
  table = subset(table,table[,2] %in% names)

  if(dim(table)[1]!=0){
    net <- as.undirected(graph_from_data_frame(table, directed=F) )
    max =  max_cliques(net,min = 2)
    for(R in 1:length(max)){
      names = induced_subgraph(net,max[[R]])
      sim = as.data.frame(vertex_attr(names))
      index = match(as.character(unlist(sim)),results[,1])
      results[index,12] = paste(results[index,12],R,sep = "|")
    }
  }
  return(results)
}
