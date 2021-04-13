#' @export
#' @import igraph
# mod1 = x
# mod2 = x2
Module_3 = function(mod1,mod2,blast,total_prob,pairwise,dict){
mod2 = mod2 %>% arrange(pVal_f) %>% select(virus) %>% unlist %>% as.character()
mod1 = mod1[[3]] %>% select(1,mod2)
df = list()
for(R1 in 2:3){
  for(R2 in (R1+1):dim(mod1)[2]){
  print(R2)
  comp = pairwise_comparator(mod1 %>% select(c(1,R1,R2)),pairwise,dict)
  print(comp[[1]])
  print(comp[[2]])

  if(comp[[1]]<.05 & comp[[2]]<.05){  #if both significant

  }
  if(comp[[1]]<.05 & comp[[2]]>.05){ # if v_i is significant
  index = which(mod1[,..R1] > 0)
  mod1[index,R2] = 0


  }
  if(comp[[1]]>.05 & comp[[2]]<.05){ # if v_j is significant
  index = which(mod1[,..R2] > 0)
  mod1[index,R1] = 0

  }
  if(comp[[1]]>.05 & comp[[2]]>.05){ # if neither is significant

  }
  }
  rest = rowSums(select(mod1,-1, -R1))
  case = mod1 %>% select(R1)
  df[[R1]] = mod1[,c(1,R1)]
  mod1 = mod1 %>% filter(!(case > 0 & rest == 0))

}
}
