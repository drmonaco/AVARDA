#' @export
#' @import igraph

Module_3 = function(mod1,mod2,blast,total_prob,pairwise,dict){
mod2 = mod2 %>% arrange(pVal_f) %>% select(virus) %>% unlist %>% as.character()
mod1 = mod1[[3]] %>% select(1,mod2)
df = list()
R1 = 2
key.R = 1
while(R1<=dim(mod1)[2]){
  print(dim(mod1)[2])
  print(R1)
  vi = colnames(mod1)[R1]
  print(vi)
  if(R1 == dim(mod1)[2]){
    df[[key.R]] = mod1 %>% select(1,2)
    mod1 = mod1 %>% select(-vi)
  }
  if(R1 <= dim(mod1)[2]){
  for(R2 in (R1+1):dim(mod1)[2]){
    comp = pairwise_comparator(mod1 %>% select(c(1,R1,R2)),pairwise,dict)
    if(comp[[1]]<.05 & comp[[2]]<.05){  #if both significant

    }
    if(comp[[1]]<.05 & comp[[2]]>.05){ # if v_i is significant
    #index = which(mod1[,..R1] > 0)
    index = which(mod1[,R1] > 0)
    mod1[index,R2] = 0
    }
    if(comp[[1]]>.05 & comp[[2]]<.05){ # if v_j is significant
    #index = which(mod1[,..R2] > 0)
    index = which(mod1[,R2] > 0)
    mod1[index,R1] = 0

    }
    if(comp[[1]]>.05 & comp[[2]]>.05){ # if neither is significant

    }
  }

  rest = rowSums(select(mod1,-1, -R1))
  case = mod1 %>% select(R1)
  df[[key.R]] = mod1 %>% select(1,2)
  mod1 = mod1 %>% filter(!(case > 0 & rest == 0)) %>% select(-vi)
  f = Module_2(x =mod1 ,dict = dict,total = total) %>% select(virus) %>% unlist() %>% as.character()
  mod1 = mod1 %>% select(c(V1,f))
  #R1 = R1+1
  }
  key.R = key.R +1
  }
return(df)
}
