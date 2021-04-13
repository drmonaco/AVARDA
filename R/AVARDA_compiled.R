#' @export
#' @import doParallel
#' @import foreach

AVARDA_compiled = function(case_data,blast,total_prob,pairwise,dict,threshold = 5){
  registerDoParallel(  registerDoParallel(detectCores()))
  zeta = foreach(R = 2:(dim(case_data)[2]),.combine=rbind) %dopar%{ #cycle through each patient column by column (goal is so be serialized)
    data.R = data2 %>% select(1,R)
    colnames(data.R)[2]
    x = Module_1(data.R,blast,total,threshold = 10)
    print("Module1_finished")
    x2 = Module_2(x =x[[3]] ,dict = dict,total = total)
    print("Module2_finished")

    x3 = Module_3(x,x2,blast,total,pairwise,dict)
    print("Module3_finished")

    x4 =sapply(x3,compiler,blast = blast,total_prob = total,dict = dict) %>% t() %>% as.data.frame() %>% mutate(pBH = p.adjust(pVal))
    x4 = x4 %>% mutate(Sample_ID = colnames(data.R)[2])
    return(x4)
    }
  return(zeta)
}
