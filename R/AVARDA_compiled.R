#' Primary function to run AVARDA modules together
#' @import doParallel
#' @import foreach
#' @import dplyr
#' @import parallel
#' @import tibble
#' @import tictoc
#' @param case_data A data.frame of PhIPseq data that identifies hits by some theshold, can be binary hits or some value
#' @param blast A data.frame listing the blastp, or tblastn, alignments of PhIP-seq peptides to organisms of interact
#' @param total_prob A data.frame listing the degree of representation each virus has in the phip-seq library
#' @param pairwise A number.
#' @param dict A number.
#' @param threshold A number.
#' @param cores A number.
#' @param mod_3 A number.
#' @return Output of AVARDA algorithm; a data frame listing the likelihood of each virus causing in infection in a given sample
#' @examples
#' \dontrun{do later}
#' @export

AVARDA_compiled = function(case_data,blast,total_prob,pairwise,dict,threshold = 5,cores,mod_3){
  if(cores == "all"){
  registerDoParallel(parallel::detectCores())
  }
  registerDoParallel(cores)
  zeta = foreach(R = 2:(dim(case_data)[2]),.combine=rbind,.packages = c("dplyr","AVARDA")) %dopar%{ #cycle through each patient column by column (goal is so be serialized)
  #zeta = foreach(R = 2:3,.combine=rbind) %dopar%{ #cycle through each patient column by column (goal is so be serialized)
    #print(R)
    data.R = case_data %>% select(1,R)
    # tic("Module1_finished")
    x = Module_1(case_data = data.R,blast = blast,total_prob = total_prob,threshold = threshold)
    # toc()
    if(dim(x[[2]])[1] != 0){
      # tic("Module2_finished")
      x2 = Module_2(x =x[[3]] ,dict = dict,total_prob = total_prob)
      # toc()
      if(dim(x2)[1]!= 0){
      #tic("Module3_finished")
      x3 = Module_3(mod1 = x,mod2 = x2,blast = blast,total_prob = total_prob,pairwise = pairwise,dict = dict,mod_3)
      # toc()

      x4 =sapply(x3,compiler,blast = blast,total_prob = total_prob,dict = dict) %>% t() %>% as.data.frame() %>% mutate(pBH = p.adjust(pVal))
      x4 = x4 %>% mutate(Sample_ID = colnames(data.R)[2])
      return(x4)
      }
      if(dim(x2)[1]== 0){
        return(NULL)
      }
      }
    if(!is.null(dim(x)[[3]])){
      return(NULL)
      }
  }
  return(zeta)
}
