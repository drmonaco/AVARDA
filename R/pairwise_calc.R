pairwise_calc  = function(total,pairwise,rank1,dict){
  ## just initializing some variables
  unique_probs =pairwise
  total_probs = total
  fullmatrix_sorted_ij = rank1
  final_matrix = data.frame(matrix(0,nrow = dim(fullmatrix_sorted_ij)[1],ncol =dim(fullmatrix_sorted_ij)[2])) #initialize the final matrix of peptide-virus alignments post reassignment
  N_rank = rownames(fullmatrix_sorted_ij)
  z = N_rank
  output = data.frame(matrix(ncol = 10, nrow =  dim(fullmatrix_sorted_ij)[2])) #initialize matrix for final data with significance data
  sim_tag = data.frame(matrix(ncol = 2)) #empty vector to fill with virus pairs that are sim-tagged
  x1 = 1
  ##

  pb <- txtProgressBar(min = 0, max = dim(fullmatrix_sorted_ij)[2], style = 3)
  for(R1 in 1: dim(fullmatrix_sorted_ij)[2]){ # test
    # for(R1 in 1: 13){ # test

    virus_i = colnames(fullmatrix_sorted_ij)[R1] # get name of virus_i
    virus_i_hits = as.data.frame(fullmatrix_sorted_ij[R1]) #get vector of virus_i peptide alignments
    v_i = subset(virus_i_hits,virus_i_hits>0) # get index of all peptide alignments >0
    v_i_xr = subset(virus_i_hits,virus_i_hits>0 & virus_i_hits<80) # get index of all peptide alignments <80 and >0 (the xr)
    if(sum(virus_i_hits>=80)>=3){
      R2 = 1
      while(R2 <= dim(fullmatrix_sorted_ij)[2]){
        if(R1 < R2){ #skip iterations of virus_i comparing to previously evaluated viruses
          binary_i_j = ifelse(fullmatrix_sorted_ij>0, 1, 0)
          shared = binary_i_j[,R1]*binary_i_j[,R2] # vector multiplication to get the peptides with alignments to virus_i and virus_j
          virus_j_hits = as.data.frame(fullmatrix_sorted_ij[R2]) #get vector of virus_j alignments to enriched peptides
          virus_j = colnames(fullmatrix_sorted_ij)[R2] #get name of virus_j
          if(sum(shared) != 0){ #if both viruses are completely unique no reassignment can occur so skip this
            probability_i_j = unique_probs[grep(paste0(virus_j,"$",collapse = ""),row.names(unique_probs)),grep(paste0(virus_i,"$",collapse = ""),row.names(unique_probs))] #probabiltiy that a random peptide will align exclusively to virus_i relative to virus_j
            probability_j_i = unique_probs[grep(paste0(virus_i,"$",collapse = ""),row.names(unique_probs)),grep(paste0(virus_j,"$",collapse = ""),row.names(unique_probs))]

            v_j = subset(virus_j_hits,virus_j_hits>0) #get all v_j alignments
            v_j_xr = subset(virus_j_hits,virus_j_hits>0 & virus_j_hits<80) # get v_j xr alignments

            p_i_j =  binom_test(rownames(v_i),rownames(v_i_xr),rownames(v_j),N_rank,probability_i_j,dict) #likelihood that virus_i has an infection that is unique relative to virus_j
            p_j_i =  binom_test(rownames(v_j),rownames(v_j_xr),rownames(v_i),N_rank,probability_j_i,dict) #likelihood that virus_j has an infection that is unique relative to virus_i

            if(p_i_j[1] <= .05 & p_j_i[1] >.05 & p_i_j[2] >=3 | p_j_i[2] <3){ #if virus_j loses (and virus_i has at least 3 aligning peptideds
              fullmatrix_sorted_ij[(shared==1),R2] = 0 #virus_j loses all alignments that are shared with virus_i (even if they are xr in v_i and evidence in v_j)

            }
            if(p_i_j[1] > .05 & p_j_i[1] <=.05 & p_j_i[2] >=3 ){ #if virus_i loses
              fullmatrix_sorted_ij[(shared==1),R1] = 0

            }
            if(p_i_j[1] > .05 & p_j_i[1] > .05 | (p_j_i[2] <3 & p_i_j[2] < 3 )){ #if both viruses fail to have enough unique evidence (but both have three peptides)
              sim_tag[x1,] = cbind(virus_i,virus_j) #note that the pair is indistinguishable
              #print(x1)
              x1 = x1+1
            }
          }
        }
        R2= R2+1
      }

      ###at this point all virus_j have been compared
      order = c() #initialize
      virus = c() #initialize
      if(R1<dim(fullmatrix_sorted_ij)[2]){
        z =N_rank
        virus_i_prob = binom_test(rownames(v_i),rownames(v_i_xr),NULL,z,total_probs[grep(paste0(colnames(fullmatrix_sorted_ij[R1]),"$",collapse = ""),unlist(total_probs[,1])),2],dict)

        if(virus_i_prob[1]<=.05 & is.null(dim(fullmatrix_sorted_ij[which(fullmatrix_sorted_ij[R1]>0),(R1+1):dim(fullmatrix_sorted_ij)[2]]))==FALSE){ ##if virus_i was significant then remove all peptides that aligned exclusively to virus_i from future calculations
          N_remove = names(which(apply(fullmatrix_sorted_ij[which(fullmatrix_sorted_ij[R1]>0),(R1+1):dim(fullmatrix_sorted_ij)[2]],1,max)==0))
          z =N_rank[!N_rank %in% N_remove]
        }

        ### Need to rerank remaining viruses after the v_i_j comparisons
        for(R in (R1+1):dim(fullmatrix_sorted_ij)[2]){ # go iteartively through viruses
          virus_i_xr2 = rownames(fullmatrix_sorted_ij)[fullmatrix_sorted_ij[,R]>0 & fullmatrix_sorted_ij[,R]<80] # get all peptides total with alignments to virus x
          v_i2 = rownames(fullmatrix_sorted_ij)[fullmatrix_sorted_ij[,R]>0] # get peptides that are evidence for virus x
          probability = total_probs[grep(paste0(colnames(fullmatrix_sorted_ij[R]),"$",collapse = ""),unlist(total_probs[,1])),2] # call the null probability for virus x

          if(length(v_i2[!v_i2 %in% virus_i_xr2])[1] >0){ #if there is at least 1 evidence peptide for virus_x reorder
            order[R-R1] = binom.test(length(v_i2[!v_i2 %in% virus_i_xr2])[1],length(z),as.numeric(probability),"greater")[[3]]
          }
          if(length(v_i2[!v_i2 %in% virus_i_xr2])[1] ==0){ #otherwise just skip
            order[R-R1] = 1
          }
        }

        fullmatrix_sorted_ij = cbind(fullmatrix_sorted_ij[1:R1],fullmatrix_sorted_ij[order(order)+R1])
        print(dim(fullmatrix_sorted_ij))
      }
      ###


    }
    ## at this point all viruses are evaluated so just need final ranking step
    ## the below data is meta data we output when investigating data
    final = fullmatrix_sorted_ij[R1]
    final_i = subset(final,final>0)
    final_xr = subset(final,final>0 & final<80)
    probability_i = total_probs[grep(paste0(colnames(fullmatrix_sorted_ij[R1]),"$",collapse = ""),unlist(total_probs[,1])),2] # call the null probability for virus x
    final_rank = binom_test(rownames(final_i),rownames(final_xr),0,N_rank,probability_i,dict)
    final_i_e = subset(final,final>=80)
    output[R1,] = c(virus_i,final_rank[1],paste(rownames(final_i_e),collapse = "|"),paste(rownames(final_xr),collapse = "|"),length(N_rank),length(rownames(final_i_e)),length(rownames(final_xr)),final_rank[2],final_rank[3],probability_i)

    N_rank = z    ## set the N_rank for the next virus_i to be reduced by those assigned to the previous virus_i
    setTxtProgressBar(pb, R1)

  }
  close(pb)
  index = which(output[,2]!=1)
  a123 = cbind(output,1)
  a123[index,11] = p.adjust(output[index,2],"BH")
  last = list(a123,sim_tag,fullmatrix_sorted_ij)
  return(last)
}
