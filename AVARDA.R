## Libraries required to run

library(tidyverse)
library(data.table)
library(igraph)
library(doParallel)
library(foreach)
#

AVARDA = function(case_path,thresh,dict_path,total_path,pairwise_path,blast_path,out_path,out_name){
  # dict =  data.frame(fread("unzip -cq ./bin/df.csv.zip",data.table = FALSE)) # read in peptide-peptide dictionary
  # total = data.frame(fread("unzip -cq ./bin/total_probability_xr.csv",data.table = FALSE)) #read in null probability table
  # pairwise = data.frame(fread("unzip -cq ./bin/unique_probabilities.csv.zip",data.table = FALSE),row.names=1) #read in null probability table
  # blast = data.frame(fread("unzip -cq ./bin/VirScan_filtered_virus_blast.csv.zip",data.table = FALSE),row.names=1) #read in virus-peptide alignment filtered table
  # case = data.frame(fread("./input/test1234.txt", data.table = FALSE)) # this is the case of interest to be analyzed - changeable

  dict =  data.frame(fread(dict_path,data.table = FALSE)) # read in peptide-peptide dictionary
  total = data.frame(fread(total_path,data.table = FALSE)) #read in null probability table
  pairwise = data.frame(fread(pairwise_path,data.table = FALSE),row.names=1) #read in null probability table
  blast = data.frame(fread(blast_path,data.table = FALSE),row.names=1) #read in virus-peptide alignment filtered table
  case = data.frame(fread(case_path, data.table = FALSE,check.names = FALSE)) # this is the case of interest to be analyzed - changeable

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

  binom_test  = function(v_i,v_xr,v_i_j,N_rank,null_prob){
    v_total = v_i[!v_i %in% v_xr] # all virus i aligning minus the xr
    v_total = v_total[!v_total %in% v_i_j] # all virus i evidence minus any shared with virus j
    v_total_f = filter(dict,v_total) # this calculates N_rank
    #if(length(N_rank)!=length(unlist(N_rank_2))){
    N_rank = N_rank[!N_rank %in% v_xr]
    N_rank = N_rank[!N_rank %in% (v_i %in% v_i_j)] # this is by default zero for total binom calculation
    N_rank_f = filter(dict,N_rank)
    #}
    if(N_rank_f == 0){
      return(NULL)
    }
    x = binom.test(v_total_f,N_rank_f,unlist(null_prob),"greater")[[3]]
    output = list(x,v_total_f,N_rank_f)
    return(output)
  }

  total_calc  = function(case,column,thresh,total,blast){ # subsetting function that also ranks peptides
    ## Processing of the blast virus matrix on a case by case basic
    total_probs = total
    enriched = data.frame(case[,c(1,column+1)],row.names = 1) # get column of from data frame of patient z-scores
    enriched = subset(enriched,enriched >=thresh) # subset on enrichment matrix peptides that are greater than >= 10
    print(colnames(enriched))
    if(dim(enriched)[1]>1600){
      enriched = enriched[order(-enriched),1,drop=FALSE]
      enriched = enriched[1:1600,1,drop = FALSE]
    }
    if(dim(enriched)[1] > 0){
      blast_subset = subset(blast,row.names(blast) %in% row.names(enriched)) %>% as.data.frame() # subset fullmatrix to only include peptides above threshold
      blast_subset=blast_subset[, colSums(ifelse(blast_subset>80, 1, 0)) > 2,drop=FALSE] # subset full matrix to only include viruses with at least three enriched peptide
      if(is.data.frame(blast_subset)==TRUE){
        if(dim(blast_subset)[2]>0){
          ###
          fullmatrix_sorted = blast_subset
          v_i_j = NULL
          order = c()
          virus = c()
          ###

          pb <- txtProgressBar(min = 0, max =dim(fullmatrix_sorted)[2], style = 3)

          for(R in 1:dim(fullmatrix_sorted)[2]){ # go iteratively through viruses
            virus_i_xr = rownames(fullmatrix_sorted)[fullmatrix_sorted[,R]>0 & fullmatrix_sorted[,R]<80] # get peptides that align crossreactively to virus_i
            v_i = rownames(fullmatrix_sorted)[fullmatrix_sorted[,R]>0] # get all peptides align to virus_i
            probability = total_probs[grep(paste0(colnames(fullmatrix_sorted[R]),"$",collapse = ""),unlist(total_probs[,1])),2] # call the null probability for virus_i
            x = binom_test(v_i,virus_i_xr,v_i_j,row.names(blast_subset),probability) # run binom function that takes an into account xr and shared peptides (which is 0 for this step)
            order[R] = x[1] # reports the p-value for a viruses likelihood of infection
            virus[R] = colnames(fullmatrix_sorted[R]) # get which virus is being compared in interation
            setTxtProgressBar(pb, R)

          }

          results = fullmatrix_sorted[,order(unlist(order)),drop = F] #reorder the subset virus-peptide blast matirx by likelihood of infection
          print("results")
          return(results)
        }
        return(NULL)
      }
      return(NULL)
    }
    return(NULL)
  }

  pairwise_calc  = function(total,pairwise,rank1){
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

              p_i_j =  binom_test(rownames(v_i),rownames(v_i_xr),rownames(v_j),N_rank,probability_i_j) #likelihood that virus_i has an infection that is unique relative to virus_j
              p_j_i =  binom_test(rownames(v_j),rownames(v_j_xr),rownames(v_i),N_rank,probability_j_i) #likelihood that virus_j has an infection that is unique relative to virus_i

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
          virus_i_prob = binom_test(rownames(v_i),rownames(v_i_xr),NULL,z,total_probs[grep(paste0(colnames(fullmatrix_sorted_ij[R1]),"$",collapse = ""),unlist(total_probs[,1])),2])

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
      final_rank = binom_test(rownames(final_i),rownames(final_xr),0,N_rank,probability_i)
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




  enrich = thresh # later implement so this can be changed??
  registerDoParallel(detectCores())
  plate = list()


  zeta = foreach(R = 1:(dim(case)[2]-1),.combine=rbind) %dopar%{ #cycle through each patient column by column (goal is so be serialized)


    #zeta = foreach(R = 35,.combine=rbind) %dopar%{ #cycle through each patient column by column (goal is so be serialized)

    # zeta = foreach(R = 1:1,.combine=rbind) %dopar%{ #cycle through each patient column by column (goal is so be serialized)
    # zeta = foreach(R = 37,.combine=rbind) %dopar%{ #cycle through each patient column by column (goal is so be serialized)
  #   zeta = for(R in 1:(dim(case)[2]-1)){ #cycle through each patient column by column (goal is so be serialized)

    rank <- total_calc(case,R,enrich,total,blast) # run the subsetting step given input of the total null probs, sample column, enrichment threshold and hte Virus blast matrix.

    if(is.null(rank) == FALSE){
      sorted_table=pairwise_calc(total,pairwise, rank) # take subset matrix from above and do all reassignments with pairwise and total null probabilities
      sorted_table_2 = simtag(sorted_table) # add indistinguishable tags
      name = colnames(case[R+1])
      output = as.data.frame(sorted_table_2)
      output[,12] =  gsub("^0\\||^0", '', output[,12])
      colnames(output) = c("Virus","P-value","Evidence Peptides","XR peptides","N-rank #","Evidence Peptide #","XR Peptide #","Filtered Evidence #","Filtered N-rank #","Null Probability","BH P-value","Indistinguishablity Groups")
      fwrite(output,file = paste0(out_path,name,".csv"))
      pool = cbind(name,output)
      return(pool)
    }
  }
  fwrite(as.data.frame(zeta[zeta[,9]>=3 & zeta[,12]<=.05,]),file = paste0(out_path,out_name,"AVARDA_compiled_full_output",".csv"))
  empty_virus_1 = colnames(blast)[which(colnames(blast)%in%zeta$Virus == FALSE)]
  empty_virus = data.frame(matrix(NA,ncol = length(unique(zeta$name))+1,nrow = length(empty_virus_1)))
  empty_virus[,1] = empty_virus_1
  fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`Evidence Peptides`) %>% spread(name,`Evidence Peptides`,fill = 0),empty_virus)),file = paste0(out_path,out_name,"AVARDA_evidence_pep",".csv"))
  fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`Evidence Peptide #`) %>% spread(name,`Evidence Peptide #`,fill = 0),empty_virus)),file = paste0(out_path,out_name,"AVARDA_unfiltered_evidence_number",".csv"))
  fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`Filtered Evidence #`) %>% spread(name,`Filtered Evidence #`,fill = 0),empty_virus)),file = paste0(out_path,out_name,"AVARDA_breadth",".csv"))
  fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`BH P-value`) %>% spread(name,`BH P-value`,fill = 1),empty_virus)),file = paste0(out_path,out_name,"AVARDA_post_p_value_BH",".csv"))
  fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`P-value`) %>% spread(name,`P-value`,fill = 1),empty_virus)),file = paste0(out_path,out_name,"AVARDA_post_p_value",".csv"))
}

