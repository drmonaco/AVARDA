total_calc  = function(case,column,thresh,total,blast,dict){ # subsetting function that also ranks peptides
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
          x = binom_test(v_i,virus_i_xr,v_i_j,row.names(blast_subset),probability,dict) # run binom function that takes an into account xr and shared peptides (which is 0 for this step)
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
