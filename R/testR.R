testR = function(case_path,thresh,dict_path,total_path,pairwise_path,blast_path,out_path,out_name){

  dict =  data.frame(fread(dict_path,data.table = FALSE)) # read in peptide-peptide dictionary
  total = data.frame(fread(total_path,data.table = FALSE)) #read in null probability table
  pairwise = data.frame(fread(pairwise_path,data.table = FALSE),row.names=1) #read in null probability table
  blast = data.frame(fread(blast_path,data.table = FALSE),row.names=1) #read in virus-peptide alignment filtered table
  case = data.frame(fread(case_path, data.table = FALSE,check.names = FALSE)) # this is the case of interest to be analyzed - changeable


  enrich = thresh # later implement so this can be changed??
  registerDoParallel(detectCores())
  plate = list()


  # zeta = foreach(R = 1:(dim(case)[2]-1),.combine=rbind) %dopar%{ #cycle through each patient column by column (goal is so be serialized)
    R = 1
    rank <- total_calc(case,R,enrich,total,blast) # run the subsetting step given input of the total null probs, sample column, enrichment threshold and hte Virus blast matrix.

    # if(is.null(rank) == FALSE){
    #   sorted_table=pairwise_calc(total,pairwise, rank) # take subset matrix from above and do all reassignments with pairwise and total null probabilities
    #   sorted_table_2 = simtag(sorted_table) # add indistinguishable tags
    #   name = colnames(case[R+1])
    #   output = as.data.frame(sorted_table_2)
    #   output[,12] =  gsub("^0\\||^0", '', output[,12])
    #   colnames(output) = c("Virus","P-value","Evidence Peptides","XR peptides","N-rank #","Evidence Peptide #","XR Peptide #","Filtered Evidence #","Filtered N-rank #","Null Probability","BH P-value","Indistinguishablity Groups")
    #   fwrite(output,file = paste0(out_path,name,".csv"))
    #   pool = cbind(name,output)
    #   return(pool)
    #  }
  # }
  # fwrite(as.data.frame(zeta[zeta[,9]>=3 & zeta[,12]<=.05,]),file = paste0(out_path,out_name,"AVARDA_compiled_full_output",".csv"))
  # empty_virus_1 = colnames(blast)[which(colnames(blast)%in%zeta$Virus == FALSE)]
  # empty_virus = data.frame(matrix(NA,ncol = length(unique(zeta$name))+1,nrow = length(empty_virus_1)))
  # empty_virus[,1] = empty_virus_1
  # fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`Evidence Peptides`) %>% spread(name,`Evidence Peptides`,fill = 0),empty_virus)),file = paste0(out_path,out_name,"AVARDA_evidence_pep",".csv"))
  # fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`Evidence Peptide #`) %>% spread(name,`Evidence Peptide #`,fill = 0),empty_virus)),file = paste0(out_path,out_name,"AVARDA_unfiltered_evidence_number",".csv"))
  # fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`Filtered Evidence #`) %>% spread(name,`Filtered Evidence #`,fill = 0),empty_virus)),file = paste0(out_path,out_name,"AVARDA_breadth",".csv"))
  # fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`BH P-value`) %>% spread(name,`BH P-value`,fill = 1),empty_virus)),file = paste0(out_path,out_name,"AVARDA_post_p_value_BH",".csv"))
  # fwrite(rbindlist(list(zeta %>% dplyr::select(name, Virus,`P-value`) %>% spread(name,`P-value`,fill = 1),empty_virus)),file = paste0(out_path,out_name,"AVARDA_post_p_value",".csv"))
}
