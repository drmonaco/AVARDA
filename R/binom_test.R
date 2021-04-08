#' @export


binom_test  = function(v_i,v_xr,v_i_j,N_rank,null_prob,dict){
  v_total = v_i[!v_i %in% v_xr] # all virus i aligning minus the xr
  v_total = v_total[!v_total %in% v_i_j] # all virus i evidence minus any shared with virus j
  v_total_f = filter.R(dict,v_total) # this calculates N_rank
  #if(length(N_rank)!=length(unlist(N_rank_2))){
  N_rank = N_rank[!N_rank %in% v_xr]
  N_rank = N_rank[!N_rank %in% (v_i %in% v_i_j)] # this is by default zero for total binom calculation
  N_rank_f = filter.R(dict,N_rank)
  #}
  if(N_rank_f == 0){
    return(NULL)
  }
  x = binom.test(v_total_f,N_rank_f,unlist(null_prob),"greater")[[3]]
  output = list(x,v_total_f,N_rank_f)
  return(output)
}
