## Custom functions required for analyses presented in
## (paper details)
## Run this script to load all necessary functions to complete
## The scripts presented in the 'scripts' folder'

## Detrending function for diel rhythmicity analysis
detrending_function<-function(row,x){
  regression<-lm(row~x)
  if(summary(regression)$coefficients[2,4]<0.05){
    return(regression$residuals)
  }else{
    return(row)
  }
}

## Calculating the circular distance between two peak times 
calculate_med_dist<-function(pair){
  difference<-abs(pair[1]-pair[2])
  return(min(difference,6-difference))
}
## Finding the mean pairwise circular distance between all peaktimes in a set
generate_mean_stat<-function(vec){
  all_pairwise<-t(combn(vec,2))
  all_differences<-apply(all_pairwise,1,calculate_med_dist)
  return(mean(all_differences))
}
## Using the number of diel transcripts for one KO, simulate a distribution of
## average pairwise circular distance between peaktimes of KOs drawn at random from all diel KOs
## We do this because the peaktimes are not uniformly distributed over the 6 peak times
generate_null_ensemble<-function(vec,pop,niter=1e4){
  number_samples<-length(vec)
  bootstrap_ensemble<-t(replicate(niter,sample(x=pop,size=number_samples,replace=TRUE)))
  print('Simulating Null Distribution of Test Statistic')
  null_dn<-apply(bootstrap_ensemble,1,generate_mean_stat)
  return(null_dn)
}
## Use the null ensemble to simulate a p-value
simulate_p<-function(vec,nd){
  sample_statistic<-generate_mean_stat(vec)
  simulated_p<-length(which(nd<=sample_statistic))/length(nd)
  return(simulated_p)
}
## A function to do all of these things and keep you updated on what's happening
full_wrapper<-function(vec,pop,niter=1e4){
  print('Reading Data')
  vec<-vec
  pop<-pop
  print('Monte Carlo Simulating Null Distribution (the slow part)...')
  nd<-generate_null_ensemble(vec,pop,niter=niter)
  print('Null distribution generated')
  print('Calculating empirical p-value')
  pval<-simulate_p(vec,nd)
  print('Test complete, gathering outputs')
  output_list<-list(pval=pval,md=generate_mean_stat(vec),null=nd)
  return(output_list)
}

