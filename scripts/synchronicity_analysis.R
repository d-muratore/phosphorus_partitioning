## Performing diel synchronicity analysis

## If not already run, please run the following scripts to prepare the necessary intermediate data files:
#source('load_libraries.R')
#source('custom_functions.R')
#source('diel_analysis.R')

## Implementing on our data
diel_long_annotated<-fread('../data/intermediate_data_files/diel_long_annotated.csv')
ko_counts<-diel_long_annotated %>%
  group_by(id) %>%
  slice_head() %>%
  ungroup() %>%
  filter(is.na(updated_function)==FALSE,!(updated_function %in% c('','Metabolite','Uncharacterized'))) %>%
  group_by(updated_function) %>%
  count() %>%
  arrange(desc(n))

## Getting one entry per transcript and renormalizing the rank times
small_anno<-diel_long_annotated %>%
group_by(id) %>%
  slice_head() %>%
  ungroup() %>%
  mutate(time_rank=mprt/400)

pop_times<-floor(small_anno$time_rank)

## Retaining only KOs with diel expression in at least 4 different taxa
keeping_kos<-ko_counts[ko_counts$n>=4,1]
output_storage<-list()
## Applying the statistical hypothesis testing framework
## EXPECT THIS ANALYSIS TO TAKE SEVERAL HOURS RUN OUT OF PARALLEL ON
## A PERSONAL MACHINE
for(i in 1:length(keeping_kos$updated_function)){
  ko_of_i<-keeping_kos$updated_function[i]
  sub_vec<-floor(small_anno$time_rank[which(small_anno$updated_function==ko_of_i)])
  print(paste0('Conducting hypothesis test for ',keeping_kos$updated_function[i]))
  output_storage[[i]]<-full_wrapper(sub_vec,pop_times)
  print(i)
}

## Familywise error correction adjusted BH procedure
output_ps<-do.call(rbind,lapply(output_storage,function(x) x$pval))
output_diffs<-do.call(rbind,lapply(output_storage, function(x) x$md))
q<-0.1
i<-1:length(output_ps)
m<-length(output_ps)
ordered_ps<-output_ps[order(output_ps,decreasing=FALSE)]
stop_reject<-min(which(ordered_ps>=i*q/m))
decisions<-rep(c('reject','fail'),c(stop_reject,length(ordered_ps)-stop_reject))
plot(1:m,ordered_ps)
lines(1:m,i*q/m,col='red',lwd=2)
abline(v=stop_reject)
pval_frame<-data.frame(pval=output_ps,
                       diffs=output_diffs,
                       ko=keeping_kos$updated_function,
                       n_tax=ko_counts$n[ko_counts$n>=4]) %>%
  arrange(pval) %>%
  mutate(rejects=decisions)

## Adding pathway assignments to output table
for(i in 1:nrow(pval_frame)){
  pval_frame$pathway[i]<-small_anno$pathway[which(small_anno$updated_function==pval_frame$ko[i])[1]]
}

## Calling out annotations for Figure 1c
phosphorus_functions<-c('K01113','K07657','K01077','K07636',
                         'K07659',
                        'K11081','K02044','K02040',
                        'K02036','K02037','K02038','K02041','K02042',
                        'K05781','K06217','K00937',
                        'K22468','K01507')
phosphorus_kos<-data.frame(code=phosphorus_functions,
                           name=c('phoD',
                                  'phoB','phoAB','phoR','ompR',
                                  'phnS','phnD','pstS','pstB','pstC',
                                  'pstA','phnC','phnE','phnK','phoH',
                                  'ppk1','ppk2','ppa'))

nitrogen_kos<-data.frame(code=c('K00265','K00266','K11959','K03320',
                                'K01915','K00284',
                                'K04751','K04752',
                                'K00265;K00284','K01915;K01949','K04751;K04752'),
                         name=c('gltB','gltD','urtA','amt',
                                'glnA','gltS',
                                'glnB','glnK','gltB/gltS','glnA','glnB/glnK'))

synch_kos<-data.frame(code=c('K03798','K00239',
                             'K02274','K04043','K02703','K02706'),
                      name=c('ftsH','sdhA','coxA','dnaK','psbA','psbD'))


## Writing outputs
write.csv(pval_frame,'../data/intermediate_data_files/synchronicity_analysis.csv',quote=FALSE)
