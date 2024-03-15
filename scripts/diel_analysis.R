## Conducting diel periodicity analysis for 
## VST-normalized transcript count tables
## this analysis will produce the necessary data and tables
## to put together Figure 1b and to run
## the analyses required for Figure 1c and Figure 2.
source('load_libraries.R')
source('custom_functions.R')

## Reading in environmental covariate data
sample_data<-data.table::fread('../data/ship_data/INVIRT19_metaT_metadata.txt')
surface_samples<-sample_data %>%
  filter(depth_ID=='SRF') %>%
  arrange(condition)

## Reading in metatranscriptome data
euk_metaT<-data.table::fread('../data/transcript_count_data/INVIRT19_timeseries_eukaryotes_new.txt') %>%
  select(-V1)
cyano_metaT<-data.table::fread('../data/transcript_count_data/INVIRT19_timeseries_cyanobacteria_new.txt') %>%
  select(-V1)
het_metaT<-data.table::fread('../data/transcript_count_data/INVIRT19_timeseries_heterotrophic_bacteria_new.txt') %>%
  select(-V1)

## Analyzing eukaryotic fraction for diel transcripts
## Reformatting table
euk_tax_table<-euk_metaT %>%
  select(c(original_group_phylum,New_Name)) %>%
  distinct() %>%
  dplyr::rename(og_name=original_group_phylum)

euk_wide<-euk_metaT[sample %in% surface_samples$sample] %>%
  mutate(sample=factor(sample,levels=surface_samples$sample)) %>%
  arrange(sample)
euk_vst<-euk_wide %>%
  distinct() %>%
  filter(KEGG_KO!='-') %>%
  mutate(updated_id=paste0(gsub(" ",'',original_group_phylum),'_',KEGG_KO)) %>%
  pivot_wider(id_cols=c('original_group_phylum','KEGG_KO','New_Name','updated_id'),
              names_from='sample',
              values_from='sum_ko_vst')
euk_vst_cleaned<-euk_vst %>%
  select(c(updated_id,starts_with('SX')))

## Detrending and setting up eukaryotic metatranscriptome data for rain analysis
euk_numeric<-as.matrix(euk_vst_cleaned[,-1])
rownames(euk_numeric)<-euk_vst_cleaned$updated_id
euk_numeric<-t(scale(t(euk_numeric)))
euk_clean<-euk_numeric[!is.na(rowSums(euk_numeric)),]
surface_sample_times<-surface_samples$condition[which(surface_samples$sample %in% colnames(euk_vst_cleaned)[-1])]
possible_days<-seq(12,17)
possible_times<-c('0000','0400','0800','1200','1600','2000')
possible_codes<-paste(rep(possible_days,each=length(possible_times)),rep(possible_times,length(possible_days)),sep='_')
possible_codes<-paste0(possible_codes,'_SRF')
possible_codes<-possible_codes[-c(1:5,34:36)]
time_matches<-sapply(possible_codes,function(x) length(which(surface_sample_times==x)))
times_from_go<-seq(0,length(time_matches)-1)*4
x_times<-rep(times_from_go,times=time_matches)
euk_detrended<-t(apply(euk_clean,1,function(y) detrending_function(row=y,x=x_times)))

## Run rain test
euk_rain_results<-rain(t(euk_detrended),
                       deltat=4,
                       period=24,
                       measure.sequence=time_matches,
                       method='independent',
                       verbose=TRUE)
write.csv(file='../data/intermediate_data_files/euk_rain_results.csv',x=euk_rain_results)

## Now we do the same for the other two tables
## Analyzing cyanobacterial portion
cyano_tax_table<-cyano_metaT %>%
  select(c(original_group_genus,New_Name)) %>%
  distinct() %>%
  dplyr::rename(og_name=original_group_genus)


cyano_wide<-cyano_metaT[sample %in% surface_samples$sample] %>%
  mutate(sample=factor(sample,levels=surface_samples$sample)) %>%
  arrange(sample)
cyano_vst<-cyano_wide %>%
  distinct() %>%
  filter(KEGG_KO!='-') %>%
  mutate(updated_id=paste0(gsub(" ",'',original_group_genus),'_',KEGG_KO)) %>%
  pivot_wider(id_cols=c('original_group_genus','KEGG_KO','New_Name','updated_id'),
              names_from='sample',
              values_from='sum_ko_vst')
cyano_vst_cleaned<-cyano_vst %>%
  select(c(updated_id,starts_with('SX')))

cyano_numeric<-as.matrix(cyano_vst_cleaned[,-1])
rownames(cyano_numeric)<-cyano_vst_cleaned$updated_id
cyano_numeric<-t(scale(t(cyano_numeric)))
cyano_clean<-cyano_numeric[!is.na(rowSums(cyano_numeric)),]
cyano_detrended<-t(apply(cyano_clean,1,function(y) detrending_function(row=y,x=x_times)))

cyano_rain_results<-rain(t(cyano_detrended),
                         deltat=4,
                         period=24,
                         measure.sequence=time_matches,
                         method='independent',
                         verbose=TRUE)
write.csv(file='../data/intermediate_data_files/cyano_rain_results.csv',x=cyano_rain_results)

## Analyzing heterotroph data
het_tax_table<-het_metaT %>%
  select(c(original_group_order,New_Name)) %>%
  distinct() %>%
  dplyr::rename(og_name=original_group_order)

het_wide<-het_metaT[sample %in% surface_samples$sample] %>%
  mutate(sample=factor(sample,levels=surface_samples$sample)) %>%
  arrange(sample)
het_vst<-het_wide %>%
  distinct() %>%
  filter(KEGG_KO!='-') %>%
  mutate(updated_id=paste0(gsub(" ",'',original_group_order),'_',KEGG_KO)) %>%
  pivot_wider(id_cols=c('original_group_order','KEGG_KO','New_Name','updated_id'),
              names_from='sample',
              values_from='sum_ko_vst')
het_vst_cleaned<-het_vst %>%
  select(c(updated_id,starts_with('SX')))

het_numeric<-as.matrix(het_vst_cleaned[,-1])
rownames(het_numeric)<-het_vst_cleaned$updated_id
het_numeric<-t(scale(t(het_numeric)))
het_clean<-het_numeric[!is.na(rowSums(het_numeric)),]
het_detrended<-t(apply(het_clean,1,function(y) detrending_function(row=y,x=x_times)))

het_rain_results<-rain(t(het_detrended),
                       deltat=4,
                       period=24,
                       measure.sequence=time_matches,
                       method='independent',
                       verbose=TRUE)
write.csv(file='../data/intermediate_data_files/het_rain_results.csv',x=het_rain_results)

## Running additional viral transcript timeseries
## for analysis to maintain the same familywise error correction
## procedure as in the initial analysis of this dataset
## viral transcripts are not discussed for the purposes of this brief
## communication
source('aux_transcripts.R')
full_tax_lookup<-rbind(euk_tax_table,cyano_tax_table,het_tax_table) %>%
  mutate(og_name=gsub('_X','',
                      gsub(" ",'',gsub('_XX','',og_name))))
write.csv(file='../data/intermediate_data_files/full_tax_lookup.csv',x=full_tax_lookup)
full_rain_results<-rbind(euk_rain_results,
                         cyano_rain_results,
                         het_rain_results,
                         virus_rain_results,
                         euk_virus_rain_results)
write.csv(file='../data/intermediate_data_files/full_transcript_rain_results.csv',x=full_rain_results)

## Adjusted-BH p-value correction procedure
full_rain_ps<-full_rain_results$pVal[order(full_rain_results$pVal)]
bh_full<-1:length(full_rain_ps)*0.1/length(full_rain_ps)
bh_pass_full<-max(which(full_rain_ps<=bh_full))
full_sig_subset<-full_rain_results[which(full_rain_results$pVal %in% full_rain_ps[1:bh_pass_full]),] %>%
  rownames_to_column('id') %>%
  mutate(KO=gsub('^.*_','',id),
         og_name=gsub('_.*$','',id)) %>%
  left_join(full_tax_lookup) %>%
  mutate(New_Name=ifelse(is.na(New_Name),'Virus',New_Name),
         og_name=ifelse(New_Name=='Virus',
                        gsub('__full.*$','',id),
                        og_name),
         KO=ifelse(New_Name=='Virus',gsub(';','',str_extract(id,'PFAM.*;')),
                   KO))

diel_metaT_dt<-rbind(virus_detrended,
                     cyano_detrended,
                     euk_detrended,
                     het_detrended,
                     ncldv_detrended)
diel_metaT_dt<-diel_metaT_dt[which(rownames(diel_metaT_dt) %in% full_sig_subset$id),] %>%
  as.data.frame() %>%
  rownames_to_column(var='id') %>%
  pivot_longer(starts_with('S'),
               names_to='sample',
               values_to='z_score') %>%
  left_join(surface_samples,by=c('sample'='sample')) %>%
  group_by(id,condition) %>%
  summarize(average_z=mean(z_score)) %>%
  pivot_wider(id_cols='id',
              names_from='condition',
              values_from='average_z')

write.csv(diel_metaT_dt,'../data/intermediate_data_files/diel_metat_dt.csv')

