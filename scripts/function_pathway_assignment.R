## Doing diel transcript annotation to prepare for synchronicity analysis
source('load_libraries.R')
source('custom_functions.R')
source('diel_analysis.R')

## Reshaping the results from the diel rain analysis
diel_metaT_dt<-data.table::fread('../data/intermediate_data_files/diel_metat_dt.csv')
full_diel<-rbind(as.matrix(diel_metaT_dt[,-c(1,2)]))
rownames(full_diel)<-c(as.character(diel_metaT_dt$id))

## Getting sample metadata
sample_lookup<-data.table::fread('../data/ship_data/INVIRT19_metaT_metadata.txt')
surface<-sample_lookup %>%
  filter(depth_m=='5') %>%
  arrange(day,time) %>%
  mutate(hours_from_go=(((day-12)*2400+time)-1600)/100)

full_tax_lookup<-fread('../data/intermediate_data_files/full_tax_lookup.csv',header=TRUE)[,-1]

ncldvs_2<-fread('../data/annotation_data/ncldv_annotations.csv')
cyano_rain_results<-fread('../data/intermediate_data_files/cyano_rain_results.csv')
het_rain_results<-fread('../data/intermediate_data_files/het_rain_results.csv')
euk_rain_results<-fread('../data/intermediate_data_files/euk_rain_results.csv')

## Doing mean peak rank time calculations
colnames(full_diel)<-c(seq(4,80,by=4),c(88,96,100,104,108,112))
diel_long<-as.data.frame(full_diel) %>%
  rownames_to_column(var='id') %>%
  pivot_longer(-id,names_to='hours_from_go',values_to='z_score') %>%
  mutate(hours_from_go=as.numeric(hours_from_go)) %>%
  left_join(surface) %>%
  select(-c(sample,cast,depth_m,depth_ID,CTD_Btl,condition)) %>%
  distinct() %>%
  group_by(id) %>%
  mutate(rank_time=rank(z_score)) %>%
  group_by(id,time) %>%
  mutate(mean_rank_time=mean(rank_time)) %>%
  ungroup() %>%
  group_by(id) %>%
  mutate(mprt=time[which.max(mean_rank_time)])

## Fixing some spelling typos in taxon names
taxa_to_change<-c('Prymnesiophyceae',
                  'Oceanospirillales',
                  'Chromatiales',
                  'Rhizobiales')

## Associating appropriate metadata with transcripts and recoding 
## some ambiguous taxonomic assignments to 'other', etc
just_mprt<-diel_long %>%
  select(c(id,mprt)) %>%
  distinct() %>%
  mutate(tax_code=gsub('_ko.*$','',gsub('__full.*$','',id)),
         function_code=gsub('^.*_','',id),
         function_code=ifelse(function_code=='DRAM',
                              gsub('^.*PFAM','PFAM',id),
                              function_code),
         tax_code=gsub('\\.','-',tax_code),
         tax_code=gsub('_XX','',tax_code),
         tax_code=gsub('_X','',tax_code)) %>%
  left_join(full_tax_lookup %>% mutate(og_name=gsub('\\.','-',og_name)),by=c('tax_code'='og_name')) %>%
  mutate(edited_name=ifelse(is.na(New_Name),
                            ifelse(str_detect(tax_code,'^T[0-9]'),
                                   'Virus',
                                   ifelse(str_detect(tax_code,'RNA'),'RNAVirus',
                                          ifelse(str_detect(tax_code,'fna'),'NCLDV',
                                                 'Metabolite'))),
                            New_Name),
         abbreviated_name=ifelse(!(edited_name %in% c('Virus','Metabolite','RNAVirus','NCLDV')),'Cellular',
                                 edited_name),
         updated_name=ifelse(tax_code %in% taxa_to_change, tax_code, edited_name),
         updated_name=ifelse(updated_name=='Flacobacteria_X',
                             'Other Flavobacteria',
                             ifelse(str_detect(updated_name,'Gammaproteo'),
                                    'Other Gammaproteobacteria',
                                    ifelse(str_detect(updated_name,'Alphaproteo'),
                                           'Other Alphaproteobacteria',
                                           updated_name))),
         updated_function=ifelse(abbreviated_name=='Cellular',
                                 gsub('\\.',';',
                                      gsub('ko.','',function_code)),
                                 ifelse(abbreviated_name=='Virus',
                                        gsub('\\.no_DRAM','',
                                             gsub('PFAM_','',function_code)),
                                        ifelse(abbreviated_name == 'RNAVirus',
                                               'RNAV_marker',
                                               ifelse(abbreviated_name=='NCLDV',
                                                      ncldvs_2$consensus_function[which(gsub('-','.',ncldvs_2$full_gene_name)==id)],
                                                      'Metabolite'))))) %>%
  mutate(updated_function=gsub('ko:','',gsub(',',';',updated_function))) %>%
  mutate(class_2=ifelse(abbreviated_name=='Cellular',
                        ifelse(id %in% gsub('-','.',gsub(':','.',gsub(',','.',rownames(cyano_rain_results)))),
                               'Cyanobacteria',
                               ifelse(id %in% gsub('-','.',gsub(':','.',gsub(',','.',rownames(het_rain_results)))),
                                      'Heterotrophic Bacteria',
                                      ifelse(id %in% gsub('-','.',gsub(':','.',gsub(',','.',rownames(euk_rain_results)))),
                                             'Eukaryotic Phytoplankton',
                                             abbreviated_name))),abbreviated_name),
         class_2=ifelse(class_2=='Virus','Bacteriophage',class_2))

diel_long_annotated<-left_join(diel_long,just_mprt)
print('completed annotation merging')

## Doing KEGG pathway assignment based on two manually
## curated lists of KO -> PATH associations
## see Muratore, Boysen, Harke et al 2022 for details
pathways<-rep(NA,nrow(diel_long_annotated))
kegg_paths<-read.csv('../data/annotation_data/kegg_pathway_1.csv')
kegg_round2<-read.csv('../data/annotation_data/kegg_pathway_2.csv') %>%
  filter(consensus_pathway!='')
## Making assignments KO by KO
for(i in 1:nrow(kegg_paths)){
  hits<-grep(kegg_paths$kegg_KO[i],diel_long_annotated$updated_function)
  pathways[hits]<-kegg_paths$pathway[i]
  print(i)
}
for(i in 1:nrow(kegg_round2)){
  hits<-grep(kegg_round2$updated_function[i],diel_long_annotated$updated_function)
  pathways[hits]<-kegg_round2$consensus_pathway[i]
  print(i)
}

diel_long_annotated<-data.frame(diel_long_annotated,pathway=pathways)
write.csv(file='../data/intermediate_data_files/diel_long_annotated.csv',x=diel_long_annotated)

