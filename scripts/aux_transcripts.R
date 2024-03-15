## Doing the diel analysis of the other transcripts we included in this
## effort but are not emphasizing in this study
## for the purpose of retaining the same familywise error correction

## Reading in data
virus_metaT<-data.table::fread('../data/transcript_count_data/INVIRT19_vOTU_vst_anno_filtered_subset.txt') %>%
  mutate(full_id=paste0(gene_ID,';',
                        'PFAM_',
                        annotation_PFAM,
                        ';',
                        ifelse(is.na(gene_id_dramv)||gene_id_dramv=='',
                               'no_DRAM',
                               gene_id_dramv)))
ncldvs_2<-data.table::fread('../data/annotation_data/ncldv_annotations.csv')
ncldvs<-data.table::fread('../data/transcript_count_data/INVIRT19_timeseries_NCLDV_vst_normalized_cts_withannotation.txt')
rnav<-data.table::fread('../data/transcript_count_data/INVIRT19_timeseries_rnavirusRDRP_vst_noramlized_cts.txt')
ssdna_marker<-data.table::fread('../data/transcript_count_data/INVIRT19_timeseries_ssDNArep_vst_noramlized_cts.txt')

## Analyzing viral fraction

virus_filtered<-virus_metaT[-grep('PFAM_-;no_DRAM',virus_metaT$full_id),]
virus_wide<-as.data.frame(virus_filtered)[,c('full_id',surface_samples$sample)]

virus_numeric<-as.matrix(virus_wide[,-1])
rownames(virus_numeric)<-virus_wide[,1]
virus_numeric<-t(scale(t(virus_numeric)))
virus_clean<-virus_numeric[!is.na(rowSums(virus_numeric)),]
virus_detrended<-t(apply(virus_clean,1,function(y) detrending_function(row=y,x=x_times)))

virus_rain_results<-rain(t(virus_detrended),
                         deltat=4,
                         period=24,
                         measure.sequence=time_matches,
                         method='independent',
                         verbose=TRUE)

ncldv_wide<-rbind(as.data.frame(ncldvs)[,c('full_gene_name',surface_samples$sample)],
                  as.data.frame(rnav)[,c('full_gene_name',surface_samples$sample)],
                  as.data.frame(ssdna_marker)[,c('full_gene_name',surface_samples$sample)])
ncldv_numeric<-as.matrix(ncldv_wide[,-1])
rownames(ncldv_numeric)<-ncldv_wide[,1]
ncldv_numeric<-t(scale(t(ncldv_numeric)))
ncldv_clean<-ncldv_numeric[!is.na(rowSums(ncldv_numeric)),]
ncldv_detrended<-t(apply(ncldv_clean,1,function(y) detrending_function(row=y,x=x_times)))
euk_virus_tax_table<-data.frame(id=c(ncldvs_2$full_gene_name,rnav$full_gene_name,ssdna_marker$full_gene_name),
                                  og_name=c(ncldvs_2$genome,rep('RNAvirus',nrow(rnav)),rep('SSDNAvirus',nrow(ssdna_marker))),
                                  KO=c(ncldvs_2$consensus_function,rep('RNAvirus_marker',nrow(rnav)),rep('SSDNAvirus_marker',nrow(ssdna_marker)))) %>%
  mutate(New_Name=ifelse(og_name %in% c('RNAvirus','SSDNAvirus'),
                         og_name,'NCLDV'))

euk_virus_rain_results<-rain(t(ncldv_detrended),
                               deltat=4,
                               period=24,
                               measure.sequence=time_matches,
                               method='independent',
                               verbose=TRUE)

