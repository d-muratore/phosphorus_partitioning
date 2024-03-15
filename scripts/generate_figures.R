## Making figures
library(tidyverse)
library(patchwork)
library(rnaturalearthdata)

## Starting with Figure 1a
nav_track<-data.table::fread('../data/ship_data/nav_track.csv') %>%
  mutate(local_time=lubridate::with_tz(iso_time,tzone='Etc/GMT+4'))
ctd_smoothed<-fread('../data/intermediate_data_files/invirt_2019_ctd.csv')
surface_ctd<-ctd_smoothed %>%
  filter(rounded_depth==-10)

world<-map_data('world')

# Plot statement
f1a<-ggplot()+
  geom_polygon(data=filter(world,region=='Bermuda'),
               aes(x=long,y=lat,group=group),
               fill='cornsilk1',
               col='darkolivegreen4')+
  geom_path(data=nav_track,
            mapping=aes(x=ship_longitude,y=ship_latitude),
            size=1,
            col='darkgrey')+
  ggsn::scalebar(data=(nav_track %>% rename('long'='ship_longitude',
                                            'lat'='ship_latitude')),
                 dist_unit='km',model='WGS84',transform=TRUE,
                 dist=10,location='bottomleft',
                 anchor=c('x'=-64.41,'y'=31.625))+
  geom_point(aes(y=31.66,
                 x=-64.16),
             size=4,col='red')+
  ggrepel::geom_label_repel(aes(y=31.66,
                                x=-64.16,
                                label='BATS'),nudge_y=-0.025)+
  ggrepel::geom_label_repel(data=filter(surface_ctd,cast==2),
                            aes(x=binned_lon,
                                y=binned_lat,
                                label='12 Oct 8pm'),
                            nudge_y=0.05)+
  ggrepel::geom_label_repel(data=filter(surface_ctd,cast==29),
                            aes(x=binned_lon,
                                y=binned_lat,
                                label='17 Oct 8am'),
                            nudge_y=0.05)+
  coord_map(xlim=c(-64.425,-64.1),
            ylim=c(31.6,31.9))+
  geom_point(data=surface_ctd,
             aes(x=binned_lon,
                 y=binned_lat),
             size=3)+
  ylab('Latitude Degrees North')+
  xlab('Longitude Degrees East')+
  theme_minimal()+
  theme(text=element_text(size=14),
        panel.background=element_rect(fill='#F1F8FF',
                                      color='grey65'),
        panel.grid.major=element_line(color='grey45'))

## Moving to Figure 1b
diel_long_annotated<-fread('../data/intermediate_data_files/diel_long_annotated.csv')
f1b<-ggplot(diel_long_annotated %>% filter(pathway=='Transport'))+
  geom_rect(data=data.frame(xmin=c(0,22,46,70,94)+4,
                            xmax=c(10,34,58,82,106)+4,
                            ymin=-Inf,ymax=Inf),
            aes(xmin=xmin,
                xmax=xmax,
                ymin=ymin,ymax=ymax),fill='grey',alpha=0.5)+
  geom_line(aes(x=hours_from_go,y=z_score,group=id),alpha=0.1,linewidth=0.5)+
  geom_smooth(aes(x=hours_from_go,y=z_score,group=mprt))+
  facet_wrap(~factor(mprt,
                     levels=c(0,400,800,1200,1600,2000),
                     labels=c('00:00 Peak',
                              '04:00 Peak',
                              '08:00 Peak',
                              '12:00 Peak',
                              '16:00 Peak',
                              '20:00 Peak')),ncol=3)+
  ylab('Relative Transcript Abundance\n(z-score)')+
  xlab('Hours from First Sample')+
  theme_bw()+
  scale_x_continuous(expand=c(0,0))+
  theme(strip.background=element_blank(),
        text=element_text(size=16))

## Figure 1c
pval_frame<-read.csv('../data/intermediate_data_files/synchronicity_analysis.csv')[,-1]

## Identifying specific KOs for plott annotations
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

f1c<-ggplot(filter(pval_frame,str_detect(ko,'K')))+
  geom_boxplot(aes(x=factor(ifelse(n_tax>=20,20,n_tax),levels=seq(4,20)),y=diffs*4),
               outlier.color=NA)+
  geom_jitter(data=filter(pval_frame,str_detect(ko,'K'),!(ko %in% c(phosphorus_functions,nitrogen_kos$code,
                                                                    synch_kos$code))),
              aes(x=factor(ifelse(n_tax>=20,20,n_tax),levels=seq(4,20)),y=diffs*4,col=rejects))+
  ggrepel::geom_label_repel(data=filter(pval_frame,ko %in% c(phosphorus_functions,nitrogen_kos$code,synch_kos$code)) %>% 
                              left_join(rbind(phosphorus_kos,nitrogen_kos,synch_kos),by=c('ko'='code')),
                            aes(x=factor(ifelse(n_tax>=20,20,n_tax),levels=seq(4,20)),y=diffs*4,label=name),force=20)+
  geom_point(data=filter(pval_frame,ko %in% c(phosphorus_functions,nitrogen_kos$code,synch_kos$code)) %>% 
               left_join(rbind(phosphorus_kos,nitrogen_kos,synch_kos),by=c('ko'='code')),
             aes(x=factor(ifelse(n_tax>=20,20,n_tax),levels=seq(4,20)),
                 y=diffs*4,fill=rejects),
             size=4,col='black',shape=24)+
  theme_bw()+
  scale_x_discrete(drop=FALSE,breaks=seq(4,20),
                   labels=c(seq(4,19),'>=20'))+
  ylab('Average Peak Time Difference [Hours]')+
  xlab('# Taxa w/Diel Expression')+
  scale_color_brewer(name='Significantly Synchronous',palette='Set2',
                     labels=c('Asynchronous','Synchronous'))+
  scale_fill_brewer(name='Significantly Synchronous',palette='Set2',
                    labels=c('Asynchronous','Synchronous'))+
  guides(fill=guide_none())+
  theme(legend.position='bottom',
        text=element_text(size=14))


## Finishing composition of Figure 1
full_comm_figure1<-(f1a+f1b)/f1c+
  plot_annotation(tag_levels=list(c('a)','b)','c)')))
ggsave('../figures/figure_1.pdf',device='pdf',
       scale=1.25)


## Generating Figure 2

phosphorus_frame<-diel_long_annotated %>%
  filter((updated_function %in% phosphorus_functions)|pathway %in% c('Phosphonate and phosphinate metabolism')) %>%
  mutate(p_function=ifelse(updated_function %in% c('K02044','K02041','K02042','K05781','K11081'),
                           'Phosphonate uptake',
                           ifelse(updated_function %in% c('K02040','K02036','K02037','K02038'),
                                  'Phosphate uptake',
                                  ifelse(updated_function %in% c('K07657','K07636','K01077'),
                                         'APase',
                                         ifelse(updated_function %in% c('K00937','K22468','K01507'),
                                                'Pi Assimilation via ATP Synthase',
                                                'C-P Lyase-Related Functions')))))

f2<-ggplot(filter(phosphorus_frame,p_function %in% c('Phosphate uptake',
                                                 'Phosphonate uptake')) %>%
         mutate(plot_name=ifelse(updated_name %in% c('Diatoms',
                                                     'Dinophyceae',
                                                     'Synechococcus',
                                                     'Prochlorococus',
                                                     'Rhodobacteriales',
                                                     'Rhodospiralles'),
                                 updated_name,
                                 'Other Heterotrophic Bacteria'),
                plot_name=ifelse(plot_name=='Prochlorococus','Prochlorococcus',
                                 plot_name)))+
  geom_rect(data=data.frame(xmin=c(0,22,46,70,94)+4,
                            xmax=c(10,34,58,82,106)+4,
                            ymin=-Inf,ymax=Inf),
            aes(xmin=xmin,
                xmax=xmax,
                ymin=ymin,ymax=ymax),fill='grey',alpha=0.5)+
  geom_vline(data=data.frame(mprt=0,x=seq(8,112,by=24)),
             aes(xintercept=x),
             linetype='dashed')+
  geom_vline(data=data.frame(mprt=0400,x=seq(12,112,by=24)),
             aes(xintercept=x),
             linetype='dashed')+
  geom_vline(data=data.frame(mprt=0800,x=seq(16,112,by=24)),
             aes(xintercept=x),
             linetype='dashed')+
  geom_vline(data=data.frame(mprt=1200,x=seq(20,112,by=24)),
             aes(xintercept=x),
             linetype='dashed')+
  geom_vline(data=data.frame(mprt=1600,x=seq(24,112,by=24)),
             aes(xintercept=x),
             linetype='dashed')+
  geom_vline(data=data.frame(mprt=2000,x=seq(4,112,by=24)),
             aes(xintercept=x),
             linetype='dashed')+
  geom_point(aes(x=hours_from_go,y=z_score,col=factor(plot_name,
                                                      levels=c('Rhodobacteriales',
                                                               'Rhodospiralles',
                                                               'Other Heterotrophic Bacteria',
                                                               'Dinophyceae',
                                                               'Diatoms',
                                                               'Prochlorococcus',
                                                               'Synechococcus'))))+
  geom_line(aes(x=hours_from_go,y=z_score,col=factor(plot_name,
                                                     levels=c('Rhodobacteriales',
                                                              'Rhodospiralles',
                                                              'Other Heterotrophic Bacteria',
                                                              'Dinophyceae',
                                                              'Diatoms',
                                                              'Prochlorococcus',
                                                              'Synechococcus'))))+
  facet_wrap(~p_function+factor(mprt,
                                levels=c(800,1200,1600,2000,0,400),
                                labels=c('08:00 Peak',
                                         '12:00 Peak',
                                         '16:00 Peak',
                                         '20:00 Peak',
                                         '00:00 Peak',
                                         '04:00 Peak')),ncol=2,dir='v',
             drop=FALSE)+
  xlab('Hours from First Sample')+
  ylab('Relative Expression Level')+
  scale_color_manual(name='Taxon',
                     values=c('#E06711',
                              '#F3B684',
                              '#C38C64',
                              '#3651E6',
                              '#7896E8',
                              '#269217',
                              '#435642'))+
  theme_bw()+
  theme(text=element_text(size=14),
        strip.background=element_blank())
ggsave('../figures/figure_2.pdf',device='pdf',scale=1.2)
