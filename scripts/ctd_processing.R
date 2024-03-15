## Process the CTD data from AE1926
# Start by loading the necessary libraries
source('load_libraries.R')
# Read in the data collected about which CTD cast associates with which samples and which times
bottle_data<-read.delim('../data/ship_data/INVIRT19_metaT_metadata.txt')
# Read in the CTD profiles (available on BCO-DMO under project ID: Invirt)
ctd_data<-read.csv('../data/ship_data/all_casts.csv') %>%
  dplyr::select(-X)

# Binning the data into 0.5m bins and applying k-smooth interpolation
ctd_smooth_ready<-ctd_data %>%
  mutate(rounded_depth=round(z/0.5)*0.5) %>%
  group_by(Cast_ID,rounded_depth,ISO_DateTime_UTC) %>%
  summarize(across(everything(),
                   list(median),
                   .names='binned_{.col}')) %>%
  mutate(time=as.POSIXct(gsub('Z','',ISO_DateTime_UTC),
                         format='\t%Y-%m-%dT%H:%M:%S',
                         tz='UTC'),
         cast=as.numeric(gsub('^.*InVirTC','',Cast_ID))) %>%
  ungroup() %>%
  group_by(cast) %>%
  mutate(lat=mean(binned_Latitude),
         lon=mean(binned_Longitude)) %>%
  dplyr::select(-c(binned_Latitude,binned_Longitude)) %>%
  mutate(time_local=lubridate::with_tz(time,tzone='Etc/GMT+4')) %>%
  dplyr::select(c(cast,time,time_local,rounded_depth,starts_with('binned')))

## Getting ready to apply smoothing to each cast in sequence
variable_columns<-grep('binned',colnames(ctd_smooth_ready))
ctd_smoothed<-split(ctd_smooth_ready,ctd_smooth_ready$cast)

## Function to apply the ksmoothing.
sample_smoothing<-function(wt,variable_columns){
  wt[,variable_columns]<-apply(wt[,variable_columns],
                               2,
                               function(x) ksmooth(wt$rounded_depth,
                                                   x,
                                                   bandwidth=5)$y)
  return(wt)
}

## Applying to each cast and identifying mixed layer depth by the 0.125 deviation in sigma-theta from 10m criterion. 
ctd_smoothed<-do.call(rbind,lapply(ctd_smoothed,function(x) sample_smoothing(wt=x,variable_columns=variable_columns)))
ctd_smoothed<-ctd_smoothed %>%
  group_by(cast) %>%
  mutate(mld=max(rounded_depth[binned_sigma_theta>=(0.125+binned_sigma_theta[(round(rounded_depth)==-10)][1])]))

# Writing the output
write.csv(ctd_smoothed,'../data/intermediate_data_files/invirt_2019_ctd.csv',quote=FALSE,row.names=FALSE)
