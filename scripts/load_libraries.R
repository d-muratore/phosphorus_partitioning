## Libraries required to run code in "Diel Partitioning in Microbial Phosphorus Acquisition in the Sargasso Sea" Muratore and Gilbert et al.
library(tidyverse) # I use v1.3.2
library(dplyr) # v1.0.10
library(patchwork) # v1.1.2
library(rain) # v1.28.0
library(data.table) # v1.14.6
library(lubridate) # v1.9.0
library(rnaturalearthdata) # v0.1.0
library(rnaturalearth) # v0.1.0
library(ggrepel) # v0.9.2
library(ggsn) # v0.5.0

## Heads up, a few of these libraries (rain, ggsn) are no longer maintained on CRAN so if you
## have to install them for the first time you will need to either devtools::install from 
## git or download the tarball from the CRAN archive and install locally.
## the version of rain I use for this is RAIN version 1.28.0 from 2015-09-01. The version of ggsn is 0.5.0.
## Recommendations for installing rain and ggsn:
#install.packages('maptools', repos='http://R-Forge.R-project.org')
#devtools::install_github('oswaldosantos/ggsn')
#BiocManager::install('rain')