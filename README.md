**

## README

**
**Introduction**
Instructions for running code to reproduce analyses in "Diel Partitioning in Microbial Phosphorus Acquisition in the Sargasso Sea." The only softwares that should be required in order to run this analysis are R v3.4 or greater and the R dependencies listed in the file load_libraries.R. These dependencies should be able to be installed via CRAN, or in the case of libraries listed in load_libraries.R, manually installed by downloading tarballs from CRAN archives (versions for those packages are listed in the comments). 

All data in the data subdirectory can be found and downloaded from the figshare connected to this git repository (figshare data DOI: 10.6084/m9.figshare.25253389)

**Overview of Subdirectories in Repository**

*Data*
Contains all data files required to run the code and generate figures for this study. Data files are divided into a few categories:

 - annotation_data: Data linking KEGG orthologue (KO) assignments of
   genes to a KEGG Pathway assignment, as well as annotation data for
   viral genes used in the overall analysis of this dataset but not
   discussed in this manuscript. KEGG pathway assignments were
   originally collected from the KEGG pathway mapper and then manually
   curated to identify the most relevant pathways for each orthologue,
   as well as to assign pathways to genes that had missing or nonsense
   pathway assignments.
		- kegg_pathway_1.csv - []
		- kegg_pathway_2.csv - []
		- ncldv_annotations.csv - []
   
 - ship_data: Sample metadata for transcriptomes including time,
   location, sample depth, CTD bottle number, etc. Also includes ship
   navigation data from the ship underway system and CTD sensor outputs
   (available on bco-dmo for cruise AE1926 under project name "INVIRT"
   https://www.bco-dmo.org/project/775816).
		   - all_casts.csv - []
		   - nav_track.csv - []
		   - INVIRT19_metaT_metadata.txt - []
   
  - transcript_count_data: Processed metatranscriptome mapping to
   assembles at the gene-level with taxonomic and functional annotations
   as described in Methods in the text. All raw sequence data are
   available under the URL
   https://genome.jgi.doe.gov/portal/Infvirtimeseries/Infvirtimeseries.info.html.
   A text file with the accession numbers for each individual
   metatranscriptome is also posted in this repository as
   metatranscriptome_accession_information.txt. 
		   - INVIRT19_timeseries_cyanobacteria_new.txt - []
		   - INVIRT19_timeseries_eukaryotes_new.txt - []
		   - INVIRT19_timeseries_heterotrophic_bacteria_new.txt - []
		   - INVIRT19_timeseries_NCLDV_vst_normalized_cts_withannotation.txt - []
		   - INVIRT19_timeseries_rnavirusRDRP_vst_noramlized_cts.txt - []
		   - INVIRT19_timeseries_ssDNArep_vst_noramlized_cts.txt - []
		   - INVIRT19_vOTU_vst_anno_filtered_subset.txt - []
   
 - intermediate_data_files: These are data files produced by the code
   contained in the 'scripts' directory that are used for figure
   generation and downstream analyses, please see the comments where
   those files are generated in the scripts for more details about each
   individual file.
		   -cyano_rain_results.csv - []
		   -diel_long_annotated.csv - []
		   -diel_metat_dt.csv - []
		   -euk_rain_results.csv - []
		   -full_tax_lookup.csv - []
		   -full_transcript_rain_results.csv - []
		   -het_rain_results.csv - []
		   -invirt_2019_ctd.csv - []
		   -synchronicity_analysis.csv - []
   
 - supplemental_data: the nicely formatted supplementary data files
   included in extended data for the manuscript and cited in the
   manuscript text.
	   -SD_1.csv: Ship navigation data with sampling locations and sample environmental data for all metatranscriptome samples.
	   -SD_2.csv: Results of RAIN diel periodicity analysis. Columns include KEGG orthologue tested, taxonomic assignment of KEGG orthologue, RAIN test p-value, and a reject/fail to reject decision on the null hypothesis that the transcript does not have a 24-hour diel oscillation.
	   -SD_3.csv: Peak time and functional and taxnomic annotation for genes with diel periodicity.
	   -SD_4.csv: Results of KEGG orthologue synchronicity analysis. Columns show the KEGG Orthologue number of the gene tested, the number of taxa in the dataset with diel expression of this KO, the pathway assignment of the KO, the average difference in peak rank time between all organisms with diel expression of the KO, the associated empirical p-value (10,000 Monte Carlo simulations) and the reject/fail to reject decision of the null hypothesis that the transcript is not synchronized.

*Figures* 
All completed, generated figures files as specified in the script generate_figures.R.

*Scripts*
All scripts required to do the analysis and generated figures and intermediate data files.
Describing each of the scripts in the order that they should be run:

 1. load_libraries.R: Loading all of the necessary packages and dependencies to run the other scripts in this subdirectory.

2. custom_functions.R - functions written for the purpose of these analyses that are used in the other scripts, in particular the functions required to conduct the synchronicity analysis in synchronicity_analysis.R

3. diel_analysis.R - using the metatranscriptomic count data to run the rain analysis (see Methods) to detect 24-hour periodicity in transcript abundances. This script will also call aux_transcripts.R, which runs the exact same analysis but for the viral transcripts that are included for the purpose of multiple hypothesis testing correction but are not discussed in this manuscript. This file will generate the intermediate data files "cyano_rain_results.csv", "het_rain_results.csv", "euk_rain_results.csv", "full_tax_lookup.csv", "full_transcript_rain_results.csv", and "diel_metat_dt.csv".

4. aux_transcripts.R - []

5. function_pathway_assignment.R - this script will take the intermediate data files generated in diel_analysis.R, merge the results with the pathway annotation data, and produce tables that will facilitate the analysis of different functional genes and the synchronicity analysis. This script produces the intermediate data file "diel_long_annotated.csv"

6. synchronicity_analysis.R - this script conducts the synchronicity analysis presented in manuscript figure 1C as coined in Muratore et al., 2022 "Complex marine microbial communities partition metabolism of scarce resources over the diel cycle". WARNING: this script take a long time to run as the number of monte carlo simulations of the test statistic for the purpose of hypothesis testing is high to assure stable p-values. This script produces the intermediate data file "synchronicity_analysis.csv"

7. ctd_processing.R - this script takes the ship navigation and CTD sensor data, along with the metatranscriptome sample metadata and generates the files necessary to produce figure 1a. This script is not required for any of the transcriptomic analysis, but is required for generate_figures.R. It produces the intermediate data file "invirt_2019_ctd.csv"

8. generate_figures.R - [].

For any questions about running the code, please feel free to contact Daniel Muratore dmuratore[at]santafe[dot]edu. 






