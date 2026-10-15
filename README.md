# Hyphosphere_microbiome_MetaAnalysis
This repository contains all the code and files needed to reproduce the analyses performed by Blanchard et al., 2027, for the hyphosphere microbiome meta-analysis.
	-NCBI SRA accession numbers are given in column Run_Browser of the Metada_meta_analyse.xlsx file
	-The mothur and cutadapat commands to process the data are listed in the files Pipeline_16S_Metadata.txt and Slurm_cutadapt_metanalyse.txt
	-The final mothur OTU table is in the filecombo.list.trim.contigs.good.cuta_16SFR_RV_RV.good.unique.denovo.vsearch.pick.dgc.0.01.abund.0.01.pick.shared and the associated OTU taxonomy in file Table_Taxonomy.xlsx
	-The R code to analyse the data and reproduce Figures 2 and S1 are given in file Code_Fig2.R
	-The produced binary tables at 4 taxonomic levels are given as excel files including a "core_30%" sheet listing the core taxa.
