############################################
# Script wrapperfor COVID-TCR paper
# AUTHOR: Elisa Rosati
# November 2021
# email: e.rosati@ikmb.uni-kiel.de
############################################
###########################################
# SETTINGS
########################################### 
local=""
output_folder="your output folder"
script_folder="folder with your scripts"
metadata_folder="folder with your metadata"
####################################################################################
#FUNCTIONS
source(paste0(script_folder,"COVID_TCR_preprocessing_functions.R"))
##########################################################################################
#create needed folders
dir.create(paste0(output_folder,"starting_seurat_objects"))
dir.create(paste0(output_folder,"marker_genes"))
dir.create(paste0(output_folder,"TCR"))
######################################################
#script 1 - create single seurat objects
input_folder_mckiff="your input folder, most likely where you store your cell ranger outputs"
TCR_folder_mckiff="same of input folder or different folder where you stored your cell ranger outputs"
input_folder_bacher="your input folder, most likely where you store your cell ranger outputs"
TCR_folder_bacher="same of input folder or different folder where you stored your cell ranger outputs"

sample_name_path_position=7 #change as needed
source(paste0(script_folder,"1_make_seurat_objects.R"))

##########################################
#script 2 - merge and batch correct seurat objects
source(paste0(script_folder,"2_load_merge_batchCorrect_seurat_object.R"))

##########################################
#script 3 - perform differential gene expression analysis and save UMAp coordinates
source(paste0(script_folder,"3_gex_analysis_cluster_naming.R"))

###############################################
#script 4 - create separate TCR files
TCR_folder_mckiff="same of input folder or different folder where you stored your cell ranger outputs"
TCR_folder_bacher="same of input folder or different folder where you stored your cell ranger outputs"
sample_name_path_position_TCR_mckiff=8 #change as needed

source(paste0(script_folder,"4_create_separate_TCRinfo.R"))

##########################################
#script 5 - merge seurat metadata, TCR and UMAP info

source(paste0(script_folder,"5_merge_seurat_metadata_TCRinfo_UMAP.R"))



