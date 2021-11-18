############################################
# Script 1 of 5: Load single-cell data and make Seurat objects
# AUTHOR: Elisa Rosati
# June 2021
# email: e.rosati@ikmb.uni-kiel.de
############################################
####################################################################################
#INTRO INFO
#GE_folder must have the following structure:
#must contain 1 folder per sample with inside the matrices files (in this case raw matrices are used)
#myname variable indicates the sample name given to the folder containing the matrices


#################################################################################
#START: LOAD DATA FROM RAW MATRICES
#load GEX data, create and filter seurat object

#LOAD MCKIFF
mydata<-system(paste0("ls -1d ",local,input_folder_mckiff, "*/outs/raw_feature_bc_matrix"),intern=T)

#select only data about SARS-CoV-2
mydata=mydata[grep("/CoVi",mydata)]
mydata=mydata[grep("_24_|_0_",mydata,invert=T)]
#read GEO annotation file
single_sample_annotation<-readRDS(paste0(metadata_folder,"201102_GEO_single_sample_annotation.rds"))

mylist=list()
for (i in 1:length(mydata)) {
  tmp=mydata[i]
  print(tmp)
  myname=strsplit(tmp,split="/")[[1]][sample_name_path_position] #select sample name
  input=tmp
  if (dir.exists(input)) {
    s=prepare_data_mckiff(myname=myname,
                                       GE_folder=input, single_sample_annotation=single_sample_annotation,
                                       with_TCR=TRUE, remove_TCR = TRUE, TCR_folder=TCR_folder_mckiff,
                                       processing_cite=FALSE,
                                       sample_name_short="covid_Pandurangan",
                                       output_folder=output_folder)
   mylist[[myname]]<-s
  }
}

#LOAD BACHER
owndata=c(1:9,11:21)
#load metadata
metadata_bacher=fread(paste0(metadata_folder,"201009_10X_COVID_T cell samples_metadata.csv"))
for (i in owndata) {
  tmp=i
  print(tmp)
  myname=paste0("covid_",i)
    s=prepare_data_bacher(myname=myname, GE_folder=input_folder_bacher, 
                          TCR_folder=TCR_folder_bacher,
                       with_TCR=TRUE, remove_TCR = TRUE, 
                       sample_name_short="AGspecific_covid",
                       metadata=metadata_bacher,
                       output_folder=output_folder)
     mylist[[myname]]<-s
}

#SAVE SEURAT OBJECTS
for (i in 1:length(mylist)) {
  my=mylist[[i]]
  myname=names(mylist)[i]
  saveRDS(my,file=paste0(output_folder,"starting_seurat_objects/start_seurat_object_x_associated_analysis_",myname,".rds"))
}
