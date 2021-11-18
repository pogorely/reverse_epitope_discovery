############################################
# Script 4 of 5: Create separate TCR tables 
# AUTHOR: Elisa Rosati
# June 2021
# email: e.rosati@ikmb.uni-kiel.de
############################################
###########################################
#INTRO
######################################################
seurat<-readRDS(paste0(output_folder,"merged_seurat_NewClusters_v4.rds"))
seurat_meta<-readRDS(paste0(output_folder,"original_merged_metadata_xTCR_info_newClustersv2.rds"))
seurat_meta$TCR_merge=paste(seurat_meta$origlib,seurat_meta$barcode,sep=".")

#################################################################################
#load TCR data, create and filter seurat object
#TCR data should be given in input in the CellRanger output structure

#MCKIFF
mydata<-system(paste0("ls -1d ",TCR_folder_mckiff, "*/outs"),intern=T)
mydata=mydata[grep("/CoVi",mydata)]
mydata=mydata[grep("_24_|_0_",mydata,invert=T)]

single_sample_annotation<-readRDS(paste0(local,"metadata/201102_GEO_single_sample_annotation.rds"))
mylist=list()
for (i in 1:length(mydata)) {
  tmp=mydata[i]
  print(tmp)
  myname=gsub("_TCR","",strsplit(tmp,split="/")[[1]][sample_name_path_position_TCR_mckiff]) #recover folder name (change number depending on folder path)
  input=tmp
  if (dir.exists(input)) {
    s=prepare_TCR_data_mckiff(myname=myname,
                       TCR_folder=TCR_folder_mckiff)
   mylist[[myname]]<-s
  }
}

#BACHER
owndata=c(1:9,11:21)
for (i in owndata) {
  tmp=i
  print(tmp)
  myname=paste0("covid_",i)
    s=prepare_TCR_data_bacher(myname=myname,
                       TCR_folder=TCR_folder_bacher)
     mylist[[myname]]<-s
}

#SAVE SINGLE TCR FILES
for (i in 1:length(mylist)) {
  my=mylist[[i]]
  myname=names(mylist)[i]
  saveRDS(my,file=paste0(output_folder,"TCR/start_TCRinfo_",myname,".rds"))
  mylist[[i]]$barcode=gsub("-1","",mylist[[i]]$barcode)
  mylist[[i]]$origin=rep(names(mylist)[i],nrow(mylist[[i]]))
  mylist[[i]]$TCR_merge=paste(mylist[[i]]$origin,mylist[[i]]$barcode,sep=".")
}

