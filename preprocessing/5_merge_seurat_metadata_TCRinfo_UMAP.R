############################################
# Script 5 of 5: Merge TCR data and UMAP positions to Seurat meta.data
# AUTHOR: Elisa Rosati
# June 2021
# email: e.rosati@ikmb.uni-kiel.de
############################################
#INTRO
######################################################
seurat_meta<-readRDS(paste0(output_folder,"original_merged_metadata_xTCR_info_newClustersv2.rds"))
seurat_meta$TCR_merge=paste(seurat_meta$origlib,seurat_meta$barcode,sep=".")
umap_x_plot<-readRDS(paste0(output_folder,"umap_x_plot_Seuratmerged.rds"))
#################################################################################
#LOAD SINGLE TCR FILES
TCRs=list.files(paste0(output_folder,"TCR"))
mylist=list()
for (i in 1:length(TCRs)) {
  mylist[[i]]<-readRDS(paste0(output_folder,"TCR/",TCRs[i]))
  mylist[[i]]$barcode=gsub("-1","",mylist[[i]]$barcode)
  myname=gsub("start_TCRinfo_","",TCRs[i])
  myname=gsub(".rds","",myname)
  mylist[[i]]$origin=rep(myname,nrow(mylist[[i]]))
  mylist[[i]]$TCR_merge=paste(mylist[[i]]$origin,mylist[[i]]$barcode,sep=".")
}

##############################################################################################
#merge GEXmetadata + TCR + UMAP
bigTCR=do.call(rbind,mylist)
seurat_metaTCR=merge(seurat_meta,bigTCR,by=c("TCR_merge"),all.x=T)
seurat_metaTCR[is.na(seurat_metaTCR)]<-"none"
saveRDS(seurat_metaTCR,file=paste0(output_folder,"seurat_metadata_withTCR.rds"))


seurat_metaTCR_UMAP=merge(seurat_metaTCR,umap_x_plot, by.x=c("barcode.x","origlib"),by.y=c("barcode","lib"))
seurat_metaTCR_UMAP$orig.hospital=gsub("none","healthy",seurat_metaTCR_UMAP$orig.hospital)
seurat_metaTCR_UMAP$orig.hospital=gsub("Yes","yes",seurat_metaTCR_UMAP$orig.hospital)
seurat_metaTCR_UMAP$orig.hospital=gsub("No","no",seurat_metaTCR_UMAP$orig.hospital)
saveRDS(seurat_metaTCR_UMAP,file=paste0(output_folder,"mergedSeurat_metadata_TCR_UMAPpositions.rds"))

##########################################################################################
