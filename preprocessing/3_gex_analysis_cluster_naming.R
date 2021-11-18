############################################
# Script 3 of 5: Find clusters and perform gene expression analysis and cluster naming
# AUTHOR: Elisa Rosati
# June 2021
# email: e.rosati@ikmb.uni-kiel.de
############################################
###########################################
#LOAD BATCH CORRECTED AND HARMONIZED SEURAT OBJECT
seurat<-readRDS(paste0(output_folder,"merged_seurat_harmonizedUMAP_v3.rds"))

#################################################################################
#DIFFERENTIAL GEX ANALYSIS

marker_gene_list=list()
for(i in 0:(length(unique(seurat$seurat_clusters))-1)){
  print(i)
  cluster.markers <- FindMarkers(seurat, group.by = "seurat_clusters",ident.1 = as.character(i), min.pct = 0.25,test.use = "MAST",only.pos = FALSE)
  write.table(cluster.markers,file=paste0(output_folder,"marker_genes/","cluster",as.character(i),"_marker_genes.tsv"),sep="\t",col.names=TRUE,row.names = TRUE)
  cluster.markers$seurat_cluster=rep(i,nrow(cluster.markers))
  number=i+1
  marker_gene_list[[number]]<-cluster.markers
  gc()
}


all_markers<-do.call(rbind,marker_gene_list)
write.table(all_markers,file=paste0(output_folder,"marker_genes/","all_clusters_markers_seurat_cluster.tsv"),sep="\t",col.names=TRUE,row.names = TRUE)
############################################################
#NAME CLUSTERS
seurat@meta.data$new_cluster_name=seurat@meta.data$seurat_clusters
seurat@meta.data$new_cluster_name=gsub("^0$","Tfh",seurat@meta.data$new_cluster_name)  
seurat@meta.data$new_cluster_name=gsub("^1$","Cytotoxic",seurat@meta.data$new_cluster_name) 
seurat@meta.data$new_cluster_name=gsub("^2$","Th17-like",seurat@meta.data$new_cluster_name)
seurat@meta.data$new_cluster_name=gsub("^3","Transitional Tcm/Tfh",seurat@meta.data$new_cluster_name) 
seurat@meta.data$new_cluster_name=gsub("^4","Th17",seurat@meta.data$new_cluster_name)
seurat@meta.data$new_cluster_name=gsub("^5","Effector memory",seurat@meta.data$new_cluster_name) 
seurat@meta.data$new_cluster_name=gsub("^6","Tcm",seurat@meta.data$new_cluster_name)
seurat@meta.data$new_cluster_name=gsub("^7","Transitional Tfh",seurat@meta.data$new_cluster_name) 
seurat@meta.data$new_cluster_name=gsub("^8","Type I IFN-signature 1",seurat@meta.data$new_cluster_name)
seurat@meta.data$new_cluster_name=gsub("^9","Th1",seurat@meta.data$new_cluster_name) #TfH 1 
seurat@meta.data$new_cluster_name=gsub("^10$","Type I IFN-signature 2",seurat@meta.data$new_cluster_name) #TfH 1 
seurat@meta.data$new_cluster_name=gsub("^11$","Cytotoxic-like",seurat@meta.data$new_cluster_name) #TfH 1 
seurat@meta.data$new_cluster_name=gsub("^12$","Cycling",seurat@meta.data$new_cluster_name) #TfH 1 

###############################################################################################################
#RENUMBER CLUSTERS
seurat@meta.data$new_cluster_number=rep(1,nrow(seurat@meta.data))
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Transitional Tfh",gsub(1,2,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Th1",gsub(1,3,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Transitional Tcm/Tfh",gsub(1,4,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Tcm",gsub(1,5,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Th17",gsub(1,6,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Th17-like",gsub(1,7,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Effector memory",gsub(1,8,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Type I IFN-signature 1",gsub(1,9,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Type I IFN-signature 2",gsub(1,10,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Cytotoxic-like",gsub(1,11,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Cytotoxic",gsub(1,12,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=ifelse(seurat@meta.data$new_cluster_name=="Cycling",gsub(1,13,seurat@meta.data$new_cluster_number),seurat@meta.data$new_cluster_number)
seurat@meta.data$new_cluster_number=factor(seurat@meta.data$new_cluster_number,levels=1:13)
###############################################################################################################
#SAVE METADATA AND SEURAT OBJECT
tt=seurat@meta.data
saveRDS(tt,file=paste0(output_folder,"original_merged_metadata_xTCR_info_newClustersv2.rds"))

#this is the final seurat object
saveRDS(seurat,file=paste0(output_folder,"merged_seurat_NewClusters_v4.rds"))
####################################################################################################
#SAVE MARKER GENES WITH CLUSTER ANNOTATION
cluster_match<-unique(seurat@meta.data[,c("seurat_clusters","new_cluster_number")])

all_markers<-merge(all_markers,cluster_match,by.x="seurat_cluster",by.y="seurat_clusters",all.x=T)

write.table(all_markers,file=paste0(output_folder,"marker_genes/","all_clusters_new_number.tsv"),sep="\t",col.names=TRUE,row.names = TRUE)

################################################################################################################
#ADD UMAP INFO AND SAVE IT
umap_x_plot= seurat@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(barcode=seurat@meta.data$barcode,
                            lib=seurat@meta.data$origlib)
saveRDS(umap_x_plot,file=paste0(output_folder,"umap_x_plot_Seuratmerged.rds"))


