############################################
# Script 2 of 5: Merge Seurat objects, filtering and correct batch effects
# AUTHOR: Elisa Rosati
# June 2021
# email: e.rosati@ikmb.uni-kiel.de
############################################
###########################################
#LOAD SINGLE OBJECTS
files=list.files(paste0(output_folder,"starting_seurat_objects/"))
files_names=gsub("\\.rds","",files)
files_names=gsub("start_seurat_object_x_associated_analysis_","",files_names)

mylist=list()
for (i in 1:length(files)) {mylist[[files_names[i]]]<-readRDS(paste0(output_folder,"starting_seurat_objects/",files[i]))}

#massage the data to make the two datasets having the same metadata columns and values
for (i in 1:length(mylist)) {
  print(i)
  if (any(grepl("dataset",colnames(mylist[[i]]@meta.data)))) {
    print("Bacher data")
    mylist[[i]]@meta.data$origlib=mylist[[i]]@meta.data$orig.donor
    mylist[[i]]@meta.data$barcode=gsub("-1","",row.names(mylist[[i]]@meta.data))
    mylist[[i]]@meta.data=mylist[[i]]@meta.data[,-which(colnames(mylist[[i]]@meta.data) %in% c("dataset"))]
    mylist[[i]]@meta.data$orig.ident=rep("covid_Bacher",nrow(mylist[[i]]@meta.data))
  } else {
    print("Mckiff data")
    mylist[[i]]@meta.data$diagnosis=rep("COVID-19",nrow( mylist[[i]]@meta.data))
    mylist[[i]]@meta.data$diagnosis=ifelse(grepl("None", mylist[[i]]@meta.data$orig.hospital),gsub("COVID-19","Healthy", mylist[[i]]@meta.data$diagnosis), mylist[[i]]@meta.data$diagnosis)
    mylist[[i]]@meta.data=mylist[[i]]@meta.data[,-which(colnames(mylist[[i]]@meta.data) %in% c("orig.peptide","orig.stim_time","orig.HT_ID.global","orig.sex"))]
    mylist[[i]]@meta.data$orig.hospital=gsub("none","healthy",mylist[[i]]@meta.data$orig.hospital)
    mylist[[i]]@meta.data$orig.hospital=gsub("Yes","yes",mylist[[i]]@meta.data$orig.hospital)
    mylist[[i]]@meta.data$orig.hospital=gsub("No","no",mylist[[i]]@meta.data$orig.hospital)
  }
}

###################################################################################################
#MERGE SEURAT OBJECTS
seurat <- merge(mylist[[1]], y = mylist[2:length(mylist)])
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat,pattern = "^MT-")
seurat <- subset(seurat, subset = nFeature_RNA > 400 & nFeature_RNA < 3500 & percent.mt < 5) 

#REMOVE TCR, MITOCHONDRIAL AND RIBOSOMAL GENES
counts <- GetAssayData(seurat, assay = "RNA")
counts <- counts[-which(grepl( "TRAC|TRBC|TRBV|TRAV|TRGV|TRDV|^TRBJ|^TRAJ|^MT|^RPS|^RPL",rownames(counts))),]
seurat <- subset(seurat, features = rownames(counts))

#REMOVE GENES PRESENT IN LESS THAN 1% OF CELLS
counts=seurat[["RNA"]]@counts
n.cells.per.gene <- Matrix::rowSums(x = counts > 0)
n.cells.per.gene<- n.cells.per.gene/ncol(seurat)*100
genes_more001=names(n.cells.per.gene[n.cells.per.gene>0.01])
seurat <- subset(seurat, features=genes_more001)

tt=seurat@meta.data
saveRDS(tt,file=paste0(output_folder,"original_merged_metadata_xTCR_info_v1.rds"))
#####################################################################
#NORMALIZE, SCALE, MAKE PCA
all.genes <- rownames(seurat)
seurat<- Seurat::NormalizeData(seurat,verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE,features = all.genes) %>% 
  RunPCA(pc.genes = seurat@var.genes, npcs = 40, verbose = FALSE)

#saveRDS(seurat,file=paste0(output_folder,date1,"_","merged_seurat_scaledPCA_v2.rds"))
################################################
#seurat<-readRDS(paste0(output_folder,date1,"_","merged_seurat_scaledPCA_v2.rds"))

#Harmony batch effect correction by project, sequencing run and processing batch
seurat <- seurat %>% 
  RunHarmony(c("orig.ident","orig.project_id","orig.run"), plot_convergence = TRUE, 
             max.iter.harmony=40, iter.max = 40)

######################################################################
#calculate umap
seurat <- seurat %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) 

#find clusters
seurat <- seurat %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

saveRDS(seurat,file=paste0(output_folder,"merged_seurat_harmonizedUMAP_v3.rds"))
#############################################################################
plot1<-DimPlot(seurat, reduction = "umap",group.by = "seurat_clusters", label=T)+ labs(title="By Seurat cluster") 

tiff(paste0(output_folder,"dimplot_bySeurat.tiff"), width=6,height=6,units="in",res=300)
print(plot1)
dev.off()