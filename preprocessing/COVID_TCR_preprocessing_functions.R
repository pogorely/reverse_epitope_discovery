############################################
# Functions: Functions for single-cell gene expression and TCR analysis
# AUTHOR: Elisa Rosati
# November 2021
# email: e.rosati@ikmb.uni-kiel.de
############################################
#this script has been tested on
#R version 4.1.2 (2021-11-01)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.3 LTS
##########################################
#INSTALL OR LOAD PACKAGES
library(data.table)
if (!require('stringdist')) install.packages('stringdist'); library('stringdist')
if (!require('Seurat')) install.packages('Seurat'); library('Seurat')
if (!require('harmony')) install.packages('harmony'); library('harmony')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('MAST')) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")}
  
  BiocManager::install("MAST"); library('MAST')
}
######################################################################
#FUNCTIONS

prepare_TCR_data<- function(myname,TCR_folder) {
  y=fread(paste0(TCR_folder,myname,"_TCR/outs/clonotypes.csv"))
  to_remove=c()
  #remove TCRs which:
  #1) are singletons with more than 1 TCR alpha and more than 1 TCR beta
  #2) have more than 2 alphas and 2 betas
  #3) have more than one beta and no alpha or more than one alpha and no beta
  
  for (i in 1:nrow(y)) {
    if (str_count(y$cdr3s_aa[i],"TRB") > 1 & str_count(y$cdr3s_aa[i],"TRA") > 1 & y$frequency[i]==1) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRB") > 2 | str_count(y$cdr3s_aa[i],"TRA")> 2) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRB") > 1 & str_count(y$cdr3s_aa[i],"TRA")==0) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRA") > 1 & str_count(y$cdr3s_aa[i],"TRB")==0) {
      to_remove=c(to_remove,i)
    }
  }
  z=y[-to_remove,]
  
  #prepare contig file by merging same clonotype rows info and join alpha and beta info together
  contig=fread(paste0(TCR_folder,myname,"_TCR/outs/all_contig_annotations.csv"))
  contig=contig[contig$raw_clonotype_id %in% z$clonotype_id,]
  contig=contig[(contig$raw_clonotype_id!="None") & (contig$cdr3!="None"),]
  contig$VJ=paste(contig$v_gene,contig$j_gene,sep="_")
  contig2=contig[,c("barcode","chain","v_gene","j_gene","VJ","cdr3","cdr3_nt","raw_clonotype_id")] %>%
    group_by(barcode) %>%
    summarise(chain = paste(chain, collapse="."),
              v_gene = paste(v_gene, collapse="."),
              #  j_gene = paste(j_gene, collapse="."),
              VJ = paste(VJ, collapse="."),
              cdr3 = paste(cdr3, collapse="."),
              cdr3_nt = paste(cdr3_nt, collapse="."),
              raw_clonotype_id=unique(raw_clonotype_id))
  contig2=merge(contig2,z[,c("clonotype_id","frequency")],by.x=c("raw_clonotype_id"),by.y=c("clonotype_id"))
  return(contig2)
}

prepare_TCR_data_mckiff<- function(myname,TCR_folder) {
  y=fread(paste0(TCR_folder,myname,"_TCR/outs/clonotypes.csv"))
  to_remove=c()
  #remove TCRs which:
  #1) are singletons with more than 1 TCR alpha and more than 1 TCR beta
  #2) have more than 2 alphas and 2 betas
  #3) have more than one beta and no alpha or more than one alpha and no beta
  
  for (i in 1:nrow(y)) {
    if (str_count(y$cdr3s_aa[i],"TRB") > 1 & str_count(y$cdr3s_aa[i],"TRA") > 1 & y$frequency[i]==1) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRB") > 2 | str_count(y$cdr3s_aa[i],"TRA")> 2) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRB") > 1 & str_count(y$cdr3s_aa[i],"TRA")==0) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRA") > 1 & str_count(y$cdr3s_aa[i],"TRB")==0) {
      to_remove=c(to_remove,i)
    }
  }
  z=y[-to_remove,]
  
  #prepare contig file by merging same clonotype rows info and join alpha and beta info together
  contig=fread(paste0(TCR_folder,myname,"_TCR/outs/all_contig_annotations.csv"))
  contig=contig[contig$raw_clonotype_id %in% z$clonotype_id,]
  contig=contig[(contig$raw_clonotype_id!="None") & (contig$cdr3!="None"),]
  contig$VJ=paste(contig$v_gene,contig$j_gene,sep="_")
  contig2=contig[,c("barcode","chain","v_gene","j_gene","VJ","cdr3","cdr3_nt","raw_clonotype_id")] %>%
    group_by(barcode) %>%
    summarise(chain = paste(chain, collapse="."),
              v_gene = paste(v_gene, collapse="."),
              #  j_gene = paste(j_gene, collapse="."),
              VJ = paste(VJ, collapse="."),
              cdr3 = paste(cdr3, collapse="."),
              cdr3_nt = paste(cdr3_nt, collapse="."),
              raw_clonotype_id=unique(raw_clonotype_id))
  contig2=merge(contig2,z[,c("clonotype_id","frequency")],by.x=c("raw_clonotype_id"),by.y=c("clonotype_id"))
  return(contig2)
}


prepare_TCR_data_bacher<- function(myname,TCR_folder) {
  y=fread(paste0(TCR_folder,myname,"_TCR/clonotypes.csv"))
  to_remove=c()
  #remove TCRs which:
  #1) are singletons with more than 1 TCR alpha and more than 1 TCR beta
  #2) have more than 2 alphas and 2 betas
  #3) have more than one beta and no alpha or more than one alpha and no beta
  
  for (i in 1:nrow(y)) {
    if (str_count(y$cdr3s_aa[i],"TRB") > 1 & str_count(y$cdr3s_aa[i],"TRA") > 1 & y$frequency[i]==1) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRB") > 2 | str_count(y$cdr3s_aa[i],"TRA")> 2) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRB") > 1 & str_count(y$cdr3s_aa[i],"TRA")==0) {
      to_remove=c(to_remove,i)
    } else if (str_count(y$cdr3s_aa[i],"TRA") > 1 & str_count(y$cdr3s_aa[i],"TRB")==0) {
      to_remove=c(to_remove,i)
    }
  }
  z=y[-to_remove,]
  
  #prepare contig file by merging same clonotype rows info and join alpha and beta info together
  contig=fread(paste0(TCR_folder,myname,"_TCR/all_contig_annotations.csv"))
  contig=contig[contig$raw_clonotype_id %in% z$clonotype_id,]
  contig=contig[(contig$raw_clonotype_id!="None") & (contig$cdr3!="None"),]
  contig$VJ=paste(contig$v_gene,contig$j_gene,sep="_")
  contig2=contig[,c("barcode","chain","v_gene","j_gene","VJ","cdr3","cdr3_nt","raw_clonotype_id")] %>%
    group_by(barcode) %>%
    summarise(chain = paste(chain, collapse="."),
              v_gene = paste(v_gene, collapse="."),
              #  j_gene = paste(j_gene, collapse="."),
              VJ = paste(VJ, collapse="."),
              cdr3 = paste(cdr3, collapse="."),
              cdr3_nt = paste(cdr3_nt, collapse="."),
              raw_clonotype_id=unique(raw_clonotype_id))
  contig2=merge(contig2,z[,c("clonotype_id","frequency")],by.x=c("raw_clonotype_id"),by.y=c("clonotype_id"))
  return(contig2)
}

#prepare gene expression data from mckiff
prepare_data_mckiff<-function(myname, metadata=NA, remove_TCR=FALSE, with_TCR=TRUE, processing_cite=FALSE, 
                              sample_name_short="default", single_sample_annotation,
                              GE_folder="", 
                              TCR_folder="",
                              output_folder="") {
  #define folders
  input_folder=GE_folder
  dirname <- paste0(input_folder)
  counts_matrix_filename = paste0(dirname)
  
  #########################################
  #filter TCR information for duplets and prepare it to add on seurat metadata
  if (with_TCR==TRUE) {contig2<-prepare_TCR_data(myname,TCR_folder=TCR_folder)}
  ####################################################
  #create seurat object
  print("loading data")
  covid_1 <- Read10X(data.dir = counts_matrix_filename)  # Seurat function to read in 10x count data
  if (length(covid_1)==2) {
    rna=covid_1$`Gene Expression`
    cite=covid_1$`Antibody Capture`
    seurat <- CreateSeuratObject(counts = rna, project = sample_name_short, min.cells = 5, min.features = 200) 
    if (processing_cite==TRUE) {
      cite=cite[,colnames(seurat)]
      seurat[["HTO"]] <- CreateAssayObject(counts =cite)
      seurat<-NormalizeData(object = seurat, assay = "HTO", normalization.method = "CLR")
    } 
    
  } else {
    seurat <- CreateSeuratObject(counts = covid_1, project = sample_name_short, min.cells = 5, min.features = 200) 
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat,pattern = "^MT-")
    seurat <- subset(seurat, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 5) 
  }
  if (myname %in% names(single_sample_annotation)) {
    my_annotation=single_sample_annotation[[myname]]
    my_annotation$cell_barcode2=gsub("[0-9]","", my_annotation$cell_barcode2)
    seurat$barcode=gsub("-1","",row.names(seurat@meta.data))
    seurat2<-subset(seurat, cells=paste0(my_annotation$cell_barcode2,"-1"))
    my_annotation=my_annotation[my_annotation$cell_barcode2 %in% seurat2@meta.data$barcode,]
    my_annotation$cell_barcode2=factor(my_annotation$cell_barcode2,levels=seurat2@meta.data$barcode)
    my_annotation=my_annotation[order(my_annotation$cell_barcode2),]
    seurat2@meta.data$origlib=my_annotation$origlib
    seurat2@meta.data$orig.project_id=my_annotation$orig.project_id
    seurat2@meta.data$orig.run=my_annotation$orig.run
    seurat2@meta.data$orig.peptide=my_annotation$orig.peptide
    seurat2@meta.data$orig.stim_time=my_annotation$orig.stim_time
    seurat2@meta.data$orig.HT_ID.global=my_annotation$orig.HT_ID.global
    seurat2@meta.data$orig.donor=my_annotation$orig.donor
    seurat2@meta.data$orig.sex=my_annotation$orig.sex
    seurat2@meta.data$orig.hospital=my_annotation$orig.hospital
    seurat2@meta.data$orig.severity=my_annotation$orig.severity
  } else {seurat2=seurat}
  
  
  #if necessary, remove TCR gene expression information
  #this may be necessary to cluster cells on pure gene expression independently from clonotype information
  if (remove_TCR==TRUE) {
    counts <- GetAssayData(seurat2, assay = "RNA")
    counts <- counts[-which(grepl( "TRAC|TRBC|TRBV|TRAV|TRGV|TRDV",rownames(counts ))),]
    seurat2 <- subset(seurat2, features = rownames(counts))
  }
  #add variables to metadata
  print("seurat object created and cleaned")
  
  #keep only cells with TCR info
  if (with_TCR==TRUE) {
    #seurat2@meta.data=merge(seurat2@meta.data,contig2,by.x=0,by.y="barcode",all.x=T)
    cell_barcodes=contig2$barcode
    seurat2<-subset(seurat2, cells= cell_barcodes)
  }
  
  #fill possible NA created where TCR info is missing
  seurat2@meta.data[is.na(seurat2@meta.data)]<-"None"
  Idents(seurat2)<-"orig.HT_ID.global"
  seurat2<-subset(seurat2, ident=c("Singlet"))
  print("attributes added")
  return(seurat2)
}

#prepare gene expression data from bacher
prepare_data_bacher<-function(myname, metadata, remove_TCR=FALSE, with_TCR=TRUE,
                              GE_folder="", 
                              local="",
                              TCR_folder="",
                              output_folder="",sample_name_short) {
  input_folder=paste0(GE_folder, myname)
  dirname <- paste0(local,input_folder)
  counts_matrix_filename = paste0(dirname)
  
  #########################################
  if (with_TCR==TRUE) {contig2<-prepare_TCR_data(myname,TCR_folder=TCR_folder)}
  ####################################################
  print("loading data")
  covid_1 <- Read10X(data.dir = counts_matrix_filename)  # Seurat function to read in 10x count data
  seurat <- CreateSeuratObject(counts = covid_1, project = sample_name_short, min.cells = 5, min.features = 200) 
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat,pattern = "^MT-")
  seurat <- subset(seurat, subset = nFeature_RNA > 400 & nFeature_RNA < 3000 & percent.mt < 5) 
  if (remove_TCR==TRUE) {
    counts <- GetAssayData(seurat, assay = "RNA")
    counts <- counts[-which(grepl( "TRAC|TRBC|TRBV|TRAV|TRGV|TRDV",rownames(counts ))),]
    seurat <- subset(seurat, features = rownames(counts))
  }
  print("seurat object created and cleaned")
  seurat@meta.data$orig.donor <- rep(myname, nrow(seurat@meta.data))
  seurat@meta.data$dataset <- rep("internal", nrow(seurat@meta.data))
  seurat@meta.data$orig.hospital=rep(metadata[metadata$patient==myname,]$Hospitalized,nrow(seurat@meta.data))
  seurat@meta.data$orig.project_id=rep(metadata[metadata$patient==myname,]$Experiment,nrow(seurat@meta.data))
  seurat@meta.data$diagnosis=rep(metadata[metadata$patient==myname,]$diagnosis,nrow(seurat@meta.data))
  seurat@meta.data$orig.run=rep(metadata[metadata$patient==myname,]$Sequencing_run,nrow(seurat@meta.data))
  
  #keep only cells with TCR info
  if (with_TCR==TRUE) {
    #seurat@meta.data=merge(seurat@meta.data,contig2,by.x=0,by.y="barcode",all.x=T)
    cell_barcodes=contig2$barcode
    seurat<-subset(seurat, cells= cell_barcodes)
  }
  seurat@meta.data[is.na(seurat@meta.data)]<-"None"
  print("attributes added")
  # saveRDS(seurat, paste0(local,output_folder,myname,".rds"))
  return(seurat)
}
########################################################################