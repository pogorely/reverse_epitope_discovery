#CD4 COVID-19 TCR analysis. 
#edit this path to point to your conga tcrdist cpp distribution
path_tcrdist="../conga/tcrdist_cpp/"

skip_lookup_generation=T # skip lookup table generation which takes a lot of time and requires full download of Emerson et al and Snyder et al cohorts and load precomputed version from file instead. 

#paths to datasets (not included!), change for your system if you want to redo the full pipeline
path_to_Emerson_cohort<-"../Data/Emerson_et_al/"
path_to_Snyder_cohort<-"../Data/Snyder_et_al/"
  

source("functions.R")

seu_exp4_n <- split_clones(readRDS("mergedSeurat_metadata_TCR_UMAPpositions.rds"))
seu_exp4_n_beta<-seu_exp4_n[seu_exp4_n$type=="TRB",]

##generation of the lookup table to quickly count number of donors in COVID and Emerson cohorts matching th data 
if (!skip_lookup_generation){
 cdr3s<-unique(c(seu_exp4_n_beta$CDR3))
 emfls<-c(list.files(path_to_Snyder_cohort,pattern = ".tsv",full.names = T),list.files(path_to_Emerson_cohort,pattern = ".tsv.gz",full.names = T))
 emfls_s<-c(list.files(path_to_Snyder_cohort,pattern = ".tsv",full.names = F),list.files(path_to_Emerson_cohort,pattern = ".tsv.gz",full.names = F))
 lookup<-list()
 for (i in 1:length(emfls))
 {
   print(i)
   tmp<-fread(emfls[i])  
   tmp$donor<-emfls_s[i]
   lookup[[i]]<-tmp[aminoAcid%in%cdr3s,.(aminoAcid,vFamilyName,donor),]
 }
lookup<-do.call(rbind,lookup)
lookup$CDR3_vfam<-paste0(lookup$aminoAcid,"_",lookup$vFamilyName)
save(lookup, file="lookup_full.rda")
}

load("lookup_full.rda")
load("covid_control.rda")
lookup<-as.data.table(lookup)
lookup_occurs<-lookup[,.(covid_donors=sum(unique(donor)%in%covid),control_donors=sum(unique(donor)%in%control)),CDR3_vfam]
setkey(lookup_occurs,CDR3_vfam)

seu_exp4_n_beta$CDR3b_vfam<-paste0(seu_exp4_n_beta$CDR3,"_",V_family_to_adaptive(gsub("*01","",seu_exp4_n_beta$bestVGene,fixed=T)))
seu_exp4_n_beta$n_covid<-lookup_occurs[seu_exp4_n_beta$CDR3b_vfam,covid_donors,]
seu_exp4_n_beta$n_control<-lookup_occurs[seu_exp4_n_beta$CDR3b_vfam,control_donors,]
seu_exp4_n_beta<-as.data.table(seu_exp4_n_beta)

#add p-value from Fisher test
seu_exp4_n_beta_full<-seu_exp4_n_beta[!is.na(n_covid),p_value:=fisher.test(matrix(c(unique(n_covid),1414-unique(n_covid),unique(n_control),786-unique(n_control)),nrow=2)),.(CDR3b_vfam)]
seu_exp4_n_beta_full[,p_val_BH:=p.adjust(p_value,method="BH"),]

#covid enriched clonotypes: MIRA assignment
mira4<-fread("mira_mixcr_CD4.tsv")
mira4pep<-fread("peptide-hits-cii.csv")

#generation of unique clonotypes from data
seu_exp4_n_beta_full[,CDR3.amino.acid.sequence:=CDR3,]
seu_exp4_n_beta_full[,betaV:=bestVGene,]
mrd_n4<-do_mira1(seu_exp4_n_beta_full[p_val_BH<0.05&(n_covid>(n_control*2)),,],mira_df =mira4,peptide_df = mira4pep) 
mrd_n4s<-split_clones(mrd_n4)
mrd_n4sa<-mrd_n4s[type=="TRA"&chain%in%c("TRA.TRB","TRB.TRA"),,]
setnames(mrd_n4sa,c("CDR3","bestVGene","CDR3.amino.acid.sequence","betaV"),c("conga_cdr3a","conga_va","conga_cdr3b","conga_vb"))
mrd_n4sa$conga_va<-gsub("DV","/DV",mrd_n4sa$conga_va,fixed = T)
#unique TCR alpha/beta clonotypes compatible with TCRdist implementation
mrd_n4sau<-filter_TCRdiste(mrd_n4sa[!duplicated(paste0(conga_cdr3b,conga_cdr3a,conga_vb,conga_va)),])

#create similarity network using TCRdist
grb<-make.sequence.graph_tcrdist(mrd_n4sau,max_dist = 120)
grb<-set.vertex.attribute(grb,"mira_match",V(grb),mrd_n4sau$mira_match)
grb<-set.vertex.attribute(grb,"pep_source",V(grb),mrd_n4sau$pep_source)
grb<-set.vertex.attribute(grb,"donor",V(grb),mrd_n4sau$orig.donor)
grb<-set.vertex.attribute(grb,"betw",V(grb),betweenness(grb)>quantile(betweenness(grb),probs = 0.99))
grb<-set.vertex.attribute(grb, 'label', V(grb), clusters(grb)$membership)

#create filtered similarity network
grb_filt<-grb
grb_filt<-delete.vertices(grb_filt,which(betweenness(grb_filt)>quantile(betweenness(grb_filt),probs = 0.99)))
grb_filt<-set.vertex.attribute(grb_filt, 'label', V(grb_filt), clusters(grb_filt)$membership)

mrd_n4sau$betw<-betweenness(grb)
mrd_n4sau_filt<-mrd_n4sau[!(betweenness(grb)>quantile(betweenness(grb),probs = 0.99)),]
mrd_n4sau_filt$clusters=clusters(grb_filt)$membership
mrd_n4sau_filt$cl=clusters(grb_filt)$membership

mrd_n4sau_filt[,va:=conga_va,]
mrd_n4sau_filt[,ja:=bestJGene,]
mrd_n4sau_filt[,vb:=conga_vb,]
mrd_n4sau_filt[,jb:=sapply(strsplit(CDR3_VJ,"_"),"[[",4),]
mrd_n4sau_filt[chain=="TRA.TRB",cdr3a_nucseq:=sapply(strsplit(cdr3_nt,split = ".",fixed = T),"[[",1),]
mrd_n4sau_filt[chain=="TRA.TRB",cdr3b_nucseq:=sapply(strsplit(cdr3_nt,split = ".",fixed = T),"[[",2),]
mrd_n4sau_filt[chain=="TRB.TRA",cdr3a_nucseq:=sapply(strsplit(cdr3_nt,split = ".",fixed = T),"[[",2),]
mrd_n4sau_filt[chain=="TRB.TRA",cdr3b_nucseq:=sapply(strsplit(cdr3_nt,split = ".",fixed = T),"[[",1),]


### MIRA specificity assignment and generation of unique clonotype table for COVID-depleted clonotypes.
mrd_no4<-do_mira1(seu_exp4_n_beta_full[p_val_BH<0.05&(n_covid<(n_control*2)),,],mira_df =mira4,peptide_df = mira4pep)

mrd_no4s<-split_clones(mrd_no4)
mrd_no4sa<-mrd_no4s[type=="TRA"&chain%in%c("TRA.TRB","TRB.TRA"),,]
setnames(mrd_no4sa,c("CDR3","bestVGene","CDR3.amino.acid.sequence","betaV"),c("conga_cdr3a","conga_va","conga_cdr3b","conga_vb"))
mrd_no4sa$conga_va<-gsub("DV","/DV",mrd_no4sa$conga_va,fixed = T)
mrd_no4sau<-filter_TCRdiste(mrd_no4sa[!duplicated(paste0(conga_cdr3b,conga_cdr3a,conga_vb,conga_va)),])
mrd_no4sau[,va:=conga_va,]
mrd_no4sau[,ja:=bestJGene,]
mrd_no4sau[,vb:=conga_vb,]
mrd_no4sau[,jb:=sapply(strsplit(CDR3_VJ,"_"),"[[",4),]
mrd_no4sau[chain=="TRA.TRB",cdr3a_nucseq:=sapply(strsplit(cdr3_nt,split = ".",fixed = T),"[[",1),]
mrd_no4sau[chain=="TRA.TRB",cdr3b_nucseq:=sapply(strsplit(cdr3_nt,split = ".",fixed = T),"[[",2),]
mrd_no4sau[chain=="TRB.TRA",cdr3a_nucseq:=sapply(strsplit(cdr3_nt,split = ".",fixed = T),"[[",2),]
mrd_no4sau[chain=="TRB.TRA",cdr3b_nucseq:=sapply(strsplit(cdr3_nt,split = ".",fixed = T),"[[",1),]

##HLA association analysis

load("joined_meta.rda")
load("respred_mat.rda")

publics<-mrd_n4sau_filt
resl<-list()
for (pub in 1:length(publics$CDR3b_vfam)){
  if (pub%%100==0)print(pub)
  resf<-numeric(ncol(respred_mat))
  for (i in 1:ncol(respred_mat)){
    resf[i]<-fisher.test(table(get_occur_for_query_short(publics$CDR3b_vfam[pub])$num>0,respred_mat[,i]>0.2),alternative = "greater")$p.value
  }
  resl[[pub]]<-resf
}

reslm<-do.call(rbind,resl)
colnames(reslm)=colnames(respred_mat)
reslmBH<-matrix(p.adjust(reslm,method = "BH"),ncol=ncol(reslm),nrow=nrow(reslm))
colnames(reslmBH)=colnames(respred_mat)

publics$best_HLA=""
publics$best_pval=0
for (pub in 1:length(publics$CDR3b_vfam)){
  publics$best_HLA[pub]<-names(which.min(reslm[pub,]))
  publics$best_pval[pub]<-reslm[pub,][which.min(reslm[pub,])]
  publics$best_pval_BH[pub]<-reslmBH[pub,][which.min(reslmBH[pub,])]
}



SI2<-publics[,.(origlib,
                  orig.donor,
                  tcrdist120_cluster_id=cl,
                  va,
                  ja,
                  clone_size=frequency,
                  cdr3a=conga_cdr3a,
                  cdr3a_nt=cdr3a_nucseq,
                  vb,
                  jb,
                  cdr3b=conga_cdr3b,
                  cdr3b_nt=cdr3b_nucseq,
                  n_covid,
                  n_control,
                  p_val_covid_assoc=p_value,
                  p_val_covid_assoc_BH=p_val_BH,
                  best_HLA_assoc=best_HLA,
                  p_val_best_HLA_assoc=best_pval,
                  p_val_best_HLA_assoc_BH=best_pval_BH, 
                  MIRA_pool=mira_match,
                  MIRA_pool_source=pep_source),]

SI3<-mrd_no4sau[,.(origlib,
              orig.donor,
              va,
              ja,
              clone_size=frequency,
              cdr3a=conga_cdr3a,
              cdr3a_nt=cdr3a_nucseq,
              vb,
              jb,
              cdr3b=conga_cdr3b,
              cdr3b_nt=cdr3b_nucseq,
              n_covid,
              n_control,
              p_val_covid_assoc=p_value,
              p_val_covid_assoc_BH=p_val_BH,                  
              MIRA_pool=mira_match,
              MIRA_pool_source=pep_source),]
write.tsv(SI2,fname="SI_table2.tsv") 
write.tsv(SI3,fname="SI_table3.tsv")            


