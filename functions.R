library(data.table)
library(igraph)
library(stringdist)
library(stringr)
library(ggplot2)

#MISC_FUNS...
most_pop<-function(vec){
  names(sort(table(vec),decreasing = T))[1]  
}

write.tsv<-function(x,fname,rows=F){
  write.table(x,sep="\t",quote = F,row.names = rows,file = fname)
}

return_occurs<-function(CDR3_vfam_seqs)
{res_stat<-list()
for (i in 1:length(CDR3_vfam_seqs)){
  if(i%%100==0)print(i)
  query=CDR3_vfam_seqs[i]
  occur_tbl<-unique(lookup[CDR3_vfam%in%query,.(freq=sum(`frequencyCount (%)`),num=.N,read=sum(`count (templates/reads)`),ranks=min(rank)),.(donor)])
  occur_tbl_full<-na_to0(merge(joined_meta,occur_tbl,by="donor",all=T))
  res_stat[[i]]<-occur_tbl_full
}
res_stat
}

get_occur_for_query<-function(CDR3_vfam_seqs){
  occur_tbl<-unique(lookup[CDR3_vfam%in%CDR3_vfam_seqs,.(freq=sum(`frequencyCount (%)`),num=.N,read=sum(`count (templates/reads)`),ranks=min(rank)),.(donor)])
  occur_tbl_full<-na_to0(merge(joined_meta,occur_tbl,by="donor",all=T))
  occur_tbl_full
}

get_occur_for_query_short<-function(CDR3_vfam_seqs){
  occur_tbl<-unique(lookup[CDR3_vfam%in%CDR3_vfam_seqs,.(num=.N),.(donor)])
  occur_tbl_full<-na_to0(merge(joined_meta,occur_tbl,by="donor",all=T))
  occur_tbl_full
}

std <- function(x) sd(x)/sqrt(length(x))

analyse_occurs<-function(occur_list){
  cov<-cbind(t(sapply(occur_list,function(x)x[num>0&ranks<10001,table(c("COVID","control",study))-1,])),t(sapply(occur_list,function(x)x[num>0,table(c("COVID","control",study))-1,])), # num donors covid/control
             t(sapply(occur_list,function(x)x[,mean(freq),study][,V1,])),
             t(sapply(occur_list,function(x)x[,std(freq),study][,V1,]))) #mean freq in covid and control
  colnames(cov)<-c("occur_control_10K","occur_COVID_10K","occur_control","occur_COVID","mean_freq_COVID","mean_freq_control","SE_freq_COVID","SE_freq_control")
  cov<-as.data.frame(cov,stringsAsFactors = F)
  #cov$fisher.test10K
  fisher_p_10K<-numeric(nrow(cov))
  fisher_p<-numeric(nrow(cov))
  for (i in 1:nrow(cov)){
    fisher_p_10K[i]<-fisher.test(matrix(c(cov$occur_control_10K[i],786-cov$occur_control_10K[i],cov$occur_COVID_10K[i],1414-cov$occur_COVID_10K[i]),ncol=2),alternative="less")$p.value
    fisher_p[i]<-fisher.test(matrix(c(cov$occur_control[i],786-cov$occur_control[i],cov$occur_COVID[i],1414-cov$occur_COVID[i]),ncol=2),alternative="less")$p.value
  }
  cov$fisher_p_val_10K<-fisher_p_10K
  cov$fisher_p_val<-fisher_p
  cov
  #1414 covid, 786 control
}

V_family_to_adaptive<-function(Vname){
  nums<-gsub("TRBV","",sapply(strsplit(Vname,split = "-"),"[[",1))
  paste0("TCRBV",str_pad(string = nums,side = "left",pad = "0",width = 2))
}

na_to0<-function (x) {
  x[is.na(x)]<-0
  x
}

make.sequence.graph<-function (.data, .name = '',max_errs=1) {
  G <- graph.empty(n = length(.data), directed=F)
  tmp<-stringdistmatrix(.data,.data,method="hamming")
  G <- add.edges(G, t(which(tmp<=max_errs,arr.ind=T)))
  G <- igraph::simplify(G)
  G <- set.vertex.attribute(G, 'label', V(G), .data)
  #print(G)
  G
}

split_clones<-function(df){
  seu_exp<-df[rep(1:nrow(df),times=sapply(strsplit(df$VJ,split=".",fixed=T),length)),] 
  seu_exp$type=unlist(strsplit(df$chain,split=".",fixed=T))
  seu_exp$CDR3VJ=paste0(unlist(strsplit(df$cdr3,split=".",fixed=T)),"_",unlist(strsplit(df$VJ,split=".",fixed=T)))
  seu_exp$CDR3<-sapply(strsplit(seu_exp$CDR3VJ,split="_",fixed=T),"[[",1)
  seu_exp$bestVGene<-sapply(strsplit(seu_exp$CDR3VJ,split="_",fixed=T),"[[",2)
  seu_exp$bestJGene<-sapply(strsplit(paste0(seu_exp$CDR3VJ,"_none"),,split="_",fixed=T),"[[",3)
  seu_exp 
}

#move to functions
all_other_letters<-function(str,ind=8){
  aa<-c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", 
        "P", "Q", "R", "S", "T", "V", "W", "Y")
  paste0(substr(str,1,ind-1),aa,substr(str,ind+1,nchar(str)))
}

#move to functions
all_other_variants_one_mismatch<-function(str){
  unique(as.vector(sapply(2:(nchar(str)-1),all_other_letters,str=str)))
}

#move to functions
do_mira1<-function(df,mira_df,peptide_df){
  df$CDR3_VJ<-paste0(df$CDR3.amino.acid.sequence,"_",df$bestVGene,"_",df$bestVGene,"_",df$bestJGene,"_",df$bestJGene)
  miraind1vj<-list()
  nei<-list()
  for (i in 1:nrow(df))
    nei[[i]]<-paste0(all_other_variants_one_mismatch(df$CDR3.amino.acid.sequence[i]),"_",df$bestVGene[i],"_",df$bestVGene[i],"_",df$bestJGene[i],"_",df$bestJGene[i])
  nei<-unlist(nei)
  miramerge1vj<-mira_df$CDR3_VJ%in%nei
  mira1_mergeVJ<-rbind(df[,.(CDR3=CDR3.amino.acid.sequence,CDR3_VJ,donor=NA,peptide=NA),],mira_df[miramerge1vj,][,.(CDR3=aaSeqCDR3,CDR3_VJ,donor=Experiment,peptide=`Amino Acids`),])
  print("search done")
  gr<-make.sequence.graph(mira1_mergeVJ[,CDR3_VJ,],max_errs = 1)
  mira1_mergeVJ$clust=clusters(gr)$membership
  end_of_ours<-nrow(df)
  egr<-get.edgelist(gr)
  pepmat<-cbind(mira1_mergeVJ[egr[egr[,1]<=end_of_ours&egr[,2]>end_of_ours,1],,],mira1_mergeVJ[egr[egr[,1]<=end_of_ours&egr[,2]>end_of_ours,2],peptide,])
  tbl<-as.matrix(table(pepmat$CDR3_VJ,pepmat$V2))
  flt<-(tbl>1)+1-1
  fltt<-flt[,colSums(flt)>0]
  flttt<-fltt[rowSums(fltt)>0,]
  tmp<-character()
  for (i in 1:nrow(flttt))
    tmp[[i]]<-paste0(names(flttt[i,])[flttt[i,]!=0],collapse="|")
  mtc<-match(df$CDR3_VJ,row.names(flttt))
  df$mira_match<-unlist(tmp)[mtc]
  peptide_df[,mira_match:=`Amino Acids`,]
  setkey(peptide_df,mira_match)
  df$pep_source=peptide_df[df$mira_match,ORF,]
  df
}


#R_TCRdist call.
path_cpp="../conga/tcrdist_cpp/bin/find_neighbors"
db_path="../conga/tcrdist_cpp/db/tcrdist_info_human.txt"
perm<-readLines("../conga/tcrdist_cpp/db/tcrdist_info_human.txt")
permv<-c(unique(gsub("VBdist ","",sapply(strsplit(perm[grepl("VBdist",perm)],"*",fixed=T),"[[",1))),unique(gsub("VAdist ","",sapply(strsplit(perm[grepl("VAdist",perm)],"*",fixed=T),"[[",1))))
out_pref="tmp"
library(igraph)

TCRdistcpp<-function(tbl_tcr,allele=T)
{
  if(allele)
  {write.tsv(tbl_tcr[,.(va=paste0(conga_va,"*01"),cdr3a=conga_cdr3a,vb=paste0(conga_vb,"*01"),cdr3b=conga_cdr3b),],fname="tmp_tcr.tsv")}
  else{
    write.tsv(tbl_tcr[,.(va=conga_va,cdr3a=conga_cdr3a,vb=conga_vb,cdr3b=conga_cdr3b),],fname="tmp_tcr.tsv")
  }
  cmd<-paste0(path_cpp," -f ","tmp_tcr.tsv"," --only_tcrdists -d ",db_path, " -o ", out_pref)
  system(cmd)
  as.matrix(fread(paste0(out_pref,"_tcrdists.txt")))
}



filter_TCRdiste<-function(df){
  df<-df[!grepl("_",conga_cdr3a,fixed=T),,]
  df<-df[!grepl("*",conga_cdr3a,fixed=T),,]
  df<-df[conga_va%in%permv,,]
  df<-df[!grepl("_",conga_cdr3b,fixed=T),,]
  df<-df[!grepl("*",conga_cdr3b,fixed=T),,]
  df<-df[conga_vb%in%permv,,]
  df
}

make.sequence.graph_tcrdist<-function (.data, .name = '',allele_id=T,max_dist=120) 
{
  G <- graph.empty(n = nrow(.data), directed=F)
  distmat<-TCRdistcpp(.data,allele = allele_id)
  G <- add.edges(G, t(which(distmat<=max_dist,arr.ind=T)))
  G <- igraph::simplify(G)
  G <- set.vertex.attribute(G, 'label', V(G), paste0(.data$conga_va,"_",.data$conga_cdr3a,"_",.data$conga_vb,"_","_",.data$conga_cdr3b))
  G
}
