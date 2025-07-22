setwd("D:/课题/细菌—真菌互作/3_R/")
rm(list=ls())

library(tidyverse)
library(MMUPHin)##批次效应校正

cli<-read.csv("2.rawData/clinical/final/filter/cli_whole_filter.csv",check.names = F,row.names = 1)
table(cli$dataset)

##分数据集计算相对丰度
rela<-function(cli,cancer,cohort){
  fungi<-read.csv(paste0("3.mirco_data/",cohort,"/Fungi.S.bracken_mpa.kreport2_genus.txt"),sep = "\t",row.names = 1,check.names = F)
  bac<-read.csv(paste0("3.mirco_data/",cohort,"/Bacteria.S.bracken_mpa.kreport2_genus.txt"),sep = "\t",row.names = 1,check.names = F)
  
  colnames(fungi)<-gsub(".S.bracken.kreport2","",colnames(fungi))
  colnames(bac)<-gsub(".S.bracken.kreport2","",colnames(bac))
  
  metadata<-cli
  metadata<-subset(metadata,cancer_type==cancer&dataset==cohort)
  fungi<-fungi[,rownames(metadata)]
  bac<-bac[,rownames(metadata)]
  
  rela_fungi<- sweep(fungi,2, colSums(fungi), FUN = "/")##行为微生物
  write.csv(fungi,paste0("5.fungi_data/",cancer,"/raw/genus/","fungi_",cancer,"_",cohort,"_raw_filter.csv"))
  write.csv(rela_fungi,paste0("5.fungi_data/",cancer,"/relabundance/genus/","fungi_",cancer,"_",cohort,"_rela_filter.csv"))
  
  rela_bac<- sweep(bac,2, colSums(bac), FUN = "/")##行为微生物
  write.csv(bac,paste0("4.bacteria_data/",cancer,"/raw/genus/","bac_",cancer,"_",cohort,"_raw_filter.csv"))
  write.csv(rela_bac,paste0("4.bacteria_data/",cancer,"/relabundance/genus/","bac_",cancer,"_",cohort,"_rela_filter.csv"))
  
}


rela(cli,"RCC","PRJEB22863")
rela(cli,"NSCLC","PRJEB22863")
rela(cli,"NSCLC","PRJNA1023797")
rela(cli,"Melanoma","PRJEB43119")
rela(cli,"Melanoma","PRJNA397906")
rela(cli,"Melanoma","PRJNA399742")
rela(cli,"Melanoma","PRJNA541981")
rela(cli,"Melanoma","PRJNA762360")
rela(cli,"Melanoma","PRJNA770295")

##文件合并
##relabundance
rela_merge<-function(text1,text2,cancer,cohort,dataout){
  dat1<-read.csv(paste0(text1,cancer,text2,cancer,"_",cohort,"_rela_filter.csv"),check.names = F)
  colnames(dat1)[1]<-"id"
  dataout<-dat1
}

##fungi
dat_rcc<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","RCC","PRJEB22863")
dat_63<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","NSCLC","PRJEB22863")
dat_97<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","NSCLC","PRJNA1023797")
dat_19<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","Melanoma","PRJEB43119")
dat_06<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","Melanoma","PRJNA397906")
dat_42<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","Melanoma","PRJNA399742")
dat_81<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","Melanoma","PRJNA541981")
dat_60<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","Melanoma","PRJNA762360")
dat_95<-rela_merge("5.fungi_data/","/relabundance/genus/fungi_","Melanoma","PRJNA770295")

dat<-merge(dat_06,dat_19,by="id",all = T)
dat<-merge(dat,dat_42,by="id",all = T)
dat<-merge(dat,dat_60,by="id",all = T)
dat<-merge(dat,dat_63,by="id",all = T)
dat<-merge(dat,dat_81,by="id",all = T)
dat<-merge(dat,dat_97,by="id",all = T)
dat<-merge(dat,dat_95,by="id",all = T)
dat<-merge(dat,dat_rcc,by="id",all = T)
rownames(dat)<-dat$id
dat<-dat[,-1]
dat[is.na(dat)]<-0
write.csv(dat,"5.fungi_data/fungi_whole_rela_genus_filter.csv")

##批次效应校正
cli$study<-as.factor(cli$study)
dat<-dat[,rownames(cli)]
fungi_adj<- adjust_batch(feature_abd = dat,
                         batch = "study",
                         data = cli,
                         control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

fungi_mmuphin <- fungi_adj$feature_abd_adj
write.csv(fungi_mmuphin,"5.fungi_data/fungi_whole_rela_genus_filter_MMUPHIN.csv")

cli_mela<-subset(cli,cancer_type=="Melanoma")
mela<-dat[,rownames(cli_mela)]
write.csv(mela,"5.fungi_data/Melanoma/relabundance/genus/fungi_mela_rela_genus_filter.csv")

##批次效应校正
cli_mela$dataset<-as.factor(cli_mela$dataset)
mela_adj<- adjust_batch(feature_abd = mela,
                         batch = "dataset",
                         data = cli_mela,
                         control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

mela_mmuphin <- mela_adj$feature_abd_adj
write.csv(mela_mmuphin,"5.fungi_data/Melanoma/relabundance/genus/MMUPHIN/fungi_mela_rela_genus_filter_MMUPHIN.csv")


cli_NSCLC<-subset(cli,cancer_type=="NSCLC")
NSCLC<-dat[,rownames(cli_NSCLC)]
write.csv(NSCLC,"5.fungi_data/NSCLC/relabundance/genus/fungi_NSCLC_rela_genus_filter.csv")

##批次效应校正
cli_NSCLC$dataset<-as.factor(cli_NSCLC$dataset)
NSCLC_adj<- adjust_batch(feature_abd = NSCLC,
                        batch = "dataset",
                        data = cli_NSCLC,
                        control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

NSCLC_mmuphin <- NSCLC_adj$feature_abd_adj
write.csv(NSCLC_mmuphin,"5.fungi_data/NSCLC/relabundance/genus/MMUPHIN/fungi_NSCLC_rela_genus_filter_MMUPHIN.csv")

RCC<-dat[,rownames(subset(cli,cancer_type=="RCC"))]
write.csv(RCC,"5.fungi_data/RCC/relabundance/genus/fungi_RCC_rela_genus_filter.csv")

##bacteria
dat_rcc<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","RCC","PRJEB22863")
dat_63<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","NSCLC","PRJEB22863")
dat_97<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","NSCLC","PRJNA1023797")
dat_19<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","Melanoma","PRJEB43119")
dat_06<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","Melanoma","PRJNA397906")
dat_42<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","Melanoma","PRJNA399742")
dat_81<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","Melanoma","PRJNA541981")
dat_60<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","Melanoma","PRJNA762360")
dat_95<-rela_merge("4.bacteria_data/","/relabundance/genus/bac_","Melanoma","PRJNA770295")

dat<-merge(dat_06,dat_19,by="id",all = T)
dat<-merge(dat,dat_42,by="id",all = T)
dat<-merge(dat,dat_60,by="id",all = T)
dat<-merge(dat,dat_63,by="id",all = T)
dat<-merge(dat,dat_81,by="id",all = T)
dat<-merge(dat,dat_97,by="id",all = T)
dat<-merge(dat,dat_95,by="id",all = T)
dat<-merge(dat,dat_rcc,by="id",all = T)
rownames(dat)<-dat$id
dat<-dat[,-1]
dat[is.na(dat)]<-0
write.csv(dat,"4.bacteria_data/bac_whole_rela_genus_filter.csv")

##批次效应校正
dat<-dat[,rownames(cli)]
bac_adj<- adjust_batch(feature_abd = dat,
                         batch = "study",
                         data = cli,
                         control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

bac_mmuphin <- bac_adj$feature_abd_adj
write.csv(bac_mmuphin,"4.bacteria_data/bac_whole_rela_genus_filter_MMUPHIN.csv")

mela<-dat[,rownames(cli_mela)]
write.csv(mela,"4.bacteria_data/melanoma/relabundance/genus/bac_mela_rela_genus_filter.csv")

##批次效应校正
mela_adj<- adjust_batch(feature_abd = mela,
                        batch = "dataset",
                        data = cli_mela,
                        control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

mela_mmuphin <- mela_adj$feature_abd_adj
write.csv(mela_mmuphin,"4.bacteria_data/melanoma/relabundance/genus/MMUPHIN/bac_mela_rela_genus_filter_MMUPHIN.csv")

NSCLC<-dat[,rownames(cli_NSCLC)]
write.csv(NSCLC,"4.bacteria_data/NSCLC/relabundance/genus/bac_NSCLC_rela_genus_filter.csv")

##批次效应校正
NSCLC_adj<- adjust_batch(feature_abd = NSCLC,
                         batch = "dataset",
                         data = cli_NSCLC,
                         control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

NSCLC_mmuphin <- NSCLC_adj$feature_abd_adj
write.csv(NSCLC_mmuphin,"4.bacteria_data/NSCLC/relabundance/genus/MMUPHIN/bac_NSCLC_rela_genus_filter_MMUPHIN.csv")

RCC<-dat[,rownames(subset(cli,cancer_type=="RCC"))]
write.csv(RCC,"4.bacteria_data/RCC/relabundance/genus/bac_RCC_rela_genus_filter.csv")


##raw
raw_merge<-function(text1,text2,cancer,cohort,dataout){
  dat1<-read.csv(paste0(text1,cancer,text2,cancer,"_",cohort,"_raw_filter.csv"),check.names = F)
  colnames(dat1)[1]<-"id"
  dataout<-dat1
}
##fungi
dat_rcc<-raw_merge("5.fungi_data/","/raw/genus/fungi_","RCC","PRJEB22863")
dat_63<-raw_merge("5.fungi_data/","/raw/genus/fungi_","NSCLC","PRJEB22863")
dat_97<-raw_merge("5.fungi_data/","/raw/genus/fungi_","NSCLC","PRJNA1023797")
dat_19<-raw_merge("5.fungi_data/","/raw/genus/fungi_","Melanoma","PRJEB43119")
dat_06<-raw_merge("5.fungi_data/","/raw/genus/fungi_","Melanoma","PRJNA397906")
dat_42<-raw_merge("5.fungi_data/","/raw/genus/fungi_","Melanoma","PRJNA399742")
dat_81<-raw_merge("5.fungi_data/","/raw/genus/fungi_","Melanoma","PRJNA541981")
dat_60<-raw_merge("5.fungi_data/","/raw/genus/fungi_","Melanoma","PRJNA762360")
dat_95<-raw_merge("5.fungi_data/","/raw/genus/fungi_","Melanoma","PRJNA770295")

dat<-merge(dat_06,dat_19,by="id",all = T)
dat<-merge(dat,dat_42,by="id",all = T)
dat<-merge(dat,dat_60,by="id",all = T)
dat<-merge(dat,dat_63,by="id",all = T)
dat<-merge(dat,dat_81,by="id",all = T)
dat<-merge(dat,dat_97,by="id",all = T)
dat<-merge(dat,dat_95,by="id",all = T)
dat<-merge(dat,dat_rcc,by="id",all = T)
rownames(dat)<-dat$id
dat<-dat[,-1]
dat[is.na(dat)]<-0
write.csv(dat,"5.fungi_data/fungi_whole_raw_genus_filter.csv")

mela<-dat[,rownames(subset(cli,cancer_type=="Melanoma"))]
write.csv(mela,"5.fungi_data/Melanoma/raw/genus/fungi_mela_raw_genus_filter.csv")

NSCLC<-dat[,rownames(subset(cli,cancer_type=="NSCLC"))]
write.csv(NSCLC,"5.fungi_data/NSCLC/raw/genus/fungi_NSCLC_raw_genus_filter.csv")

RCC<-dat[,rownames(subset(cli,cancer_type=="RCC"))]
write.csv(RCC,"5.fungi_data/RCC/raw/genus/fungi_RCC_raw_genus_filter.csv")

##bacteria
dat_rcc<-raw_merge("4.bacteria_data/","/raw/genus/bac_","RCC","PRJEB22863")
dat_63<-raw_merge("4.bacteria_data/","/raw/genus/bac_","NSCLC","PRJEB22863")
dat_97<-raw_merge("4.bacteria_data/","/raw/genus/bac_","NSCLC","PRJNA1023797")
dat_19<-raw_merge("4.bacteria_data/","/raw/genus/bac_","Melanoma","PRJEB43119")
dat_06<-raw_merge("4.bacteria_data/","/raw/genus/bac_","Melanoma","PRJNA397906")
dat_42<-raw_merge("4.bacteria_data/","/raw/genus/bac_","Melanoma","PRJNA399742")
dat_81<-raw_merge("4.bacteria_data/","/raw/genus/bac_","Melanoma","PRJNA541981")
dat_60<-raw_merge("4.bacteria_data/","/raw/genus/bac_","Melanoma","PRJNA762360")
dat_95<-raw_merge("4.bacteria_data/","/raw/genus/bac_","melanoma","PRJNA770295")

dat<-merge(dat_06,dat_19,by="id",all = T)
dat<-merge(dat,dat_42,by="id",all = T)
dat<-merge(dat,dat_60,by="id",all = T)
dat<-merge(dat,dat_63,by="id",all = T)
dat<-merge(dat,dat_81,by="id",all = T)
dat<-merge(dat,dat_97,by="id",all = T)
dat<-merge(dat,dat_95,by="id",all = T)
dat<-merge(dat,dat_rcc,by="id",all = T)
rownames(dat)<-dat$id
dat<-dat[,-1]
dat[is.na(dat)]<-0
write.csv(dat,"4.bacteria_data/bac_whole_raw_genus_filter.csv")

mela<-dat[,rownames(subset(cli,cancer_type=="Melanoma"))]
write.csv(mela,"4.bacteria_data/Melanoma/raw/genus/bac_mela_raw_genus_filter.csv")

NSCLC<-dat[,rownames(subset(cli,cancer_type=="NSCLC"))]
write.csv(NSCLC,"4.bacteria_data/NSCLC/raw/genus/bac_NSCLC_raw_genus_filter.csv")

RCC<-dat[,rownames(subset(cli,cancer_type=="RCC"))]
write.csv(RCC,"4.bacteria_data/RCC/raw/genus/bac_RCC_raw_genus_filter.csv")
