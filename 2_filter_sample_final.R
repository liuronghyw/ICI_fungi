##filter samples 

setwd("D:/课题/细菌—真菌互作/3_R/")
rm(list=ls())

library(tidyverse)
library(ggplot2)
library(MMUPHin)##批次效应校正

##提取细菌、真菌、病毒、古菌的reads####
raw_reads<-function(cohort,micro){
  archaea<-read.csv(paste0("3.mirco_data/",cohort,"/Archaea.S.bracken_mpa.kreport2"),check.names = F,row.names = 1,sep="\t")
  bacteria<-read.csv(paste0("3.mirco_data/",cohort,"/Bacteria.S.bracken_mpa.kreport2"),check.names = F,row.names = 1,sep="\t")
  fungi<-read.csv(paste0("3.mirco_data/",cohort,"/Fungi.S.bracken_mpa.kreport2"),check.names = F,row.names = 1,sep="\t")
  Viruses<-read.csv(paste0("3.mirco_data/",cohort,"/Viruses.S.bracken_mpa.kreport2"),check.names = F,row.names = 1,sep="\t")
  dat<-rbind(archaea,bacteria,fungi,Viruses)
  dat<-dat[c("k__Archaea","k__Bacteria","k__Eukaryota","k__Viruses"),]
  colnames(dat)<-gsub(".S.bracken.kreport2","",colnames(dat))
  rownames(dat)<-gsub("k_","reads_",rownames(dat))
  micro<-dat
}


##分数据集
cli<-read.csv("2.rawData/clinical/final/cli_whole_samples_new.csv",check.names = F,row.names = 1)
table(cli$dataset)
micr_63<-raw_reads("PRJEB22863")
micr_19<-raw_reads("PRJEB43119")
micr_06<-raw_reads("PRJNA397906")
micr_42<-raw_reads("PRJNA399742")
micr_81<-raw_reads("PRJNA541981")
micr_97<-raw_reads("PRJNA1023797")
micr_60<-raw_reads("PRJNA762360")
micr_95<-raw_reads("PRJNA770295")

microb<-cbind(micr_06,micr_19,micr_42,micr_60,micr_63,micr_81,micr_97,micr_95)
microb<-data.frame(t(microb))
microb<-microb[rownames(cli),]
microb$all_reads<-apply(microb,1,sum)
cli_sub<-cbind(cli,microb)
cli_sub$fungi_ratio<-cli_sub$reads__Eukaryota/cli_sub$all_reads
cli_sub$bac_ratio<-cli_sub$reads__Bacteria/cli_sub$all_reads
cli_sub$log_fungi<-log10(cli_sub$reads__Eukaryota)
cli_sub$log_bac<-log10(cli_sub$reads__Bacteria)
cli_sub$log_fungi_bac<-log10(cli_sub$reads__Eukaryota/cli_sub$reads__Bacteria)
cli_sub$log_fungi_all<-log10(cli_sub$reads__Eukaryota/cli_sub$all_reads)

write.csv(cli_sub,"2.rawData/clinical/final/cli_whole_samples_ratio.csv")

##过滤标准0：reads大于1M###
cli_sub<-subset(cli_sub,cli_sub$all_reads>=1000000)##过滤30个样本


##过滤标准1：真菌reads数在1%——0.005%####
cli_sub<-subset(cli_sub,cli_sub$fungi_ratio>=0.0001&cli_sub$fungi_ratio<=0.01)##过滤155个样本（0.01-0.0001标准共118个样本）
#cli_sub$fungi_ratio>=0.0001&
##过滤标准2：比对fungi数大于5%(或真菌数要达到20个）####
##melanoma####
num_fungi<-function(cohort,num){
  fungi<-read.csv(paste0("3.mirco_data/",cohort,"/Fungi.S.bracken_mpa.kreport2_species.txt"),check.names = F,sep="\t")
  colnames(fungi)<-gsub(".S.bracken.kreport2","",colnames(fungi))
  num<-fungi
}

fungi_63<-num_fungi("PRJEB22863")
fungi_19<-num_fungi("PRJEB43119")
fungi_06<-num_fungi("PRJNA397906")
fungi_42<-num_fungi("PRJNA399742")
fungi_81<-num_fungi("PRJNA541981")
fungi_97<-num_fungi("PRJNA1023797")
fungi_60<-num_fungi("PRJNA762360")
fungi_95<-num_fungi("PRJNA770295")


comb_fungi<-merge(fungi_19,fungi_95,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_06,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_42,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_81,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_60,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_97,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_63,by="#Classification",all = T)

rownames(comb_fungi)<-comb_fungi$`#Classification`
comb_fungi<-comb_fungi[,-1]

fungi_count<-comb_fungi[,rownames(cli)]###后面计算每个菌多少样本中比对到
fungi_count[is.na(fungi_count)]<-0
fungi_count<-fungi_count[,colSums(fungi_count)!=0]##有两个样本未比对到真菌

comb_fungi<-comb_fungi[,rownames(cli_sub)]

comb_fungi[is.na(comb_fungi)]<-0

comb_fungi<-subset(comb_fungi,rowSums(comb_fungi>0)!=0)
comb_fungi<-as.data.frame(t(comb_fungi))
#comb_fungi_mela$ratio<-rowSums(comb_fungi_mela>0)/917##917种真菌
comb_fungi$ratio<-rowSums(comb_fungi>0)
quantile(comb_fungi$ratio,probs = 0.05)

comb_fungi<-subset(comb_fungi,comb_fungi$ratio>=0)##过滤72个样本

##bacteria文件
num_bac<-function(cohort,num){
  bac<-read.csv(paste0("3.mirco_data/",cohort,"/Bacteria.S.bracken_mpa.kreport2_species.txt"),check.names = F,sep="\t")
  colnames(bac)<-gsub(".S.bracken.kreport2","",colnames(bac))
  num<-bac
}


bac_63<-num_bac("PRJEB22863")
bac_19<-num_bac("PRJEB43119")
bac_06<-num_bac("PRJNA397906")
bac_42<-num_bac("PRJNA399742")
bac_81<-num_bac("PRJNA541981")
bac_97<-num_bac("PRJNA1023797")
bac_60<-num_bac("PRJNA762360")
bac_95<-num_bac("PRJNA770295")


comb_bac<-merge(bac_19,bac_95,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_06,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_42,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_81,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_60,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_97,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_63,by="#Classification",all = T)

rownames(comb_bac)<-comb_bac$`#Classification`
comb_bac<-comb_bac[,-1]

bac_count<-comb_bac[,rownames(cli)]###后面计算每个菌多少样本中比对到
bac_count[is.na(bac_count)]<-0
bac_count<-bac_count[,colnames(fungi_count)]

comb_bac<-comb_bac[,rownames(comb_fungi)]
comb_bac[is.na(comb_bac)]<-0
comb_bac<-subset(comb_bac,rowSums(comb_bac>0)!=0)##8784种细菌

comb_bac<-as.data.frame(t(comb_bac))
#comb_bac_mela$ratio<-rowSums(comb_bac_mela>0)/8769
comb_bac$ratio<-rowSums(comb_bac>0)
quantile(comb_bac$ratio,probs = 0.05)
comb_bac<-subset(comb_bac,comb_bac$ratio>=0)##1个样本

##过滤标准3：若该样本最大相对丰度大于70%，视为疑似污染,暂定不需要####
##fungi
comb_fungi<-comb_fungi[rownames(comb_bac),][,-ncol(comb_fungi)]
fungi_rela<- sweep(comb_fungi, 1, rowSums(comb_fungi), FUN = "/")
fungi_rela$max<-apply(fungi_rela,1,max)

fungi_rela$filter<-ifelse(fungi_rela$max>=1,"filtered","retained")
fungi_rela<-subset(fungi_rela,fungi_rela$filter=="retained")#过滤12个样本
fungi_rela<-fungi_rela[,1:(ncol(fungi_rela)-2)]##去除max和filter列

##bacteria
comb_bac<-comb_bac[rownames(fungi_rela),][,-ncol(comb_bac)]
bac_rela<- sweep(comb_bac, 1, rowSums(comb_bac), FUN = "/")
bac_rela$max<-apply(bac_rela,1,max)

bac_rela$filter<-ifelse(bac_rela$max>=1,"filtered","retained")

bac_rela<-subset(bac_rela,bac_rela$filter=="retained")#过滤1个样本
bac_rela<-bac_rela[,1:(ncol(bac_rela)-2)]

fungi_rela<-fungi_rela[rownames(bac_rela),]

##输出过滤后临床信息####
cli_final<-cli[rownames(fungi_rela),]
cli_final_ratio<-cli_sub[rownames(fungi_rela),]
cli_mela<-subset(cli_final,cancer_type=="Melanoma")
cli_NSCLC<-subset(cli_final,cancer_type=="NSCLC")
cli_RCC<-subset(cli_final,cancer_type=="RCC")

write.csv(cli_final,"2.rawData/clinical/final/filter/cli_whole_filter.csv",row.names = T)
write.csv(cli_final_ratio,"2.rawData/clinical/final/filter/cli_whole_filter_ratio.csv",row.names = T)
write.csv(cli_mela,"2.rawData/clinical/final/filter/cli_mela_filter.csv",row.names = T)
write.csv(cli_NSCLC,"2.rawData/clinical/final/filter/cli_NSCLC_filter.csv",row.names = T)
write.csv(cli_RCC,"2.rawData/clinical/final/filter/cli_RCC_filter.csv",row.names = T)

###输出fungi文件####
table(cli_final$dataset)
fungi<-function(cancer,cohort){
  fungi<-t(fungi_rela[rownames(subset(cli_final,cancer_type==cancer&dataset==cohort)),])
  write.csv(fungi,paste0("5.fungi_data/",cancer,"/relabundance/species/","fungi_",cancer,"_",cohort,"_rela_filter.csv"))
}

fungi("NSCLC","PRJEB22863")
fungi("NSCLC","PRJNA1023797")
fungi("RCC","PRJEB22863")
fungi("Melanoma","PRJEB43119")
fungi("Melanoma","PRJNA397906")
fungi("Melanoma","PRJNA399742")
fungi("Melanoma","PRJNA541981")
fungi("Melanoma","PRJNA762360")
fungi("Melanoma","PRJNA770295")


fungi_mela<-t(fungi_rela[rownames(cli_mela),])
write.csv(fungi_mela,"5.fungi_data/Melanoma/relabundance/species/fungi_mela_whole_rela_filter.csv")

fungi_NSCLC<-t(fungi_rela[rownames(cli_NSCLC),])
write.csv(fungi_NSCLC,"5.fungi_data/NSCLC/relabundance/species/fungi_NSCLC_whole_rela_filter.csv")

fungi_RCC<-t(fungi_rela[rownames(cli_RCC),])
write.csv(fungi_RCC,"5.fungi_data/RCC/relabundance/species/fungi_RCC_whole_rela_filter.csv")

write.csv(fungi_rela,"5.fungi_data/fungi_whole_rela_species_filter.csv")
comb_fungi<-comb_fungi[rownames(fungi_rela),]
write.csv(comb_fungi,"5.fungi_data/fungi_whole_raw_species_filter.csv")

##原始丰度
fungi_raw<-function(cancer,cohort){
  fungi<-t(comb_fungi[rownames(subset(cli_final,cancer_type==cancer&dataset==cohort)),])
  write.csv(fungi,paste0("5.fungi_data/",cancer,"/raw/species/","fungi_",cancer,"_",cohort,"_raw_filter.csv"))
}

fungi_raw("NSCLC","PRJEB22863")
fungi_raw("NSCLC","PRJNA1023797")
fungi_raw("RCC","PRJEB22863")
fungi_raw("Melanoma","PRJEB43119")
fungi_raw("Melanoma","PRJNA397906")
fungi_raw("Melanoma","PRJNA399742")
fungi_raw("Melanoma","PRJNA541981")
fungi_raw("Melanoma","PRJNA762360")
fungi_raw("Melanoma","PRJNA770295")

###输出bacteria文件####
table(cli_final$dataset)
bac<-function(cancer,cohort){
  bac<-t(bac_rela[rownames(subset(cli_final,cancer_type==cancer&dataset==cohort)),])
  write.csv(bac,paste0("4.Bacteria_data/",cancer,"/relabundance/species/","bac_",cancer,"_",cohort,"_rela_filter.csv"))
}

bac("NSCLC","PRJEB22863")
bac("NSCLC","PRJNA1023797")
bac("RCC","PRJEB22863")
bac("Melanoma","PRJEB43119")
bac("Melanoma","PRJNA397906")
bac("Melanoma","PRJNA399742")
bac("Melanoma","PRJNA541981")
bac("Melanoma","PRJNA762360")
bac("Melanoma","PRJNA770295")

bac_mela<-t(bac_rela[rownames(cli_mela),])
write.csv(bac_mela,"4.Bacteria_data/Melanoma/relabundance/species/bac_mela_whole_rela_filter.csv")

bac_NSCLC<-t(bac_rela[rownames(cli_NSCLC),])
write.csv(bac_NSCLC,"4.Bacteria_data/NSCLC/relabundance/species/bac_NSCLC_whole_rela_filter.csv")

bac_RCC<-t(bac_rela[rownames(cli_RCC),])
write.csv(bac_RCC,"4.Bacteria_data/RCC/relabundance/species/bac_RCC_whole_rela_filter.csv")

write.csv(bac_rela,"4.Bacteria_data/bac_whole_rela_species_filter.csv")
comb_bac<-comb_bac[rownames(bac_rela),]
write.csv(comb_bac,"4.Bacteria_data/bac_whole_raw_species_filter.csv")

##原始丰度
bac_raw<-function(cancer,cohort){
  bac<-t(comb_bac[rownames(subset(cli_final,cancer_type==cancer&dataset==cohort)),])
  write.csv(bac,paste0("4.Bacteria_data/",cancer,"/raw/species/","bac_",cancer,"_",cohort,"_raw_filter.csv"))
}

bac_raw("NSCLC","PRJEB22863")
bac_raw("NSCLC","PRJNA1023797")
bac_raw("RCC","PRJEB22863")
bac_raw("Melanoma","PRJEB43119")
bac_raw("Melanoma","PRJNA397906")
bac_raw("Melanoma","PRJNA399742")
bac_raw("Melanoma","PRJNA541981")
bac_raw("Melanoma","PRJNA762360")
bac_raw("Melanoma","PRJNA770295")

##微生物批次效应校正####
##fungi
cli_mela$study<-as.factor(cli_mela$study)
fungi_mela_adj<- adjust_batch(feature_abd = fungi_mela,
                                 batch = "study",
                                 data = cli_mela,
                                 control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

fungi_mela_mmuphin <- fungi_mela_adj$feature_abd_adj
write.csv(fungi_mela_mmuphin,"5.fungi_data/Melanoma/relabundance/species/MMUPHIN/fungi_mela_whole_rela_filter_MMUPHIN.csv")

cli_NSCLC$study<-as.factor(cli_NSCLC$study)
fungi_NSCLC_adj<- adjust_batch(feature_abd = fungi_NSCLC,
                              batch = "study",
                              data = cli_NSCLC,
                              control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

fungi_NSCLC_mmuphin <- fungi_NSCLC_adj$feature_abd_adj
write.csv(fungi_NSCLC_mmuphin,"5.fungi_data/NSCLC/relabundance/species/MMUPHIN/fungi_NSCLC_whole_rela_filter_MMUPHIN.csv")

cli_final$study<-as.factor(cli_final$study)
fungi_rela<-t(fungi_rela)
fungi_adj<- adjust_batch(feature_abd = fungi_rela,
                               batch = "study",
                               data = cli_final,
                               control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

fungi_mmuphin <- fungi_adj$feature_abd_adj
write.csv(fungi_mmuphin,"5.fungi_data/fungi_whole_rela_species_filter_MMUPHIN.csv")

##bacteria
bac_mela_adj<- adjust_batch(feature_abd = bac_mela,
                              batch = "study",
                              data = cli_mela,
                              control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

bac_mela_mmuphin <- bac_mela_adj$feature_abd_adj
write.csv(bac_mela_mmuphin,"4.bacteria_data/Melanoma/relabundance/species/MMUPHIN/bac_mela_whole_rela_filter_MMUPHIN.csv")


bac_NSCLC_adj<- adjust_batch(feature_abd = bac_NSCLC,
                               batch = "study",
                               data = cli_NSCLC,
                               control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

bac_NSCLC_mmuphin <- bac_NSCLC_adj$feature_abd_adj
write.csv(bac_NSCLC_mmuphin,"4.bacteria_data/NSCLC/relabundance/species/MMUPHIN/bac_NSCLC_whole_rela_filter_MMUPHIN.csv")

bac_rela<-t(bac_rela)
bac_adj<- adjust_batch(feature_abd = bac_rela,
                         batch = "study",
                         data = cli_final,
                         control = list(verbose = FALSE))##丰度文件行为微生物；临床文件行为样本

bac_mmuphin <- bac_adj$feature_abd_adj
write.csv(bac_mmuphin,"4.bacteria_data/bac_whole_rela_species_filter_MMUPHIN.csv")

