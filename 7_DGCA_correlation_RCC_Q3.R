
setwd("D:/课题/细菌—真菌互作/3_R/")
rm(list = ls())

##加载包####
require(DT) # for databale
require(reshape2) # for dcast
require(corrplot)
require(DGCA)
require(ggthemes)  # for theme_calc()
library(tidyverse)
library(patchwork)
library(plyr)##ddply函数
library(aplot)##分组条形图拼图
library(ggpubr)##检验

# 微生物名属名缩写
abbreviate_microbial_name <- function(name) {
  # 使用strsplit拆分名称
  parts <- strsplit(name, " ")[[1]]
  
  # 处理第一个单词，保留首字母缩写
  abbreviated_first <- paste0(substr(parts[1], 1, 1), ".")
  
  # 组合首字母缩写和其余部分
  abbreviated_name <- paste(abbreviated_first, paste(parts[-1], collapse = " "))
  
  return(abbreviated_name)
}

##计算有显著性关系中某个菌出现的频次
freq<-function(data,df_node,dataout){
  a<-as.data.frame(table(data$Source))
  b<-as.data.frame(table(data$Target))
  c<-merge(a,b,by="Var1",all = T)
  c[is.na(c)]<-0
  c$Freq<-c$Freq.x+c$Freq.y
  dataout<-merge(df_node,c[,c(1,4)],by.x = "id",by.y="Var1")
}

##RCC
#仅采用response meta分析的结果来分析
##fungi
core_fungi_RCC<-read.csv("7.featureseletions/logistic_response/fungi/RCC/Q3/Logistic_PRJEB22863_RCC_response_fungi_pvalue.csv",header = T,check.names = F)
colnames(core_fungi_RCC)[1]<-"id"
core_fungi_RCC$id<-gsub("_"," ",core_fungi_RCC$id)
core_fungi_RCC<-subset(core_fungi_RCC,core_fungi_RCC$pvalue<0.05)
fungi_rela<-read.csv("5.fungi_data/RCC/relabundance/species/fungi_RCC_PRJEB22863_rela_filter.csv",header = T,check.names = F,row.names = 1)
rownames(fungi_rela)<-gsub("_"," ",rownames(fungi_rela))
core_fungi_rela<-fungi_rela[core_fungi_RCC$id,]
core_fungi_rela$type<-"Fungi"
core_fungi_rela$color<-"#A7BE68"

##bacteria
core_bac_RCC<-read.csv("7.featureseletions/logistic_response/bacteria/RCC/Q3/Logistic_PRJEB22863_RCC_response_bac_pvalue.csv",header = T,check.names = F)
colnames(core_bac_RCC)[1]<-"id"
core_bac_RCC$id<-gsub("_"," ",core_bac_RCC$id)
core_bac_RCC<-subset(core_bac_RCC,core_bac_RCC$pvalue<0.05)
bac_rela<-read.csv("4.bacteria_data/RCC/relabundance/species/bac_RCC_PRJEB22863_rela_filter.csv",header = T,check.names = F,row.names = 1)
rownames(bac_rela)<-gsub("_"," ",rownames(bac_rela))
core_bac_rela<-bac_rela[core_bac_RCC$id,]
core_bac_rela$type<-"Bacteria"
core_bac_rela$color<-"#EEB977"

##metadata
meta_RCC<-read.csv("2.rawData/clinical/final/filter/cli_RCC_filter.csv",header = T,row.names = 1)
RCC_fungi_bac<-rbind(core_fungi_rela,core_bac_rela)
table(meta_RCC$dataset)
RCC_fungi_bac_rela<-RCC_fungi_bac[,rownames(meta_RCC)]
rownames(RCC_fungi_bac_rela)<-gsub("\\[|\\]", "",rownames(RCC_fungi_bac_rela))
write.csv(RCC_fungi_bac_rela,"8.correlation/RCC/Q3/RCC_fungi_bac_core_rela_Q3.csv")

rownames(RCC_fungi_bac)<-gsub("\\[|\\]", "",rownames(RCC_fungi_bac))
RCC_fungi_bac$id<-rownames(RCC_fungi_bac)
RCC_fungi_bac$Label<-sapply(RCC_fungi_bac$id, abbreviate_microbial_name)
node<-RCC_fungi_bac[,c(ncol(RCC_fungi_bac)-1,ncol(RCC_fungi_bac)-3,ncol(RCC_fungi_bac)-2,ncol(RCC_fungi_bac))]
write.csv(node,"8.correlation/RCC/Q3/RCC_fungi_bac_node_Q3.csv",row.names = F)


##Comparing pairwise correlations across conditions####
##Runs the full discovery of differential correlation (ddcor) section 
#for comparing pairwise correlations across conditions in the Differential 
#Gene Correlation Analysis (DGCA) package.
##不同组别间相关性系数差异分析
response_design_mat <- lapply(meta_RCC[colnames(RCC_fungi_bac_rela), 'response_code'], function(x){
  if (x == "1") {
    return(c(1,0))
  }else{
    return(c(0,1))
  }
})%>%as.data.frame()%>%t()
colnames(response_design_mat)<-c("Response","Non-response")

##spearman                   
set.seed(123)
ddcor_res <- ddcorAll(inputMat = RCC_fungi_bac_rela, 
                      corrType = 'spearman', 
                      design = response_design_mat,#design必须是数值型的
                      nPerms = 20,
                      compare = c("Response","Non-response"))

##将两物种间分组
group_ddcor_index <- apply(ddcor_res[,c(1,2)], 1, function(x){
  if (sum(x %in% rownames(core_fungi_rela)) == 1) {
    return('intersect')
  }else if (sum(x %in% rownames(core_fungi_rela)) == 2){
    return('fungi')
  }else{
    return('bacteria')
  }
})


all_dgca <- ddcor_res
all_dgca$Response_pVal_adj <-p.adjust(all_dgca$Response_pVal, method = "BH")
all_dgca$'Non-response_pVal_adj' <-p.adjust(all_dgca$`Non-response_pVal`, method = "BH")
all_dgca$type <- group_ddcor_index
all_dgca$abs_zScoreDiff <- abs(all_dgca$zScoreDiff)#abs()取绝对值
all_dgca<-all_dgca[,c(1:4,12,5:6,13,7:11,14:15)]##列排下序
write.csv(all_dgca,"8.correlation/RCC/Q3/DCGA_spearman_RCC_Q3.csv",row.names=F)

all_dgca_R<-all_dgca[all_dgca$Response_pVal_adj<0.05,]
all_dgca_R<-all_dgca_R[,c(1:3,14)]
colnames(all_dgca_R)<-c("Source",	"Target",	"r","type")
all_dgca_R$sig<-ifelse(all_dgca_R$r>0,"P","N")
all_dgca_R$color<-ifelse(all_dgca_R$sig=="P","#BE68A7","#67A7BE")
all_dgca_R$abs_cor<-abs(all_dgca_R$r)
all_dgca_R<-all_dgca_R[all_dgca_R$abs_cor>=0.3,]
node_R<-node[union(all_dgca_R$Source,all_dgca_R$Target),]
node_R<-freq(all_dgca_R,node_R)

all_dgca_NR<-all_dgca[all_dgca$`Non-response_pVal_adj`<0.05,]
all_dgca_NR<-all_dgca_NR[,c(1:2,6,14)]
colnames(all_dgca_NR)<-c("Source",	"Target",	"r","type")
all_dgca_NR$sig<-ifelse(all_dgca_NR$r>0,"P","N")
all_dgca_NR$color<-ifelse(all_dgca_NR$sig=="P","#BE68A7","#67A7BE")
all_dgca_NR$abs_cor<-abs(all_dgca_NR$r)
all_dgca_NR<-all_dgca_NR[all_dgca_NR$abs_cor>=0.3,]
node_NR<-node[union(all_dgca_NR$Source,all_dgca_NR$Target),]
node_NR<-freq(all_dgca_NR,node_NR)

write.csv(all_dgca_R,"8.correlation/RCC/Q3/DCGA_spearman_RCC_R_Q3.csv",row.names=F)
write.csv(all_dgca_NR,"8.correlation/RCC/Q3/DCGA_spearman_RCC_NR_Q3.csv",row.names=F)

write.csv(node_R,"8.correlation/RCC/Q3/DCGA_spearman_RCC_R_node_Q3.csv",row.names=F)
write.csv(node_NR,"8.correlation/RCC/Q3/DCGA_spearman_RCC_NR__node_Q3.csv",row.names=F)

##pearson
set.seed(123)
ddcor_res_pearson <- ddcorAll(inputMat = RCC_fungi_bac_rela, 
                              corrType = 'pearson', 
                              design = response_design_mat,#design必须是数值型的
                              nPerms = 20,
                              compare = c("Response","Non-response"))

##将两物种间分组
group_ddcor_index_pearson <- apply(ddcor_res_pearson[,c(1,2)], 1, function(x){
  if (sum(x %in% rownames(core_fungi_rela)) == 1) {
    return('intersect')
  }else if (sum(x %in% rownames(core_fungi_rela)) == 2){
    return('fungi')
  }else{
    return('bacteria')
  }
})


all_dgca_pearson <- ddcor_res_pearson
all_dgca_pearson$Response_pVal_adj <-p.adjust(all_dgca_pearson$Response_pVal, method = "bonferroni")
all_dgca_pearson$'Non-response_pVal_adj' <-p.adjust(all_dgca_pearson$`Non-response_pVal`, method = "bonferroni")
all_dgca_pearson$type <- group_ddcor_index_pearson
all_dgca_pearson$abs_zScoreDiff <- abs(all_dgca_pearson$zScoreDiff)#abs()取绝对值
all_dgca_pearson<-all_dgca_pearson[,c(1:4,12,5:6,13,7:11,14:15)]##列排下序
write.csv(all_dgca_pearson,"8.correlation/RCC/Q3/DCGA_pearson_RCC_Q3.csv",row.names=F)

all_dgca_pearson_R<-all_dgca_pearson[all_dgca_pearson$Response_pVal_adj<0.05,]
all_dgca_pearson_R<-all_dgca_pearson_R[,c(1:3,14)]
colnames(all_dgca_pearson_R)<-c("Source",	"Target",	"r","type")
all_dgca_pearson_R$sig<-ifelse(all_dgca_pearson_R$r>0,"P","N")
all_dgca_pearson_R$color<-ifelse(all_dgca_pearson_R$sig=="P","#BE68A7","#67A7BE")
all_dgca_pearson_R$abs_cor<-abs(all_dgca_pearson_R$r)
all_dgca_pearson_R<-all_dgca_pearson_R[all_dgca_pearson_R$abs_cor>=0.3,]
node_pearson_R<-node[union(all_dgca_pearson_R$Source,all_dgca_pearson_R$Target),]
node_pearson_R<-freq(all_dgca_pearson_R,node_pearson_R)

all_dgca_pearson_NR<-all_dgca_pearson[all_dgca_pearson$`Non-response_pVal_adj`<0.05,]
all_dgca_pearson_NR<-all_dgca_pearson_NR[,c(1:2,6,14)]
colnames(all_dgca_pearson_NR)<-c("Source",	"Target",	"r","type")
all_dgca_pearson_NR$sig<-ifelse(all_dgca_pearson_NR$r>0,"P","N")
all_dgca_pearson_NR$color<-ifelse(all_dgca_pearson_NR$sig=="P","#BE68A7","#67A7BE")
all_dgca_pearson_NR$abs_cor<-abs(all_dgca_pearson_NR$r)
all_dgca_pearson_NR<-all_dgca_pearson_NR[all_dgca_pearson_NR$abs_cor>=0.3,]
node_pearson_NR<-node[union(all_dgca_pearson_NR$Source,all_dgca_pearson_NR$Target),]
node_pearson_NR<-freq(all_dgca_pearson_NR,node_pearson_NR)

write.csv(all_dgca_pearson_R,"8.correlation/RCC/Q3/DCGA_pearson_RCC_R_Q3.csv",row.names=F)
write.csv(all_dgca_pearson_NR,"8.correlation/RCC/Q3/DCGA_pearson_RCC_NR_Q3.csv",row.names=F)

write.csv(node_pearson_R,"8.correlation/RCC/Q3/DCGA_pearson_RCC_R_node_Q3.csv",row.names=F)
write.csv(node_pearson_NR,"8.correlation/RCC/Q3/DCGA_pearson_RCC_NR__node_Q3.csv",row.names=F)

