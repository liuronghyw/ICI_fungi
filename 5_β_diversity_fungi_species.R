setwd("D:/课题/细菌—真菌互作/3_R/")
rm(list=ls())


library(vegan)
library(ggplot2)
library(ggrepel)
library(ape)
library(ggforce)


#读取数据(特征表)
βDiversity<-function(raw,cli,datout){
  ##导入相对丰度及临床信息文件
  df<- read.csv(raw, header = T, check.names = F,row.names = 1)
  df<-as.data.frame(t(df))##行为样本
  group<-read.csv(cli, row.names = 1, header = TRUE, check.names = FALSE)
  #group<-subset(group,!is.na(response_code))
  df<-df[rownames(group),]##行为样本
  #计算距离
  distance<-vegdist(df,method='bray')
  
  pcoa<- cmdscale(distance,eig=TRUE)
  
  #提取前两个分类解释
  plot_data<-data.frame({pcoa$point})[1:2]
  
  #前两个分类解释命名
  names(plot_data)[1:2]<-c('PCoA1','PCoA2') 
  plot_data$eig=pcoa$eig
  
  data<-plot_data[match(rownames(group),rownames(plot_data)),]
  data<-data.frame(group,plot_data)
  datout<-data
}


##melanoma####
mela<-βDiversity("5.fungi_data/melanoma/relabundance/species/MMUPHIN/fungi_mela_whole_rela_filter_MMUPHIN.csv",
                 "2.rawData/clinical/final/filter/cli_mela_filter.csv")
mela$response_code<-gsub("0","NR",mela$response_code)
mela$response_code<-gsub("1","R",mela$response_code)



##response
mela_res<-subset(mela,!is.na(response_code)&!is.na(gender)&!is.na(BMI_group))
mela_raw<-as.data.frame(t(read.csv("5.fungi_data/melanoma/relabundance/species/MMUPHIN/fungi_mela_whole_rela_filter_MMUPHIN.csv", header = T, check.names = F,row.names = 1)))
mela_raw<-mela_raw[rownames(mela_res),]
mela_res$response_code<-as.factor(mela_res$response_code)


set.seed(123)
dune.div<-adonis2(mela_raw~response_code+gender+study+BMI_group,data=mela_res,permutations=999,method="bray")
dune_adonis<-paste0("PERMANOVA P = ",round(dune.div$`Pr(>F)`,3))

p1<-ggplot(mela_res,aes(x=PCoA1,y=PCoA2,shape=response_code,color=response_code))+
  geom_point(alpha=1,size=1)+
  stat_ellipse(level=0.95,size=0.5)+
  labs(x=paste("PCoA1(",format(100*mela$eig[1]/sum(mela$eig),digits = 4),"%)",sep=""),
       y=paste("PCoA2(",format(100*mela$eig[2]/sum(mela$eig),digits = 4),"%)",sep=""),
       title="Melanoma",
       subtitle = dune_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  scale_color_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  theme(panel.background = element_rect(fill='white',colour = 'black'),
        axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y = element_text(colour = 'black',size=12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.8,0.13),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
        plot.title = element_text(size=14,hjust = 0.5,colour = 'black'),
        plot.subtitle =element_text(size=13,hjust = 0.5,face = "italic",colour = 'black') )
p1

ggsave("10.pics/supplemental figures/β-diversity/species/mela_species_β_diversity_pcoa_res.tiff",plot=p1,width = 100,height = 100,units = "mm",dpi=300)
##协变量对pcoa的影响
##pcoa1
##cohort
p2<-ggplot(mela,aes(study,PCoA1,fill=study))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Cohort")+
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=2.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p2

##age

mela_age<-subset(mela,!is.na(age_group))
p3<-ggplot(mela_age,aes(age_group,PCoA1,fill=age_group))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Age")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)),  # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.4,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p3

##BMI
for(i in 1:nrow(mela)){
  if(is.na(mela$BMI[i])){mela$BMI_group[i]=NA}
  else if(mela$BMI[i]<25){mela$BMI_group[i]="lean"}
  else if(mela$BMI[i]>=25&mela$BMI[i]<=30){mela$BMI_group[i]="overweight"}
  else if(mela$BMI[i]>30){mela$BMI_group[i]="obese"}
}
mela_BMI<-subset(mela,!is.na(BMI_group))
mela_BMI$BMI_group<-factor(mela_BMI$BMI_group,levels = c("lean","overweight","obese"))
p4<-ggplot(mela_BMI,aes(BMI_group,PCoA1,fill=BMI_group))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="BMI")+
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.9,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle =60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p4

##country
p5<-ggplot(mela,aes(country,PCoA1,fill=country))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Country")+
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P  ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=1.3,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p5
##gender
mela_gender<-subset(mela,!is.na(gender))
p6<-ggplot(mela_gender,aes(gender,PCoA1,fill=gender))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Gender")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=0.4,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p6

p7<-p2+p3+p4+p5+p6+plot_layout(ncol = 5) 
p7

##pcoa2
##cohort
p8<-ggplot(mela,aes(study,PCoA2,fill=study))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Cohort")+
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=1.8,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p8

##age
p9<-ggplot(mela_age,aes(age_group,PCoA2,fill=age_group))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Age")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)),  # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.4,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p9

##BMI

p10<-ggplot(mela_BMI,aes(BMI_group,PCoA2,fill=BMI_group))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="BMI")+
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.9,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p10

##country
p11<-ggplot(mela,aes(country,PCoA2,fill=country))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Country")+
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=1.3,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p11
##gender
p12<-ggplot(mela_gender,aes(gender,PCoA2,fill=gender))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Gender")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=0.4,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p12

p13<-p8+p9+p10+p11+p12+plot_layout(ncol = 5) 
p13

p13_sub<-p2+p3+p4+p6+p8+p9+p10+p12+plot_layout(ncol = 4,nrow = 2) 
p13_sub
ggsave("10.pics/supplemental figures/β-diversity/species/mela_species_pcoa1_confounder_boxplot.tiff",plot=p7,width = 200,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/β-diversity/species/mela_species_pcoa2_confounder_boxplot.tiff",plot=p13,width = 200,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/β-diversity/species/mela_species_pcoa_confounder_boxplot_patchwork.tiff",plot=p13_sub,width = 150,height =250,units = "mm",dpi=300)

##NSCLC####
NSCLC<-βDiversity("5.fungi_data/NSCLC/relabundance/species/MMUPHIN/fungi_NSCLC_whole_rela_filter_MMUPHIN.csv",
                  "2.rawData/clinical/final/filter/cli_NSCLC_filter.csv")
NSCLC$response_code<-ifelse(NSCLC$response_code==0,"NR","R")
NSCLC<-subset(NSCLC,!is.na(response_code)&!is.na(age_group)&!is.na(gender))
##response
NSCLC_raw<-as.data.frame(t(read.csv("5.fungi_data/NSCLC/relabundance/species/MMUPHIN/fungi_NSCLC_whole_rela_filter_MMUPHIN.csv", header = T, check.names = F,row.names = 1)))
NSCLC$response_code<-as.factor(NSCLC$response_code)
NSCLC_raw<-NSCLC_raw[rownames(NSCLC),]

set.seed(123)
dune.div<-adonis2(NSCLC_raw~response_code+study+gender+age_group,data=NSCLC,permutations=999,method="bray")
dune_adonis<-paste0("PERMANOVA P = ",round(dune.div$`Pr(>F)`,2))

p14<-ggplot(NSCLC,aes(x=PCoA1,y=PCoA2,shape=response_code,color=response_code))+
  geom_point(alpha=1,size=1)+
  stat_ellipse(level=0.95,size=0.5)+
  labs(x=paste("PCoA1(",format(100*NSCLC$eig[1]/sum(NSCLC$eig),digits = 4),"%)",sep=""),
       y=paste("PCoA2(",format(100*NSCLC$eig[2]/sum(NSCLC$eig),digits = 4),"%)",sep=""),
       title="NSCLC",
       subtitle = dune_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  scale_color_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  theme(panel.background = element_rect(fill='white',colour = 'black'),
        axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y = element_text(colour = 'black',size=12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.15,0.2),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
        plot.title = element_text(size=14,hjust = 0.5),
        plot.subtitle =element_text(size=13,hjust = 0.5,face = "italic",colour = 'black') )
p14
ggsave("10.pics/supplemental figures/β-diversity/species/NSCLC_species_β_diversity_pcoa_res.tiff",plot=p14,width = 100,height = 100,units = "mm",dpi=300)
##协变量对pcoa的影响
##pcoa1
##cohort
p15<-ggplot(NSCLC,aes(study,PCoA1,fill=study))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Cohort")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p15

##age
for(i in 1:nrow(NSCLC)){
  if(is.na(NSCLC$age_bin[i])){NSCLC$age_group[i]=NA}
  else if(NSCLC$age_bin[i]=="15-20"|NSCLC$age_bin[i]=="20-25"|
          NSCLC$age_bin[i]=="25-30"|NSCLC$age_bin[i]=="30-35"|
          NSCLC$age_bin[i]=="35-40"|NSCLC$age_bin[i]=="40-45"|
          NSCLC$age_bin[i]=="45-50"|NSCLC$age_bin[i]=="50-55"|
          NSCLC$age_bin[i]=="55-60")
  {NSCLC$age_group[i]="< 60"}
  else if(NSCLC$age_bin[i]=="60-65"|NSCLC$age_bin[i]=="65-70"|
          NSCLC$age_bin[i]=="70-75"|NSCLC$age_bin[i]=="75-80"|
          NSCLC$age_bin[i]=="80-85"|NSCLC$age_bin[i]=="85-90"|
          NSCLC$age_bin[i]=="90-95")
  {NSCLC$age_group[i]="≥ 60"}
}

NSCLC_age<-subset(NSCLC,!is.na(age_group))
p16<-ggplot(NSCLC_age,aes(age_group,PCoA1,fill=age_group))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Age")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)),  # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p16

##gender
NSCLC_gender<-subset(NSCLC,!is.na(gender))
p17<-ggplot(NSCLC_gender,aes(gender,PCoA1,fill=gender))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Gender")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p17

p18<-p15+p16+p17+plot_layout(ncol = 3) 
p18

##pcoa2
##cohort
p19<-ggplot(NSCLC,aes(study,PCoA2,fill=study))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Cohort")+
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle =60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p19

##age
p20<-ggplot(NSCLC_age,aes(age_group,PCoA2,fill=age_group))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Age")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)),  # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p20

##gender
p21<-ggplot(NSCLC_gender,aes(gender,PCoA2,fill=gender))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Gender")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p21

p22<-p19+p20+p21+plot_layout(ncol = 3) 
p22

p22_sub<-p15+p16+p17+p19+p20+p21+plot_layout(ncol = 3,nrow =2)
p22_sub
ggsave("10.pics/supplemental figures/β-diversity/species/NSCLC_species_pcoa1_confounder_boxplot.tiff",plot=p18,width = 120,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/β-diversity/species/NSCLC_species_pcoa2_confounder_boxplot.tiff",plot=p22,width = 120,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/β-diversity/species/NSCLC_species_pcoa_confounder_boxplot_patchwork.tiff",plot=p22_sub,width = 100,height =250,units = "mm",dpi=300)

##RCC####
RCC<-βDiversity("5.fungi_data/RCC/relabundance/species/fungi_RCC_PRJEB22863_rela_filter.csv",
                "2.rawData/clinical/final/filter/cli_RCC_filter.csv")
RCC$response_code<-ifelse(RCC$response_code==0,"NR","R")

##response
RCC_raw<-as.data.frame(t(read.csv("5.fungi_data/RCC/relabundance/species/fungi_RCC_whole_rela_filter.csv", header = T, check.names = F,row.names = 1)))
RCC$response_code<-as.factor(RCC$response_code)

set.seed(123)
dune.div<-adonis2(RCC_raw~response_code+gender+age_group,data=RCC,permutations=999,method="bray")
dune_adonis<-paste0("PERMANOVA P = ",round(dune.div$`Pr(>F)`,2))

p23<-ggplot(RCC,aes(x=PCoA1,y=PCoA2,shape=response_code,color=response_code))+
  geom_point(alpha=1,size=1)+
  stat_ellipse(level=0.95,size=0.5)+
  labs(x=paste("PCoA1(",format(100*RCC$eig[1]/sum(RCC$eig),digits = 4),"%)",sep=""),
       y=paste("PCoA2(",format(100*RCC$eig[2]/sum(RCC$eig),digits = 4),"%)",sep=""),
       title="RCC",
       subtitle = dune_adonis)+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  scale_color_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  theme(panel.background = element_rect(fill='white',colour = 'black'),
        axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y = element_text(colour = 'black',size=12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.85,0.85),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
        plot.title = element_text(size=14,hjust = 0.5),
        plot.subtitle =element_text(size=13,hjust = 0.5,face = "italic",colour = 'black') )
p23
ggsave("10.pics/supplemental figures/β-diversity/species/RCC_species_β_diversity_pcoa_res.tiff",plot=p23,width = 100,height = 100,units = "mm",dpi=300)
##协变量对pcoa的影响
##pcoa1
##age

RCC_age<-subset(RCC,!is.na(age_group))
p24<-ggplot(RCC_age,aes(age_group,PCoA1,fill=age_group))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Age")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)),  # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p24

##gender
RCC_gender<-subset(RCC,!is.na(gender))
p25<-ggplot(RCC_gender,aes(gender,PCoA1,fill=gender))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Gender")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p25

p26<-p24+p25+plot_layout(ncol = 2) 
p26

##pcoa2
##age
p27<-ggplot(RCC_age,aes(age_group,PCoA2,fill=age_group))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Age")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)),  # 显示p值
                     size = 4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))

p27

##gender
p28<-ggplot(RCC_gender,aes(gender,PCoA2,fill=gender))+
  geom_boxplot(width=0.3,size=0.5,outlier.color =NA,alpha=0.8)+
  theme_bw()+
  labs(title ="Gender")+
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =4,          # 修改字体大小
                     position = position_nudge(x=0.2,y = 0))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 60,hjust = 1),
        legend.position ="none" ,
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(hjust=0.5,size=12),
        legend.title = element_text(size=12))
p28

p29<-p27+p28+plot_layout(ncol =2) 
p29
p29_sub<-p24+p25+p27+p28+plot_layout(ncol =2,nrow = 2) 
p29_sub

ggsave("10.pics/supplemental figures/β-diversity/species/RCC_species_pcoa1_confounder_boxplot.tiff",plot=p26,width = 80,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/β-diversity/species/RCC_species_pcoa2_confounder_boxplot.tiff",plot=p29,width = 80,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/β-diversity/species/RCC_species_pcoa_confounder_boxplot_patchwork.tiff",plot=p29_sub,width = 80,height = 250,units = "mm",dpi=300)

write.csv(mela,"6.microbial_diversity/β_diversity/mela_fungi_species_β_diversity.csv")
write.csv(NSCLC,"6.microbial_diversity/β_diversity/NSCLC_fungi_species_β_diversity.csv")
write.csv(RCC,"6.microbial_diversity/β_diversity/RCC_fungi_species_β_diversity.csv")

p30<-p1+p14+p23
p30
ggsave("10.pics/supplemental figures/β-diversity/species/species_mela_NSCLC_β_diversity.tiff",plot=p30,width = 300,height = 100,units = "mm",dpi=300)
