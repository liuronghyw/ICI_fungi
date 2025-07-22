rm(list=ls())
setwd("D:/课题/细菌—真菌互作/3_R/")

library(tidyverse)
require(otuSummary)
library(ggvenn)
library(stringr)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(patchwork)
library(gghalves) # Compose Half-Half Plots Using Your Favourite Geoms


cli<-read.csv("2.rawData/clinical/final/cli_whole_samples_ratio.csv",row.names = 1,check.names = F)

##α多样性
Alpha<-function(datain,dataout){
  datain<-read.csv(datain,header=T,check.names = F,row.names = 1)
  alpha_res <- alphaDiversity(datain, siteInCol = F,
                              taxhead = "taxonomy", threshold = 1, 
                              percent = FALSE, write = FALSE)
  allBio<-as.data.frame(alpha_res["allBio"])
  colnames(allBio)<-gsub("allBio.","",colnames(allBio))
  dataout<-allBio
}

##species水平####
##fungi
fungi_species<-Alpha("5.fungi_data/fungi_whole_raw_species_filter.csv")
cli<-select(cli,c("age_bin","gender","study","cancer_type","country","BMI_group","response_code","age_group"))
table(cli$age_bin)
cli<-cli[rownames(fungi_species),]

fungi_species_cli<-cbind(fungi_species,cli)


for(i in 1:nrow(fungi_species_cli)){
  if(is.na(fungi_species_cli$response_code[i])){fungi_species_cli$response_code[i]=NA}
  else if(fungi_species_cli$response_code[i]==1){fungi_species_cli$response_code[i]="R"}
  else if(fungi_species_cli$response_code[i]==0){fungi_species_cli$response_code[i]="NR"}
}



table(fungi_species_cli$study,fungi_species_cli$cancer_type)
for(i in 1:nrow(fungi_species_cli)){
  if(fungi_species_cli$cancer_type[i]=="RCC"){fungi_species_cli$study_sub[i]="RoutyB_2018(RCC)"}
  else if(fungi_species_cli$cancer_type[i]=="NSCLC"&fungi_species_cli$study[i]=="RoutyB_2018"){fungi_species_cli$study_sub[i]="RoutyB_2018(NSCLC)"}
  else if(fungi_species_cli$cancer_type[i]=="Melanoma"){fungi_species_cli$study_sub[i]=paste0(fungi_species_cli$study[i],"(Melanoma)")}
  else if(fungi_species_cli$cancer_type[i]=="NSCLC"&fungi_species_cli$study[i]=="DerosaL_2024"){fungi_species_cli$study_sub[i]="DerosaL_2024(NSCLC)"}
}

write.csv(fungi_species_cli,"6.microbial_diversity/α_diversity/fungi_α_diversity_species.csv")
##response####
##分数据集
res<-subset(fungi_species_cli,!is.na(response_code))
res$study_sub<-factor(res$study_sub,levels = c("FrankelAE_2017(Melanoma)","LeeKA_2022(Melanoma)","MatsonV_2018(Melanoma)",
                                               "McCullochJA_2022(Melanoma)","SpencerCN_2021(Melanoma)","DerosaL_2024(NSCLC)",
                                               "RoutyB_2018(NSCLC)","RoutyB_2018(RCC)"))

p1<-ggplot(res,aes(response_code,shannon,fill=response_code))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=response_code),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  scale_fill_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  facet_wrap(~ study_sub, nrow = 2)+
  geom_signif(comparisons = list(c("NR","R")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = wilcox.test, ##计算方法
              y_position = 6.5,#图中横线位置 设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.5,color="black")+
  theme_bw()+
  labs(fill="Response")+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(
        axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
        legend.position = "right",
        #legend.key.size = unit(0.5, "cm"),
        strip.text = element_text(size = 8,colour = "black"),
        legend.title = element_text(size=12,colour = "black"),
        legend.text = element_text(size=12,colour = "black"))

p1
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_cohort_response.tiff",plot = p1,width = 220,height = 180,units = "mm",dpi = 300)

##合并数据集
p2<-ggplot(res,aes(response_code,shannon,fill=response_code))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=response_code),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  scale_fill_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  geom_signif(comparisons = list(c("NR","R")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = wilcox.test, ##计算方法
              y_position = 6.5,#图中横线位置 设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.5,color="black")+
  theme_bw()+
  labs(fill="Response")+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(
        axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
        legend.position = "right",
        #legend.key.size = unit(0.5, "cm"),
        legend.title = element_text(size=12,colour = "black"),
        legend.text = element_text(size=12,colour = "black"))
p2
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_all_response.tiff",plot = p2,width = 200,height = 180,units = "mm",dpi = 300)

##不同癌种
p3<-ggplot(res,aes(response_code,shannon,fill=response_code))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=response_code),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  scale_fill_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  facet_wrap(~ cancer_type)+
  geom_signif(comparisons = list(c("NR","R")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = wilcox.test, ##计算方法
              y_position = 6.5,#图中横线位置 设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.5,color="black")+
  theme_bw()+
  labs(fill="Response")+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(
        axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
        legend.position = "right",
        #legend.key.size = unit(0.5, "cm"),
        strip.text = element_text(size = 12,colour = "black"),
        legend.title = element_text(size=12,colour = "black"),
        legend.text = element_text(size=12,colour = "black"))
p3
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_cancer_type_response.tiff",plot = p3,width = 200,height = 180,units = "mm",dpi = 300)

##癌种之间####
p4<-ggplot(fungi_species_cli,aes(cancer_type,shannon,fill=cancer_type))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=cancer_type),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P ", ..p.format..)), # 显示p值
                     size =3.8,          # 修改字体大小
                     position = position_nudge(x=1,y = -0.4))+
  theme_bw()+
  labs(fill="Cancer",x=NULL)+
  scale_fill_manual(values = c("#B2DF8A","#f2b56f","#b8aeeb"))+
  scale_y_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        #legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))

p4
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_cancer_type.tiff",plot = p4,width = 200,height = 180,units = "mm",dpi = 300)

p4_3<-p4+p3+plot_layout(guides = 'collect', widths = c(1, 2))
p4_3
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_cancer_type_sub.tiff",plot = p4_3,width = 200,height = 100,units = "mm",dpi = 300)

##age####
age<-subset(fungi_species_cli,!is.na(age_group))
p5<-ggplot(age,aes(age_group,shannon,fill=age_group))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=age_group),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)),  # 显示p值
                     size = 3.8,          # 修改字体大小
                     position = position_nudge(x=0.4,y = -0.4))+
  theme_bw()+
  labs(fill="Age",x=NULL)+
  scale_fill_manual(values = c("#B2DF8A","#f2b56f"))+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        #legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
p5
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_age.tiff",plot = p5,width = 200,height = 180,units = "mm",dpi = 300)

##gender####
gender<-subset(fungi_species_cli,!is.na(gender))
p6<-ggplot(gender,aes(gender,shannon,fill=gender))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=gender),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  stat_compare_means(method = "wilcox.test", 
                     aes(label = paste0("P = ", ..p.format..)),  # 显示p值
                     size = 3.8,          # 修改字体大小
                     position = position_nudge(x=0.4,y = -0.4))+
  theme_bw()+
  labs(fill="Gender",x=NULL)+
  scale_fill_manual(values = c("#B2DF8A","#f2b56f"))+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        #legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
p6
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_gender.tiff",plot = p6,width = 200,height = 180,units = "mm",dpi = 300)

##BMI####
BMI<-subset(fungi_species_cli,!is.na(BMI_group))
BMI$BMI_group<-gsub("Obese","obese",BMI$BMI_group)
BMI$BMI_group<-factor(BMI$BMI_group,levels = c("lean","overweight","obese"))
p7<-ggplot(BMI,aes(BMI_group,shannon,fill=BMI_group))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=BMI_group),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P = ", ..p.format..)), # 显示p值
                     size =3.8,          # 修改字体大小
                     position = position_nudge(x=1,y = -0.4))+
  theme_bw()+
  labs(fill="BMI",x=NULL)+
  scale_fill_manual(values = c("#B2DF8A","#f2b56f","#b8aeeb"))+
  scale_y_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        #legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))
p7
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_BMI.tiff",plot = p7,width = 200,height = 180,units = "mm",dpi = 300)

##cohort####
p8<-ggplot(fungi_species_cli,aes(study,shannon,fill=study))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=BMI_group),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P ", ..p.format..)), # 显示p值
                     size =3.8,          # 修改字体大小
                     position = position_nudge(x=3.5,y = -0.4))+
  theme_bw()+
  labs(fill="Cohort",x=NULL)+
  scale_fill_manual(values = c("#B2DF8A","#f2b56f","#b8aeeb","#88d8db",
                               "#fae69e","#f2a7da" ,"#bbde7a","#68A7BE"))+
  scale_y_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        #legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10))
  #guides(fill = guide_legend(nrow = 4, byrow = TRUE))  # 设置图例为三行，每行三个

p8
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_cohort.tiff",plot = p8,width = 200,height = 180,units = "mm",dpi = 300)

p8_1<-p8+p1+plot_layout(guides = 'collect', widths = c(1,3))
p8_1
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_cohort_sub.tiff",plot = p8_1,width = 300,height = 120,units = "mm",dpi = 300)

##country####
p9<-ggplot(fungi_species_cli,aes(country,shannon,fill=country))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=BMI_group),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  stat_compare_means(method = "kruskal.test", 
                     aes(label = paste0("P ", ..p.format..)), # 显示p值
                     size =3.8,          # 修改字体大小
                     position = position_nudge(x=2,y = -0.4))+
  theme_bw()+
  labs(fill="Country",x=NULL)+
  scale_fill_manual(values = c("#B2DF8A","#f2b56f","#b8aeeb","#88d8db","#fae69e","#f2a7da"))+
  scale_y_continuous(limits = c(0,8),breaks = seq(0,8,1))+
  theme(axis.text = element_text(color = 'black',size=12),
        axis.title = element_text(color = 'black',size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        #legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12))+
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))  # 设置图例为三行，每行三个

p9
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_country.tiff",plot = p9,width = 200,height = 180,units = "mm",dpi = 300)

p10<-(p8+p6+p4)/(p9+p5+p7)
p10
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_cancer_age_gender_BMI.tiff",plot = p10,width = 300,height = 250,units = "mm",dpi = 300)

##不同BMI间
BMI<-subset(BMI,!is.na(response_code))
p11<-ggplot(BMI,aes(response_code,shannon,fill=response_code))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  geom_jitter(aes(fill=response_code),shape=21,size=1.5,width=0.2,alpha=0.8)+#添加抖动点
  scale_fill_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  facet_wrap(~ BMI_group)+
  geom_signif(comparisons = list(c("NR","R")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = wilcox.test, ##计算方法
              y_position = 6.5,#图中横线位置 设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.5,color="black")+
  theme_bw()+
  labs(fill="Response")+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(
    axis.text = element_text(color = 'black',size=12),
    axis.title = element_text(color = 'black',size=12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
    legend.position = "top",
    #legend.key.size = unit(0.5, "cm"),
    strip.text = element_text(size = 12,colour = "black"),
    legend.title = element_text(size=12,colour = "black"),
    legend.text = element_text(size=12,colour = "black"))
p11
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_BMI_response.tiff",plot = p11,width = 200,height = 180,units = "mm",dpi = 300)

##不同Gender间
gender<-subset(gender,!is.na(response_code))
p12<-ggplot(gender,aes(response_code,shannon,fill=response_code))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  geom_jitter(aes(fill=response_code),shape=21,size=1.5,width=0.2,alpha=0.8)+#添加抖动点
  scale_fill_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  facet_wrap(~ gender)+
  geom_signif(comparisons = list(c("NR","R")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = wilcox.test, ##计算方法
              y_position = 6.5,#图中横线位置 设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.5,color="black")+
  theme_bw()+
  labs(fill="Response")+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(
    axis.text = element_text(color = 'black',size=12),
    axis.title = element_text(color = 'black',size=12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
    legend.position = "top",
    #legend.key.size = unit(0.5, "cm"),
    strip.text = element_text(size = 12,colour = "black"),
    legend.title = element_text(size=12,colour = "black"),
    legend.text = element_text(size=12,colour = "black"))
p12
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_gender_response.tiff",plot = p12,width = 200,height = 180,units = "mm",dpi = 300)

##不同age间
age<-subset(age,!is.na(response_code))
p13<-ggplot(age,aes(response_code,shannon,fill=response_code))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=response_code),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  scale_fill_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  facet_wrap(~ age_group)+
  geom_signif(comparisons = list(c("NR","R")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = wilcox.test, ##计算方法
              y_position = 6.5,#图中横线位置 设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.5,color="black")+
  theme_bw()+
  labs(fill="Response")+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(
    axis.text = element_text(color = 'black',size=12),
    axis.title = element_text(color = 'black',size=12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
    legend.position = "right",
    #legend.key.size = unit(0.5, "cm"),
    strip.text = element_text(size = 12,colour = "black"),
    legend.title = element_text(size=12,colour = "black"),
    legend.text = element_text(size=12,colour = "black"))
p13
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_age_response.tiff",plot = p13,width = 200,height = 180,units = "mm",dpi = 300)
p13_5<-p5+p13+plot_layout(guides = 'collect', widths = c(1, 2))
p13_5
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_age_sub.tiff",plot = p13_5,width = 200,height = 100,units = "mm",dpi = 300)

##不同gender间
gender<-subset(gender,!is.na(response_code))
p14<-ggplot(gender,aes(response_code,shannon,fill=response_code))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=response_code),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  scale_fill_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  facet_wrap(~ gender)+
  geom_signif(comparisons = list(c("NR","R")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = wilcox.test, ##计算方法
              y_position = 6.5,#图中横线位置 设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.5,color="black")+
  theme_bw()+
  labs(fill="Response")+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(
    axis.text = element_text(color = 'black',size=12),
    axis.title = element_text(color = 'black',size=12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
    legend.position = "right",
    #legend.key.size = unit(0.5, "cm"),
    strip.text = element_text(size = 12,colour = "black"),
    legend.title = element_text(size=12,colour = "black"),
    legend.text = element_text(size=12,colour = "black"))
p14
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_gender_response.tiff",plot = p14,width = 200,height = 180,units = "mm",dpi = 300)
p14_6<-p6+p14+plot_layout(guides = 'collect', widths = c(1, 2))
p14_6
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_gender_sub.tiff",plot = p14_6,width = 200,height = 100,units = "mm",dpi = 300)

##不同BMI间
BMI<-subset(BMI,!is.na(response_code))
p15<-ggplot(BMI,aes(response_code,shannon,fill=response_code))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.6,color=NA,alpha=0.8)+
  geom_boxplot(width=0.3,size=0.8,outlier.color =NA,alpha=0.8)+
  #geom_jitter(aes(fill=response_code),shape=21,size=1,width=0.2,alpha=0.8)+#添加抖动点
  scale_fill_manual(values = c("R" = "#68A7BE","NR" = "#EE7E77"))+
  facet_wrap(~ BMI_group)+
  geom_signif(comparisons = list(c("NR","R")),# 设置需要比较的组
              map_signif_level = F, #是否使用星号显示
              test = wilcox.test, ##计算方法
              y_position = 6.5,#图中横线位置 设置
              tip_length = c(0.02,0.02),#横线下方的竖线设置
              size=0.5,color="black")+
  theme_bw()+
  labs(fill="Response")+
  scale_y_continuous(limits = c(0,7),breaks = seq(0,7,1))+
  theme(
    axis.text = element_text(color = 'black',size=12),
    axis.title = element_text(color = 'black',size=12),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
    legend.position = "right",
    #legend.key.size = unit(0.5, "cm"),
    strip.text = element_text(size = 12,colour = "black"),
    legend.title = element_text(size=12,colour = "black"),
    legend.text = element_text(size=12,colour = "black"))
p15
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_BMI_response.tiff",plot = p15,width = 200,height = 180,units = "mm",dpi = 300)
p15_7<-p7+p15+plot_layout(guides = 'collect', widths = c(1, 2))
p15_7
ggsave("10.pics/supplemental figures/α-diversity/species/fungi_species_α_diversity_BMI_sub.tiff",plot = p15_7,width = 200,height = 100,units = "mm",dpi = 300)
