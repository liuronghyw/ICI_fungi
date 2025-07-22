setwd("D:/课题/细菌—真菌互作/3_R/")
rm(list=ls())

library(tidyverse)
library(ggplot2)
library(patchwork)
library(plyr)##ddply函数
library(lemon)##拼图共享一个图例

cli_whole<-read.csv("2.rawData/clinical/final/cli_whole_samples_new.csv")
colnames(cli_whole)[1]<-"id"

##过滤标准1：样本最大相对丰度分布####
##fungi####
fungi_19<-read.csv("3.mirco_data/PRJEB43119/Fungi.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(fungi_19)<-gsub(".S.bracken.kreport2","",colnames(fungi_19))

fungi_06<-read.csv("3.mirco_data/PRJNA397906/Fungi.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(fungi_06)<-gsub(".S.bracken.kreport2","",colnames(fungi_06))

fungi_42<-read.csv("3.mirco_data/PRJNA399742/Fungi.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(fungi_42)<-gsub(".S.bracken.kreport2","",colnames(fungi_42))

fungi_81<-read.csv("3.mirco_data/PRJNA541981/Fungi.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(fungi_81)<-gsub(".S.bracken.kreport2","",colnames(fungi_81))

fungi_60<-read.csv("3.mirco_data/PRJNA762360/Fungi.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(fungi_60)<-gsub(".S.bracken.kreport2","",colnames(fungi_60))

fungi_95<-read.csv("3.mirco_data/PRJNA770295/Fungi.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(fungi_95)<-gsub(".S.bracken.kreport2","",colnames(fungi_95))

fungi_97<-read.csv("3.mirco_data/PRJNA1023797/Fungi.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(fungi_97)<-gsub(".S.bracken.kreport2","",colnames(fungi_97))

fungi_63<-read.csv("3.mirco_data/PRJEB22863/Fungi.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(fungi_63)<-gsub(".S.bracken.kreport2","",colnames(fungi_63))

comb_fungi<-merge(fungi_19,fungi_95,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_06,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_42,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_81,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_60,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_97,by="#Classification",all = T)
comb_fungi<-merge(comb_fungi,fungi_63,by="#Classification",all = T)

rownames(comb_fungi)<-comb_fungi$`#Classification`
comb_fungi<-comb_fungi[,-1]
comb_fungi<-comb_fungi[,cli_whole$id]

comb_fungi[is.na(comb_fungi)]<-0

comb_fungi<-as.data.frame(t(comb_fungi))
comb_fungi<-subset(comb_fungi,rowSums(comb_fungi)!=0)##去除未比对到真菌样本
write.csv(comb_fungi,"3.mirco_data/fungi/fungi_all_cohorts_raw.csv")
comb_fungi_rela <- sweep(comb_fungi, 1, rowSums(comb_fungi), FUN = "/")##计算相对丰度

comb_fungi_rela$max<-apply(comb_fungi_rela,1,max)
comb_fungi_rela$max<-as.numeric(format(comb_fungi_rela$max,digits=2))
median(comb_fungi_rela$max)

##作图
p1<-ggplot(comb_fungi_rela,aes(x=max,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,1,0.1),binwidth =0.1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3.5,,breaks=seq(0,1,0.1))+
  labs(x="Maximum relative abundance of fungi in per sample",y="Count of sample")+
  geom_vline(aes(xintercept=median(max)), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 0.13, y = 520, label = "Median = 0.13",size=4)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,1,0.1))+
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size =12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p1


##bacteria####
bac_19<-read.csv("3.mirco_data/PRJEB43119/Bacteria.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(bac_19)<-gsub(".S.bracken.kreport2","",colnames(bac_19))

bac_06<-read.csv("3.mirco_data/PRJNA397906/Bacteria.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(bac_06)<-gsub(".S.bracken.kreport2","",colnames(bac_06))

bac_42<-read.csv("3.mirco_data/PRJNA399742/Bacteria.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(bac_42)<-gsub(".S.bracken.kreport2","",colnames(bac_42))

bac_81<-read.csv("3.mirco_data/PRJNA541981/Bacteria.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(bac_81)<-gsub(".S.bracken.kreport2","",colnames(bac_81))

bac_60<-read.csv("3.mirco_data/PRJNA762360/Bacteria.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(bac_60)<-gsub(".S.bracken.kreport2","",colnames(bac_60))

bac_95<-read.csv("3.mirco_data/PRJNA770295/Bacteria.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(bac_95)<-gsub(".S.bracken.kreport2","",colnames(bac_95))

bac_63<-read.csv("3.mirco_data/PRJEB22863/Bacteria.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(bac_63)<-gsub(".S.bracken.kreport2","",colnames(bac_63))

bac_97<-read.csv("3.mirco_data/PRJNA1023797/Bacteria.S.bracken_mpa.kreport2_species.txt",sep="\t",check.names = F)
colnames(bac_97)<-gsub(".S.bracken.kreport2","",colnames(bac_97))

comb_bac<-merge(bac_19,bac_95,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_06,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_42,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_81,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_60,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_97,by="#Classification",all = T)
comb_bac<-merge(comb_bac,bac_63,by="#Classification",all = T)

rownames(comb_bac)<-comb_bac$`#Classification`
comb_bac<-comb_bac[,-1]
comb_bac<-comb_bac[,rownames(comb_fungi)]

comb_bac[is.na(comb_bac)]<-0

comb_bac<-as.data.frame(t(comb_bac))

write.csv(comb_bac,"3.mirco_data/bacteria/bac_all_cohorts_raw.csv")
comb_bac_rela <- sweep(comb_bac, 1, rowSums(comb_bac), FUN = "/")##计算相对丰度

comb_bac_rela$max<-apply(comb_bac_rela,1,max)
comb_bac_rela$max<-as.numeric(format(comb_bac_rela$max,digits=2))
median(comb_bac_rela$max)
##作图
p2<-ggplot(comb_bac_rela,aes(x=max,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,1,0.1),binwidth =0.1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3.5,,breaks=seq(0,1,0.1))+
  labs(x="Maximum relative abundance of bacteria in per sample",y="Count of sample")+
  geom_vline(aes(xintercept=median(max)), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 0.22, y = 500, label = "Median = 0.22",size=4)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,1,0.1))+
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p2

##过滤标准2:比对到微生物数量####
##fungi####
comb_fungi$count<-rowSums(comb_fungi>0)
median(comb_fungi$count)
nfungi<-ncol(comb_fungi)

for(i in 1:nrow(comb_fungi)){
  if(comb_fungi[i,nfungi]<=10){comb_fungi$group[i]<-"0"}
  else if(10<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=20){comb_fungi$group[i]<-"1.1"}
  else if(20<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=30){comb_fungi$group[i]<-"2.1"}
  else if(30<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=40){comb_fungi$group[i]<-"3.1"}
  else if(40<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=50){comb_fungi$group[i]<-"4.1"}
  else if(50<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=100){comb_fungi$group[i]<-"5.1"}
  else if(100<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=200){comb_fungi$group[i]<-"6.1"}
  else if(200<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=300){comb_fungi$group[i]<-"7.1"}
  else if(300<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=400){comb_fungi$group[i]<-"8.1"}
  else if(400<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=500){comb_fungi$group[i]<-"9.1"}
  else if(500<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=600){comb_fungi$group[i]<-"10.1"}
  else if(600<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=700){comb_fungi$group[i]<-"11.1"}
  else if(700<comb_fungi[i,nfungi]&comb_fungi[i,nfungi]<=800){comb_fungi$group[i]<-"12.1"}
  else if(comb_fungi[i,nfungi]>800){comb_fungi$group[i]<-"13.1"}
}
comb_fungi$group<-as.numeric(comb_fungi$group)
median(comb_fungi$group)

p3<-ggplot(comb_fungi,aes(x=group,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,14,1),binwidth =1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3.5,breaks=seq(0,14,1))+
  labs(x="The number of fungi mapped in per sample",y="Count of sample")+
  geom_vline(aes(xintercept=6.2), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 6.1, y = 550, label = "Median = 124",size=4)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,14,1),labels = c("0","10","20","30","40","50","100","200","300","400","500","600","700","800","900"))+
  scale_y_continuous(breaks = seq(0,500,100))+
  theme(axis.text = element_text(size =12), 
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p3

comb_fungi$ratio<-comb_fungi$count/(nfungi-1)##所有数据集交集共比对918个真菌
median(comb_fungi$ratio)
p4<-ggplot(comb_fungi,aes(x=ratio,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,1,0.05),binwidth =0.05,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3.5,breaks=seq(0,1,0.05))+
  labs(x="The ratio of fungi mapped in per sample",y="Count of sample")+
  geom_vline(aes(xintercept=0.133), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 0.13, y = 360, label = "Median = 0.13",size=4)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,1,0.05))+
  theme(axis.text.x = element_text(size =12,angle = 45,hjust=1), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p4

##bacteria####
comb_bac$count<-rowSums(comb_bac>0)
median(comb_bac$count)
nbac<-ncol(comb_bac)

for(i in 1:nrow(comb_bac)){
  if(comb_bac[i,nbac]<=200){comb_bac$group[i]<-"0"}
  else if(200<comb_bac[i,nbac]&comb_bac[i,nbac]<=400){comb_bac$group[i]<-"1.1"}##.1的目的是避免包括界限数字，不然画直方图时0与1组会一个柱形图
  else if(400<comb_bac[i,nbac]&comb_bac[i,nbac]<=600){comb_bac$group[i]<-"2.1"}
  else if(600<comb_bac[i,nbac]&comb_bac[i,nbac]<=800){comb_bac$group[i]<-"3.1"}
  else if(800<comb_bac[i,nbac]&comb_bac[i,nbac]<=1000){comb_bac$group[i]<-"4.1"}
  else if(1000<comb_bac[i,nbac]&comb_bac[i,nbac]<=2000){comb_bac$group[i]<-"5.1"}
  else if(2000<comb_bac[i,nbac]&comb_bac[i,nbac]<=3000){comb_bac$group[i]<-"6.1"}
  else if(3000<comb_bac[i,nbac]&comb_bac[i,nbac]<=4000){comb_bac$group[i]<-"7.1"}
  else if(4000<comb_bac[i,nbac]&comb_bac[i,nbac]<=5000){comb_bac$group[i]<-"8.1"}
  else if(5000<comb_bac[i,nbac]&comb_bac[i,nbac]<=6000){comb_bac$group[i]<-"9.1"}
  else if(6000<comb_bac[i,nbac]&comb_bac[i,nbac]<=7000){comb_bac$group[i]<-"10.1"}
  else if(7000<comb_bac[i,nbac]&comb_bac[i,nbac]<=8000){comb_bac$group[i]<-"11.1"}
}

comb_bac$group<-as.numeric(comb_bac$group)
median(comb_bac$group)

p5<-ggplot(comb_bac,aes(x=group,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,12,1),binwidth =1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3.5,breaks=seq(0,12,1))+
  labs(x="The number of bacteria mapped in per sample",y="Count of sample")+
  geom_vline(aes(xintercept=5.43), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 5.4, y = 630, label = "Median = 1415",size=4)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,12,1),labels = c("0","200","400","600","800","1000","2000","3000","4000","5000","6000","7000","8000"))+
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size =12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p5

comb_bac$ratio<-comb_bac$count/(nbac-1)##所有数据集交集共比对8789个细菌
median(comb_bac$ratio)
p6<-ggplot(comb_bac,aes(x=ratio,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,1,0.05),binwidth =0.05,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3.5,breaks=seq(0,1,0.05))+
  labs(x="The ratio of bacteria mapped in per sample",y="Count of sample")+
  geom_vline(aes(xintercept=0.16), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 0.16, y = 300, label = "Median = 0.16",size=4)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,1,0.05))+
  theme(axis.text.x = element_text(size =12,angle = 45,hjust=1), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size =12),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p6

p7<-(p1+p2)/(p3+p5)
p7
ggsave("10.pics/supplemental figures/filter_samples/Figure1_filter_samples_standards.tiff",plot = p7,width = 240,height = 180,units = "mm",dpi=300)

p8<-(p1+p2)/(p4+p6)
p8
ggsave("10.pics/supplemental figures/filter_samples/Figure1_filter_samples_standards_ratio.tiff",plot = p8,width = 260,height = 180,units = "mm",dpi=300)

ggsave("10.pics/supplemental figures/filter_samples/Figure1_filter_samples_standards_fungi_rela.tiff",plot = p1,width = 150,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/filter_samples/Figure1_filter_samples_standards_bac_rela.tiff",plot = p2,width = 150,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/filter_samples/Figure1_filter_samples_standards_fungi_count.tiff",plot = p3,width = 150,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/filter_samples/Figure1_filter_samples_standards_bac_count.tiff",plot = p5,width = 150,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/filter_samples/Figure1_filter_samples_standards_fungi_ratio.tiff",plot = p4,width = 150,height = 100,units = "mm",dpi=300)
ggsave("10.pics/supplemental figures/filter_samples/Figure1_filter_samples_standards_bac_ratio.tiff",plot = p6,width = 150,height = 100,units = "mm",dpi=300)


##过滤标准3：真菌reads数在0.005%-1%#####
cli_ratio<-read.csv("2.rawData/clinical/final/cli_whole_samples_ratio.csv",check.names = F,row.names = 1)
cli_ratio<-cli_ratio[rownames(comb_fungi),]
median(cli_ratio$fungi_ratio)

for(i in 1:nrow(cli_ratio)){
  if(cli_ratio$fungi_ratio[i]<=0.00001){cli_ratio$group[i]<-"0"}
  else if(0.00001<cli_ratio$fungi_ratio[i]&cli_ratio$fungi_ratio[i]<=0.00005){cli_ratio$group[i]<-"1.1"}##.1的目的是避免包括界限数字，不然画直方图时0与1组会一个柱形图
  else if(0.00005<cli_ratio$fungi_ratio[i]&cli_ratio$fungi_ratio[i]<=0.0001){cli_ratio$group[i]<-"2.1"}
  else if(0.0001<cli_ratio$fungi_ratio[i]&cli_ratio$fungi_ratio[i]<=0.0005){cli_ratio$group[i]<-"3.1"}
  else if(0.0005<cli_ratio$fungi_ratio[i]&cli_ratio$fungi_ratio[i]<=0.001){cli_ratio$group[i]<-"4.1"}
  else if(0.001<cli_ratio$fungi_ratio[i]&cli_ratio$fungi_ratio[i]<=0.005){cli_ratio$group[i]<-"5.1"}
  else if(0.005<cli_ratio$fungi_ratio[i]&cli_ratio$fungi_ratio[i]<=0.01){cli_ratio$group[i]<-"6.1"}
  else if(cli_ratio$fungi_ratio[i]>0.01){cli_ratio$group[i]<-"7.1"}
}

cli_ratio$group<-as.numeric(cli_ratio$group)
median(cli_ratio$group)

p9<-ggplot(cli_ratio,aes(x=group,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,8,1),binwidth =1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,breaks=seq(0,8,1))+
  labs(x="The reads ratio of fungi mapped",y="Count of sample")+
  geom_vline(aes(xintercept=3.45), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 3.5, y = 850, label = "Median = 0.024%",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,8,1),labels = c("0","0.001%","0.005%","0.01%","0.05%","0.1%","0.5%","1%","100%"))+
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size =10),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p9


median(cli_ratio$log_fungi_all)

p9_all<-ggplot(cli_ratio,aes(x=log_fungi_all,y=after_stat(count))) +
  geom_histogram(breaks=seq(-5.6,-0.4,0.3),binwidth =0.1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,breaks=seq(-5.6,-0.4,0.3))+
  labs(x=expression(log[10](Fungi/Microbe)),y="Count of sample")+
  geom_vline(aes(xintercept=-3.61), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = -3.61, y = 450, label = "Median = -3.61",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(-5.6,-0.4,0.3))+
  theme(axis.text.x = element_text(size = 10,angle = 45,hjust=1), 
        axis.text.y = element_text(size =10),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p9_all

ggsave("10.pics/Figure1/figure1_filter_samples_fungi_reads_ratio.tiff",plot = p9,width = 200,height = 100,units = "mm",dpi=300)
ggsave("10.pics/Figure1/figure1_filter_samples_fungi_reads_ratio_sub.tiff",plot = p9_all,width = 200,height = 100,units = "mm",dpi=300)

table(cli_ratio$cancer_type) 
cli_mela<-subset(cli_ratio,cli_ratio$cancer_type=="Melanoma")
median(cli_mela$log_fungi_all)
p9_mela<-ggplot(cli_mela,aes(x=log_fungi_all,y=after_stat(count))) +
  geom_histogram(breaks=seq(-5.6,-0.4,0.3),binwidth =0.1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,breaks=seq(-5.6,-0.4,0.3))+
  labs(x=expression(log[10](Fungi/Microbe)),y="Count of sample")+
  geom_vline(aes(xintercept=-3.40), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = -3.40, y = 150, label = "Median = -3.39",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(-5.6,-0.4,0.3))+
  theme(axis.text.x = element_text(size = 10,angle = 45,hjust=1), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p9_mela
ggsave("10.pics/Figure1/figure1_filter_samples_fungi_reads_ratio_mela.tiff",plot = p9_mela,width = 200,height = 100,units = "mm",dpi=300)

cli_NSCLC<-subset(cli_ratio,cli_ratio$cancer_type=="NSCLC")
median(cli_NSCLC$log_fungi_all)
p9_NSCLC<-ggplot(cli_NSCLC,aes(x=log_fungi_all,y=after_stat(count))) +
  geom_histogram(breaks=seq(-5.6,-0.4,0.3),binwidth =0.1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,breaks=seq(-5.6,-0.4,0.3))+
  labs(x=expression(log[10](Fungi/Microbe)),y="Count")+
  geom_vline(aes(xintercept=-3.7), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = -3.7, y = 300, label = "Median = -3.70",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(-5.6,-0.4,0.3))+
  theme(axis.text.x = element_text(size =10,angle = 45,hjust=1), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size =10),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p9_NSCLC
ggsave("10.pics/Figure1/figure1_filter_samples_fungi_reads_ratio_NSCLC.tiff",plot = p9_NSCLC,width = 200,height = 100,units = "mm",dpi=300)

cli_RCC<-subset(cli_ratio,cli_ratio$cancer_type=="RCC")
median(cli_RCC$log_fungi_all)
p9_RCC<-ggplot(cli_RCC,aes(x=log_fungi_all,y=after_stat(count))) +
  geom_histogram(breaks=seq(-5.6,-0.4,0.3),binwidth =0.1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,breaks=seq(-5.6,-0.4,0.3))+
  labs(x=expression(log[10](Fungi/Microbe)),y="Count of sample")+
  geom_vline(aes(xintercept=-3.74), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = -3.74, y = 50, label = "Median = -3.74",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(-5.6,-0.4,0.3))+
  theme(axis.text.x = element_text(size = 10,angle = 45,hjust=1), 
        axis.text.y = element_text(size =10),
        axis.title = element_text(size = 10),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p9_RCC
ggsave("10.pics/Figure1/figure1_filter_samples_fungi_reads_ratio_RCC.tiff",plot = p9_RCC,width = 200,height = 100,units = "mm",dpi=300) 

###每种菌的平均相对丰度分布####
##fungi
fungi<-as.data.frame(t(comb_fungi_rela[,1:(ncol(comb_fungi_rela)-1)]))##有列max
fungi$mean<-apply(fungi,1,mean)
quantile(fungi$mean,probs = 0.5)

for(i in 1:nrow(fungi)){
  if(fungi$mean[i]<=0.000001){fungi$group[i]<-"0"}
  else if(0.000001<fungi$mean[i]&fungi$mean[i]<=0.000005){fungi$group[i]<-"1.1"}##.1的目的是避免包括界限数字，不然画直方图时0与1组会一个柱形图
  else if(0.000005<fungi$mean[i]&fungi$mean[i]<=0.00001){fungi$group[i]<-"2.1"}
  else if(0.00001<fungi$mean[i]&fungi$mean[i]<=0.00005){fungi$group[i]<-"3.1"}
  else if(0.00005<fungi$mean[i]&fungi$mean[i]<=0.0001){fungi$group[i]<-"4.1"}
  else if(0.0001<fungi$mean[i]&fungi$mean[i]<=0.0005){fungi$group[i]<-"5.1"}
  else if(0.0005<fungi$mean[i]&fungi$mean[i]<=0.001){fungi$group[i]<-"6.1"}
  else if(0.001<fungi$mean[i]&fungi$mean[i]<=0.005){fungi$group[i]<-"7.1"}
  else if(0.005<fungi$mean[i]&fungi$mean[i]<=0.01){fungi$group[i]<-"8.1"}
  else if(fungi$mean[i]>0.01){fungi$group[i]<-"9.1"}
}

fungi$group<-as.numeric(fungi$group)
median(fungi$group)

p10<-ggplot(fungi,aes(x=group,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,10,1),binwidth =1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,breaks=seq(0,10,1))+
  labs(x="The average relative abundance of fungi",y="Count of fungi")+
  geom_vline(aes(xintercept=6.1), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 6.1, y = 420, label = "Median = 0.054%",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,10,1),labels = c("0","0.0001%","0.0005%","0.001%","0.005%","0.01%","0.05%","0.1%","0.5%","1%","10%"))+
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size =10),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p10

##bacteria
bac<-as.data.frame(t(comb_bac_rela[,1:(ncol(comb_bac_rela)-1)]))##有列max
bac$mean<-apply(bac,1,mean)
quantile(bac$mean,probs = 0.95)

for(i in 1:nrow(bac)){
  if(bac$mean[i]<=0.0000001){bac$group[i]<-"0"}
  else if(0.0000001<bac$mean[i]&bac$mean[i]<=0.0000005){bac$group[i]<-"1.1"}##.1的目的是避免包括界限数字，不然画直方图时0与1组会一个柱形图
  else if(0.0000005<bac$mean[i]&bac$mean[i]<=0.000001){bac$group[i]<-"2.1"}
  else if(0.000001<bac$mean[i]&bac$mean[i]<=0.000005){bac$group[i]<-"3.1"}
  else if(0.000005<bac$mean[i]&bac$mean[i]<=0.00001){bac$group[i]<-"4.1"}
  else if(0.00001<bac$mean[i]&bac$mean[i]<=0.00005){bac$group[i]<-"5.1"}
  else if(0.00005<bac$mean[i]&bac$mean[i]<=0.0001){bac$group[i]<-"6.1"}
  else if(0.0001<bac$mean[i]&bac$mean[i]<=0.0005){bac$group[i]<-"7.1"}
  else if(0.0005<bac$mean[i]&bac$mean[i]<=0.001){bac$group[i]<-"8.1"}
  else if(0.001<bac$mean[i]&bac$mean[i]<=0.005){bac$group[i]<-"9.1"}
  else if(0.005<bac$mean[i]&bac$mean[i]<=0.01){bac$group[i]<-"10.1"}
  else if(bac$mean[i]>0.01){bac$group[i]<-"11.1"}
}

bac$group<-as.numeric(bac$group)
quantile(bac$group,probs = 0.95)

p11<-ggplot(bac,aes(x=group,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,12,1),binwidth =1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,breaks=seq(0,12,1))+
  labs(x="The average relative abundance of bacteria",y="Count of bacteria")+
  geom_vline(aes(xintercept=5.1), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x =5.1, y = 4200, label = "95%: 0.0011%",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,12,1),labels = c("0","0.00001%","0.00005%","0.0001%","0.0005%","0.001%","0.005%","0.01%","0.05%","0.1%","0.5%","1%","15%"))+
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size =10),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p11

p12<-p10+p11
p12
ggsave("10.pics/Figure1/figure1_selection_condition.tiff",plot = p12,width = 250,height = 100,units = "mm",dpi=300) 


##每种菌在样本中检测的比例####
###melanoma
mela<-read.csv("5.fungi_data/Melanoma/relabundance/species/fungi_mela_whole_rela_filter.csv",check.names = F,row.names = 1)
mela$ratio<-rowSums(mela>0)/ncol(mela)
median(mela$ratio)

p13<-ggplot(mela,aes(x=ratio,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,1,0.05),binwidth =0.1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,,breaks=seq(0,1,0.05))+
  labs(x="Ratio of fungi mapped in samples",y="Count of fungi",title = "Melanoma")+
  geom_vline(aes(xintercept=median(ratio)), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 0.49, y = 200, label = "Median = 0.496",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,1,0.05))+
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size =10),
        plot.title = element_text(size = 10,hjust=0.5),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p13

###NSCLC
NSCLC<-read.csv("5.fungi_data/NSCLC/relabundance/species/fungi_NSCLC_whole_rela_filter.csv",check.names = F,row.names = 1)
NSCLC$ratio<-rowSums(NSCLC>0)/ncol(NSCLC)
median(NSCLC$ratio)

p14<-ggplot(NSCLC,aes(x=ratio,y=after_stat(count))) +
  geom_histogram(breaks= seq(0,1,0.05),binwidth =0.1,color="#508AB2",fill="grey")+  ##binwidth代表间隔宽度
  theme_classic()+
  geom_text(aes(label=gsub("^0$","",as.character(after_stat(count)))),stat = "bin",vjust=-0.5,size=3,,breaks=seq(0,1,0.05))+
  labs(x="Ratio of fungi mapped in samples",y="Count of fungi",title = "NSCLC")+
  geom_vline(aes(xintercept=median(ratio)), color="#EB4A59", linetype="dashed", linewidth=0.5)+
  annotate("text", x = 0.10, y = 300, label = "Median = 0.11",size=3)+        ##添加文本注释
  scale_x_continuous(breaks = seq(0,1,0.05))+
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size =10),
        axis.text.x = element_text(angle = 45,hjust = 1),
        plot.title = element_text(size = 10,hjust=0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"))
p14

p15<-p13+p14
p15
ggsave("10.pics/Figure1/figure1_selection_condition_fungi_mapped_ratio.tiff",plot = p15,width = 250,height = 100,units = "mm",dpi=300) 

