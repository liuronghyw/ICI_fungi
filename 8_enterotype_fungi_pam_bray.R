rm(list=ls())
setwd("D:/课题/细菌—真菌互作/3_R/")

library(plyr)
library(reshape2)
library(magrittr)
library(vegan)
library(ade4)
library(cluster)
library(sparcl)
library(factoextra)
library(ggplot2)
library(patchwork)
library(tidyr)

############
###1.数据预处理
##读入真菌——属丰度信息
rela_fungi_adj<-read.csv("5.fungi_data/fungi_whole_rela_genus_filter_MMUPHIN.csv",header=T,check.names = F,row.names = 1)

##读入临床文件用于匹配
clin<-read.csv("2.rawData/clinical/final/filter/cli_whole_filter.csv",check.names = F,row.names = 1)

rela_fungi_adj <- subset(rela_fungi_adj,(rowSums(rela_fungi_adj>=0.001)/ncol(rela_fungi_adj))>=0.4)#至少40%的样本相对丰度为0.1%
#rela_fungi_adj <- subset(rela_fungi_adj,rowMeans(rela_fungi_adj)>=0.001&(rowSums(rela_fungi_adj>0)/ncol(rela_fungi_adj))>=0.5)
rela_fungi_adj<-t(rela_fungi_adj)

rela_fungi_scale<-scale(rela_fungi_adj,center = F, scale = apply(rela_fungi_adj, 2, sd, na.rm = TRUE))

pam.clustering=function(x,k) {
  # x is a distance matrix and k the number of clusters\
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

##选择最佳聚类数
tiff("10.pics/supplemental figures/Enterotype/fungi/fungi_clusters_silhouette.tiff",width = 100,height =100,units = "mm",res = 300)
fviz_nbclust(rela_fungi_adj,pam,method = "silhouette")+# 总平方和方法silhouette (Within Sum of Squares)
  ggtitle("Optimal number of fungal enterotypes") +
  theme(plot.title = element_text(hjust = 0.5))  # 居中标题
dev.off()
##bray距离，采用pam聚类
data.dist = vegdist(rela_fungi_scale,method='bray')
data.cluster = pam.clustering(data.dist,k=3)

cluster_result = as.data.frame(data.cluster)
cluster_result$id <- rownames(rela_fungi_scale)
colnames(cluster_result)<-c("cluster_fungi","id")
rownames(cluster_result)<-cluster_result$id
clin<-cbind(clin,cluster_result)

## PCoA plot of fungal enterotype clustering results

pcoa<- cmdscale(data.dist,eig=TRUE)
#提取前两个分类解释
plot_data<-data.frame({pcoa$point})[1:2]
head(plot_data)
#前两个分类解释命名
names(plot_data)[1:2]<-c('PCoA1','PCoA2') 
eig=pcoa$eig
plot_data$id<-rownames(plot_data)
samples<-merge(clin,plot_data,by="id")

#samples$cluster_fungi<-paste0("Cluster ",samples$cluster_fungi)
samples$cluster_fungi<-gsub("1","Asper_high",samples$cluster_fungi)
samples$cluster_fungi<-gsub("2","Sacc_type",samples$cluster_fungi)
samples$cluster_fungi<-gsub("3","Asper_low",samples$cluster_fungi)
#samples$cluster_fungi<-gsub("4","Sacc_type",samples$cluster_fungi)
samples$cluster_fungi<-factor(samples$cluster_fungi,levels = c("Asper_high","Sacc_type","Asper_low"))
table(samples$cluster_fungi)

#作图
rela_fungi_scale<-rela_fungi_scale[samples$id,]

set.seed(123)
dune.div<-adonis2(rela_fungi_scale~cluster_fungi,data=samples,permutations=999,method="bray")
dune_adonis<-paste0("P = ",round(dune.div$`Pr(>F)`,3))

p1<-ggplot(samples,aes(x=PCoA1,y=PCoA2,shape=cluster_fungi,color=cluster_fungi))+
  geom_point(size=1)+
  stat_ellipse(level=0.9,size=0.5)+
  theme_bw()+
  labs(subtitle = dune_adonis,title = "PCoA of fungal enterotype")+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits = 3),"%)",sep=""),
       y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits = 3),"%)",sep=""))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  scale_color_manual(values = c("Sacc_type"="#84c3b7","Asper_low"="#b8aeeb","Asper_high"="#f2b56f"))+#
  theme(axis.title.x=element_text(colour = 'black',size=12),
        axis.title.y = element_text(colour = 'black',size=12),
        legend.text = element_text(face="italic",size = 12,colour = "black"),
        axis.text = element_text(size = 12,colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(size=12,hjust = 0.5,colour = "black"),
        plot.subtitle = element_text(size=12,hjust = 0.5,face="italic",colour = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm"),
        legend.position ="bottom")

p1
ggsave("10.pics/Figure5/fungi_3_enterotype_clusters.tiff",width = 120,height =120,units = "mm",dpi = 300,plot = p1)

write.csv(samples,"9.enterotype/fungi_cluster_pcoa_clin.csv",row.names = F)

##不同cluster优势属
rela_fungi_adj<-as.data.frame(rela_fungi_adj)
var_fungi= colnames(rela_fungi_adj)
rela_fungi_adj$id<-rownames(rela_fungi_adj)
X <- merge(rela_fungi_adj,clin,by="id") # 您的实际数据
rownames(X)<-X$id
X <-X[,c(2:ncol(rela_fungi_adj),ncol(X))]
#计算每个聚类中的属的丰度

cluster_attribute_votes <- matrix(0, nrow=max(X$cluster_fungi), ncol=length(var_fungi))
for (k in 1:max(X$cluster_fungi)) {
  for (j in 1:length(var_fungi)) {
    cluster_attribute_votes[k, j] <- sum(X[X$cluster_fungi == k, j])
  }
}
colnames(cluster_attribute_votes) <- colnames(X[,1:length(var_fungi)]) 
#标准化每个聚类的属的丰度

cluster_size <- tabulate(X$cluster_fungi)
cluster_attribute_abundance <- cluster_attribute_votes / cluster_size

# 找出cluster1的最大值
max_value_row_1 <- max(cluster_attribute_abundance[1, ])
max_col_row_1 <- colnames(cluster_attribute_abundance)[which.max(cluster_attribute_abundance[1, ])]

# 找出cluster2的最大值
max_value_row_2 <- max(cluster_attribute_abundance[2, ])
max_col_row_2 <- colnames(cluster_attribute_abundance)[which.max(cluster_attribute_abundance[2, ])]

# 找出cluster3的最大值
max_value_row_3 <- max(cluster_attribute_abundance[3, ])
max_col_row_3 <- colnames(cluster_attribute_abundance)[which.max(cluster_attribute_abundance[3, ])]



# 找出平均相对丰度top10的属
abundance<-as.data.frame(t(rela_fungi_adj[,-ncol(rela_fungi_adj)]))##去除id列
abundance$mean<-rowMeans(abundance)

df_main_genus<-as.data.frame(t(abundance[,-ncol(abundance)]))
df_main_genus_cli<-cbind(df_main_genus,clin[,c("study","cluster_fungi","cancer_type")])
#df_main_genus_cli$cluster_fungi<-paste0("Cluster ",df_main_genus_cli$cluster_fungi)
df_main_genus_cli$cluster_fungi<-gsub("1","Asper_high",df_main_genus_cli$cluster_fungi)
df_main_genus_cli$cluster_fungi<-gsub("2","Sacc_type",df_main_genus_cli$cluster_fungi)
df_main_genus_cli$cluster_fungi<-gsub("3","Asper_low",df_main_genus_cli$cluster_fungi)
df_main_genus_cli$cluster_fungi<-factor(df_main_genus_cli$cluster_fungi,levels = c("Asper_high","Sacc_type","Asper_low"))

##不同cluster的genus相对丰度比较
library(ggsignif)

p2<-ggplot(df_main_genus_cli,aes(cluster_fungi,Aspergillus,fill=cluster_fungi))+
  geom_boxplot(outlier.shape = NA,width=0.7)+
  geom_jitter(aes(fill=cluster_fungi),shape=21,position = position_jitter(0.2),
              size=2,alpha=0.5)+
  geom_signif(comparisons = list(c("Asper_high","Sacc_type"),
                                 c("Asper_low","Sacc_type"),
                                 c("Asper_high","Asper_low")),
              map_signif_level = T, 
              test = "wilcox.test",
              textsize = 4,
              y_position = c(0.76,0.80,0.88),
              tip_length = c(0,0,0),
              size=0.5,color="black")+
  theme_bw()+
  theme(axis.title = element_text(color="black",size = 12),
        #panel.grid = element_blank(),
        axis.text = element_text(color="black",size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1,face="italic"),
        legend.position = "none",
        plot.title = element_text(size = 12,hjust = 0.5,face="italic"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm") )+
  labs(y="Relative abundance",title = "Aspergillus",x=NULL)+
  scale_y_continuous(limits = c(-0.01,1.01),
                     breaks = seq(0,1,0.1),
                     labels =  seq(0,1,0.1))+
  scale_fill_manual(values =c("Sacc_type"="#84c3b7","Asper_low"="#b8aeeb","Asper_high"="#f2b56f"))
p2


p3<-ggplot(df_main_genus_cli,aes(cluster_fungi,Saccharomyces,fill=cluster_fungi))+
  geom_boxplot(outlier.shape = NA,width=0.7)+
  geom_jitter(aes(fill=cluster_fungi),shape=21,position = position_jitter(0.2),
              size=2,alpha=0.5)+
  geom_signif(comparisons = list(c("Asper_high","Sacc_type"),
                                 c("Asper_low","Sacc_type"),
                                 c("Asper_high","Asper_low")),
              map_signif_level = T, 
              test = "wilcox.test",
              textsize = 4,
              y_position = c(0.76,0.80,0.88),
              tip_length = c(0,0,0),
              size=0.5,color="black")+
  theme_bw()+
  theme(axis.title = element_text(color="black",size = 12),
        #panel.grid = element_blank(),
        axis.text = element_text(color="black",size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1,face="italic"),
        legend.position = "none",
        plot.title = element_text(size = 12,hjust = 0.5,face="italic"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm") )+
  labs(y="Relative abundance",title = "Saccharomyces",x=NULL)+
  scale_y_continuous(limits = c(-0.01,1.01),
                     breaks = seq(0,1,0.1),
                     labels =  seq(0,1,0.1))+
  scale_fill_manual(values = c("Sacc_type"="#84c3b7","Asper_low"="#b8aeeb","Asper_high"="#f2b56f"))
p3

p4<-ggplot(df_main_genus_cli,aes(cluster_fungi,Penicillium,fill=cluster_fungi))+
  geom_boxplot(outlier.shape = NA,width=0.7)+
  geom_jitter(aes(fill=cluster_fungi),shape=21,position = position_jitter(0.2),
              size=2,alpha=0.5)+
  geom_signif(comparisons = list(c("Asper_high","Sacc_type"),
                                 c("Asper_low","Sacc_type"),
                                 c("Asper_high","Asper_low")),
              map_signif_level = T, 
              test = "wilcox.test",
              textsize = 3.5,
              y_position = c(0.76,0.80,0.88),
              tip_length = c(0,0,0),
              size=0.5,color="black")+
  theme_bw()+
  theme(axis.title = element_text(color="black",size = 12),
        #panel.grid = element_blank(),
        axis.text = element_text(color="black",size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1,face="italic"),
        legend.position = "none",
        plot.title = element_text(size = 12,hjust = 0.5,face="italic"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm") )+
  labs(y="Relative abundance",title = "Penicillium",x=NULL)+
  scale_y_continuous(limits = c(-0.01,1.01),
                     breaks = seq(0,1,0.1),
                     labels =  seq(0,1,0.1))+
  scale_fill_manual(values = c("Sacc_type"="#84c3b7","Asper_low"="#b8aeeb","Asper_high"="#f2b56f"))
p4

p5<-p2+p4+p3
p5
ggsave("10.pics/Figure5/fungi_3_enterotype_max_genus_rela.tiff",width = 250,height =120,units = "mm",dpi = 300,plot = p5)

p5_sub<-p2+p3
p5_sub
ggsave("10.pics/Figure5/fungi_3_enterotype_max_genus_rela_df.tiff",width = 150,height =100,units = "mm",dpi = 300,plot = p5_sub)

##不同肿瘤中丰度差异
df_main_genus_cli_Asper<-subset(df_main_genus_cli,df_main_genus_cli$cluster_fungi=="Asper_type")
p2_sub<-ggplot(df_main_genus_cli_Asper,aes(cancer_type,Aspergillus,fill=cancer_type))+
  geom_boxplot(outlier.shape = NA,width=0.7)+
  geom_jitter(aes(fill=cancer_type),shape=21,position = position_jitter(0.2),
              size=2.5,alpha=0.5)+
  geom_signif(comparisons = list(c("Melanoma","NSCLC")),
              map_signif_level = T, 
              test = "wilcox.test",
              textsize = 4,
              y_position = c(0.76,0.80,0.88),
              tip_length = c(0,0,0),
              size=0.8,color="black")+
  theme_bw()+
  theme(axis.title = element_text(color="black",size = 12),
        #panel.grid = element_blank(),
        axis.text = element_text(color="black",size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 12,hjust = 0.5,face="italic"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm") )+
  labs(y="Relative abundance",title = "Aspergillus",x=NULL)+
  scale_y_continuous(limits = c(-0.01,1.01),
                     breaks = seq(0,1,0.1),
                     labels =  seq(0,1,0.1))
  #scale_fill_manual(values =c("Sacc_type"="#84c3b7","Pen_type"="#b8aeeb","Asper_type"="#f2b56f"))
p2_sub

df_main_genus_cli_Sacc<-subset(df_main_genus_cli,df_main_genus_cli$cluster_fungi=="Sacc_type")
p3_sub<-ggplot(df_main_genus_cli_Sacc,aes(cancer_type,Saccharomyces,fill=cancer_type))+
  geom_boxplot(outlier.shape = NA,width=0.7)+
  geom_jitter(aes(fill=cancer_type),shape=21,position = position_jitter(0.2),
              size=2.5,alpha=0.5)+
  geom_signif(comparisons = list(c("Melanoma","NSCLC")),
              map_signif_level = T, 
              test = "wilcox.test",
              textsize = 4,
              y_position = c(0.76,0.80,0.88),
              tip_length = c(0,0,0),
              size=0.8,color="black")+
  theme_bw()+
  theme(axis.title = element_text(color="black",size = 12),
        #panel.grid = element_blank(),
        axis.text = element_text(color="black",size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1,face="italic"),
        legend.position = "none",
        plot.title = element_text(size = 12,hjust = 0.5,face="italic"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm") )+
  labs(y="Relative abundance",title = "Saccharomyces",x=NULL)+
  scale_y_continuous(limits = c(-0.01,1.01),
                     breaks = seq(0,1,0.1),
                     labels =  seq(0,1,0.1))
  #scale_fill_manual(values = c("Sacc_type"="#84c3b7","Pen_type"="#b8aeeb","Asper_type"="#f2b56f"))
p3_sub

p3_sub2<-ggplot(df_main_genus_cli_Sacc,aes(cancer_type,Aspergillus,fill=cancer_type))+
  geom_boxplot(outlier.shape = NA,width=0.7)+
  geom_jitter(aes(fill=cancer_type),shape=21,position = position_jitter(0.2),
              size=2.5,alpha=0.5)+
  geom_signif(comparisons = list(c("Melanoma","NSCLC")),
              map_signif_level = T, 
              test = "wilcox.test",
              textsize = 4,
              y_position = c(0.76,0.80,0.88),
              tip_length = c(0,0,0),
              size=0.8,color="black")+
  theme_bw()+
  theme(axis.title = element_text(color="black",size = 12),
        #panel.grid = element_blank(),
        axis.text = element_text(color="black",size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 12,hjust = 0.5,face="italic"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm") )+
  labs(y="Relative abundance",title = "Aspergillus",x=NULL)+
  scale_y_continuous(limits = c(-0.01,1.01),
                     breaks = seq(0,1,0.1),
                     labels =  seq(0,1,0.1))
#scale_fill_manual(values = c("Sacc_type"="#84c3b7","Pen_type"="#b8aeeb","Asper_type"="#f2b56f"))
p3_sub2


df_main_genus_cli_Pen<-subset(df_main_genus_cli,df_main_genus_cli$cluster_fungi=="Pen_type")
p4_sub<-ggplot(df_main_genus_cli_Pen,aes(cancer_type,Penicillium,fill=cancer_type))+
  geom_boxplot(outlier.shape = NA,width=0.7)+
  geom_jitter(aes(fill=cancer_type),shape=21,position = position_jitter(0.2),
              size=2.5,alpha=0.5)+
  geom_signif(comparisons = list(c("Melanoma","NSCLC")),
              map_signif_level = T, 
              test = "wilcox.test",
              textsize = 3.5,
              y_position = c(0.76,0.80,0.88),
              tip_length = c(0,0,0),
              size=0.8,color="black")+
  theme_bw()+
  theme(axis.title = element_text(color="black",size = 12),
        #panel.grid = element_blank(),
        axis.text = element_text(color="black",size = 12),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 12,hjust = 0.5,face="italic"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),units = "cm") )+
  labs(y="Relative abundance",title = "Penicillium",x=NULL)+
  scale_y_continuous(limits = c(-0.01,1.01),
                     breaks = seq(0,1,0.1),
                     labels =  seq(0,1,0.1))
  #scale_fill_manual(values = c("Sacc_type"="#84c3b7","Pen_type"="#b8aeeb","Asper_type"="#f2b56f"))
p4_sub




###cluter中国家分布
f1<-as.data.frame(table(clin$cluster,clin$country))
colnames(f1)<-c("cluster","country","con")
f1$cluster<-paste0("cluster ",f1$cluster)
clu_count<-ddply(f1,'cluster',transform,ratio=con/sum(con)*100)

p6<-ggplot(clu_count,aes(x=cluster,y=ratio,fill=country))+geom_bar(stat = 'identity',width = 0.8)+
  theme( panel.grid = element_blank(),panel.background = element_blank(), strip.text = element_text(size =5)) +
  labs(x="",fill="Country",y="Ratio",title ="Country" )+
  scale_fill_manual(values = c("#f2b56f","#84c3b7","#b8aeeb","#b0a875","#f2a7da","#EE7E77"))+
  scale_x_discrete(labels = c("Cluster 1(n=333)","Cluster 2(n=310)","Cluster 3(n=277)"))+
  theme(axis.text = element_text(size =10), axis.title = element_text(size = 10), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size=unit(0.4,'cm'),
        plot.title = element_text(size=10,hjust=0.5),
        plot.subtitle = element_text(size=10,hjust=0.5,face="italic"))+
  scale_y_continuous(breaks = seq(0,100,25),labels=c('0.00','0.25','0.50','0.75','1.00'))+
  geom_text(aes(label=format(ratio/100,digits=1)),position=position_stack(vjust=0.5),size=3,color="white")

p6
table(clin$cluster)

##不同洲/大陆
table(clin$cluster,clin$country)
group<-c("Cluster1","Cluster2","Cluster3")
Europe<-c(246,247,186)
NorthAmerica<-c(87,63,91)
f1<-data.frame(group,Europe,NorthAmerica)
f2<-f1%>%pivot_longer(cols = c(Europe:NorthAmerica),names_to = 'Continent',values_to = 'con')
clu_cont<-ddply(f2,'group',transform,ratio=con/sum(con)*100)

p7<-ggplot(clu_cont,aes(x=group,y=ratio,fill=Continent))+geom_bar(stat = 'identity',width = 0.8)+
  theme( panel.grid = element_blank(),panel.background = element_blank(), strip.text = element_text(size =5),
         axis.line = element_blank()) +
  labs(x="",fill="Continent",y="Ratio",title = "Continent")+
  scale_fill_manual(values = c("#f2b56f","#84c3b7"))+
  scale_x_discrete(labels = c("Cluster 1(n=333)","Cluster 2(n=310)","Cluster 3(n=277)"))+
  theme(axis.text = element_text(size =10), axis.title = element_text(size = 10), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size=unit(0.4,'cm'),
        plot.title = element_text(size=10,hjust=0.5),
        plot.subtitle = element_text(size=10,hjust=0.5,face="italic"))+
  scale_y_continuous(breaks = seq(0,100,25),labels=c('0.00','0.25','0.50','0.75','1.00'))+
  geom_text(aes(label=format(ratio/100,digits=2)),position=position_stack(vjust=0.5),size=3,color="white")

p7

###cluter不同肿瘤分布
f1<-as.data.frame(table(clin$cluster,clin$cancer_type))
colnames(f1)<-c("cluster","cancer_type","con")
f1$cluster<-paste0("Cluster ",f1$cluster)
clu_cancer<-ddply(f1,'cluster',transform,ratio=con/sum(con)*100)

p8<-ggplot(clu_cancer,aes(x=cluster,y=ratio,fill=cancer_type))+geom_bar(stat = 'identity',width = 0.8)+
  theme( panel.grid = element_blank(),panel.background = element_blank(), strip.text = element_text(size =5)) +
  labs(x="",fill="Cancer type",y="Ratio",title = "Cancer type")+
  scale_fill_manual(values = c("#f2b56f","#84c3b7","#b8aeeb","#b0a875","#f2a7da"))+
  #scale_x_discrete(labels = c("Cluster 1(n=333)","Cluster 2(n=310)","Cluster 3(n=277)"))+
  theme(axis.text = element_text(size =10), axis.title = element_text(size = 10), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position="right",
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size=unit(0.4,'cm'),
        plot.title = element_text(size=10,hjust=0.5),
        plot.subtitle = element_text(size=10,hjust=0.5,face="italic"))+
  scale_y_continuous(breaks = seq(0,100,25),labels=c('0.00','0.25','0.50','0.75','1.00'))+
  geom_text(aes(label=format(ratio/100,digits=1)),position=position_stack(vjust=0.5),size=3,color="white")

p8
p9<-p6+p7+p8+plot_layout(guides = 'collect')
p9
ggsave("10.pics/supplemental figures/Enterotype/fungi/fungi_enterotype_country_continent_cancer.jpeg",width =200,height =150,units = "mm",dpi = 300,plot = p9)
