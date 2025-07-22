### clinical information
### Liurong
### 2023-09-27

##install.packages("readr")
library("readr")

setwd("D:/课题/细菌—真菌互作/3_R/2.rawData/")
rm(list=ls())

###################################################################################
###  PRJNA397906
clin_397906<-read.csv("clinical/metadata/PRJNA397906_clinical.csv",header=T)
id_397906<-read.csv("clinical/metadata/PRJNA397906_id.csv",header=T)

colnames(clin_397906)[1]<-"Patient"
clin_397906_whole<-merge(clin_397906,id_397906,by.x="Patient",by.y="Patient")

names(clin_397906_whole)[names(clin_397906_whole) =="Run"] <- "id"
names(clin_397906_whole)[names(clin_397906_whole) =="Age"] <- "age"
names(clin_397906_whole)[names(clin_397906_whole) =="RECIST.Category"] <- "recist"
names(clin_397906_whole)[names(clin_397906_whole) =="Change.in.Tumor.Size...."] <- "recist_ratio"

nr<-nrow(clin_397906_whole)
for(i in 1:nr){
  if(is.na(clin_397906_whole$Sex[i])==T){clin_397906_whole$gender[i]=NA}
  if(is.na(clin_397906_whole$Sex[i])==F){
    if (clin_397906_whole$Sex[i]=="F"){clin_397906_whole$gender[i]="Female"}
    else if (clin_397906_whole$Sex[i]=="M"){clin_397906_whole$gender[i]="Male"}
  }
}


nr<-nrow(clin_397906_whole)
for(i in 1:nr){
  if(is.na(clin_397906_whole$recist[i])==T){clin_397906_whole$response_code[i]=NA}
  if(is.na(clin_397906_whole$recist[i])==F){
    if (clin_397906_whole$recist[i]=="Response"){clin_397906_whole$response_code[i]=1}
    else if (clin_397906_whole$recist[i]=="Progression"){clin_397906_whole$response_code[i]=0}
    else if (clin_397906_whole$recist[i]=="Stable"){clin_397906_whole$response_code[i]=0}
  }
}

clin_397906_whole$os<-NA
clin_397906_whole$os_event<-NA
clin_397906_whole$pfs<-NA
clin_397906_whole$pfs_event<-NA
clin_397906_whole$dataset<-"PRJNA397906"
clin_397906_whole$cancer_type<-"Melanoma"
clin_397906_whole$LDH<-NA
clin_397906_whole$BMI<-NA
clin_397906_whole$pfs_12_months<-NA
clin_397906_whole$pfs_6_months<-NA
clin_397906_whole$irAEs<-NA
clin_397906_whole$irAEs_grade<-NA
clin_397906_whole$study<-"FrankelAE_2017"
clin_397906_whole$country<-"USA"
clin_397906_whole$Instrument<-"Illumina HiSeq 2000"

clin_397906_Count<-read.csv("read_number/PRJNA397906.csv",header=T)
colnames(clin_397906_Count)<-c("id","read_count")
clin_397906_whole<-merge(clin_397906_whole,clin_397906_Count,by.x="id",by.y="id")

clin_397906_sub<-clin_397906_whole[,c("id","age","gender","os","os_event","pfs","pfs_event","drug","recist",
                                      "study","dataset","cancer_type","recist_ratio","LDH","BMI","response_code","target","read_count","pfs_6_months",
                                      "pfs_12_months","irAEs","irAEs_grade","country","Instrument")]


################################################################################################################
### PRJNA399742

clin_399742<-read.csv("clinical/metadata/human_16S.sampleinfo.csv",header=T)
id_399742<-read.csv("clinical/metadata/PRJNA399742_Run_table.csv",header=T)

colnames(clin_399742)[1]<-"subject_id"
clin_399742_whole<-merge(clin_399742,id_399742,by.x="subject_id",by.y="subject_id")

names(clin_399742_whole)[names(clin_399742_whole) =="Run"] <- "id"
names(clin_399742_whole)[names(clin_399742_whole) =="Response"] <- "response"
names(clin_399742_whole)[names(clin_399742_whole) =="BORR"] <- "recist"
names(clin_399742_whole)[names(clin_399742_whole) =="RECIST."] <- "recist_ratio"
names(clin_399742_whole)[names(clin_399742_whole) =="Sex"] <- "gender"
names(clin_399742_whole)[names(clin_399742_whole) =="PFS.6.months..yes.no."] <- "pfs_6_months"

nr<-nrow(clin_399742_whole)
for(i in 1:nr){
  if(is.na(clin_399742_whole$recist[i])==T){clin_399742_whole$response_code[i]=NA}
  if(is.na(clin_399742_whole$recist[i])==F){
    if (clin_399742_whole$recist[i]=="Complete Response"){clin_399742_whole$response_code[i]=1}
    else if (clin_399742_whole$recist[i]=="Partial Response"){clin_399742_whole$response_code[i]=1}
    else if (clin_399742_whole$recist[i]=="Progressive Disease"){clin_399742_whole$response_code[i]=0}
    else if (clin_399742_whole$recist[i]=="Stable Disease"){clin_399742_whole$response_code[i]=0}
  }
}


for(i in 1:nr){
  if(is.na(clin_399742_whole$pfs_6_months[i])==T){clin_399742_whole$pfs_6_months[i]=NA}
  if(is.na(clin_399742_whole$pfs_6_months[i])==F){
    if (clin_399742_whole$pfs_6_months[i]=="Yes"){clin_399742_whole$pfs_6_months[i]=1}
    else if (clin_399742_whole$pfs_6_months[i]=="No"){clin_399742_whole$pfs_6_months[i]=0}
    else if (clin_399742_whole$pfs_6_months[i]=="No, Still on TX"){clin_399742_whole$pfs_6_months[i]=NA}
  }
}


clin_399742_whole$os<-NA
clin_399742_whole$os_event<-NA
clin_399742_whole$pfs<-NA
clin_399742_whole$pfs_event<-NA
clin_399742_whole$dataset<-"PRJNA399742"
clin_399742_whole$study<-"MatsonV_2018"
clin_399742_whole$cancer_type<-"Melanoma"
clin_399742_whole$age<-NA
clin_399742_whole$LDH<-NA
clin_399742_whole$BMI<-NA
clin_399742_whole$pfs_12_months<-NA
clin_399742_whole$irAEs<-NA
clin_399742_whole$irAEs_grade<-NA
clin_399742_whole$country<-"USA"
clin_399742_whole$Instrument<-paste("Illumina",clin_399742_whole$Instrument)

clin_399742_Count<-read.csv("read_number/PRJNA399742.csv",header=T)
colnames(clin_399742_Count)<-c("id","read_count")
clin_399742_whole<-merge(clin_399742_whole,clin_399742_Count,by.x="id",by.y="id")

clin_399742_sub<-clin_399742_whole[,c("id","age","gender","os","os_event","pfs","pfs_event","drug","recist",
                                      "study","dataset","cancer_type","recist_ratio","LDH","BMI","response_code","target","read_count","pfs_6_months",
                                      "pfs_12_months","irAEs","irAEs_grade","country","Instrument")]


####################################################################################################
## PRJNA541981
clin_541981<-read.csv("clinical/metadata/PRJNA541981_clinical.csv",header=T)

names(clin_541981)[names(clin_541981) =="AGE"] <- "age"
names(clin_541981)[names(clin_541981) =="dead"] <- "os_event"
names(clin_541981)[names(clin_541981) =="any_prog"] <- "pfs_event"
names(clin_541981)[names(clin_541981) =="immunotherapy"] <- "target"
names(clin_541981)[names(clin_541981) =="sex"] <- "gender"

clin_541981$gender<-ifelse(clin_541981$gender=="female","Female","Male")

clin_541981$os<-clin_541981$os_months/12
clin_541981$pfs<-clin_541981$pfs_months/12

nr<-nrow(clin_541981)
for(i in 1:nr){
  if(is.na(clin_541981$pfs[i])==T | is.na(clin_541981$pfs_event[i])==T){clin_541981$pfs_6_months[i]=NA}
  if(is.na(clin_541981$pfs[i])==F & is.na(clin_541981$pfs_event[i])==F){
    if (clin_541981$pfs[i]>=0.5){clin_541981$pfs_6_months[i]=1}
    else if (clin_541981$pfs[i]<0.5 & clin_541981$pfs_event[i]==1){clin_541981$pfs_6_months[i]=0}
    else if (clin_541981$pfs[i]<0.5 & clin_541981$pfs_event[i]==0){clin_541981$pfs_6_months[i]=NA}
  }
}

nr<-nrow(clin_541981)
for(i in 1:nr){
  if(is.na(clin_541981$pfs[i])==T | is.na(clin_541981$pfs_event[i])==T){clin_541981$pfs_12_months[i]=NA}
  if(is.na(clin_541981$pfs[i])==F & is.na(clin_541981$pfs_event[i])==F){
    if (clin_541981$pfs[i]>=1){clin_541981$pfs_12_months[i]=1}
    else if (clin_541981$pfs[i]<1 & clin_541981$pfs_event[i]==1){clin_541981$pfs_12_months[i]=0}
    else if (clin_541981$pfs[i]<1 & clin_541981$pfs_event[i]==0){clin_541981$pfs_12_months[i]=NA}
  }
}

clin_541981$dataset<-"PRJNA541981"
clin_541981$study<-"PetersBA_2019"
clin_541981$cancer_type<-"Melanoma"
clin_541981$drug<-NA
clin_541981$recist_ratio<-NA
clin_541981$recist<-NA
clin_541981$response_code<-NA
clin_541981$irAEs<-NA
clin_541981$irAEs_grade<-NA
clin_541981$country<-"USA"

clin_541981_Count<-read.csv("read_number/PRJNA541981.csv",header=T)
colnames(clin_541981_Count)<-c("id","read_count")
clin_541981<-merge(clin_541981,clin_541981_Count,by.x="id",by.y="id")

clin_541981<-subset(clin_541981,time=="Baseline",drop=T)

clin_541981_sub<-clin_541981[,c("id","age","gender","os","os_event","pfs","pfs_event","drug","recist",
                                "study","dataset","cancer_type","recist_ratio","LDH","BMI","response_code","target","read_count","pfs_6_months",
                                "pfs_12_months","irAEs","irAEs_grade","country","Instrument")]


########################################################################################################
## PRJEB22863

clin_22863<-read.csv("clinical/metadata/PRJEB22863_clinical.csv",header=T)
id_22863<-read.csv("clinical/metadata/PRJEB22863_SraRunInfo.csv",header=T)
colnames(clin_22863)[1]<-"BioSample"

clin_22863_whole<-merge(id_22863,clin_22863,by.x="BioSample",by.y="BioSample")

names(clin_22863_whole)[names(clin_22863_whole) =="Run"] <- "id"
names(clin_22863_whole)[names(clin_22863_whole) =="host_age"] <- "age"
names(clin_22863_whole)[names(clin_22863_whole) =="host_disease_status"] <- "cancer_type"
names(clin_22863_whole)[names(clin_22863_whole) =="host_sex"] <- "gender"
names(clin_22863_whole)[names(clin_22863_whole) =="response"] <- "recist"
clin_22863_whole$gender<-ifelse(clin_22863_whole$gender=="female","Female","Male")

nr<-nrow(clin_22863_whole)
for(i in 1:nr){
  if(is.na(clin_22863_whole$recist[i])==T){clin_22863_whole$response_code[i]=NA}
  if(is.na(clin_22863_whole$recist[i])==F){
    if (clin_22863_whole$recist[i]=="Partialresponse"){clin_22863_whole$response_code[i]=1}
    else if (clin_22863_whole$recist[i]=="Dead"){clin_22863_whole$response_code[i]=0}
    else if (clin_22863_whole$recist[i]=="Progression"){clin_22863_whole$response_code[i]=0}
    else if (clin_22863_whole$recist[i]=="Stable"){clin_22863_whole$response_code[i]=0}
  }
}

nr<-nrow(clin_22863_whole)
for(i in 1:nr){
  if(is.na(clin_22863_whole$PFS_6_month[i])==T){clin_22863_whole$pfs_6_months[i]=NA}
  if(is.na(clin_22863_whole$PFS_6_month[i])==F){
    if (clin_22863_whole$PFS_6_month[i]=="Above"){clin_22863_whole$pfs_6_months[i]=1}
    else if (clin_22863_whole$PFS_6_month[i]=="Below"){clin_22863_whole$pfs_6_months[i]=0}
    else if (clin_22863_whole$PFS_6_month[i]=="Notdetermined"){clin_22863_whole$pfs_6_months[i]=NA}
  }
}

clin_22863_whole$dataset<-"PRJEB22863"
clin_22863_whole$study<-"RoutyB_2018"
clin_22863_whole$drug<-NA
clin_22863_whole$recist_ratio<-NA
clin_22863_whole$recist<-NA
clin_22863_whole$LDH<-NA
clin_22863_whole$BMI<-NA
clin_22863_whole$os<-NA
clin_22863_whole$os_event<-NA
clin_22863_whole$pfs<-NA
clin_22863_whole$pfs_event<-NA
clin_22863_whole$pfs_12_months<-NA
clin_22863_whole$irAEs<-NA
clin_22863_whole$irAEs_grade<-NA
clin_22863_whole$country<-"France"
clin_22863_whole$Instrument<-"Ion Torrent Proton"


clin_22863_Count<-read.csv("read_number/PRJEB22863.csv",header=T)
colnames(clin_22863_Count)<-c("id","read_count")
clin_22863_whole<-merge(clin_22863_whole,clin_22863_Count,by.x="id",by.y="id")

clin_22863_sub<-clin_22863_whole[,c("id","age","gender","os","os_event","pfs","pfs_event","drug","recist",
                                    "study","dataset","cancer_type","recist_ratio","LDH","BMI","response_code","target","read_count","pfs_6_months",
                                    "pfs_12_months","irAEs","irAEs_grade","country","Instrument")]




###########################################################################################################
### PRJNA770295

clin_770295<-read.csv("clinical/metadata/PRJNA770295_clinical.csv",header=T)

names(clin_770295)[names(clin_770295) =="Run"] <- "id"
names(clin_770295)[names(clin_770295) =="Age"] <- "age"
names(clin_770295)[names(clin_770295) =="Sex"] <- "gender"
clin_770295$gender<-ifelse(clin_770295$gender=="female","Female","Male")

clin_770295$pfs<-clin_770295$pfs_d/365
clin_770295$pfs_event<-clin_770295$pfsevent
clin_770295$response_code<-clin_770295$response

nr<-nrow(clin_770295)
for(i in 1:nr){
  if(is.na(clin_770295$pfs[i])==T | is.na(clin_770295$pfs_event[i])==T){clin_770295$pfs_6_months[i]=NA}
  if(is.na(clin_770295$pfs[i])==F & is.na(clin_770295$pfs_event[i])==F){
    if (clin_770295$pfs[i]>=0.5){clin_770295$pfs_6_months[i]=1}
    else if (clin_770295$pfs[i]<0.5 & clin_770295$pfs_event[i]==1){clin_770295$pfs_6_months[i]=0}
    else if (clin_770295$pfs[i]<0.5 & clin_770295$pfs_event[i]==0){clin_770295$pfs_6_months[i]=NA}
  }
}

nr<-nrow(clin_770295)
for(i in 1:nr){
  if(is.na(clin_770295$pfs[i])==T | is.na(clin_770295$pfs_event[i])==T){clin_770295$pfs_12_months[i]=NA}
  if(is.na(clin_770295$pfs[i])==F & is.na(clin_770295$pfs_event[i])==F){
    if (clin_770295$pfs[i]>=1){clin_770295$pfs_12_months[i]=1}
    else if (clin_770295$pfs[i]<1 & clin_770295$pfs_event[i]==1){clin_770295$pfs_12_months[i]=0}
    else if (clin_770295$pfs[i]<1 & clin_770295$pfs_event[i]==0){clin_770295$pfs_12_months[i]=NA}
  }
}

clin_770295$dataset<-"PRJNA770295"
clin_770295$study<-"SpencerCN_2021"
clin_770295$drug<-NA
clin_770295$recist_ratio<-NA
clin_770295$recist<-NA
clin_770295$os<-NA
clin_770295$os_event<-NA
clin_770295$cancer_type<-"Melanoma"
clin_770295$irAEs<-NA
clin_770295$irAEs_grade<-NA
clin_770295$country<-"USA" 

clin_770295_Count<-read.csv("read_number/PRJNA770295.csv",header=T)
colnames(clin_770295_Count)<-c("id","read_count")
clin_770295<-merge(clin_770295,clin_770295_Count,by.x="id",by.y="id")

clin_770295<-subset(clin_770295,treatment!="other systemic",drop=T)

clin_770295_sub<-clin_770295[,c("id","age","gender","os","os_event","pfs","pfs_event","drug","recist",
                                "study","dataset","cancer_type","recist_ratio","LDH","BMI","response_code","target","read_count","pfs_6_months",
                                "pfs_12_months","irAEs","irAEs_grade","country","Instrument")]



###########################################################################################################
### PRJNA1023797

clin_1023797<-read.csv("clinical/metadata/PRJNA1023797_clinical.csv",header=T)
id_1023797<-read.csv("clinical/metadata/PRJNA1023797_id.csv",header=T)

colnames(id_1023797)<-c("id","SampleID")
colnames(clin_1023797)[2]<-"SampleID"
clin_1023797<-merge(id_1023797,clin_1023797,by.x="SampleID",by.y="SampleID")
table(clin_1023797$target)
clin_1023797<-subset(clin_1023797,target=="PD-1"|target=="PD-L1"|target=="CTLA-4")

#

names(clin_1023797)[names(clin_1023797) =="Age"] <- "age"
names(clin_1023797)[names(clin_1023797) =="Gender"] <- "gender"
names(clin_1023797)[names(clin_1023797) =="Death"] <- "os_event"
names(clin_1023797)[names(clin_1023797) =="PD"] <- "pfs_event"
names(clin_1023797)[names(clin_1023797) =="Center"] <- "country"

nr<-nrow(clin_1023797)
table(clin_1023797$Best_response)
for(i in 1:nr){
  if(is.na(clin_1023797$Best_response[i])==T){clin_1023797$response_code[i]=NA}
  if(is.na(clin_1023797$Best_response[i])==F){
    if (clin_1023797$Best_response[i]=="CR"){clin_1023797$response_code[i]=1}
    else if (clin_1023797$Best_response[i]=="PR"|clin_1023797$Best_response[i]=="PR/dissociated response"){clin_1023797$response_code[i]=1}
    else if (clin_1023797$Best_response[i]=="PD"|
             clin_1023797$Best_response[i]=="PD/Death"){clin_1023797$response_code[i]=0}
    else if (clin_1023797$Best_response[i]=="SD"|
             clin_1023797$Best_response[i]=="Death"|
             clin_1023797$Best_response[i]=="SD/dissociated response"|
             clin_1023797$Best_response[i]=="SD/PD"){clin_1023797$response_code[i]=0}
  }
}

clin_1023797$recist<-clin_1023797$Best_response

clin_1023797$os<-clin_1023797$OS/12
clin_1023797$pfs<-clin_1023797$PFS/12

for(i in 1:nr){
  if(is.na(clin_1023797$pfs[i])==T | is.na(clin_1023797$pfs_event[i])==T){clin_1023797$pfs_12_months[i]=NA}
  if(is.na(clin_1023797$pfs[i])==F & is.na(clin_1023797$pfs_event[i])==F){
    if (clin_1023797$pfs[i]>=1){clin_1023797$pfs_12_months[i]=1}
    else if (clin_1023797$pfs[i]<1 & clin_1023797$pfs_event[i]==1){clin_1023797$pfs_12_months[i]=0}
    else if (clin_1023797$pfs[i]<1 & clin_1023797$pfs_event[i]==0){clin_1023797$pfs_12_months[i]=NA}
  }
}

for(i in 1:nr){
  if(is.na(clin_1023797$pfs[i])==T | is.na(clin_1023797$pfs_event[i])==T){clin_1023797$pfs_6_months[i]=NA}
  if(is.na(clin_1023797$pfs[i])==F & is.na(clin_1023797$pfs_event[i])==F){
    if (clin_1023797$pfs[i]>=0.5){clin_1023797$pfs_6_months[i]=1}
    else if (clin_1023797$pfs[i]<0.5 & clin_1023797$pfs_event[i]==1){clin_1023797$pfs_6_months[i]=0}
    else if (clin_1023797$pfs[i]<0.5 & clin_1023797$pfs_event[i]==0){clin_1023797$pfs_6_months[i]=NA}
  }
}

clin_1023797$study<-"DerosaL_2024"
clin_1023797$dataset<-"PRJNA1023797"
clin_1023797$drug<-NA
clin_1023797$recist_ratio<-NA
clin_1023797$cancer_type<-"NSCLC"
clin_1023797$LDH<-NA
clin_1023797$irAEs<-NA
clin_1023797$irAEs_grade<-NA
clin_1023797$Instrument<-"Ion Torrent Proton"

clin_1023797_Count<-read.csv("read_number/PRJNA1023797_sub.csv",header=T)
colnames(clin_1023797_Count)<-c("id","read_count")
clin_1023797<-merge(clin_1023797,clin_1023797_Count,by.x="id",by.y="id")

clin_1023797_sub<-clin_1023797[,c("id","age","gender","os","os_event","pfs","pfs_event","drug","recist",
                                "study","dataset","cancer_type","recist_ratio","LDH","BMI","response_code","target","read_count","pfs_6_months",
                                "pfs_12_months","irAEs","irAEs_grade","country","Instrument")]


###############################################################################################
## PRJNA762360
clin_762360<-read.csv("clinical/metadata/PRJNA762360_clinical.csv",header=T)

clin_762360<-subset(clin_762360,Cohort=="Pittsburgh_Early",drop=T)
##clin_762360<-subset(clin_762360,timepoint_days<=120,drop=T)

names(clin_762360)[names(clin_762360) =="Run"] <- "id"
names(clin_762360)[names(clin_762360) =="AGE"] <- "age"
names(clin_762360)[names(clin_762360) =="Treatment_Response_3_mth"] <- "response"
#names(clin_762360)[names(clin_762360) =="Study_Clin_Response"] <- "response"
names(clin_762360)[names(clin_762360) =="Death_Event"] <- "os_event"
names(clin_762360)[names(clin_762360) =="Progression_Event"] <- "pfs_event"
names(clin_762360)[names(clin_762360) =="Pre_Treatment_Disease_Stage"] <- "stage"
names(clin_762360)[names(clin_762360) =="Immunotherapy_Drug"] <- "drug"
names(clin_762360)[names(clin_762360) =="sex"] <- "gender"

clin_762360$gender<-ifelse(clin_762360$gender=="female","Female","Male")

nr<-nrow(clin_762360)
for(i in 1:nr){
  if(is.na(clin_762360$response[i])==T){clin_762360$response_code[i]=NA}
  if(is.na(clin_762360$response[i])==F){
    if (clin_762360$response[i]=="CR"|clin_762360$response[i]=="PR"){clin_762360$response_code[i]=1}
    else if (clin_762360$response[i]=="PD"|clin_762360$response[i]=="SD"){clin_762360$response_code[i]=0}
  }
}

clin_762360$os<-clin_762360$overall_survival_days/365
clin_762360$pfs<-clin_762360$PFS_days/365

nr<-nrow(clin_762360)
for(i in 1:nr){
  if(is.na(clin_762360$pfs[i])==T | is.na(clin_762360$pfs_event[i])==T){clin_762360$pfs_6_months[i]=NA}
  if(is.na(clin_762360$pfs[i])==F & is.na(clin_762360$pfs_event[i])==F){
    if (clin_762360$pfs[i]>=0.5){clin_762360$pfs_6_months[i]=1}
    else if (clin_762360$pfs[i]<0.5 & clin_762360$pfs_event[i]==1){clin_762360$pfs_6_months[i]=0}
    else if (clin_762360$pfs[i]<0.5 & clin_762360$pfs_event[i]==0){clin_762360$pfs_6_months[i]=NA}
  }
}

nr<-nrow(clin_762360)
for(i in 1:nr){
  if(is.na(clin_762360$pfs[i])==T | is.na(clin_762360$pfs_event[i])==T){clin_762360$pfs_12_months[i]=NA}
  if(is.na(clin_762360$pfs[i])==F & is.na(clin_762360$pfs_event[i])==F){
    if (clin_762360$pfs[i]>=1){clin_762360$pfs_12_months[i]=1}
    else if (clin_762360$pfs[i]<1 & clin_762360$pfs_event[i]==1){clin_762360$pfs_12_months[i]=0}
    else if (clin_762360$pfs[i]<1 & clin_762360$pfs_event[i]==0){clin_762360$pfs_12_months[i]=NA}
  }
}

clin_762360$study<-"McCullochJA_2022"
clin_762360$dataset<-"PRJNA762360"
clin_762360$cancer_type<-"Melanoma"
clin_762360$drug<-NA
clin_762360$recist_ratio<-NA
clin_762360$recist<-NA
clin_762360$country<-"USA"

clin_762360_Count<-read.csv("read_number/PRJNA762360.csv",header=T)
colnames(clin_762360_Count)<-c("id","read_count")
clin_762360<-merge(clin_762360,clin_762360_Count,by.x="id",by.y="id")

table(clin_762360$Most_Severe_IrAE)

nr<-nrow(clin_762360)
for(i in 1:nr){
  if(is.na(clin_762360$Most_Severe_IrAE[i])==T){clin_762360$irAEs[i]=NA}
  if(is.na(clin_762360$Most_Severe_IrAE[i])==F){
    if (clin_762360$Most_Severe_IrAE[i]=="Grades_1-2"){clin_762360$irAEs[i]=1
    clin_762360$irAEs_grade[i]="G1-G2" }
    else if (clin_762360$Most_Severe_IrAE[i]=="Grades_3-4"){clin_762360$irAEs[i]=1
    clin_762360$irAEs_grade[i]="G3-G4"}
    else if (clin_762360$Most_Severe_IrAE[i]=="No_AE_reported"){clin_762360$irAEs[i]=0
    clin_762360$irAEs_grade[i]=NA}
  }
}

clin_762360_sub<-clin_762360[,c("id","age","gender","os","os_event","pfs","pfs_event","drug","recist",
                                "study","dataset","cancer_type","recist_ratio","LDH","BMI","response_code","target","read_count","pfs_6_months",
                                "pfs_12_months","irAEs","irAEs_grade","country","Instrument")]



###############################################################################################
## PRJEB43119
clin_43119<-read.csv("clinical/metadata/PRJEB43119_clinical.csv",header=T)

names(clin_43119)[names(clin_43119) =="Run"] <- "id"
names(clin_43119)[names(clin_43119) =="age"] <- "age_bin"
names(clin_43119)[names(clin_43119)=="geo_loc_name_country"] <- "country"
names(clin_43119)[names(clin_43119)=="host_body_mass_index"] <- "BMI"
names(clin_43119)[names(clin_43119)=="host_sex"] <- "gender"
names(clin_43119)[names(clin_43119)=="host_age_5yr_bin"] <- "age_bin"
clin_43119$gender<-ifelse(clin_43119$gender=="female","Female","Male")


nr<-nrow(clin_43119)
for(i in 1:nr){
  if(is.na(clin_43119$ORR[i])==T){clin_43119$response_code[i]=NA}
  if(is.na(clin_43119$ORR[i])==F){
    if (clin_43119$ORR[i]=="CR"){clin_43119$response_code[i]=1}
    if (clin_43119$ORR[i]=="PR"){clin_43119$response_code[i]=1}
    else if (clin_43119$ORR[i]=="SD"){clin_43119$response_code[i]=0}
    else if (clin_43119$ORR[i]=="PD"){clin_43119$response_code[i]=0}
  }
}

nr<-nrow(clin_43119)
for(i in 1:nr){
  if(is.na(clin_43119$PFS12[i])==T){clin_43119$pfs_12_months[i]=NA}
  if(is.na(clin_43119$PFS12[i])==F){
    if (clin_43119$PFS12[i]=="yes"){clin_43119$pfs_12_months[i]=1}
    else if (clin_43119$PFS12[i]=="no"){clin_43119$pfs_12_months[i]=0}
  }
}

clin_43119$study<-"LeeKA_2022"
clin_43119$dataset<-"PRJEB43119"
clin_43119$cancer_type<-"Melanoma"
clin_43119$drug<-NA
clin_43119$recist_ratio<-NA
clin_43119$recist<-NA
clin_43119$age<-NA
clin_43119$os<-NA
clin_43119$os_event<-NA
clin_43119$pfs<-NA
clin_43119$pfs_event<-NA
clin_43119$LDH<-NA
clin_43119$target<-NA
clin_43119$pfs_6_months<-NA
clin_43119$irAEs<-NA
clin_43119$irAEs_grade<-NA
clin_43119$Instrument<-"Illumina NovaSeq"

clin_43119_Count<-read.csv("read_number/PRJEB43119.csv",header=T)
colnames(clin_43119_Count)<-c("id","read_count")
clin_43119<-merge(clin_43119,clin_43119_Count,by.x="id",by.y="id")

clin_43119_sub<-clin_43119[,c("id","age","gender","os","os_event","pfs","pfs_event","drug","recist",
                              "study","dataset","cancer_type","recist_ratio","LDH","BMI","response_code","target","read_count","pfs_6_months",
                              "pfs_12_months","irAEs","irAEs_grade","age_bin","country","Instrument")]


#############################################################################################################
### 合并数据集

clinical_sub<-rbind(clin_770295_sub,clin_22863_sub,clin_762360_sub,
                    clin_397906_sub,clin_1023797_sub,clin_541981_sub,clin_399742_sub)

min(clinical_sub$age,na.rm=T)
nr<-nrow(clinical_sub)

for(i in 1:nr){
  if(is.na(clinical_sub$age[i])==T){clinical_sub$age_bin[i]=NA}
  else if (is.na(clinical_sub$age[i])==F){
    if (clinical_sub$age[i]>=20 & clinical_sub$age[i]<25){clinical_sub$age_bin[i]="20-25"}
    else if (clinical_sub$age[i]>=25 & clinical_sub$age[i]<30){clinical_sub$age_bin[i]="25-30"}
    else if (clinical_sub$age[i]>=30 & clinical_sub$age[i]<35){clinical_sub$age_bin[i]="30-35"}
    else if (clinical_sub$age[i]>=35 & clinical_sub$age[i]<40){clinical_sub$age_bin[i]="35-40"}
    else if (clinical_sub$age[i]>=40 & clinical_sub$age[i]<45){clinical_sub$age_bin[i]="40-45"}
    else if (clinical_sub$age[i]>=45 & clinical_sub$age[i]<50){clinical_sub$age_bin[i]="45-50"}
    else if (clinical_sub$age[i]>=50 & clinical_sub$age[i]<55){clinical_sub$age_bin[i]="50-55"}
    else if (clinical_sub$age[i]>=55 & clinical_sub$age[i]<60){clinical_sub$age_bin[i]="55-60"}
    else if (clinical_sub$age[i]>=60 & clinical_sub$age[i]<65){clinical_sub$age_bin[i]="60-65"}
    else if (clinical_sub$age[i]>=65 & clinical_sub$age[i]<70){clinical_sub$age_bin[i]="65-70"}
    else if (clinical_sub$age[i]>=70 & clinical_sub$age[i]<75){clinical_sub$age_bin[i]="70-75"}
    else if (clinical_sub$age[i]>=75 & clinical_sub$age[i]<80){clinical_sub$age_bin[i]="75-80"}
    else if (clinical_sub$age[i]>=80 & clinical_sub$age[i]<85){clinical_sub$age_bin[i]="80-85"}
    else if (clinical_sub$age[i]>=85 & clinical_sub$age[i]<90){clinical_sub$age_bin[i]="85-90"}
    else if (clinical_sub$age[i]>=90 & clinical_sub$age[i]<95){clinical_sub$age_bin[i]="90-95"}
  }
}

clinical_whole<-rbind(clinical_sub,clin_43119_sub)

table(clinical_whole$dataset,clinical_whole$os_event)

median(clin_1023797_sub$os)   ### 0.92
median(clin_762360_sub$os)   ### 2.47

nr<-nrow(clinical_whole)

for(i in 1:nr){
  if(is.na(clinical_whole$os[i])==T | is.na(clinical_whole$os_event[i])==T){clinical_whole$os_1_year[i]=NA}
  if(is.na(clinical_whole$os[i])==F & is.na(clinical_whole$os_event[i])==F){
    if (clinical_whole$os[i]>=1){clinical_whole$os_1_year[i]=1}
    else if (clinical_whole$os[i]<1 & clinical_whole$os_event[i]==1){clinical_whole$os_1_year[i]=0}
    else if (clinical_whole$os[i]<1 & clinical_whole$os_event[i]==0){clinical_whole$os_1_year[i]=NA}
  }
}

for(i in 1:nr){
  if(is.na(clinical_whole$os[i])==T | is.na(clinical_whole$os_event[i])==T){clinical_whole$os_2_year[i]=NA}
  if(is.na(clinical_whole$os[i])==F & is.na(clinical_whole$os_event[i])==F){
    if (clinical_whole$os[i]>=2){clinical_whole$os_2_year[i]=1}
    else if (clinical_whole$os[i]<2 & clinical_whole$os_event[i]==1){clinical_whole$os_2_year[i]=0}
    else if (clinical_whole$os[i]<2 & clinical_whole$os_event[i]==0){clinical_whole$os_2_year[i]=NA}
  }
}

##BMI分组
for(i in 1:nr){
  if(is.na(clinical_whole$BMI[i])==T ){clinical_whole$BMI_group[i]=NA}
  if(is.na(clinical_whole$BMI[i])==F ){
    if (clinical_whole$BMI[i]<25){clinical_whole$BMI_group[i]="lean"}
    else if (clinical_whole$BMI[i]>=25 &clinical_whole$BMI[i]<=30){clinical_whole$BMI_group[i]="overweight"}
    else if (clinical_whole$BMI[i]>30){clinical_whole$BMI_group[i]="Obese"}
  }
}

##age分组
table(clinical_whole$age_bin)
for(i in 1:nr){
  if(is.na(clinical_whole$age_bin[i])==T ){clinical_whole$age_group[i]=NA}
  if(is.na(clinical_whole$age_bin[i])==F ){
    if (clinical_whole$age_bin[i]=="15-20"|clinical_whole$age_bin[i]=="20-25"|
        clinical_whole$age_bin[i]=="25-30"|clinical_whole$age_bin[i]=="30-35"|
        clinical_whole$age_bin[i]=="35-40"|clinical_whole$age_bin[i]=="40-45"|
        clinical_whole$age_bin[i]=="45-50"|clinical_whole$age_bin[i]=="50-55"|
        clinical_whole$age_bin[i]=="55-60"){clinical_whole$age_group[i]="<60"}
    else if (clinical_whole$age_bin[i]=="60-65"|clinical_whole$age_bin[i]=="65-70"|
             clinical_whole$age_bin[i]=="70-75"|clinical_whole$age_bin[i]=="75-80"|
             clinical_whole$age_bin[i]=="80-85"|clinical_whole$age_bin[i]=="85-90"|
             clinical_whole$age_bin[i]=="90-95"){clinical_whole$age_group[i]="≥60"}
  }
}

write.csv(clinical_whole,"clinical/final/cli_whole_samples_new.csv",row.names=F)

##########################################
## response 

nrow(subset(clinical_whole,is.na(response_code)==T,drop=T))


clin_resp<-subset(clinical_whole,is.na(response_code)==F,drop=T)
table(clin_resp$cancer_type)

table(clin_resp$study)

melanoma<-subset(clin_resp,cancer_type=="Melanoma",drop=T)
table(melanoma$study)

NSCLC<-subset(clin_resp,cancer_type=="NSCLC",drop=T)
table(NSCLC$study)

RCC<-subset(clin_resp,cancer_type=="RCC",drop=T)
table(RCC$study)

clin_pfs<-subset(clinical_whole,is.na(pfs_12_months)==F,drop=T)
table(clin_pfs$cancer_type)

nrow(subset(clinical_whole,is.na(pfs_12_months)==T,drop=T))


melanoma<-subset(clin_pfs,cancer_type=="Melanoma",drop=T)
table(melanoma$study)

nrow(subset(clinical_whole,is.na(os)==T,drop=T))
nrow(subset(clinical_whole,is.na(os_event)==T,drop=T))

clin_os<-subset(clinical_whole,is.na(os)==F & is.na(os_event)==F ,drop=T)
table(clin_os$cancer_type)

melanoma<-subset(clin_os,cancer_type=="Melanoma",drop=T)
table(melanoma$study)

NSCLC<-subset(clin_os,cancer_type=="NSCLC",drop=T)
table(NSCLC$study)


clin_irAEs<-subset(clinical_whole,is.na(irAEs)==F,drop=T)
table(clin_irAEs$cancer_type)

melanoma<-subset(clin_irAEs,cancer_type=="Melanoma",drop=T)
table(melanoma$study)




#############################################################################
#### table 1B(basic characteristics)

clinical<-clinical_whole

table(clinical$drug)
table(clinical$target)
table(clinical$recist)
table(clinical$cancer_type)
table(clinical$study)

table(clinical$recist,clinical$cancer_type)


characteristic<-function(datasets){
  #clin<-subset(clinical,study=="TCGA_PAAD",drop=T)
  clin<-subset(clinical,study==datasets,drop=T)
  
  nr<-nrow(clin)
  print(nr)
  
  clin$agebin50<-0
  for(i in 1:nr){
    if(is.na(clin$age[i])==T){clin$agebin50[i]="NA"}
    else if(clin$age[i]<50){clin$agebin50[i]="<50"}
    else if(clin$age[i]>=50){clin$agebin50[i]=">=50"}
  }
  
  attach(clin)
  
  print("age")
  print(mean(clin$age,na.rm=T))
  print(sd(clin$age,na.rm=T))
  
  #print("agebin")
  #print(table(age_bin,useNA="always"))
  #table(age_bin,useNA="always")/nr
  
  print("agebin50")
  print(table(agebin50,useNA="always"))
  table(agebin50,useNA="always")/nr
  
  print("gender")
  print(table(gender,useNA="always"))
  table(gender,useNA="always")/nr
  
  print("os_event")
  print(table(os_event,useNA="always"))
  
  print("os")
  print(mean(clin$os,na.rm=T))
  print(sd(clin$os,na.rm=T))
  
  print("pfs_12_months")
  print(table(pfs_12_months,useNA="always"))
  
  print("target")
  print(table(target,useNA="always"))
  table(target,useNA="always")/nr
  
  print("response_recist")
  print(table(response_code,useNA="always"))
  table(response_code,useNA="always")/nr
  
  print("irAEs")
  print(table(irAEs,useNA="always"))
  table(irAEs,useNA="always")/nr
  
  print("irAEs_grade")
  print(table(irAEs_grade,useNA="always"))
  table(irAEs_grade,useNA="always")/nr
  
}

table(clinical$study)

characteristic("FrankelAE_2017")
characteristic("McCullochJA_2022")
characteristic("RoutyB_2018")
characteristic("SpencerCN_2021")
characteristic("DerosaL_2024")
characteristic("LeeKA_2022")
characteristic("PetersBA_2019")







