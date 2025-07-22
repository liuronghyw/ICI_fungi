rm(list=ls())
setwd("D:/课题/细菌—真菌互作/3_R/")

library(tidyverse)
library(survival)

###Melanoma#####
cli_mela<-read.csv("2.rawData/clinical/final/filter/cli_mela_filter.csv",row.names = 1, fileEncoding = "UTF-8")
cli_mela$gender<-gsub("Male",1,cli_mela$gender)
cli_mela$gender<-gsub("Female",0,cli_mela$gender)
cli_mela$BMI_group<-gsub("lean",0,cli_mela$BMI_group)
cli_mela$BMI_group<-gsub("overweight",1,cli_mela$BMI_group)
cli_mela$BMI_group<-gsub("Obese",2,cli_mela$BMI_group)
cli_mela$age_group<-gsub("<60",0,cli_mela$age_group)
cli_mela$age_group<-gsub("≥60",1,cli_mela$age_group)
cli_mela$BMI_group<-as.numeric(cli_mela$BMI_group)

##gender
####response
cli_mela_gender<-subset(cli_mela,!is.na(gender)&!is.na(response_code))
a<-glm(response_code~gender, family=binomial(logit), data=cli_mela_gender)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS6
cli_mela_gender<-subset(cli_mela,!is.na(gender)&!is.na(pfs_6_months))
a<-glm(pfs_6_months~gender, family=binomial(logit), data=cli_mela_gender)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS12
cli_mela_gender<-subset(cli_mela,!is.na(gender)&!is.na(pfs_12_months))
a<-glm(pfs_12_months~gender, family=binomial(logit), data=cli_mela_gender)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4


####PFS
cli_mela_gender<-subset(cli_mela,!is.na(gender)&!is.na(pfs_event))
a<-summary(coxph(Surv(pfs,pfs_event)~gender,cli_mela_gender,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

####OS
cli_mela_gender<-subset(cli_mela,!is.na(gender)&!is.na(os_event))
a<-summary(coxph(Surv(os,os_event)~gender,cli_mela_gender,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

##BMI
##response
cli_mela_BMI<-subset(cli_mela,!is.na(BMI_group)&!is.na(response_code))
a<-glm(response_code~BMI_group, family=binomial(logit), data=cli_mela_BMI)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS6
cli_mela_BMI<-subset(cli_mela,!is.na(BMI_group)&!is.na(pfs_6_months))
a<-glm(pfs_6_months~BMI_group, family=binomial(logit), data=cli_mela_BMI)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS12
cli_mela_BMI<-subset(cli_mela,!is.na(BMI_group)&!is.na(pfs_12_months))
a<-glm(pfs_12_months~BMI_group, family=binomial(logit), data=cli_mela_BMI)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS
cli_mela_BMI<-subset(cli_mela,!is.na(BMI_group)&!is.na(pfs_event))
a<-summary(coxph(Surv(pfs,pfs_event)~BMI_group,cli_mela_BMI,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

####OS
cli_mela_BMI<-subset(cli_mela,!is.na(BMI_group)&!is.na(os_event))
a<-summary(coxph(Surv(os,os_event)~BMI_group,cli_mela_BMI,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3


##AGE
##response
cli_mela_age<-subset(cli_mela,!is.na(age_group)&!is.na(response_code))
a<-glm(response_code~age_group, family=binomial(logit), data=cli_mela_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS6
cli_mela_age<-subset(cli_mela,!is.na(age_group)&!is.na(pfs_6_months))
a<-glm(pfs_6_months~age_group, family=binomial(logit), data=cli_mela_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS12
cli_mela_age<-subset(cli_mela,!is.na(age_group)&!is.na(pfs_12_months))
a<-glm(pfs_12_months~age_group, family=binomial(logit), data=cli_mela_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4


####PFS
cli_mela_age<-subset(cli_mela,!is.na(age_group)&!is.na(pfs_event))
a<-summary(coxph(Surv(pfs,pfs_event)~age_group,cli_mela_age,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

####OS
cli_mela_age<-subset(cli_mela,!is.na(age_group)&!is.na(os_event))
a<-summary(coxph(Surv(os,os_event)~age_group,cli_mela_age,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3


##NSCLC####
cli_NSCLC<-read.csv("2.rawData/clinical/final/filter/cli_NSCLC_filter.csv",row.names = 1)
cli_NSCLC$gender<-gsub("Male",1,cli_NSCLC$gender)
cli_NSCLC$gender<-gsub("Female",0,cli_NSCLC$gender)
cli_NSCLC$BMI_group<-gsub("lean",0,cli_NSCLC$BMI_group)
cli_NSCLC$BMI_group<-gsub("overweight",1,cli_NSCLC$BMI_group)
cli_NSCLC$BMI_group<-gsub("Obese",2,cli_NSCLC$BMI_group)
cli_NSCLC$age_group<-gsub("<60",0,cli_NSCLC$age_group)
cli_NSCLC$age_group<-gsub("≥60",1,cli_NSCLC$age_group)
cli_NSCLC$BMI_group<-as.numeric(cli_NSCLC$BMI_group)

##gender
####response
cli_NSCLC_gender<-subset(cli_NSCLC,!is.na(gender)&!is.na(response_code))
a<-glm(response_code~gender, family=binomial(logit), data=cli_NSCLC_gender)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS6
cli_NSCLC_gender<-subset(cli_NSCLC,!is.na(gender)&!is.na(pfs_6_months))
a<-glm(pfs_6_months~gender, family=binomial(logit), data=cli_NSCLC_gender)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS12
cli_NSCLC_gender<-subset(cli_NSCLC,!is.na(gender)&!is.na(pfs_12_months))
a<-glm(pfs_12_months~gender, family=binomial(logit), data=cli_NSCLC_gender)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4



####PFS
cli_NSCLC_gender<-subset(cli_NSCLC,!is.na(gender)&!is.na(pfs_event))
a<-summary(coxph(Surv(pfs,pfs_event)~gender,cli_NSCLC_gender,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

####OS
cli_NSCLC_gender<-subset(cli_NSCLC,!is.na(gender)&!is.na(os_event))
a<-summary(coxph(Surv(os,os_event)~gender,cli_NSCLC_gender,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

##BMI
##response
cli_NSCLC_BMI<-subset(cli_NSCLC,!is.na(BMI_group)&!is.na(response_code))
a<-glm(response_code~BMI_group, family=binomial(logit), data=cli_NSCLC_BMI)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS6
cli_NSCLC_BMI<-subset(cli_NSCLC,!is.na(BMI_group)&!is.na(pfs_6_months))
a<-glm(pfs_6_months~BMI_group, family=binomial(logit), data=cli_NSCLC_BMI)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS12
cli_NSCLC_BMI<-subset(cli_NSCLC,!is.na(BMI_group)&!is.na(pfs_12_months))
a<-glm(pfs_12_months~BMI_group, family=binomial(logit), data=cli_NSCLC_BMI)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4


####PFS
cli_NSCLC_BMI<-subset(cli_NSCLC,!is.na(BMI_group)&!is.na(pfs_event))
a<-summary(coxph(Surv(pfs,pfs_event)~BMI_group,cli_NSCLC_BMI,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

####OS
cli_NSCLC_BMI<-subset(cli_NSCLC,!is.na(BMI_group)&!is.na(os_event))
a<-summary(coxph(Surv(os,os_event)~BMI_group,cli_NSCLC_BMI,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3


##AGE
##response
cli_NSCLC_age<-subset(cli_NSCLC,!is.na(age_group)&!is.na(response_code))
a<-glm(response_code~age_group, family=binomial(logit), data=cli_NSCLC_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS6
cli_NSCLC_age<-subset(cli_NSCLC,!is.na(age_group)&!is.na(pfs_6_months))
a<-glm(pfs_6_months~age_group, family=binomial(logit), data=cli_NSCLC_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS12
cli_NSCLC_age<-subset(cli_NSCLC,!is.na(age_group)&!is.na(pfs_12_months))
a<-glm(pfs_12_months~age_group, family=binomial(logit), data=cli_NSCLC_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4



####PFS
cli_NSCLC_age<-subset(cli_NSCLC,!is.na(age_group)&!is.na(pfs_event))
a<-summary(coxph(Surv(pfs,pfs_event)~age_group,cli_NSCLC_age,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

####OS
cli_NSCLC_age<-subset(cli_NSCLC,!is.na(age_group)&!is.na(os_event))
a<-summary(coxph(Surv(os,os_event)~age_group,cli_NSCLC_age,na.action=na.exclude))
a<-as.vector(unlist(a))
pvalue<-a$coefficients5
HR<-a$coefficients2
HR_left<-a$conf.int3
HR_right<-a$conf.int4
se<-a$coefficients3

##RCC####
cli_RCC<-read.csv("2.rawData/clinical/final/filter/cli_RCC_filter.csv",row.names = 1)
table(cli_RCC$gender)
table(cli_RCC$age_group)
table(cli_RCC$BMI_group)
cli_RCC$gender<-gsub("Male",1,cli_RCC$gender)
cli_RCC$gender<-gsub("Female",0,cli_RCC$gender)
cli_RCC$age_group<-gsub("<60",0,cli_RCC$age_group)
cli_RCC$age_group<-gsub("≥60",1,cli_RCC$age_group)


##gender
####response
cli_RCC_gender<-subset(cli_RCC,!is.na(gender)&!is.na(response_code))
a<-glm(response_code~gender, family=binomial(logit), data=cli_RCC_gender)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS6
cli_RCC_gender<-subset(cli_RCC,!is.na(gender)&!is.na(pfs_6_months))
a<-glm(pfs_6_months~gender, family=binomial(logit), data=cli_RCC_gender)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4


##AGE
##response
cli_RCC_age<-subset(cli_RCC,!is.na(age_group)&!is.na(response_code))
a<-glm(response_code~age_group, family=binomial(logit), data=cli_RCC_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS6
cli_RCC_age<-subset(cli_RCC,!is.na(age_group)&!is.na(pfs_6_months))
a<-glm(pfs_6_months~age_group, family=binomial(logit), data=cli_RCC_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4

####PFS12
cli_RCC_age<-subset(cli_RCC,!is.na(age_group)&!is.na(pfs_12_months))
a<-glm(pfs_12_months~age_group, family=binomial(logit), data=cli_RCC_age)
ci<-confint(a)
ci<-unlist(ci)

b<-unlist(summary(a))
pvalue<-b$coefficients8
or<-exp(b$coefficients2)

orl<-exp(ci[2])
orr<-exp(ci[4])
se<-b$coefficients4


