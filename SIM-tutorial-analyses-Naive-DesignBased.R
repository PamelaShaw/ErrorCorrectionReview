### This code implements the naive analysis (that ignores measurment error),
### IPW analysis, and generalized raking analysis (GR).
### GR is done two ways: one with naive influence functions and
### once with influence functions constructed with multiple imputation.
### Original code date 6/27/2026

rm(list=ls())

#version.string R version 4.4.2 (2024-10-31)


CIfunc<-function(beta,SE){ cbind(Est=beta, LL=beta-1.96*SE,UL=beta + 1.96*SE)}
statsfunc<-function(fit){ 
	beta<-summary(fit)$coef[,1]
	SE<- summary(fit)$coef[,2]
	pval<- summary(fit)$coef[,4]
	
return(cbind(Est=beta, SE=SE, LL=(beta-1.96*SE), UL=(beta + 1.96*SE), pvalue=pval))
}

statsfuncCox<-function(fit){ 
	beta<-summary(fit)$coef[,1]
	SE<- summary(fit)$coef[,3]
	pval<- summary(fit)$coef[,5]
	
return(cbind(Est=beta, SE=SE, LL=(beta-1.96*SE), UL=(beta + 1.96*SE), pvalue=pval))
}


library(survival) ### coxph
library(survey) ### calibrate, svy_glm()
library(lava)  #### influence function calculation for glm models
library(mice)  ### MI using mice
library(xtable)  ### Convert to latex table
library(Hmisc)


##### Code should work for everyone below here.

d<-read.csv("synthetic-mother-child-data-2025-03-18.csv") 

## No errors in mat_age_delivery.ph1, so this is now suppressed and we simply use 
## mat_age_delivery everywhere

d$mat_age_delivery<-d$mat_age_delivery.ph1
d$mat_age_delivery.ph1<-NULL


d$mat3race<-ifelse(d$mat4race=="Black" | d$mat4race=="White",d$mat4race,"Other")
d$mat3race<-relevel(factor(d$mat3race),ref="White")

d$mat3race.ph1<-ifelse(d$mat4race.ph1=="Black" | d$mat4race.ph1=="White",d$mat4race.ph1,"Other")
d$mat3race.ph1<-relevel(factor(d$mat3race.ph1),ref="White")


d$mat4race.ph1<-relevel(factor(d$mat4race.ph1),ref="White")
d$mat4race<-relevel(factor(d$mat4race),ref="White")

table(d$mat4race,d$mat3race)
table(d$mat4race.ph1,d$mat3race.ph1)

d$mat_diabetes3.ph1<-relevel(factor(d$mat_diabetes3.ph1),ref="None")
d$mat_diabetes3<-relevel(factor(d$mat_diabetes3),ref="None")

### get rid of slash in category
d$insurancer<-ifelse(d$insurancer=="Public/Other","PublicOther","Private")
d$insurancer.ph1<-ifelse(d$insurancer.ph1=="Public/Other","PublicOther","Private")

### Code as a factor
d$insurancer<-relevel(factor(d$insurancer),ref="Private")
d$insurancer.ph1<-relevel(factor(d$insurancer.ph1),ref="Private")

### for raking also important to code as factor
d$wave4Strata<-factor(d$wave4Strata)

#### make new variables to have regression coef as per 5 units of BMI and per 10 years of age
d$est_bmi5_preg_mother<- d$est_bmi_preg_mother/5
d$est_bmi5_preg_mother.ph1 <- d$est_bmi_preg_mother.ph1/5 
d$mat_age10_delivery<- d$mat_age_delivery/10

# maternal weight gain ~ maternal BMI + age + race + ethnicity + depression + insurance + smoking; 
## exclude those that are not singleton

d.singleton<- d[d$singleton.ph1==1 & (is.na(d$singleton) | d$singleton==1),]  ## Limiting weight gain analysis to mothers giving birth to only one child



#### IPW design and fit
IPWdesignC<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d.singleton,method="simple")

ipwfitC<-svyglm(estWtChangePerWk ~ est_bmi5_preg_mother + mat_age10_delivery + mat4race +
                     hispanic + depressionr + insurancer + tobacPregr, family=gaussian, design=IPWdesignC)

summary(ipwfitC)


### Raking with naive influence functions
naivefit<-glm(estWtChangePerWk.ph1 ~ est_bmi5_preg_mother.ph1  + mat_age10_delivery  + mat4race.ph1  +
                     hispanic.ph1  + depressionr.ph1  + insurancer.ph1  + tobacPregr.ph1 ,data=d.singleton,family=gaussian)
infMat<-IC(naivefit)
coefnum<-length(coef(naivefit))

#### compare naive with ipwfit
cbind(naive=coef(naivefit),ipw=coef(ipwfitC))

#### Default, generalized raking with all influence functions from regression
#### rakeformula = ~inf1+...+infk, where k = # of coef fit by regression
rakeformula<-"~inf1"
for(i in 1:coefnum){
	varname<-paste0("inf",i)
	d.singleton$inf <- infMat[,i]
	names(d.singleton)[names(d.singleton)=="inf"]<-varname
	if(i>1) rakeformula = paste0(rakeformula,"+",varname)
}
rakeformula = paste0(rakeformula,"+","wave4Strata")


### Update design since created new variables
IPWdesignC<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d.singleton,method="simple")
rakeDesignC<-calibrate(IPWdesignC,formula=formula(rakeformula),phase=2,calfun="raking")

rakefitC<-svyglm(estWtChangePerWk ~ est_bmi5_preg_mother + mat_age10_delivery + mat4race +
                     hispanic + depressionr + insurancer + tobacPregr, family=gaussian, design=rakeDesignC)
                     
## Compare fit
#### Note can see a less similar to real data that race in naive was incorrectly showing protecting against weight gain
#### Nice to see good gain in efficeincy by raking.

### naive
Nans<-statsfunc(naivefit)
row.names(Nans) <- c("Intercept","BMI5","Age10","Asian","Black","Other","Hispanic","Depression","Public Insurance","Smoking")
round(Nans,3)
#xtable(round(Nans,3),digits=3)

#### Naive IPW with incorrect SE
##IPW
IPWans<-statsfunc(ipwfitC)
row.names(IPWans) <- c("Intercept","BMI5","Age10","Asian","Black","Other","Hispanic","Depression","Public Insurance","Smoking")
round(IPWans,3)
#xtable(round(IPWans,3),digits=3)



## Raking
GRans<-statsfunc(rakefitC)
row.names(GRans) <- c("Intercept","BMI5","Age10","Asian","Black","Other","Hispanic","Depression","Public Insurance","Smoking")
round(GRans,3)
#xtable(round(GRans,3),digits=3)

#### Alternate Raking algorithm for Maternal weight
#### Maybe more efficient: use multiple imputation to generate raking variable
#define a dataset with formula covariates and the outcome variable and impute missing data
set.seed(3243)
formula.out<-formula(ipwfitC)
#data.mi<-model.frame(formula.out,data=d.singleton,na.action=NULL)
data.mi<-d.singleton[,c("estWtChangePerWk.ph1","est_bmi5_preg_mother.ph1",
                     "mat4race.ph1","hispanic.ph1","depressionr.ph1","insurancer.ph1","tobacPregr.ph1",
                     "estWtChangePerWk","est_bmi5_preg_mother","mat_age10_delivery",
                     "mat4race","hispanic","depressionr","insurancer","tobacPregr")]


phase_2_indx <- (d.singleton$R == 1)
miss_phase2 <- data.mi[d.singleton$R == 1, ]
miss_cols <- which(colSums(is.na(data.mi)) > 0)
miss_phase2[, miss_cols] <- NA
data.mi2 <- rbind(data.mi, miss_phase2)
fake_phase_2 <- ((1:nrow(data.mi2)) %in% (nrow(data.mi) + 1):nrow(data.mi2))
data.mi <- data.mi2

# Create NimpRaking datasets where W are imputed for all N individuals.
NimpRaking<-10
init <- mice::mice(data.mi, maxit = 0)
pred.matrix <- init$predictorMatrix
data.imputed <- mice::mice(data.mi, predictorMatrix = pred.matrix, m = NimpRaking, maxit = 20, print = FALSE)

# Estimate value of influence functions by fitting a logistic regression model within each imputed dataset, cal
#culating the resulting influence functions, and then averaging values of influence functions across imputed datas
#ets.
### now have all variables so set coefnum
coefnum<-length(coef(ipwfitC))

infMat_all <- array(data = 0, dim = c(nrow(d.singleton), coefnum, NimpRaking))
for (iter in 1:NimpRaking) {
# Limiting the dataset to imputed data, i.e. removing rows where W was originally observed.
  
# explore correlation with imputed and observed data
	imp_init <- mice::complete(data.imputed, iter)
	impData_i <- imp_init[1:nrow(d.singleton), ]
	
	impData_i[phase_2_indx, ] <- imp_init[fake_phase_2, ]
	mifit <- glm(formula.out, family = gaussian, data = impData_i)
	infMat_all[, , iter] <- IC(mifit)
}
infMat <- rowMeans(infMat_all, dims = 2) ### average each element across the imputations
#### Default, generalized raking with all influence functions from regression
#### rakeformula = ~inf1+...+infk, where k = # of coef fit by regression
rakeformula<-"~inf1"
for(i in 1:coefnum){
	varname<-paste0("inf",i)
	d.singleton$inf <- infMat[,i]
	names(d.singleton)[names(d.singleton)=="inf"]<-varname
	if(i>1) rakeformula = paste0(rakeformula,"+",varname)
}
rakeformula = paste0(rakeformula,"+","wave4Strata")

IPWdesignC<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d.singleton,method="simple")
rakeDesignCv2<-calibrate(IPWdesignC,formula=as.formula(rakeformula),phase=2,calfun="raking")

rakefitCv2<-svyglm(formula.out, family=gaussian, design=rakeDesignCv2)
                     
## Raking
round(statsfunc(rakefitC),8)
## Raking v2
round(statsfunc(rakefitCv2),8)
                     
############### ############### ############### ############### ############### ############### 
############### ############### Binary Outcome ############### ############### ###############
###############     ############### ############### ############### ############### ###############   


## Child asthma status ~ everything in original manuscript, with alteration to drop a few covariates due to the fewer degrees of freedom
### dropping  hispanic + mat_diabetes3 + cesarean and changing to a 3 category race

d.A<-d[d$in.Aframe==1,]  ## Limiting analysis to eligible children

dim(d.A)

mean(d.A$asthma.ph1)
table(d.A$asthma[!is.na(d.A$asthma)]) ### 107 phase 2 cases
table(validated=d.A$asthma[!is.na(d.A$asthma)],naive=d.A$asthma.ph1[!is.na(d.A$asthma)]) ### naive vs true


IPWdesignB<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d.A,method="simple")
ipwfitB<-svyglm(asthma ~ estWtChangePerWk + est_bmi5_preg_mother + mat_age10_delivery + 
                      mat3race  + male + insurancer + tobacPregr + mat_asthma + egaWk, family=quasibinomial, design=IPWdesignB)

### Raking with naive influence functions
naivefit<-glm(asthma.ph1 ~ estWtChangePerWk.ph1 + est_bmi5_preg_mother.ph1  + mat_age10_delivery  + 
                      mat3race.ph1 + male.ph1  + insurancer.ph1  + tobacPregr.ph1  + mat_asthma.ph1,
                      data=d.A,family=binomial)
infMat<-IC(naivefit)
coefnum<-length(coef(naivefit))
#### Default, generalized raking with all influence functions from regression
#### rakeformula = ~inf1+...+infk, where k = # of coef fit by regression
rakeformula<-"~inf1"
for(i in 1:coefnum){
	varname<-paste0("inf",i)
	d.A$inf <- infMat[,i]
	names(d.A)[names(d.A)=="inf"]<-varname
	if(i>1) rakeformula = paste0(rakeformula,"+",varname)
}
rakeformula = paste0(rakeformula,"+","wave4Strata")

IPWdesignB<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d.A,method="simple")
rakeDesignB<-calibrate(IPWdesignB,formula=as.formula(rakeformula),phase=2,calfun="raking")

rakefitB<-svyglm(asthma ~ estWtChangePerWk + est_bmi5_preg_mother + mat_age10_delivery + 
                      mat3race + male + insurancer + tobacPregr + mat_asthma + egaWk, 
                      family=quasibinomial, design=rakeDesignB)
                     
## Compare fit
### naive fit
round(statsfunc(naivefit),3)
##IPW
round(statsfunc(ipwfitB),3)
## Raking
round(statsfunc(rakefitB),3)

#### Alternate Raking algorithm for asthma
#### Since missing one of the key variables, use multiple imputation to generate raking variable
#define a dataset with formula covariates and the outcome variable and impute missing data
set.seed(3243)
formula.out<-formula("asthma ~ estWtChangePerWk + est_bmi5_preg_mother + mat_age10_delivery + 
                      mat3race + male + insurancer + tobacPregr + mat_asthma + egaWk")
#data.mi<-model.frame(formula.out,data=d.A,na.action=NULL)
data.mi<-d.A[,c("estWtChangePerWk.ph1","est_bmi5_preg_mother.ph1",
               "mat3race.ph1","male.ph1","mat_asthma.ph1","insurancer.ph1","tobacPregr.ph1",
               "asthma.ph1",
               "estWtChangePerWk","est_bmi5_preg_mother","mat_age10_delivery",
               "mat3race","male","mat_asthma","insurancer","tobacPregr","egaWk",
               "asthma")]
phase_2_indx <- (d.A$R == 1)
miss_phase2 <- data.mi[d.A$R == 1, ]
miss_cols <- which(colSums(is.na(data.mi)) > 0)
miss_phase2[, miss_cols] <- NA
data.mi2 <- rbind(data.mi, miss_phase2)
fake_phase_2 <- ((1:nrow(data.mi2)) %in% (nrow(data.mi) + 1):nrow(data.mi2))
data.mi <- data.mi2

# Create NimpRaking datasets where W are imputed for all N individuals.
NimpRaking<-10
init <- mice::mice(data.mi, maxit = 0)
pred.matrix <- init$predictorMatrix
data.imputed <- mice::mice(data.mi, predictorMatrix = pred.matrix, m = NimpRaking, maxit = 20, print = FALSE)


# Estimate value of influence functions by fitting a logistic regression model within each imputed dataset, cal
#culating the resulting influence functions, and then averaging values of influence functions across imputed datas
#ets.
### now have all variables so set coefnum
coefnum<-length(coef(ipwfitB))

infMat_all <- array(data = 0, dim = c(nrow(d.A), coefnum, NimpRaking))
for (iter in 1:NimpRaking) {
# Limiting the dataset to imputed data, i.e. removing rows where W was originally observed.
	imp_init <- mice::complete(data.imputed, iter)
	impData_i <- imp_init[1:nrow(d.A), ]
	impData_i[phase_2_indx, ] <- imp_init[fake_phase_2, ]
	mifit <- glm(formula.out, family = binomial, data = impData_i)
	infMat_all[, , iter] <- IC(mifit)
}
infMat <- rowMeans(infMat_all, dims = 2) ### average each element across the imputations
#### Default, generalized raking with all influence functions from regression
#### rakeformula = ~inf1+...+infk, where k = # of coef fit by regression
rakeformula<-"~inf1"
for(i in 1:coefnum){
	varname<-paste0("inf",i)
	d.A$inf <- infMat[,i]
	names(d.A)[names(d.A)=="inf"]<-varname
	if(i>1) rakeformula = paste0(rakeformula,"+",varname)
}
rakeformula = paste0(rakeformula,"+","wave4Strata")

IPWdesignB<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d.A,method="simple")
rakeDesignBv2<-calibrate(IPWdesignB,formula=as.formula(rakeformula),phase=2,calfun="raking")

rakefitBv2<-svyglm(asthma ~ estWtChangePerWk + est_bmi5_preg_mother + mat_age10_delivery + 
                      mat3race + male + insurancer + tobacPregr + mat_asthma + egaWk, 
                      family=quasibinomial, design=rakeDesignBv2)
                     
## Compare fit
### naive fit
round(statsfunc(naivefit),3)
##IPW
round(statsfunc(ipwfitB),3)
## Raking
round(statsfunc(rakefitB),3)
## Raking v2
round(statsfunc(rakefitBv2),3)


############### ############### ############### ############### ############### ############### 
############### ############### Survival Outcome ############### ############### ###############
###############  ############### ############### ############### ############### ###############   

#### there are no inclusion restrictions

	## Childhood obesity ~ everything in original manuscript
	IPWdesignS<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d,method="simple")
	ipwfitS<-svycoxph(Surv(ttobesity,obesity) ~ estWtChangePerWk + est_bmi5_preg_mother + 
	               mat_age10_delivery + mat4race + hispanic + mat_diabetes3 + cesarean + 
	               male + depressionr + insurancer + singleton + tobacPregr + married + 
	               number_children + egaWk, design=IPWdesignS)
	summary(ipwfitS)$coef
	
	### Raking with naive influence functions
	naivefit<- coxph(Surv(ttobesity.ph1,obesity.ph1) ~ estWtChangePerWk.ph1 + est_bmi5_preg_mother.ph1 + 
	                      mat_age10_delivery + mat4race.ph1 + hispanic.ph1 + mat_diabetes3.ph1 + cesarean.ph1 + 
	                      male.ph1 + depressionr.ph1 + insurancer.ph1 + singleton.ph1 + tobacPregr.ph1, data=d)
	summary(naivefit)$coef
	
	infMat<-resid(naivefit,type="dfbeta",weighted=FALSE)
	coefnum<-length(coef(naivefit))
	#### Default, generalized raking with all influence functions from regression
	#### rakeformula = ~inf1+...+infk, where k = # of coef fit by regression
	rakeformula<-"~inf1"
	for(i in 1:coefnum){
		varname<-paste0("inf",i)
		d$inf <- infMat[,i]
		names(d)[names(d)=="inf"]<-varname
		if(i>1) rakeformula = paste0(rakeformula,"+",varname)
	}
	rakeformula = paste0(rakeformula,"+","wave4Strata")
	
	IPWdesignS<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d,method="simple")
	rakeDesignS<-calibrate(IPWdesignS,formula=formula(rakeformula),phase=2,calfun="raking")
	formula.out<-formula(ipwfitS)
	rakefitS<-svycoxph(formula.out, design=rakeDesignS)
	                     
	## Compare fit
	
	##Naive fit
	round(statsfuncCox(naivefit),3)
	##IPW
	round(statsfuncCox(ipwfitS),3)
	## Raking
	round(statsfuncCox(rakefitS),3)
	
	
### Not much improvement in SE over IPW. Can check correlation of esitmaeted vs true influence function.
trueinf<-resid(rakefitS,type="dfbeta",weighted=FALSE)
naiveinf<-infMat[d$R==1,]
for(i in 1:ncol(trueinf)) print(cor(trueinf[,i],naiveinf[,i]))

#### Try fewer raking variables, just primary exposure of weight gain as in Shepherd et al 2023
#### SE did reduce, especially for the maternal weight gain coeficient. As expected, no worse than IPW
	rakeDesignS.simple1<-calibrate(IPWdesignS,~inf1+wave4Strata,phase=2,calfun="raking")
	formula.out<-formula(ipwfitS)
	rakefitS.alt1<-svycoxph(formula.out, design=rakeDesignS.simple1)
		round(statsfuncCox(rakefitS.alt1),3)

rakeDesignS.simple2<-calibrate(IPWdesignS,~inf1,phase=2,calfun="raking")
	formula.out<-formula(ipwfitS)
	rakefitS.alt2<-svycoxph(formula.out, design=rakeDesignS.simple2)
		round(statsfuncCox(rakefitS.alt2),3)



#### Alternate Raking algorithm for Survival
#### Maybe more efficient: use multiple imputation to generate raking variable
#define a dataset with formula covariates and the outcome variable and impute missing data

set.seed(3243)
data.mi<-d[,c("estWtChangePerWk.ph1","est_bmi5_preg_mother.ph1",
               "mat4race.ph1","hispanic.ph1","mat_diabetes3.ph1","cesarean.ph1",
               "male.ph1","depressionr.ph1","insurancer.ph1","singleton.ph1",
               "tobacPregr.ph1",
               "ttobesity.ph1","obesity.ph1",
               "estWtChangePerWk","est_bmi5_preg_mother","mat_age10_delivery",
               "mat4race","hispanic","mat_diabetes3","cesarean",
               "male","depressionr","insurancer","singleton",
               "tobacPregr","married","number_children","egaWk",
               "ttobesity","obesity")]

#data.mi$interaction<-data.mi$ttobesity*data.mi$obesity
#data.mi$interaction.ph1<-data.mi$ttobesity.ph1*data.mi$obesity.ph1

phase_2_indx <- (d$R == 1)
miss_phase2 <- data.mi[d$R == 1, ]
miss_cols <- which(colSums(is.na(data.mi)) > 0)
miss_phase2[, miss_cols] <- NA
data.mi2 <- rbind(data.mi, miss_phase2)
fake_phase_2 <- ((1:nrow(data.mi2)) %in% (nrow(data.mi) + 1):nrow(data.mi2))
data.mi <- data.mi2

# Create NimpRaking datasets where W are imputed for all N individuals.
NimpRaking<-10
init <- mice::mice(data.mi, maxit = 0)
pred.matrix <- init$predictorMatrix
data.imputed <- mice::mice(data.mi, predictorMatrix = pred.matrix, m = NimpRaking, maxit = 20, print = FALSE)


# Estimate value of influence functions by fitting a logistic regression model within each imputed dataset, cal
#culating the resulting influence functions, and then averaging values of influence functions across imputed datas
#ets.
### now have all variables so set coefnum
coefnum<-length(coef(ipwfitS))

infMat_all <- array(data = 0, dim = c(nrow(d), coefnum, NimpRaking))
for (iter in 1:NimpRaking) {
# Limiting the dataset to imputed data, i.e. removing rows where W was originally observed.
	imp_init <- mice::complete(data.imputed, iter)
	impData_i <- imp_init[1:nrow(d), ]
	impData_i[phase_2_indx, ] <- imp_init[fake_phase_2, ]
	mifit <- coxph(formula.out, data = impData_i)
	infMat_all[, , iter] <- resid(mifit,type="dfbeta",weighted=FALSE)

}
infMat <- rowMeans(infMat_all, dims = 2) ### average each element across the imputations
#### Default, generalized raking with all influence functions from regression
#### rakeformula = ~inf1+...+infk, where k = # of coef fit by regression
rakeformula<-"~inf1"
for(i in 1:coefnum){
	varname<-paste0("inf",i)
	d$inf <- infMat[,i]
	names(d)[names(d)=="inf"]<-varname
	if(i>1) rakeformula = paste0(rakeformula,"+",varname)
}
rakeformula = paste0(rakeformula,"+","wave4Strata")

IPWdesignS<-twophase(id=list(~1,~1),strata=list(NULL,~wave4Strata),subset=~I(R==1),data=d,method="simple")
rakeDesignSv2<-calibrate(IPWdesignS,formula=as.formula(rakeformula),phase=2,calfun="raking")

rakefitSv2<-svycoxph(formula.out, design=rakeDesignSv2)
                     
## Raking
round(statsfuncCox(rakefitS),4)
## Raking v2
round(statsfuncCox(rakefitSv2),4)
                     


              
