
########################################################################################
####  These analyses perform multiple imputation (MI) to account for errors in variables
####   in the synthetic mother-child weight study.
########################################################################################


rm(list=ls())

#############################################
####      Reading in synthetic data
setwd("~/Dropbox/pam/Shaw Shepherd Research/synthetic-data-PCORI")
d<-read.csv("synthetic-mother-child-data-2025-03-18.csv")

head(d)

#############################################
####      Doing a little data management
d$mat3race.ph1<-ifelse(d$mat4race.ph1=="Asian","Other",d$mat4race.ph1)
d$mat3race<-ifelse(d$mat4race=="Asian","Other",d$mat4race)

d$mat3race.ph1<-relevel(factor(d$mat3race.ph1),ref="White")
d$mat3race<-relevel(factor(d$mat3race),ref="White")

d$mat4race.ph1<-relevel(factor(d$mat4race.ph1),ref="White")
d$mat4race<-relevel(factor(d$mat4race),ref="White")


d$mat_diabetes3.ph1<-relevel(factor(d$mat_diabetes3.ph1),ref="None")
d$mat_diabetes3<-relevel(factor(d$mat_diabetes3),ref="None")

## No errors in mat_age_delivery.ph1, so this is now suppressed and we simply use 
## mat_age_delivery everywhere
d$mat_age_delivery<-d$mat_age_delivery.ph1
d$mat_age_delivery.ph1<-NULL

## Making binary variables factor variables, as required by mice
d$hispanic<-as.factor(d$hispanic)
d$hispanic.ph1<-as.factor(d$hispanic.ph1)
d$depressionr<-as.factor(d$depressionr)
d$depressionr.ph1<-as.factor(d$depressionr.ph1)
d$insurancer<-as.factor(d$insurancer)
d$insurancer.ph1<-as.factor(d$insurancer.ph1)
d$tobacPregr<-as.factor(d$tobacPregr)
d$tobacPregr.ph1<-as.factor(d$tobacPregr.ph1)
d$mat_asthma.ph1<-as.factor(d$mat_asthma.ph1)
d$mat_asthma<-as.factor(d$mat_asthma)
d$asthma.ph1<-as.factor(d$asthma.ph1)
d$asthma<-as.factor(d$asthma)
d$cesarean.ph1<-as.factor(d$cesarean.ph1)
d$cesarean<-as.factor(d$cesarean)
d$obesity.ph1<-as.factor(d$obesity.ph1)
d$obesity<-as.factor(d$obesity)
d$married<-as.factor(d$married)


#####################################################
####  Maternal weight gain analysis - linear model
####  Our first model is the following:
####      maternal weight gain ~ maternal BMI + age + race + ethnicity + 
####                             depression + insurance + smoking 
####  And we want to exclude those that are not singleton

library(mice)

####  Limiting weight gain analysis to mothers giving birth to only one child
d.singleton<-d[d$singleton.ph1==1 & (is.na(d$singleton) | d$singleton==1),]
####  Creating a reduced dataset with only the necessary variables
d.mi<-d.singleton[,c("estWtChangePerWk","est_bmi_preg_mother","mat_age_delivery",
                     "mat4race","hispanic","depressionr","insurancer","tobacPregr",
                     "estWtChangePerWk.ph1","est_bmi_preg_mother.ph1",
                     "mat4race.ph1","hispanic.ph1","depressionr.ph1","insurancer.ph1",
                     "tobacPregr.ph1")]

####  These next few lines of code are to make sure that mice is properly including variables
####   and is using reasonable imputation models
init<-mice(d.mi, maxit=0)
pred.matrix <- init$predictorMatrix
pred.matrix

my.methods<-init$method
my.methods

####  Multiply imputing missing phase-2 data using mice
set.seed(1)
mice_data <- mice(d.mi, predictorMatrix=pred.matrix, method=my.methods, 
                  maxit=20, m=25, print=FALSE)
mice_data$loggedEvents

####  Fitting our model of interest using the the multiply imputed datasets
mod.mi <- with(mice_data, 
               lm(estWtChangePerWk ~ est_bmi_preg_mother + mat_age_delivery + mat4race +
                     hispanic + depressionr + insurancer + tobacPregr))

####  Pooling estimates across imputations using Rubin's rules
pooled.est<-pool(mod.mi)
mi.results<-summary(pooled.est)
mi.results


####  Displaying results and 95% confidence intervals, for MI estimates in Table 1 of manuscript 
ests.mi1<-cbind(mi.results$estimate,
                mi.results$std.error,
                mi.results$estimate-1.96*mi.results$std.error,
                mi.results$estimate+1.96*mi.results$std.error)
ests.mi<-ests.mi1
ests.mi[2,]<-ests.mi1[2,]*5 # changing BMI coefficient estimate to be per 5 kg/m2
ests.mi[3,]<-ests.mi1[3,]*10 # changing age coefficient estimate to be per 10 years
round(ests.mi,3)



####  Because of possible concerns over uncongeniality of imputation and analysis models,
####   we are bootstrapping the multiple imputation process. To do this, we use the 
####   package bootImpute.

library(bootImpute)

#### With 200 bootstrap replications and 2 imputations within each, this takes 7-8 minutes
####  to run on my Macbook Air.
imps<-bootMice(d.mi, nBoot=200, nImp=2, seed=123)
analysisImp <- function(inputData) {
   coef(lm(estWtChangePerWk ~ est_bmi_preg_mother + mat_age_delivery + mat4race +
     hispanic + depressionr + insurancer + tobacPregr, data=inputData))
}
ans<-bootImputeAnalyse(imps=imps, analysisfun=analysisImp)
ans

####  Showing results of Bootstrap-MI. Very similar to those using Rubin's rules.
####   These are presented in Supplementary Table 1.
ests.boot<-cbind(ans$ests,sqrt(ans$var),ans$ci)
ests.boot[2,]<-ests.boot[2,]*5  ## changing BMI coefficient to per 5 kg/m2
ests.boot[3,]<-ests.boot[3,]*10 ## changing age coefficient to per 10 years
round(ests.boot,3)



####  Calculating the naive estimate that only analyzes the error prone data
mod.naive<- lm(estWtChangePerWk.ph1 ~ est_bmi_preg_mother.ph1 + mat_age_delivery + mat4race.ph1 +
                 hispanic.ph1 + depressionr.ph1 + insurancer.ph1 + tobacPregr.ph1, data=d.singleton)
summary(mod.naive)
####  Displaying results from Naive analysis that are presented in Table 1 of manuscript
ests.naive1<-cbind(mod.naive$coefficients,
                   summary(mod.naive)$coefficients[,2],
                   mod.naive$coefficients-1.96*summary(mod.naive)$coefficients[,2],
                   mod.naive$coefficients+1.96*summary(mod.naive)$coefficients[,2])
ests.naive<-ests.naive1
ests.naive[2,]<-ests.naive1[2,]*5 # changing BMI coefficient to per 5 kg/m2
ests.naive[3,]<-ests.naive1[3,]*10 # changing age coefficient to per 10 years
round(ests.naive, 3)
                  

####  Figure showing naive, MI, and boot MI results.
n.covs<-length(mod.naive$coefficients)-1
gap<-0.03
factor<-10
plot(c(-4,2),c(0,1),type="n",xlab="Effect (Change in kg/wk *10)",ylab="", axes=FALSE)
axis(1)
abline(v=0,col=gray(.9))
points(ests.naive[-1,1]*factor,c(1:n.covs)/n.covs, pch=3)
for (i in 1:n.covs){
  lines(c(ests.naive[i+1,3],ests.naive[i+1,4])*factor,rep(i/n.covs,each=2))
}
points(ests.mi[-1,1]*factor,c(1:n.covs)/n.covs-gap, col=3, pch=3)
for (i in 1:n.covs){
  lines(c(ests.mi[i+1,3],ests.mi[i+1,4])*factor,rep(i/n.covs,each=2)-gap,col=3)
}
points(ests.boot[-1,1]*factor,c(1:n.covs)/n.covs-2*gap, col=2, pch=3)
for (i in 1:n.covs){
  lines(c(ests.boot[i+1,3],ests.boot[i+1,4])*factor,rep(i/n.covs-2*gap,each=2), col=2)
}
text(rep(-3.5,n.covs),c(1:n.covs)/n.covs,
      c("BMI (per 5 kg/m2)","Age (per 10 years)","Black","Asian","Other","Hispanic","Depression","No Insurance",
        "Smoking"), pos=4)
legend(x="bottomright", legend=c("naive","MI","Boot-MI"), lty=c(1,1,1), col=c(1,3,2), bty="n")






#########################################################################
####  Childhood asthma analyses - Binary outcome
####  Our model is a logistic regression model of the following form:
####      asthma ~ estWtChangePerWk + est_bmi_preg_mother + mat_age_delivery + 
####                      mat3race + maleA + insurancer + tobacPregr + mat_asthma + egaWk
####  And we are going to limit the analyses to eligible children.

####  Limiting analysis to eligible children
d.A<-d[d$in.Aframe==1,]  

dim(d.A)
summary(d.A$asthma.ph1)
table(d.A$asthma[!is.na(d.A$asthma)]) ### 107 phase 2 cases
table(validated=d.A$asthma[!is.na(d.A$asthma)],naive=d.A$asthma.ph1[!is.na(d.A$asthma)]) ### naive vs true

####  Simplifying the dataset for multiple imputation to only include relevant variables
d.mi.A<-d.A[,c("estWtChangePerWk","est_bmi_preg_mother","mat_age_delivery",
                     "mat3race","male","mat_asthma","insurancer","tobacPregr","egaWk",
                     "asthma","estWtChangePerWk.ph1","est_bmi_preg_mother.ph1",
                     "mat3race.ph1","male.ph1","mat_asthma.ph1","insurancer.ph1","tobacPregr.ph1",
                     "asthma.ph1")]

with(d.mi.A, table(male.ph1,male))
####  There were very few errors in sex of baby. 
####   mice() has issues with the male variable, presumably because
####   male.ph1 is such a good predictor of it. To avoid this problem, we are going
####   to always use the phase-1 sex, unless a specific error was found in the 
####   audit.
d.mi.A$maleA<-with(d.mi.A, ifelse(is.na(male), male.ph1, male))
d.mi.A$male<-d.mi.A$male.ph1<-NULL

####  Making sure that mice is including the correct variables and using correct models
init<-mice(d.mi.A, maxit=0)
pred.matrix <- init$predictorMatrix
my.method<-init$method

####  Multiply imputing missing data
mice_data.A <- mice(d.mi.A, predictorMatrix=pred.matrix, method=my.method, maxit=20, m=25, print=TRUE)
mice_data.A$loggedEvents

####  Fitting the model to the multiply imputed data
mod.mi.A <- with(mice_data.A, 
                glm(asthma ~ estWtChangePerWk + est_bmi_preg_mother + mat_age_delivery + 
                      mat3race + maleA + insurancer + tobacPregr + mat_asthma + egaWk, family="binomial"))

####  Pooling results across MI replications
pooled.est.A<-pool(mod.mi.A)
mi.results.A<-summary(pooled.est.A)
mi.results.A

####  Printing out results; these are in Table 2 in the manuscript
ests.mi1.A<-cbind(mi.results.A$estimate,
                  mi.results.A$std.error,
                  mi.results.A$estimate-1.96*mi.results.A$std.error,
                  mi.results.A$estimate+1.96*mi.results.A$std.error)
ests.mi.A<-ests.mi1.A
ests.mi.A[3,]<-ests.mi1.A[3,]*5 # changing per 5 BMI units
ests.mi.A[4,]<-ests.mi1.A[4,]*10 # changing per 10 years
OR.mi.A<-exp(ests.mi.A[,c(1,3,4)])
round(ests.mi.A,3)


####  Because of potential incongeniality issues, we now do a Bootstrap MI procedure
#library(bootImpute)

imps.A<-bootMice(d.mi.A, nBoot=200, nImp=2, seed=123)
analysisImp.A <- function(inputData) {
  coef(glm(asthma ~ estWtChangePerWk + est_bmi_preg_mother + mat_age_delivery + 
             mat3race + maleA + insurancer + tobacPregr + mat_asthma + egaWk, family="binomial", data=inputData))
}
ans.A<-bootImputeAnalyse(imps=imps.A, analysisfun=analysisImp.A)
ans.A
####  Formatting for presentation in Table S2 of manuscript.
ests.boot.A1<-cbind(ans.A$ests,sqrt(ans.A$var),ans.A$ci)
ests.boot.A<-ests.boot.A1
ests.boot.A[3,]<-ests.boot.A1[3,]*5    # changing BMI coefficient to per 5 kg/m2
ests.boot.A[4,]<-ests.boot.A1[4,]*10   # changing age coefficient to per 10 years
OR.boot.A<-exp(ests.boot.A[,c(1,3,4)])
round(ests.boot.A,3)


####  Naive estimate using unvalidated data. Shown in Table 2 of manuscript.
mod.naive.A<- glm(asthma.ph1 ~ estWtChangePerWk.ph1 + est_bmi_preg_mother.ph1 + mat_age_delivery + 
                    mat3race.ph1 + male.ph1 + insurancer.ph1 + tobacPregr.ph1 + mat_asthma.ph1, family="binomial", data=d.A)
summary(mod.naive.A)
ests.naive1.A<-cbind(mod.naive.A$coefficients,
                     summary(mod.naive.A)$coefficients[,2],
                     mod.naive.A$coefficients-1.96*summary(mod.naive.A)$coefficients[,2],
                     mod.naive.A$coefficients+1.96*summary(mod.naive.A)$coefficients[,2])
ests.naive.A<-ests.naive1.A
ests.naive.A["est_bmi_preg_mother.ph1",]<-ests.naive1.A["est_bmi_preg_mother.ph1",]*5 # changing per 5 BMI units
ests.naive.A["mat_age_delivery",]<-ests.naive1.A["mat_age_delivery",]*10              # changing per 10 years
OR.naive.A<-exp(ests.naive.A[,c(1,3,4)])
round(ests.naive.A,3)




#### Figure comparing naive, MI, and Boot MI estimates for logistic regression / asthma analysis
n.covs<-length(mod.naive.A$coefficients)
gap<-0.02
plot(c(-4,2),c(0,1),type="n",xlab="Asthma Odds Ratio",ylab="", axes=FALSE)
axis(1, at=log(c(0.2,0.5,1,2,5)), lab=c(0.2,0.5,1,2,5))
abline(v=0,col=gray(.9))
points(ests.naive.A[-1,1],c(1:(n.covs-1)/n.covs), pch=3, cex=0.5)
for (i in 1:(n.covs-1)){
  lines(c(ests.naive.A[i+1,3],ests.naive.A[i+1,4]),rep(i/n.covs,each=2))
}

points(ests.mi.A[-1,1],c(1:n.covs)/n.covs-gap, col=3, pch=3,cex=0.5)
for (i in 1:n.covs){
  lines(c(ests.mi.A[i+1,3],ests.mi.A[i+1,4]),rep(i/n.covs,each=2)-gap,col=3)
}

points(ests.boot.A[-1,1],c(1:n.covs)/n.covs-2*gap, col=2, pch=3,cex=0.5)
for (i in 1:n.covs){
  lines(c(ests.boot.A[i+1,3],ests.boot.A[i+1,4]),rep(i/n.covs-2*gap,each=2), col=2)
}
text(rep(-3.5,n.covs),c(1:n.covs)/n.covs,
      c("Wgt Change (per kg/wk)", "BMI (per 5 kg/m2)","Age (per 10 years)","Black","Other Race","Male Sex",
        "No Insurance","Smoking","Maternal Asthma","EGA (weeks)"), pos=4)






#########################################################################
####  Childhood obesity analyses - Time-to-event outcome
####  Our model is a Cox regression model of the following form:
####  coxph(Surv(ttobesityC,obesityC) ~ estWtChangePerWk + est_bmi_preg_mother + 
####                         mat_age_delivery + mat4race + hispanic + mat_diabetes3 + cesarean + 
####                         maleA + depressionr + insurancer + singleton + tobacPregr + married + 
####                         number_children + egaWk))
####  There are no inclusion restrictions


set.seed(20250822)

library(survival)
library(Hmisc)

####  Reducing the dataset in preparation for MI
d.mi.o<-d[,c("estWtChangePerWk","est_bmi_preg_mother",
               "mat4race","hispanic","mat_diabetes3","cesarean",
               "male","depressionr","insurancer","singleton",
               "tobacPregr","married","number_children","egaWk",
               "ttobesity","obesity",
               "estWtChangePerWk.ph1","est_bmi_preg_mother.ph1","mat_age_delivery",
               "mat4race.ph1","hispanic.ph1","mat_diabetes3.ph1","cesarean.ph1",
               "male.ph1","depressionr.ph1","insurancer.ph1","singleton.ph1",
               "tobacPregr.ph1",
               "ttobesity.ph1","obesity.ph1")]

####  Assuming that there are no errors in the sex variable
d.mi.o$maleA<-with(d.mi.o, ifelse(is.na(male), male.ph1, male))
d.mi.o$male<-d.mi.o$male.ph1<-NULL

####  mice has problems with highly correlated variables like obesity.ph1 and obesity
####   and ttobesity.ph1 and ttobesity. We'll do this imputation by hand.

####  Estimating positive predictive values (PPV) and negative predictive values (NPV) for obesity
####   based on phase-1 data. 
table(d.mi.o$obesity.ph1,d.mi.o$obesity)
ests.ppv<-log(1-binconf(table(d.mi.o$obesity.ph1,d.mi.o$obesity)[2,2],sum(table(d.mi.o$obesity.ph1,d.mi.o$obesity)[2,]),method="wilson"))
ests.ppv
sd.ests.ppv<-((abs(ests.ppv[1]-ests.ppv[2])+abs(ests.ppv[1]-ests.ppv[3]))/2)/2
ests.npv<-log(1-binconf(table(d.mi.o$obesity.ph1,d.mi.o$obesity)[1,1],sum(table(d.mi.o$obesity.ph1,d.mi.o$obesity)[1,]),method="wilson"))
ests.npv
sd.ests.npv<-((abs(ests.npv[1]-ests.npv[2])+abs(ests.npv[1]-ests.npv[3]))/2)/2

nImp<-25
est.coeff<-est.vars<-matrix(NA, nImp, 18)

for (j in 1:nImp){

  #### sampling phase-2 obesity status
  ppv.samp<-1-exp(rnorm(1,ests.ppv[1],sd.ests.ppv))
  npv.samp<-1-exp(rnorm(1,ests.npv[1],sd.ests.npv))
  obesity.samp.pos<-rbinom(length(d.mi.o$obesity),1,ppv.samp)
  obesity.samp.neg<-1-rbinom(length(d.mi.o$obesity),1,npv.samp)
  d.mi.o$obesityC<-with(d.mi.o, ifelse(!is.na(obesity), as.numeric(obesity)-1,
                                       ifelse(obesity.ph1==1, obesity.samp.pos, 
                                              ifelse(obesity.ph1==0, obesity.samp.neg, NA))))
  #####  Understanding the relationship between phase-1 and phase-2 time-to-obesity
  cor(d.mi.o$ttobesity,d.mi.o$ttobesity.ph1,use="complete.obs")
  plot(d.mi.o$ttobesity,d.mi.o$ttobesity.ph1)
  diff<-d.mi.o$ttobesity-d.mi.o$ttobesity.ph1
  summary(diff[d.mi.o$obesity==1])
  summary(diff[d.mi.o$obesity==0])
  mean(diff[d.mi.o$obesity==1]==0,na.rm=TRUE)
  mean(diff[d.mi.o$obesity==0]==0,na.rm=TRUE)
  p.tto.err0<-mean(diff[d.mi.o$obesity==0]==0,na.rm=TRUE)
  p.tto.err1<-mean(diff[d.mi.o$obesity==1]==0,na.rm=TRUE)

  sd.p.tto.err0<-sqrt(p.tto.err0*(1-p.tto.err0)/sum(d.mi.o$obesity==0,na.rm=TRUE))
  sd.p.tto.err1<-sqrt(p.tto.err1*(1-p.tto.err1)/sum(d.mi.o$obesity==1,na.rm=TRUE))

  p.err0<-rnorm(1,p.tto.err0,sd.p.tto.err0)
  p.err1<-rnorm(1,p.tto.err1,sd.p.tto.err1)

  ####  Imputing the time-to-obesity or censoring based on relationships seen in phase-2 data
  err.samp0<-rbinom(length(d.mi.o$obesity),1,p.err0)
  err.samp1<-rbinom(length(d.mi.o$obesity),1,p.err1)
  resids<-diff[diff!=0 & !is.na(diff)]
  resids.samp<-sample(resids, length(d.mi.o$obesity), replace=TRUE)
  ttobesityC1<-with(d.mi.o,
                          ifelse(!is.na(ttobesity),ttobesity,
                          ifelse(obesityC==0 & err.samp0==1, ttobesity.ph1,
                          ifelse(obesityC==0 & err.samp0==0, ttobesity.ph1+resids.samp,
                          ifelse(obesityC==1 & err.samp1==1, ttobesity.ph1,
                          ifelse(obesityC==1 & err.samp1==0, ttobesity.ph1+resids.samp,NA
                                 ))))))
  d.mi.o$ttobesityC<-pmin(ttobesityC1,6)

  ####  Best practice is to include baseline cumulative hazard function in MI model
  ####   for imputing covariates (White and Royston, 2009), so we are computing the 
  ####   baseline hazard and having it replace the time to obesity variables in the MI 
  ####   model.
  cumHaz<-survfit(Surv(ttobesityC, obesityC)~1, ctype=1, data=d.mi.o)
  cbind(cumHaz$time,cumHaz$cumhaz)
  cum.hazC<-NULL
  for (i in 1:length(d.mi.o$ttobesity.ph1)) {
    cum.hazC[i]<-max(cumHaz$cumhaz[cumHaz$time<=d.mi.o$ttobesityC[i]])
  }  
  d.mi.o$cum.hazC<-cum.hazC

  d.mi.o1<-d.mi.o
  d.mi.o1$obesity<-d.mi.o1$ttobesity<-d.mi.o1$obesity.ph1<-d.mi.o1$ttobesity.ph1<-NULL

  ####  Making sure our imputation model is set up OK. We do not want to impute with ttobesityC.
  init<-mice(d.mi.o1, maxit=0)
  pred.matrix <- init$predictorMatrix
  pred.matrix[,"ttobesityC"]<-0
  pred.matrix
  init$method

  ####  Running mice with 1 imputation for each of 25 iterations of the loop.
  mice_data.o <- mice(d.mi.o1, predictorMatrix=pred.matrix, maxit=20, m=1, print=TRUE)
  mice_data.o$loggedEvents

  mice_data.o1<-complete(mice_data.o)

  mod.mi.o <- with(mice_data.o1, 
                   coxph(Surv(ttobesityC,obesityC) ~ estWtChangePerWk + est_bmi_preg_mother + 
                         mat_age_delivery + mat4race + hispanic + mat_diabetes3 + cesarean + 
                         maleA + depressionr + insurancer + singleton + tobacPregr + married + 
                         number_children + egaWk))

  est.coeff[j,]<-mod.mi.o$coefficients
  est.vars[j,]<-diag(mod.mi.o$var)
}
est.coeff
est.vars

####  Combining using Rubin's rules
estimate<-colMeans(est.coeff)
std.error <- sqrt(colMeans(est.vars) + (1+1/nImp)*apply(est.coeff, 2, FUN=var))

####  Our estimates, reported in Table 3 of the manuscript
mi.results.o<-data.frame(estimate, std.error)
ests.mi.o<-cbind(mi.results.o$estimate,
                 mi.results.o$std.error,
                  mi.results.o$estimate-1.96*mi.results.o$std.error,
                  mi.results.o$estimate+1.96*mi.results.o$std.error)
ests.mi.o[2,]<-ests.mi.o[2,]*5 # changing per 5 BMI units
ests.mi.o[3,]<-ests.mi.o[3,]*10 # changing per 10 years
round(ests.mi.o,3)



####  I had issues using bootMice for the time-to-event imputation, so I essentially
####   repeated the above process for MI across bootstrap replications. Within each 
####   bootstrap, I am only doing 1 imputation

set.seed(20250822)
nBoot<-200
nImp.b<-1
est.coeff.b<-matrix(NA, nBoot, 18)

for (jj in 1:nBoot){
  bootsamp<-sample(1:dim(d.mi.o)[1],dim(d.mi.o)[1],replace=TRUE)
  d.boot<-d.mi.o[bootsamp,]  
  ppv.samp.b<-table(d.boot$obesity.ph1,d.boot$obesity)[2,2]/sum(table(d.boot$obesity.ph1,d.boot$obesity)[2,]) #
  npv.samp.b<-table(d.boot$obesity.ph1,d.boot$obesity)[1,1]/sum(table(d.boot$obesity.ph1,d.boot$obesity)[1,]) #
  obesity.samp.pos.b<-rbinom(length(d.boot$obesity),1,ppv.samp.b) 
  obesity.samp.neg.b<-1-rbinom(length(d.boot$obesity),1,npv.samp.b) 
  d.boot$obesityC<-with(d.boot, ifelse(!is.na(obesity), as.numeric(obesity)-1,
                                     ifelse(obesity.ph1==1, obesity.samp.pos.b, 
                                            ifelse(obesity.ph1==0, obesity.samp.neg.b, NA))))
  diff.b<-d.boot$ttobesity-d.boot$ttobesity.ph1
  p.tto.err0.b<-mean(diff.b[d.boot$obesity==0]==0,na.rm=TRUE)
  p.tto.err1.b<-mean(diff.b[d.boot$obesity==1]==0,na.rm=TRUE)
  err.samp0.b<-rbinom(length(d.boot$obesity),1,p.tto.err0.b) #p.err0.b
  err.samp1.b<-rbinom(length(d.boot$obesity),1,p.tto.err1.b) #p.err1.b
  resids.b<-diff.b[diff.b!=0 & !is.na(diff.b)]
  resids.samp.b<-sample(resids.b, length(d.boot$obesity), replace=TRUE)
  ttobesityC1.b<-with(d.boot,
                        ifelse(!is.na(ttobesity),ttobesity,
                        ifelse(obesityC==0 & err.samp0.b==1, ttobesity.ph1,
                        ifelse(obesityC==0 & err.samp0.b==0, ttobesity.ph1+resids.samp.b,
                        ifelse(obesityC==1 & err.samp1.b==1, ttobesity.ph1,
                        ifelse(obesityC==1 & err.samp1.b==0, ttobesity.ph1+resids.samp.b,NA
                               ))))))
  d.boot$ttobesityC<-pmin(ttobesityC1.b,6)

  cumHaz.b<-survfit(Surv(ttobesityC, obesityC)~1, ctype=1, data=d.boot)
  cbind(cumHaz.b$time,cumHaz.b$cumhaz)
  cum.hazC.b<-NULL
  for (i in 1:length(d.boot$ttobesity.ph1)) {
    cum.hazC.b[i]<-max(cumHaz.b$cumhaz[cumHaz.b$time<=d.boot$ttobesityC[i]])
  }  
  d.boot$cum.hazC<-cum.hazC.b

  d.boot1<-d.boot
  d.boot1$obesity<-d.boot1$ttobesity<-d.boot1$obesity.ph1<-d.boot1$ttobesity.ph1<-NULL

  init.b<-mice(d.boot1, maxit=0)
  pred.matrix.b <- init.b$predictorMatrix
  pred.matrix.b[,"ttobesityC"]<-0

  mice_data.o.b <- mice(d.boot1, predictorMatrix=pred.matrix.b, maxit=20, m=nImp.b, print=TRUE)
  mice_data.o.b$loggedEvents

  mice_data.o1.b<-complete(mice_data.o.b)

  mod.mi.o.b <- with(mice_data.o1.b, 
                 coxph(Surv(ttobesityC,obesityC) ~ estWtChangePerWk + est_bmi_preg_mother + 
                         mat_age_delivery + mat4race + hispanic + mat_diabetes3 + cesarean + 
                         maleA + depressionr + insurancer + singleton + tobacPregr + married + 
                         number_children + egaWk))

  est.coeff.b[jj,]<-mod.mi.o.b$coefficients
}

####  Displaying results of Boot MI. These are somewhat different from MI with Rubin's rules,
####   including estimates obtainedby taking the mean of the bootstrap replications. The standard
####   error estimates are currently included in Table 3 of the manuscript.
ests.boot.o<-cbind(apply(est.coeff.b, 2, mean), 
                   apply(est.coeff.b, 2, sd),
                   apply(est.coeff.b, 2, quantile, 0.025),
                   apply(est.coeff.b, 2, quantile, 0.975))
ests.boot.o[2,]<-ests.boot.o[2,]*5
ests.boot.o[3,]<-ests.boot.o[3,]*10
round(ests.boot.o,3)



####  Computing the naive estimate using Phase-1 data only. These are reported in Table 3 of manuscript.
mod.naive.o<- coxph(Surv(ttobesity.ph1,as.numeric(obesity.ph1)) ~ estWtChangePerWk.ph1 + est_bmi_preg_mother.ph1 + 
                      mat_age_delivery + mat4race.ph1 + hispanic.ph1 + mat_diabetes3.ph1 + cesarean.ph1 + 
                      male.ph1 + depressionr.ph1 + insurancer.ph1 + singleton.ph1 + tobacPregr.ph1, data=d)
summary(mod.naive.o)
ests.naive1.o<-cbind(mod.naive.o$coefficients,
                     summary(mod.naive.o)$coefficients[,3],
                     mod.naive.o$coefficients-1.96*summary(mod.naive.o)$coefficients[,3],
                     mod.naive.o$coefficients+1.96*summary(mod.naive.o)$coefficients[,3])
ests.naive.o<-ests.naive1.o
ests.naive.o["est_bmi_preg_mother.ph1",]<-ests.naive1.o["est_bmi_preg_mother.ph1",]*5 # changing per 5 BMI units
ests.naive.o["mat_age_delivery",]<-ests.naive1.o["mat_age_delivery",]*10 # changing per 10 years
HR.naive.o<-exp(ests.naive.o)              
round(ests.naive.o,3)


####  Figure comparing hazard ratios and confidence interval estimates between
####   MI with Rubin's rules and Bootstrap MI.
n.covs<-length(ests.mi.o[,1])
gap<-0.01
plot(c(-4,2),c(0,1),type="n",xlab="Obesity Hazard Ratio",ylab="", axes=FALSE)
axis(1, at=log(c(0.2,0.5,1,2,5)), lab=c(0.2,0.5,1,2,5))
abline(v=0,col=gray(.9))

points(ests.mi.o[,1],c(1:n.covs)/n.covs-gap, col=3, pch=3,cex=0.5)
for (i in 1:n.covs){
  lines(c(ests.mi.o[i,3],ests.mi.o[i,4]),rep(i/n.covs,each=2)-gap,col=3)
}

points(ests.boot.o[,1],c(1:n.covs)/n.covs-2*gap, col=2, pch=3,cex=0.5)
for (i in 1:n.covs){
  lines(c(ests.boot.o[i,3],ests.boot.o[i,4]),rep(i/n.covs-2*gap,each=2), col=2)
}
text(rep(-3.5,n.covs),c(1:n.covs)/n.covs,
      c("Wgt Change (per kg/wk)", "BMI (per 5 kg/m2)","Age (per 10 years)","Asian","Black",
        "Other","Hispanic","Gestational Diabetes","Type 1/2 Diabetes", "Cesarean","Male Sex",
        "Depression","No Insurance","Singleton","Smoking","Married","Number of Children","EGA (weeks)"), pos=4)

