

########################################################################################
####  These analyses perform semiparametric maximumu likelihood estimation to account 
####   for errors in variables in the synthetic mother-child weight study.
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

####  Limiting weight gain analysis to mothers giving birth to only one child
d.singleton<-d[d$singleton.ph1==1 & (is.na(d$singleton) | d$singleton==1),]
####  Creating a reduced dataset with only the necessary variables
d.mi<-d.singleton[,c("estWtChangePerWk","est_bmi_preg_mother","mat_age_delivery",
                     "mat4race","hispanic","depressionr","insurancer","tobacPregr",
                     "estWtChangePerWk.ph1","est_bmi_preg_mother.ph1",
                     "mat4race.ph1","hispanic.ph1","depressionr.ph1","insurancer.ph1",
                     "tobacPregr.ph1")]



#######################################################
########## Installing sleev package
#######################################################

#install.packages("sleev")
library("sleev")

########################################################################################
## sleev requires that we re-code factor variables as dummy variables. This is just the
##  opposite of mice (at least using its defaults). So basically, I'm undoing what I did
##  earlier in this code. I could have just skipped all of this, but I am trying to work 
##  with the same original dataset for both MI and SMLE analyses. Hence, I'm recoding 
##  variables here.

d.mi.num<-d.mi
d.mi.num$mat4race.black.ph1<-ifelse(d.mi.num$mat4race.ph1=="Black",1,0)
d.mi.num$mat4race.asian.ph1<-ifelse(d.mi.num$mat4race.ph1=="Asian",1,0)
d.mi.num$mat4race.other.ph1<-ifelse(d.mi.num$mat4race.ph1=="Other",1,0)
d.mi.num$mat4race.black<-ifelse(d.mi.num$mat4race=="Black",1,0)
d.mi.num$mat4race.asian<-ifelse(d.mi.num$mat4race=="Asian",1,0)
d.mi.num$mat4race.other<-ifelse(d.mi.num$mat4race=="Other",1,0)
d.mi.num$hispanic<-ifelse(d.mi.num$hispanic==1,1,0)
d.mi.num$hispanic.ph1<-ifelse(d.mi.num$hispanic.ph1==1,1,0)
d.mi.num$depressionr<-ifelse(d.mi.num$depressionr==1,1,0)
d.mi.num$depressionr.ph1<-ifelse(d.mi.num$depressionr.ph1==1,1,0)
d.mi.num$insurancer<-ifelse(d.mi.num$insurancer=="Public/Other",1,0)
d.mi.num$insurancer.ph1<-ifelse(d.mi.num$insurancer.ph1=="Public/Other",1,0)
d.mi.num$tobacPregr<-ifelse(d.mi.num$tobacPregr==1,1,0)
d.mi.num$tobacPregr.ph1<-ifelse(d.mi.num$tobacPregr.ph1==1,1,0)


### Because (X*,Z) has so many dimensions, we need to collapse the data in some manner.
###  I do this using factor analysis of mixed data, using the FactoMineR library. This 
###  library handles factor variables, so I am using the factor variables (not the dummy)
###  variables to collapse the data. Then I grab the first two factors and create B-splines
###  from them. This is slightly worrisome, because the first two factors contain only a 
###  relatively small amount of the variation in (X*,Z). But it's a necessary simplification.

#install.packages("FactoMineR")
library(FactoMineR)

phase1.vars<-d.mi[,c("est_bmi_preg_mother.ph1","mat_age_delivery","mat4race.ph1",
                     "hispanic.ph1","depressionr.ph1","insurancer.ph1","tobacPregr.ph1")]
famd<-FAMD(base=phase1.vars)
famd$eig
famd$var
famd2<-famd$ind$coord[,c(1:2)]

####  Adding the first two factors into the analysis dataset
d.mi.smle1<-d.mi.num
d.mi.smle1$famd1<-famd2[,1]
d.mi.smle1$famd2<-famd2[,2]

####  Creating B-splines using the first two factors
data.famd <- spline2ph(x=c("famd1","famd2"), size=4, degree=3, data=d.mi.smle1)

####  Fitting the SMLE. This takes hours to run.
####    Note that hn_scale has been set to 1/100 to get SEs. Results are similar with
####     hn_scale=1/200, but SEs do not converge with hn_scale=1/2 or 1 (the default).
####    Syntax is hopefully self-explanatory and it is detailed in the help file for the 
####     function. y_unval and x_unval are the error-prone unvalidated values of the 
####     outcome, y, and covariates, x. z are those covariates that have no errors.
start.time <- Sys.time()
mod.linear.smle2 <- linear2ph(y = "estWtChangePerWk",
                              y_unval = "estWtChangePerWk.ph1",
                              x = c("est_bmi_preg_mother","mat4race.asian","mat4race.black","mat4race.other",
                                    "hispanic",
                                    "depressionr","insurancer","tobacPregr"), 
                              x_unval = c("est_bmi_preg_mother.ph1","mat4race.asian.ph1","mat4race.black.ph1","mat4race.other.ph1",
                                          "hispanic.ph1",
                                          "depressionr.ph1","insurancer.ph1","tobacPregr.ph1"),
                              z = "mat_age_delivery",
                              data = data.famd, hn_scale = 1/100,
                              se = TRUE, tol = 1e-04, max_iter = 1000,
                              verbose = TRUE)
paste0("Run time: ", round(difftime(Sys.time(), start.time,
                                    units = "secs"), 3), " sec")
round(mod.linear.smle2$coefficients,3)
summary(mod.linear.smle2)





##############################################################
### Binary asthma outcome with sleev
##############################################################

#########################################################################
####  Childhood asthma analyses - Binary outcome
####  Our model is a logistic regression model of the following form:
####      asthma ~ estWtChangePerWk + est_bmi_preg_mother + mat_age_delivery + 
####                      mat3race + maleA + insurancer + tobacPregr + mat_asthma + egaWk
####  And we are going to limit the analyses to eligible children.

####  I start by re-creating the analysis dataset used for the MI analyses, for consistency
####   across code.

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

##################################################################################################
#### Now I am doing a few other data prep / manipulation to fit the SMLE:
####  sleev doesn't like factor variables, so I'm changing all factor variables to dummy 0-1 variables
####  sleev needs to have the same number of x and x_unval variables. egaWk was not collected
####   during phase-1, so I artificially create egaWk.ph1 and assign everyone to the value of 40.
d.mi.A$egaWk.ph1<-40 
d.mi.A.num<-d.mi.A
d.mi.A.num$mat3race.black.ph1<-ifelse(d.mi.A.num$mat3race.ph1=="Black",1,0)
d.mi.A.num$mat3race.other.ph1<-ifelse(d.mi.A.num$mat3race.ph1=="Other",1,0)
d.mi.A.num$mat3race.black<-ifelse(d.mi.A.num$mat3race=="Black",1,0)
d.mi.A.num$mat3race.other<-ifelse(d.mi.A.num$mat3race=="Other",1,0)
d.mi.A.num$insurancer<-ifelse(d.mi.A.num$insurancer=="Public/Other",1,0)
d.mi.A.num$insurancer.ph1<-ifelse(d.mi.A.num$insurancer.ph1=="Public/Other",1,0)
d.mi.A.num$tobacPregr<-ifelse(d.mi.A.num$tobacPregr==1,1,0)
d.mi.A.num$tobacPregr.ph1<-ifelse(d.mi.A.num$tobacPregr.ph1==1,1,0)
d.mi.A.num$asthma.ph1<-as.numeric(d.mi.A.num$asthma.ph1)
d.mi.A.num$mat3race.ph1<-NULL
d.mi.A.num$mat3race<-NULL

#########################################################################################
## I have a lot of X* and Z variables. So I am collapsing them down to 2 variables
##  using the FactoMineR package (factor analysis). This of course is a simplification,
##  but a necessary one. FactoMineR seems to like factor variables, so I'm using the 
##  original dataset with factor variables.
library(FactoMineR)
phase1.vars.A<-d.mi.A[,c("estWtChangePerWk.ph1","est_bmi_preg_mother.ph1","mat_age_delivery","mat3race.ph1",
                         "maleA","insurancer.ph1","tobacPregr.ph1","mat_asthma.ph1","egaWk.ph1")]
famd.A<-FAMD(base=phase1.vars.A)
famd.A$eig
famd.A$var
famd2.A<-famd.A$ind$coord[,c(1:2)]

####  Now I incorporate these factor variables into the analysis dataset
d.mi.smle1.A<-d.mi.A.num
d.mi.smle1.A$famd1<-famd2.A[,1]
d.mi.smle1.A$famd2<-famd2.A[,2]

####  Creating B-splines using the factor variables
data.famd.A <- spline2ph(x=c("famd1","famd2"), size=4, degree=3, data=d.mi.smle1.A)


####  Fitting the SMLE. A few things to note: 
####     I set hn_scale=0.1 instead of its default, so that I could estimate SEs.
####     This takes several hours to run.    
start.time <- Sys.time()
mod.logistic.smle2 <- logistic2ph(y = "asthma",
                                  y_unval = "asthma.ph1",
                                  x = c("estWtChangePerWk","est_bmi_preg_mother",
                                        "mat3race.black","mat3race.other",
                                        "insurancer","tobacPregr","mat_asthma",
                                        "egaWk"),
                                  x_unval = c("estWtChangePerWk.ph1","est_bmi_preg_mother.ph1",
                                              "mat3race.black.ph1","mat3race.other.ph1",
                                              "insurancer.ph1","tobacPregr.ph1","mat_asthma.ph1"),
                                  z = c("mat_age_delivery","maleA"),
                                  data = data.famd.A, hn_scale = 0.1,
                                  se = TRUE, tol = 1e-04, max_iter = 1000,
                                  verbose = TRUE)
paste0("Run time: ", round(difftime(Sys.time(), start.time,
                                    units = "secs"), 3), " sec")
mod.logistic.smle2



#### Naive estimate using only error-prone phase-1 data (without egaWk) for comparison
mod.naive.A<- glm(asthma.ph1 ~ estWtChangePerWk.ph1 + est_bmi_preg_mother.ph1 + mat_age_delivery +
                    mat3race.ph1 + male.ph1 + insurancer.ph1 + tobacPregr.ph1 + mat_asthma.ph1, family="binomial", data=d.A)
summary(mod.naive.A)

