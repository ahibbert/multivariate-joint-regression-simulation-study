#source("common_functions.R")

###Functions - run first
set.seed(1000);options(scipen=999);

create_longitudinal_dataset <- function(response,covariates,labels=NA) {
  num_time_points=ncol(response)
  if(num_time_points <=1) {print('Not enough time points')}
  
  dataset<-matrix(data=NA,ncol=2+length(covariates),nrow=0) 
  subject<-as.factor(seq(1:nrow(response)))
  
  for (t in 1:ncol(response)) {
    
    dataset_temp<-cbind(subject,t,response[,t])
    
    for (i in 1:length(covariates)) {
      if (ncol(covariates[[i]]) == 1 ) {
        covariate_for_time=covariates[[i]]
      } else {
        covariate_for_time=covariates[[i]][,t]
      }
      dataset_temp<-cbind(dataset_temp,covariate_for_time)  
    }
    
    ###Add dataset temp to full table
    dataset <- rbind(dataset,dataset_temp)
  }
  
  if(!all(is.na(labels))) {
    colnames(dataset) <- labels  
  }
  
  dataset=dataset[order(dataset$subject,dataset$subject),]
  
  return(dataset)
}

find_best_marginal_fits <- function(margins,type) {
  require(gamlss)
  fits<-list()
  for (i in 1:ncol(margins)) {
    fits[[i]]<-fitDist(margins[,i],type = type)$fits[1:6]
  }
  return(fits)
}

fit_margins <- function(dataset,mu_formula,family) {
  #i=1
  require(gamlss)
  fits<-list()
  for (i in 1:max(dataset$time)) {
    fits[[i]]<-gamlss(mu_formula,family=family,data=dataset[dataset$time==i,])
  }
  return(fits)
}

fit_copulas = function(margin_fits) {
  require(gamlss)
  require(VineCopula)
  
  num_margins=length(margin_fits)
  fit_unif=matrix(0,nrow=length(margin_fits[[1]]$residuals),ncol=num_margins)
  for (i in 1:num_margins) {
    fit_unif[,i]=pNO(margin_fits[[i]]$residuals)    
  }
  
  vinefit<-RVineStructureSelect(fit_unif)
  summary(vinefit)
  print(vinefit$pair.logLik)
  
  return(vinefit)
}


########### 0. Only needed once: Extract needed columns from rand_HRS dataset and save in rand_mvt ######
library(haven)
randFULL <- read_sas("Data/randhrs1992_2020v2.sas7bdat")
docnames<-(grep("R[0-9]+DOCTIM", colnames(randFULL),value=TRUE))
docnames_sorted<-docnames[order(as.numeric(regmatches(docnames,regexpr("[0-9]+",docnames))))]

ages<-(grep("R[0-9]+AGEY_E", colnames(randFULL),value=TRUE))
ages_sorted<-ages[order(as.numeric(regmatches(ages,regexpr("[0-9]+",ages))))]

interview_year<-(grep("R[0-9]+IWENDY", colnames(randFULL),value=TRUE))
interview_year_sorted<-interview_year[order(as.numeric(regmatches(interview_year,regexpr("[0-9]+",interview_year))))]

rand_allyearsdocvisits<-randFULL[,c("HHIDPN","RABYEAR", "RAGENDER",docnames_sorted,ages_sorted,interview_year_sorted)]

#1775 individuals with no missing data
rand_mvt<-as.data.frame(na.omit(rand_allyearsdocvisits))
save(rand_mvt,file="Data/rand_mvt.rds")

###### 1. Load RAND data subset and transform to standard longitudinal dataset #######

load("Data/rand_mvt.rds")
head(rand_mvt)

###Basic data setup
response = rand_mvt[,4:(4+4)]#[,4:18]
covariates=list()
covariates[[1]] = as.data.frame(rand_mvt[,19]) #Age 19:33 - changed to age at start to avoid correlation with time
covariates[[2]] = as.data.frame(rand_mvt[,34:(34+4)]) #Time 34:48
covariates[[3]] = as.data.frame(rand_mvt[,3]) #Gender

####Setup data as longitudinal file
dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
head(dataset)

######## 2. Get initial estimates - marginal fits for GAMLSS, VineCopula ######
##Actually we can use this as a benchmark for all future fits anyway so that's good.
##Maybe incorporate all the other models here too as other benchmarks? Maybe later actually

#margin_dist<- find_best_marginal_fits(response,type="counts")
mu_formula=formula(response~cs(age)+as.factor(gender))
margin_fits<- fit_margins(dataset,mu_formula,family="ZISICHEL") ##CHANGE THIS TO TAKE IN DATASET ABOVE
copula_fits<- fit_copulas(margin_fits)

likelihood=copula_fits$logLik
for (i in 1:length(margin_fits)) {
  likelihood=likelihood+logLik(margin_fits[[i]])
}
print(paste("Overall LogLik:", as.character(round(likelihood,2))))
print(paste("Global Deviance:", as.character(round(likelihood,2)*-2)))

df=0
for ( i in 1:length(margin_fits)) {
  df=df+margin_fits[[i]]$df.fit
}
df = df+(length(margin_fits)*(length(margin_fits)-1))

print(paste("Overall LogLik:", as.character(round(likelihood,2))))
print(paste("Global Deviance:", as.character(round(likelihood,2)*-2)))
print(paste("Degrees of Freedom:", as.character(round(df,2))))

####Global deviance: 44771.36

####TEMP FOR INVESTIGATING FITS

summary(margin_fits[[5]])
plot(margin_fits[[5]])
term.plot(margin_fits[[5]])
contour(copula_fits)

######## 3. Optimise parameters jointly somehow ######



######## 4. Comparison fits

  ###### Regular ZISICHEL Global Deviance: 47647.01
margin_model_formula=formula(response~cs(age)+as.factor(gender)+as.factor(time))
gamlss_model <- gamlss(formula = margin_model_formula,family="ZISICHEL",data=dataset)
summary(baseline_model)
plot(baseline_model)
term.plot(baseline_model)
  
  ###### GLMM ZISICHEL Global Deviance: 47647.01
margin_model_formula=formula(response~cs(age)+as.factor(gender)+as.factor(time)+random(as.factor(subject)))
gamlss_glmm_model <- gamlss(formula = margin_model_formula,family="ZISICHEL",data=dataset)
summary(gamlss_glmm_model)
plot(gamlss_glmm_model)
term.plot(gamlss_glmm_model)

#####################################Copula fits

plot.new()
par(mfrow=c(3,5))
for (i in 1:length(docnames)) {
  histDist(dataset[dataset[,"time"]==i,"random_variable"],family=ZISICHEL,main=paste("Time",i),xlim=c(0,100),ylim=c(0,.2))
}
plot(1:15,colMeans(rand_mvt)[3:17],main="Mean no. doctor visits in year",xlab="Time",ylab="Mean doctor visits")

par(mfrow=c(1,1))
boxplot(dataset[,"random_variable"]~dataset[,"time"],ylim=c(0,20),ylab="doctor visits",xlab="time")
boxplot(rand_mvt$R15DOCTIM~rand_mvt$RABYEAR,ylim=c(0,40))
  
gamlss_fits <- list()
family_spec=ZISICHEL
fit_no<-matrix(NA,ncol=0,nrow=nrow(dataset[dataset[,"time"]==1,]))
fit_unif<-matrix(NA,ncol=0,nrow=nrow(dataset[dataset[,"time"]==1,]))
for (i in 1:length(docnames)) {
  gamlss_fits[[i]]<-gamlss(dataset[dataset[,"time"]==i,"random_variable"]~1,family=family_spec,method=RS(100))
  fit_no<-cbind(fit_no,(gamlss_fits[[i]]$residuals))
  fit_unif<-cbind(fit_unif,pNO(gamlss_fits[[i]]$residuals))
}
colnames(fit_no)<-colnames(fit_unif)<-docnames_sorted

library(ggplot2)
library(ggpubr)

#pairs(fit_no)
#pairs(fit_unif)

copFits<-list()
plots<-list()
library(VineCopula)
for (i in 1:(length(docnames)-1)) {
  copFits[[i]]<-BiCopSelect(fit_unif[,i],fit_unif[,i+1])  
  data<-as.data.frame(fit_unif[,c(i,i+1)])
  colnames(data)<-c("V1","V2")
  plots[[i]]<-ggplot(data) +
    geom_density_2d(aes(x=V1,y=V2))
}

ggarrange(plotlist=plots)

library(VineCopula)
vinefit<-RVineStructureSelect(fit_unif[,c(11,12,13,14,15)])
summary(vinefit)
contour(vinefit)

vinefit$pair.logLik

#library(network);plot(vinefit)


#########Random effect fit

gamlss_glm_fit_nbi<-gamlss(random_variable~time                                                 ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=NBI,method=RS(100))
gamlss_glm_fit_zis<-gamlss(random_variable~time                                                 ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=ZISICHEL,method=RS(100))
gamlss_glm_fit_zis_timecat<-gamlss(random_variable~as.factor(time)                              ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=ZISICHEL,method=RS(100))
gamlss_re_fit_zis<-gamlss(random_variable~time+random(as.factor(patient))                       ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=ZISICHEL,method=RS(100))
gamlss_re_fit_zis_timecat<-gamlss(random_variable~as.factor(time)+random(as.factor(patient))    ,data=as.data.frame(dataset[dataset[,"time"]>=11,]),family=ZISICHEL,method=RS(100))

term.plot(gamlss_re_fit_zis_timecat,ylim="free")

results<-as.data.frame(rbind(
c(logLik(gamlss_glm_fit_nbi), "GLM NBI + time", gamlss_glm_fit_nbi$df.fit)
,c(logLik(gamlss_glm_fit_zis), "GLM ZIS  + time", gamlss_glm_fit_zis$df.fit)
,c(logLik(gamlss_glm_fit_zis_timecat), "GLM ZIS + factor(time)",gamlss_glm_fit_zis_timecat$df.fit)
,c(logLik(gamlss_re_fit_zis), "GLMM ZIS + time",gamlss_re_fit_zis$df.fit)
,c(logLik(gamlss_re_fit_zis_timecat), "GLMM ZIS + factor(time)",gamlss_re_fit_zis_timecat$df.fit)
,c(logLik(gamlss_fits[[11]])+logLik(gamlss_fits[[12]])+logLik(gamlss_fits[[13]])+logLik(gamlss_fits[[14]])+logLik(gamlss_fits[[15]])+vinefit$logLik,"GJRM ~ factor(time)",4*5+10*2)
))

colnames(results)<-c("LogLik","Model","EDF")
results<-results[c("Model","LogLik","EDF")]
results$LogLik<-round(as.numeric(results$LogLik))
results$EDF<-round(as.numeric(results$EDF))

results$AIC2 <- results$LogLik*-2+results$EDF*2
results$AIC4 <- results$LogLik*-2+results$EDF*4
results$BIC <- round(results$LogLik*-2+results$EDF*log(nrow(rand_mvt)))

#sum of marginal fits






##Calculating likelihoods


source("common_functions.R")
data<-generateMvtDist("NO",c(1,2,3),c(1,2,3),matrix(c(0,.5,.1,.5,0,.9,.1,.9,0),nrow=3))
plotDist(data,"NO")












