#The goal of this code is to develop a prototype for a 'GJRM-L' function which allows a user to specify a model
#in the form of a univariate longitudinal model e.g. lme4, gamlss with random, and have it estimated as a 
#copula-based joint regression model.
#GJRM-L(mu.formula, sigma.formula,nu.formula,zeta.forula, tau.formula, other.formula,
#margin.family=list(), copula.family=list()) #All optional except mu.formula

#### 0. Starting parameters ####
set.seed(1000);options(scipen=999);
source("common_functions.R")
library(copula)
library(gamlss)

#### 0. Functions ####

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

fit_margins <- function(mu.formula = formula(response~1), sigma.formula = ~1, 
                        nu.formula = ~1, tau.formula = ~1,dataset,family) {
  #i=1
  require(gamlss)
  fits<-list()
  for (i in 1:max(dataset$time)) {
    fits[[i]]<-gamlss(formula=mu.formula,family=family,data=dataset[dataset$time==i,])
  }
  return(fits)
}

find_best_marginal_fits <- function(dataset,type) {
  require(gamlss)
  fits<-list()
  for (i in 1:max(dataset$time)) {
    fits[[i]]<-fitDist(dataset[dataset$time==i,"response"],type = type)$fits[1:6]
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

likelihood_fitted_margin_copula <- function(fitted_margins,fitted_copulas) { 
  
  likelihood=fitted_copulas$logLik
  for (i in 1:length(fitted_margins)) {
    likelihood=likelihood+logLik(fitted_margins[[i]])
  }
  
  df=0
  for ( i in 1:length(fitted_margins)) {
    df=df+fitted_margins[[i]]$df.fit
  }
  df = df+(length(fitted_margins)*(length(fitted_margins)-1))
  
  return(
    paste("Overall LogLik:", as.character(round(likelihood,2)),
          "\nGlobal Deviance:", as.character(round(likelihood,2)*-2),
          "\nDegrees of Freedom:", as.character(round(df,2)))
  )
  
} ###List of gamlss objects, VineCopula object as input

#### 1. Load RAND data subset and transform to standard longitudinal dataset #######

load("Data/rand_mvt.rds")
head(rand_mvt)

# Basic data setup
response = rand_mvt[,4:(4+4)]#[,4:18] ####Currently limiting to just 5 margins for simplicity
covariates=list()
covariates[[1]] = as.data.frame(rand_mvt[,19]) #Age 19:33 - changed to age at start to avoid correlation with time
covariates[[2]] = as.data.frame(rand_mvt[,34:(34+4)]) #Time 34:48
covariates[[3]] = as.data.frame(rand_mvt[,3]) #Gender

# Setup data as longitudinal file
dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
head(dataset)


#### 2. Benchmark models ####

###### Regular ZISICHEL Global Deviance: 47647.01
margin_model_formula=formula(response~(age)+as.factor(gender)+as.factor(time))
gamlss_model <- gamlss(formula = margin_model_formula,family="ZISICHEL",data=dataset)
summary(gamlss_model)
plot(gamlss_model)
term.plot(gamlss_model)

###### GLMM ZISICHEL Global Deviance: 47647.01
margin_model_formula=formula(response~(age)+as.factor(gender)+as.factor(time)+random(as.factor(subject)))
gamlss_glmm_model <- gamlss(formula = margin_model_formula,family="ZISICHEL",data=dataset)
summary(gamlss_glmm_model)
plot(gamlss_glmm_model)
term.plot(gamlss_glmm_model)

#### 3. Starting parameter models ####

#Goal

# 1. Starting fit for margins and copula
#best_marginal_fits <- find_best_marginal_fits(dataset,type="counts") #Find best fit margins
fitted_margins<-fit_margins(mu.formula=formula(response~age+as.factor(gender)),dataset=dataset,family=ZISICHEL())
  #AIC(fitted_margins[[1]])
  #plot(fitted_margins[[1]])
  #term.plot(fitted_margins[[2]])

fitted_copulas<-fit_copulas(fitted_margins)
  #contour(fitted_copulas)

cat(likelihood_fitted_margin_copula(fitted_margins,fitted_copulas))

#### 4. 

#1. Calculate likelihood for given set of parameters
#2. Be able to iterate parameters in direction of optimal

#### Appendix ###












