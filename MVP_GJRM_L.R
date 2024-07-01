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
response = rand_mvt[,4:(4+2)]#[,4:18] ####Currently limiting to just 5 margins for simplicity
covariates=list()
covariates[[1]] = as.data.frame(rand_mvt[,19]) #Age 19:33 - changed to age at start to avoid correlation with time
covariates[[2]] = as.data.frame(rand_mvt[,34:(34+4)]) #Time 34:48
covariates[[3]] = as.data.frame(rand_mvt[,3]) #Gender

# Setup data as longitudinal file
dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
head(dataset)

plotDist(dataset,"ZISICHEL")

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

#1. Calculate likelihood for given set of parameters - see if we can put it into optim...

extract_parameters_from_fitted <- function(fitted_margins,fitted_copulas) {
  mu=list()
  sigma=list()
  nu=list()
  tau=list()
  for (i in 1:length(fitted_margins)) {
    mu[[i]]=coefAll(fitted_margins[[i]])[1]
    sigma[[i]]=coefAll(fitted_margins[[i]])[2]
    nu[[i]]=coefAll(fitted_margins[[i]])[3]
    tau[[i]]=coefAll(fitted_margins[[i]])[4]
  }  
  
  cop_par1=fitted_copulas$par[lower.tri(fitted_copulas$par)]
  cop_par2=fitted_copulas$par2[lower.tri(fitted_copulas$par2)]
  
  return(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau),cop_par1,cop_par2))
}

start_par<-extract_parameters_from_fitted(fitted_margins,fitted_copulas)

#### Get linear predictors and their linked value -> **Creates eta, eta_link_inv** #### 
eta = eta_linkinv = list()

num_margins=length(fitted_margins)
observation_count=nrow(dataset[dataset$time==1,])
margin_dist=ZISICHEL(mu.link = "log", sigma.link = "log", nu.link = "identity", tau.link = "logit")

for (margin_number in 1:num_margins) {
  model_matrices=matrix(nrow=observation_count,ncol=5)
  model_matrices_inv=matrix(nrow=observation_count,ncol=5)
  
  counter=1
  for (parameter in c("mu","sigma","nu","tau")) {
    model_matrices[,counter] = model.matrix(fitted_margins[[margin_number]],what=parameter)    %*% (as.matrix(coef(fitted_margins[[margin_number]],what=parameter))) #These are the linear predictors
    counter=counter+1
  }
  model_matrices[,5] = dataset[dataset$time==margin_number,"response"]
  colnames(model_matrices)=c("mu","sigma","nu","tau","x")
  eta[[margin_number]]=model_matrices
  
  model_matrices_inv[,1] = margin_dist$mu.linkinv(model_matrices[,1])
  model_matrices_inv[,2] = margin_dist$sigma.linkinv(model_matrices[,2])
  model_matrices_inv[,3] = margin_dist$nu.linkinv(model_matrices[,3])
  model_matrices_inv[,4] = margin_dist$tau.linkinv(model_matrices[,4])
  model_matrices_inv[,5] = dataset[dataset$time==margin_number,"response"]
  
  colnames(model_matrices_inv)=c("mu","sigma","nu","tau","x")
  
  eta_linkinv[[margin_number]] = model_matrices_inv

}

#### Get density and probability functions for the values from model -> ** Creates margin_p and margin_d ####

margin_d=margin_p=matrix(nrow=observation_count,ncol=num_margins)

for (margin_number in 1:num_margins) {

  margin_d[,margin_number] = dZISICHEL(eta[[margin_number]][,5]
            ,mu=    eta_linkinv[[margin_number]][,"mu"]
            ,sigma= eta_linkinv[[margin_number]][,"sigma"]
            ,nu=    eta_linkinv[[margin_number]][,"nu"]
            ,tau=   eta_linkinv[[margin_number]][,"tau"])
  
  margin_p[,margin_number] = pZISICHEL(eta[[margin_number]][,5]
                                ,mu=    eta_linkinv[[margin_number]][,"mu"]
                                ,sigma= eta_linkinv[[margin_number]][,"sigma"]
                                ,nu=    eta_linkinv[[margin_number]][,"nu"]
                                ,tau=   eta_linkinv[[margin_number]][,"tau"])

} 

##CHECK - First run should equal original gamlss

for (margin_number in 1:num_margins) {
  if(sum(log(margin_d[,margin_number]))-logLik(fitted_margins[[margin_number]])>1) {print("ERROR: margin_d does not agree with original gamlss LogLik - something has gone wrong!!!")}
}

#num_margins=length(fitted_margins)
#observation_count=nrow(dataset[dataset$time==1,])

#### Get density of the copula function for each pair of margins ####

copest_temp =BiCopEst(margin_p[,1],margin_p[,2],family = BiCopName("N")) #method="itau"? Should be faster if it's just starting fit
parest=c(copest_temp$par,copest_temp$par2)

copula_d=matrix(nrow=observation_count,ncol=num_margins-1)
copula_pdf=tCopula
for (margin_number in 1:(num_margins-1)) {
  copula_d[,margin_number]<-BiCopPDF(margin_p[,margin_number],margin_p[,margin_number+1],family = 1,par=parest[1],par2=parest[2])
}
colnames(copula_d)=paste(1:(num_margins-1),",",2:(num_margins),sep="")

i=1; if(colSums(log(copula_d))[1]-copest_temp$logLik>1) {print("ERROR: copula_d does not agree with original VineCopula LogLik - something has gone wrong!!!")}

#### Return Log Likelihood ####
log_lik_results<-c(colSums(log(margin_d)),colSums(log(copula_d)),sum(colSums(log(margin_d)),colSums(log(copula_d))))
names(log_lik_results)[1:num_margins]=1:num_margins
names(log_lik_results)[length(log_lik_results)]="Total"
print(log_lik_results)


#### Getting first and second derivatives 


margin_dldm = matrix(0,ncol=num_margins,nrow=observation_count)

sigma=(sigma_vector[i]);mu=mu_vector[i];y=dataset_matrix[j,i]

all_derivatives=c(names(margin_dist)[grepl("dl", names(margin_dist))],names(margin_dist)[grepl("d2l", names(margin_dist))])

i=1
derivatives_calculated=matrix(nrow=observation_count,ncol=length((all_derivatives)))
for (d_fun in all_derivatives) {
  print(d_fun)
  derivatives_calculated[,i]=margin_dist[[d_fun]](y=eta_linkinv[[1]][,"x"],mu=eta_linkinv[[1]][,"mu"],sigma=eta_linkinv[[1]][,"sigma"],nu=eta_linkinv[[1]][,"nu"],tau=eta_linkinv[[1]][,"tau"])
  i=i+1
}
colnames(derivatives_calculated)=all_derivatives
print(derivatives_calculated)





#### Extract score function and hessian for the given margin + copula
#### Get all the first order derivatives





#i=1
#(fitted_margins[[i]]$mu.x)%*%as.vector(fitted_margins[[i]]$mu.coefficients)

#Input: Parameters for marginal fits and copula fits
#Output: Likelihood and maybe hessian?


#2. Be able to iterate parameters in direction of optimal

#### Appendix ###












