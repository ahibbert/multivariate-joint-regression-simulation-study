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
                        nu.formula = ~1, tau.formula = ~1,dataset,margin_dist) {
  #i=1
  require(gamlss)
  fits<-list()
  for (i in 1:max(dataset$time)) {
    fits[[i]]<-gamlss(formula=mu.formula,family=margin_dist,data=dataset[dataset$time==i,])
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

fit_copulas = function(margin_fits,copula_dist="N",method="vine") {
  require(gamlss)
  require(VineCopula)
  num_margins=length(margin_fits)
  fit_unif=matrix(0,nrow=length(margin_fits[[1]]$residuals),ncol=num_margins)
  for (i in 1:num_margins) {
    fit_unif[,i]=pNO(margin_fits[[i]]$residuals)    
  }
  
  if (method=="novine") {
    cop_fit=list()  
    for (i in 1:(num_margins-1)){
      cop_fit[[i]]=BiCopEst(fit_unif[,i],fit_unif[,i+1],family=as.numeric(BiCopName(copula_dist)))
    }
    names(cop_fit)=paste(1:(num_margins-1),(1:(num_margins-1))+1,sep=",")
    return(cop_fit)
  } else {
    vinefit<-RVineStructureSelect(fit_unif)
    return(vinefit)
  }
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

###### Regular ZISICHEL logLik(gamlss_model) = -13591.09 (df=8) AIC: 27198.18 | 27250.86
margin_model_formula=formula(response~(age)+as.factor(gender)+as.factor(time))
gamlss_model <- gamlss(formula = margin_model_formula,family="ZISICHEL",data=dataset)
summary(gamlss_model)
plot(gamlss_model)
term.plot(gamlss_model)

###### GLMM ZISICHEL logLik(gamlss_glmm_model) = -11187.34 (df=1042.17) AIC: 24459.02 | BIC: 31322.54
margin_model_formula_glmm=formula(response~(age)+as.factor(gender)+as.factor(time)+random(as.factor(subject)))
gamlss_glmm_model <- gamlss(formula = margin_model_formula_glmm,family="ZISICHEL",data=dataset)
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

fitted_copulas<-fit_copulas(fitted_margins,copula_dist="N",method="novine")
#fitted_copulas<-fit_copulas(fitted_margins,copula_dist="N",method="vine")
#likelihood_fitted_margin_copula(fitted_margins,fitted_copulas) ###Not much better results with vine versus optimising non-vine
#contour(fitted_copulas)


log_lik_list=vector()
num_margins=length(fitted_margins)
for (i in 1:num_margins) {
  log_lik_list=c(log_lik_list, logLik(fitted_margins[[i]]))
}
for (i in 1:(num_margins-1)) {
  log_lik_list=c(log_lik_list, fitted_copulas[[i]]$logLik)
}

log_lik_list=c(log_lik_list,sum(log_lik_list))

names(log_lik_list)<- c(1:num_margins,paste(1:(num_margins-1),(1:(num_margins-1))+1,sep=","),"Total")


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
    names(mu[[i]])=paste(names(mu[[i]]),i,sep=".")
    names(sigma[[i]])=paste(names(sigma[[i]]),i,sep=".")
    names(nu[[i]])=paste(names(nu[[i]]),i,sep=".")
    names(tau[[i]])=paste(names(tau[[i]]),i,sep=".")
  }  
  
  theta1=theta2=vector()
  for (i in 1:length(fitted_copulas)) {
    theta1[i]=fitted_copulas[[i]]$par
    theta2[i]=fitted_copulas[[i]]$par2
    names(theta1)[i]=paste("theta1",names(fitted_copulas)[i])
    names(theta2)[i]=paste("theta2",names(fitted_copulas)[i])
  }
  
  return(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau),theta1,theta2))
}

start_par<-extract_parameters_from_fitted(fitted_margins,fitted_copulas)

####Calculate Log Likelihood ####

calc_joint_likelihood <- function(input_par,start_fitted_margins,margin_dist,copula_dist,return_option="list")  {
  #input_par=start_par
  tryCatch(
  {
  
  if(return_option=="list") {
    print("Starting parameters:")
    print(as.matrix(start_par))  
  }

  #### Get linear predictors and their linked value -> **Creates eta, eta_link_inv** ####
  eta = eta_linkinv = list()
  num_margins = length(start_fitted_margins)
  observation_count = length(start_fitted_margins[[1]]$y)
  
  for (margin_number in 1:num_margins) {
    model_matrices = matrix(nrow = observation_count, ncol = 5)
    model_matrices_inv = matrix(nrow = observation_count, ncol = 5)
    
    counter = 1
    for (parameter in c("mu", "sigma", "nu", "tau")) {
      par_temp <- input_par[grepl(parameter, names(input_par))]
      par_temp <- par_temp[((margin_number - 1) * (length(par_temp) / num_margins) +
                              1):(margin_number * (length(par_temp) / num_margins))]
      model_matrices[, counter] = model.matrix(start_fitted_margins[[margin_number]], what =
                                                 parameter)    %*% par_temp #These are the linear predictors
      counter = counter + 1
    }
    model_matrices[, 5] = dataset[dataset$time == margin_number, "response"]
    colnames(model_matrices) = c("mu", "sigma", "nu", "tau", "x")
    eta[[margin_number]] = model_matrices
    
    model_matrices_inv[, 1] = margin_dist$mu.linkinv(model_matrices[, 1])
    model_matrices_inv[, 2] = margin_dist$sigma.linkinv(model_matrices[, 2])
    model_matrices_inv[, 3] = margin_dist$nu.linkinv(model_matrices[, 3])
    model_matrices_inv[, 4] = margin_dist$tau.linkinv(model_matrices[, 4])
    model_matrices_inv[, 5] = dataset[dataset$time == margin_number, "response"]
    
    colnames(model_matrices_inv) = c("mu", "sigma", "nu", "tau", "x")
    
    eta_linkinv[[margin_number]] = model_matrices_inv
    
  }
  
  #### Get density and probability functions for the values from model -> ** Creates margin_p and margin_d ####
  
  margin_d=margin_p=matrix(nrow=observation_count,ncol=num_margins)
  
  for (margin_number in 1:num_margins) {
    
    margin_d[,margin_number] = dZISICHEL(eta_linkinv[[margin_number]][,5]
                                         ,mu=    eta_linkinv[[margin_number]][,"mu"]
                                         ,sigma= eta_linkinv[[margin_number]][,"sigma"]
                                         ,nu=    eta_linkinv[[margin_number]][,"nu"]
                                         ,tau=   eta_linkinv[[margin_number]][,"tau"])
    
    margin_p[,margin_number] = pZISICHEL(eta_linkinv[[margin_number]][,5]
                                         ,mu=    eta_linkinv[[margin_number]][,"mu"]
                                         ,sigma= eta_linkinv[[margin_number]][,"sigma"]
                                         ,nu=    eta_linkinv[[margin_number]][,"nu"]
                                         ,tau=   eta_linkinv[[margin_number]][,"tau"])
    
  }
  
  
  for (margin_number in 1:num_margins) {
    if(sum(log(margin_d[,margin_number]))-logLik(fitted_margins[[margin_number]])>1) {print("ERROR: margin_d does not agree with original gamlss LogLik - something has gone wrong!!!")}
  }
  
  #num_margins=length(fitted_margins)
  #observation_count=nrow(dataset[dataset$time==1,])
  
  #### Get density of the copula function for each pair of margins ####
  theta1=input_par[grepl("theta1", names(input_par))]
  theta2=input_par[grepl("theta2", names(input_par))]
  
  copula_d=dldth=d2ldth=matrix(nrow=observation_count,ncol=num_margins-1)
  for (margin_number in 1:(num_margins-1)) {
    copula_d[,margin_number]<-BiCopPDF( margin_p[,margin_number],margin_p[,margin_number+1],family = BiCopName(copula_dist),par=theta1[i],par2=theta2[i])
    dldth[,margin_number]=BiCopDeriv(   margin_p[,margin_number],margin_p[,margin_number+1],family = BiCopName(copula_dist),par=theta1[i],par2=theta2[i],deriv="par",log=TRUE)
    d2ldth[,margin_number]=BiCopDeriv2( margin_p[,margin_number],margin_p[,margin_number+1],family = BiCopName(copula_dist),par=theta1[i],par2=theta2[i],deriv="par")
  }
  colnames(copula_d)=colnames(dldth)=colnames(d2ldth)=paste(1:(num_margins-1),",",2:(num_margins),sep="")
  
  #### Return Log Likelihood ####
  
  log_lik_results<-c(colSums(log(margin_d)),colSums(log(copula_d)),sum(colSums(log(margin_d)),colSums(log(copula_d))))
  names(log_lik_results)[1:num_margins]=1:num_margins
  names(log_lik_results)[length(log_lik_results)]="Total"
  
  
  margin_dldm = matrix(0,ncol=num_margins,nrow=observation_count)
  all_derivatives=c(names(margin_dist)[grepl("dl", names(margin_dist))],names(margin_dist)[grepl("d2l", names(margin_dist))])
  
  derivatives_calculated_all_margins=list()
  
  for (margin_number in 1:(num_margins)) {
    i=1
    derivatives_calculated=matrix(nrow=observation_count,ncol=length((all_derivatives)))
    for (d_fun in all_derivatives) {
      derivatives_calculated[,i]=margin_dist[[d_fun]](y=eta_linkinv[[margin_number]][,"x"],mu=eta_linkinv[[margin_number]][,"mu"],sigma=eta_linkinv[[margin_number]][,"sigma"],nu=eta_linkinv[[margin_number]][,"nu"],tau=eta_linkinv[[margin_number]][,"tau"])
      i=i+1
    }
    colnames(derivatives_calculated)=all_derivatives
    derivatives_calculated_all_margins[[margin_number]]=derivatives_calculated
  }
  
  return_list=list(log_lik_results,margin_d,margin_p,copula_d,dldth,d2ldth,margin_dist,eta,eta_linkinv,derivatives_calculated_all_margins)
  names(return_list)=c("log_lik_results","margin_d","margin_p","copula_d","dldth","d2ldth","margin_dist","eta","eta_linkinv","derivatives_calculated_all_margins")
  
  if (return_option == "list") {
    print("Log Likelihoods:")
    print(log_lik_results)
    return(return_list)}
  if (return_option == "log_lik") {return(as.numeric(log_lik_results["Total"]))}
  }, 
  error = function(err) {
    print("Function error")
    return(NA)
  }
  )
}
  
results= calc_joint_likelihood(input_par = start_par
                      ,start_fitted_margins = fitted_margins
                      ,margin_dist = ZISICHEL(
                        mu.link = "log",
                        sigma.link = "log",
                        nu.link = "identity",
                        tau.link = "logit"
                      ),
                      copula_dist="N",
                      return_option="list"
                      )

results$log_lik_results # LogLik=12980.0415
-results$log_lik_results["Total"]*2+2*length(start_par) #AIC= 26,004 | 

optim_results=optim(start_par
      ,calc_joint_likelihood
      ,start_fitted_margins = fitted_margins
      ,margin_dist = ZISICHEL(
        mu.link = "log",
        sigma.link = "log",
        nu.link = "identity",
        tau.link = "logit"
      )
      , copula_dist="N"
      , return_option="log_lik"
      , control=list(fnscale=-1,trace=3)
      , hessian=TRUE
      )

#2. Be able to iterate parameters in direction of optimal

#### Appendix ###












