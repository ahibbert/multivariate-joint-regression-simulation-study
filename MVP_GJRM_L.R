#The goal of this code is to develop a prototype for a 'GJRM-L' function which allows a user to specify a model
#in the form of a univariate longitudinal model e.g. lme4, gamlss with random, and have it estimated as a 
#copula-based joint regression model.
#GJRM-L(mu.formula, sigma.formula,nu.formula,zeta.forula, tau.formula, other.formula,
#margin.family=list(), copula.family=list()) #All optional except mu.formula

#### 0. Starting parameters ####
set.seed(1000);options(scipen=999);
source("common_functions.R")
library(gamlss)
library(VineCopula)
#library(copula)

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

fit_copulas = function(margin_fits,copula_dist="N",method="vine",dataset) {
  require(gamlss)
  require(VineCopula)
  
  fitted_margins_length=length(margin_fits)
  
  num_margins=length(unique(dataset$time))
  observation_count=nrow(dataset[dataset$time==unique(dataset$time)[1],])
  
  fit_unif=matrix(0,nrow=observation_count,ncol=num_margins)
  
  if(fitted_margins_length==1) {
    #margin_fits=fitted_margins #for testing
    for (i in unique(dataset$time)) {
      print(i)
      fit_unif[,i]=pZISICHEL(fitted_margins[[1]]$y[dataset$time==i]
                ,mu=    margin_fits[[1]]$mu.fv[dataset$time==i]
                ,sigma= margin_fits[[1]]$sigma.fv[dataset$time==i]
                ,nu=    margin_fits[[1]]$nu.fv[dataset$time==i]
                ,tau=   margin_fits[[1]]$tau.fv[dataset$time==i])
          }
  } else {
    
    for (i in 1:num_margins) {
      fit_unif[,i]=pNO(margin_fits[[i]]$residuals)    
    }
    
  }
  
  if (method=="novine") {
    cop_fit=list()  
    for (i in 1:(num_margins-1)){
      cop_fit[[i]]=BiCopEst(fit_unif[,i],fit_unif[,i+1],family=as.numeric(BiCopName(copula_dist)),se=TRUE)
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

logit <- function(x) {
  return(log(x/(1-x)))
}

logit_inv <- function(x) {
  return(
    exp(x)/(1+exp(x))
  )
}

log_2plus <- function(x) {
  return(
    log(x-2)
  )
}

log_2plus_inv <- function(x) {
  y=exp(x)+2
  #Adjust for error when close to two.
  y[y==2]=y[y==2]+0.00001
  return(y)
}

dlogit <- function(x) {
  return(1/(x-(x^2)))
}

dlog_2plus <- function(x) {
  return(1/(x-2))
}

calc_joint_likelihood <- function(input_par,start_fitted_margins,margin_dist,copula_dist,return_option="list",copula_link,dataset,verbose=TRUE)  {
  
  #return_option="list";start_fitted_margins=fitted_margins;input_par=start_par;copula_dist="t";margin_dist = ZISICHEL(mu.link = "log",sigma.link = "log",nu.link = "identity",tau.link = "logit");
  
  #tryCatch(
  #  {
      
      if(return_option=="list" & verbose==TRUE) {
        print("Starting parameters:")
        print(as.matrix(input_par))  
      }
      
      #### Get linear predictors and their linked value -> **Creates eta, eta_link_inv** ####
      eta = eta_linkinv = list()
      length_fitted_margins=length(start_fitted_margins)
      margin_names=unique(dataset$time)
      num_margins = length(margin_names)
      observation_count = nrow(dataset[dataset$time==margin_names[1],])
      
      mm_all=list()
      
      for (margin_number in 1:length_fitted_margins) {
        
        mm_all[[margin_number]]=list()  
        #########THIS IS WHERE I'M UP TO 
        
        if (length_fitted_margins==1) {
          model_matrices = matrix(nrow = nrow(dataset), ncol = 5)
          model_matrices_inv = matrix(nrow = nrow(dataset), ncol = 5)
        } else {
          model_matrices = matrix(nrow = observation_count, ncol = 5)
          model_matrices_inv = matrix(nrow = observation_count, ncol = 5)
        }
        colnames(model_matrices)=colnames(model_matrices_inv)= c("mu", "sigma", "nu", "tau", "response")
        
        for (parameter in c("mu", "sigma", "nu", "tau")) {
          par_temp <- input_par[grepl(paste(parameter,margin_number,sep="."), names(input_par))]
          #par_temp <- par_temp[((margin_number - 1) * (length(par_temp) / num_margins) +
          #                        1):(margin_number * (length(par_temp) / num_margins))]
          
          mm=model.matrix(start_fitted_margins[[margin_number]], what = parameter)
          mm_all[[margin_number]][[parameter]]=mm
          
          model_matrices[, parameter] = rowSums(mm * matrix(rep(par_temp,each=nrow(mm)),ncol=length(par_temp),dimnames=list(NULL,c(names(par_temp))))) #These are the linear predictors
          
        }
        model_matrices[, "response"] = dataset[dataset$time == margin_number, "response"]
        
        eta[[margin_number]] = model_matrices
        
        model_matrices_inv[, "mu"] = margin_dist$mu.linkinv(model_matrices[, "mu"])
        model_matrices_inv[, "sigma"] = margin_dist$sigma.linkinv(model_matrices[, "sigma"])
        model_matrices_inv[, "nu"] = margin_dist$nu.linkinv(model_matrices[, "nu"])
        model_matrices_inv[, "tau"] = margin_dist$tau.linkinv(model_matrices[, "tau"])
        
        if (length_fitted_margins==1) {
          model_matrices_inv[, "response"] = dataset[, "response"]
        } else 
        {
          model_matrices_inv[, "response"] = dataset[dataset$time == margin_number, "response"]
        }
        
        eta_linkinv[[margin_number]] = model_matrices_inv
        
      }
      
      eta_nomargin=eta
      eta_linkinv_nomargin=eta_linkinv
      
      eta_temp=eta_linkinv_temp=list()
      if(length_fitted_margins==1) {
        for (i in margin_names) {
          eta_temp[[i]]         =eta[[1]][dataset$time==i,]
          eta_linkinv_temp[[i]] =eta_linkinv[[1]][dataset$time==i,]
        }
      }
      eta=eta_temp
      eta_linkinv=eta_linkinv_temp
      
      #### Get density and probability functions for the values from model -> ** Creates margin_p and margin_d ####
      
      margin_d=margin_p=matrix(nrow=observation_count,ncol=num_margins)
      
      for (margin_number in 1:num_margins) {
        
        margin_d[,margin_number] = dZISICHEL(eta_linkinv[[margin_number]][,"response"]
                                             ,mu=    eta_linkinv[[margin_number]][,"mu"]
                                             ,sigma= eta_linkinv[[margin_number]][,"sigma"]
                                             ,nu=    eta_linkinv[[margin_number]][,"nu"]
                                             ,tau=   eta_linkinv[[margin_number]][,"tau"])
        
        margin_p[,margin_number] = pZISICHEL(eta_linkinv[[margin_number]][,"response"]
                                             ,mu=    eta_linkinv[[margin_number]][,"mu"]
                                             ,sigma= eta_linkinv[[margin_number]][,"sigma"]
                                             ,nu=    eta_linkinv[[margin_number]][,"nu"]
                                             ,tau=   eta_linkinv[[margin_number]][,"tau"])
        
      }
      
      if(length_fitted_margins>1) {
        for (margin_number in 1:num_margins) {
          if(sum(log(margin_d[,margin_number]))-logLik(start_fitted_margins[[margin_number]])>1) {print("ERROR: margin_d does not agree with original gamlss LogLik - something has gone wrong!!!")}
        }
      } else {
        if(sum(log(margin_d))-logLik(start_fitted_margins[[1]])>1) {print("ERROR: margin_d does not agree with original gamlss LogLik - something has gone wrong!!!")}
      }
      
      #### Get density of the copula function for each pair of margins ####
      theta=copula_link$theta.linkinv(as.vector(input_par[grepl("theta", names(input_par))]))
      zeta=copula_link$zeta.linkinv(as.vector(input_par[grepl("zeta", names(input_par))])) ###These are correct
      
      copula_d=dcdu1=dcdu2=dldth=d2ldth=dldth2=d2ldth2=dthdeta=dzdeta=matrix(0,nrow=observation_count,ncol=num_margins-1)
      for (margin_number in 1:(num_margins-1)) {
        eta[[margin_number]]=cbind(cbind(eta[[margin_number]],rep(copula_link$theta.linkfun(theta[[margin_number]]),nrow(eta[[margin_number]]))),rep(copula_link$zeta.linkfun(zeta[[margin_number]]),nrow(eta[[margin_number]])))  
        colnames(eta[[margin_number]])[(length(colnames(eta[[margin_number]]))-1):length(colnames(eta[[margin_number]]))]=c("theta","zeta")
        
        copula_d[,margin_number]=BiCopPDF( margin_p[,margin_number],margin_p[,margin_number+1],family = as.vector(BiCopName(copula_dist)),par=theta[margin_number],par2=zeta[margin_number])
        dldth[,margin_number]=BiCopDeriv(   margin_p[,margin_number],margin_p[,margin_number+1],family = as.vector(BiCopName(copula_dist)),par=theta[margin_number],par2=zeta[margin_number],deriv="par",log=TRUE)
        dcdu1[,margin_number]=BiCopDeriv(   margin_p[,margin_number],margin_p[,margin_number+1],family = as.vector(BiCopName(copula_dist)),par=theta[margin_number],par2=zeta[margin_number],deriv="u1",log=FALSE)
        dcdu2[,margin_number]=BiCopDeriv(   margin_p[,margin_number],margin_p[,margin_number+1],family = as.vector(BiCopName(copula_dist)),par=theta[margin_number],par2=zeta[margin_number],deriv="u2",log=FALSE)
        d2ldth[,margin_number]=BiCopDeriv2( margin_p[,margin_number],margin_p[,margin_number+1],family = as.vector(BiCopName(copula_dist)),par=theta[margin_number],par2=zeta[margin_number],deriv="par")
  
        dthdeta[,margin_number]=copula_link$dthdeta(eta[[margin_number]][,"theta"])
        dzdeta[,margin_number]=copula_link$dzdeta(eta[[margin_number]][,"zeta"])
        
        if (copula_dist=="t") {
          dldth2[,margin_number]=BiCopDeriv(   margin_p[,margin_number],margin_p[,margin_number+1],family = as.vector(BiCopName(copula_dist)),par=theta[margin_number],par2=zeta[margin_number],deriv="par2",log=TRUE)
          d2ldth2[,margin_number]=BiCopDeriv2(   margin_p[,margin_number],margin_p[,margin_number+1],family = as.vector(BiCopName(copula_dist)),par=theta[margin_number],par2=zeta[margin_number],deriv="par2")
        }
        
        
      }
      colnames(copula_d)=colnames(dthdeta)=colnames(dzdeta)=colnames(dldth)=colnames(d2ldth)=colnames(dldth2)=colnames(d2ldth2)=paste(1:(num_margins-1),",",2:(num_margins),sep="")
      
      #### Return Log Likelihood ####
      
      log_lik_results<-c(colSums(log(margin_d)),colSums(log(copula_d)),sum(colSums(log(margin_d)),colSums(log(copula_d))))
      names(log_lik_results)[1:num_margins]=1:num_margins
      names(log_lik_results)[length(log_lik_results)]="Total"
      
      #### Extract derivatives of margins ####
      all_derivatives=c(names(margin_dist)[grepl("dl", names(margin_dist))],names(margin_dist)[grepl("d2l", names(margin_dist))],names(margin_dist)[grepl("dr", names(margin_dist))])
      derivatives_calculated_all_margins=list()
      
      #Margin derivatives
      
      for (margin_number in 1:(length_fitted_margins)) {
        i=1 #DON'T DELETE THIS - it's an internal counter
        derivatives_calculated=matrix(nrow=(if(length_fitted_margins==1){nrow(dataset)}else{observation_count}),ncol=length((all_derivatives)))
        for (d_fun in all_derivatives) { 
          if (grepl("dr", d_fun) == TRUE) {
            derivatives_calculated[,i]=margin_dist[[d_fun]](eta_nomargin[[margin_number]][,grepl(sub(".dr","",d_fun), colnames(eta_nomargin[[margin_number]]))])
          } else {
            derivatives_calculated[,i]=margin_dist[[d_fun]](y=eta_linkinv_nomargin[[margin_number]][,"response"],mu=eta_linkinv_nomargin[[margin_number]][,"mu"],sigma=eta_linkinv_nomargin[[margin_number]][,"sigma"],nu=eta_linkinv_nomargin[[margin_number]][,"nu"],tau=eta_linkinv_nomargin[[margin_number]][,"tau"])
          }
          i=i+1
        }
        colnames(derivatives_calculated)=all_derivatives
        derivatives_calculated_all_margins[[margin_number]]=derivatives_calculated
      }
      
      return_list=list(log_lik_results,margin_d,margin_p,copula_d,dldth,dldth2,d2ldth,d2ldth2,margin_dist,copula_dist,eta,eta_linkinv,derivatives_calculated_all_margins,eta_nomargin,eta_linkinv_nomargin,dthdeta,dzdeta,mm_all,start_par,dcdu1,dcdu2)
      names(return_list)=c("log_lik_results","margin_d","margin_p","copula_d","dldth","dldth2","d2ldth","d2ldth2","margin_dist","copula_dist","eta","eta_linkinv","derivatives_calculated_all_margins","eta_nomargin","eta_linkinv_nomargin","dthdeta","dzdeta","mm_all","start_par","dcdu1","dcdu2")
      
      if (return_option == "list") {
        if (verbose==TRUE) {
          print("Log Likelihoods:")
          print(log_lik_results)
        }
        return(return_list)}
      if (return_option == "log_lik") {return(as.numeric(log_lik_results["Total"]))}
      if (return_option == "log_lik_cop") {return(sum(log_lik_results[grepl(",",names(log_lik_results))]))}
    #}, 
    #error = function(err) {
    #  print("Function error")
    #  return(NA)
    #}
  #)
}

extract_parameters_from_fitted <- function(fitted_margins,fitted_copulas,copula_link) {
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
  
  theta=zeta=vector()
  for (i in 1:length(fitted_copulas)) {
    theta[i]=copula_link$theta.linkfun(fitted_copulas[[i]]$par)
    zeta[i]=copula_link$zeta.linkfun(fitted_copulas[[i]]$par2)
    names(theta)[i]=paste("theta",names(fitted_copulas)[i])
    names(zeta)[i]=paste("zeta",names(fitted_copulas)[i])
  }
  
  return(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau),theta,zeta))
}

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

gamlss_model <- gamlss(formula = response ~ as.factor(time)+as.factor(gender)+age
                       , sigma.formula = ~ as.factor(time)+age
                       , nu.formula = ~ as.factor(time)+as.factor(gender)
                       , tau.formula = ~ age
                       , family="ZISICHEL",data=dataset,method=RS(50))
summary(gamlss_model)
plot(gamlss_model)
#term.plot(gamlss_model)
#term.plot(gamlss_model,what="sigma")
#term.plot(gamlss_model,what="nu")
#term.plot(gamlss_model,what="tau")

###### GLMM ZISICHEL logLik(gamlss_glmm_model) = -11187.34 (df=1042.17) AIC: 24459.02 | BIC: 31322.54
#margin_model_formula_glmm=formula(response~(age)+as.factor(gender)+as.factor(time)+random(as.factor(subject)))
#gamlss_glmm_model <- gamlss(formula = margin_model_formula_glmm,family="ZISICHEL",data=dataset) #So this doesn't run automatically...
#summary(gamlss_glmm_model)
#plot(gamlss_glmm_model)
#term.plot(gamlss_glmm_model)

#### 3. Separate margin and copula models ####

#Goal

fitted_margins=list()
fitted_margins[[1]]=gamlss_model

# Starting fit for margins and copula
#best_marginal_fits <- find_best_marginal_fits(dataset,type="counts") #Find best fit margins
#fitted_margins<-fit_margins(mu.formula=formula(response~age+as.factor(gender)),dataset=dataset,margin_dist=ZISICHEL())
  #AIC(fitted_margins[[1]])
  #plot(fitted_margins[[1]])
  #term.plot(fitted_margins[[2]])

fitted_copulas<-fit_copulas(fitted_margins,copula_dist="t",method="novine",dataset=dataset)
#Extract parameters from starting fits
copula_link=list(logit,logit_inv,log_2plus,log_2plus_inv,dlogit,dlog_2plus)
names(copula_link)=c("theta.linkfun","theta.linkinv","zeta.linkfun","zeta.linkinv","dthdeta","dzdeta")

#fitted_copulas<-fit_copulas(fitted_margins,copula_dist="N",method="vine")
#contour(fitted_copulas)
#summary(fitted_copulas)
#fitted_copulas$pair.logLik
#likelihood_fitted_margin_copula(fitted_margins,fitted_copulas) ###Not much better results with vine versus optimising non-vine
#contour(fitted_copulas)

log_lik_list=vector()
num_margins=length(unique(dataset$time))
length_fitted_margins=length(fitted_margins)
if(length_fitted_margins==1) {log_lik_list=logLik(fitted_margins[[1]])} else {
  for (i in 1:num_margins) {
    log_lik_list=c(log_lik_list, logLik(fitted_margins[[i]]))
  }
}
for (i in 1:(num_margins-1)) {
  log_lik_list=c(log_lik_list, fitted_copulas[[i]]$logLik)
}

log_lik_list=c(log_lik_list,sum(log_lik_list))

if(length_fitted_margins==1) {
  names(log_lik_list)<- c("gamlss_all_margin",paste(1:(num_margins-1),(1:(num_margins-1))+1,sep=","),"Total")
} else { 
  names(log_lik_list)<- c(1:num_margins,paste(1:(num_margins-1),(1:(num_margins-1))+1,sep=","),"Total")
}
print(log_lik_list)


#### 4. Log likelihood calculation + optimisation with optim ####

results= calc_joint_likelihood(input_par =end_par_matrix[nrow(end_par_matrix),]
                      ,start_fitted_margins = fitted_margins
                      ,margin_dist = ZISICHEL(
                        mu.link = "log",
                        sigma.link = "log",
                        nu.link = "identity",
                        tau.link = "logit"
                      )
                      , copula_dist="t"
                      , copula_link=copula_link
                      , return_option="list"
                      , dataset=dataset
                      )

observation_count=nrow(dataset[dataset$time==1,])
df=length(start_par)+start_par["zeta 1,2"]+start_par["zeta 2,3"]-2
results_df=c(results$log_lik_results,-results$log_lik_results["Total"]*2+2*df,-results$log_lik_results["Total"]*2+log(observation_count*num_margins)*df) #AIC= 26,004 | 
names(results_df)=c(names(results_df[1:(length(results_df)-3)]),"LogLik","AIC","BIC")
print(results_df)

optim_results=optim(start_par
      ,calc_joint_likelihood
      ,start_fitted_margins = fitted_margins
      ,margin_dist = ZISICHEL(
        mu.link = "log",
        sigma.link = "log",
        nu.link = "identity",
        tau.link = "logit"
      )
      , copula_dist="t"
      , copula_link=copula_link
      , dataset=dataset
      , return_option="log_lik"
      , control=list(fnscale=-1,trace=3)
      , hessian=TRUE
      )

library(numDeriv)
grad_l<-grad(calc_joint_likelihood,
              , x=start_par
             ,method = "simple"
               ,start_fitted_margins = fitted_margins
               ,margin_dist = ZISICHEL(
                 mu.link = "log",
                 sigma.link = "log",
                 nu.link = "identity",
                 tau.link = "logit"
               )
               , copula_dist="t"
             , copula_link=copula_link
               , return_option="log_lik_cop"
             ,dataset=dataset
)
grad_l

help(grad)

hessian_l<-hessian(calc_joint_likelihood,
             , x=start_par
             ,start_fitted_margins = fitted_margins
             ,margin_dist = ZISICHEL(
               mu.link = "log",
               sigma.link = "log",
               nu.link = "identity",
               tau.link = "logit"
             )
             , copula_dist="t"
             , copula_link=copula_link
             , return_option="log_lik"
)
hessian_l
sqrt(diag(solve(-hessian)))

#### Extract parameters and SE from separate and jointly optimised datasets -> all_out_combined ####
library("MASS")
optim_par_results<-cbind(data.frame(start_par),data.frame(optim_results$par),data.frame(diag(sqrt(ginv(-optim_results$hessian)))))
colnames(optim_par_results)<-c("separate","joint","joint se")
print(optim_par_results)

all_out=matrix(ncol=2,nrow=0)
for (i in 1:length_fitted_margins) {
  out=summary(fitted_margins[[i]])
  rownames(out)=sub(".", paste(".",i,".",sep=""), labels(unlist(coefAll(fitted_margins[[i]]))),fixed=TRUE)
  all_out=rbind(all_out,out[,c(1,2)])
}
all_out

all_out_sorted=matrix(ncol=2,nrow=0)
for (parameter in c("mu","sigma","nu","tau")) {
  all_out_sorted=rbind(all_out_sorted,all_out[grepl(parameter, rownames(all_out)),])
}
all_out_sorted

for (j in 1:2) {
  for (i in 1:length(fitted_copulas)) {
    out=(fitted_copulas[[i]])
    all_out_sorted=rbind(all_out_sorted,c(as.numeric(out[paste("par",if(j==1){""}else{j},sep="")]),as.numeric(out[paste("se",if(j==1){""}else{j},sep="")])))
    }
}

all_out_combined<-cbind(optim_par_results,all_out_sorted)[,c(4,5,2,3)]
z_stat=all_out_combined[,1]/all_out_combined[,2]
p_vals1=round(1-pNO(abs(z_stat)),5)
z_stat=all_out_combined[,3]/all_out_combined[,4]
p_vals2=round(1-pNO(abs(z_stat)),5)
all_out_combined=cbind(all_out_combined[,c(1,2)],p_vals1,all_out_combined[,c(3,4)],p_vals2)
colnames(all_out_combined)<-c("Sep. Est","Sep. SE","p-value","Joint Est.","Joint SE","p-value")
round(all_out_combined,3)

#### 5. Newton Raphson optimisation ####

#### 5.1 Inner iteration ####

####TAKES results as input

#Backfitting

######STARTING PARAMETERS 
score_function <- function(eta,dldpar,dpardeta,dlcopdpar,phi=1,verbose=TRUE) {
  #u_k
  dldpar=dldpar#+dlcopdpar commented out for now until we come back and fix it
  
  u_k=dldeta = dldpar * dpardeta
  
  #f_k and w_k
  #f_k_quasi_newton=(-(dldpar*dldpar))
  f_k=-dldpar*dldpar
  
  w_k=-f_k*(dpardeta*dpardeta) ###OK so my dpardeta is right.... and my dldpar....
  
  #w_k=matrix(1,ncol=ncol(f_k),nrow=nrow(f_k)) ####Weights of 1
  
  z_k=(1-phi)*eta+phi*(eta+(u_k/w_k))
  
  if(verbose==TRUE) {
    steps_mean=round(rbind(colMeans(eta)
                           ,colMeans(dlcopdpar)
                           ,colMeans(dpardeta)
                           ,colMeans(dpardeta*dpardeta)
                           ,colMeans(f_k)
                           ,colMeans(w_k)
                           ,colMeans(u_k)
                           ,colMeans((1/w_k)*u_k)
                           ,colMeans(z_k)
    ),8)
    rownames(steps_mean)=c("eta","dlcopdpar","dpardeta","dpardeta2","f_k","w_k","u_k","(1/w_k)*u_k","z_k")
    print(steps_mean)  
  }
  return_list=list(colMeans(z_k),u_k,f_k,w_k,z_k)
  names(return_list)=c("par","u_k","f_k","w_k","z_k")
  return(return_list)
}

newton_raphson_iteration=function(results,input_par,phi=1,verbose=c(FALSE,FALSE)) {
  margin_dist=results$margin_dist
  copula_dist=results$copula_dist
  
  steps_all=list()
  par_out=list()
  
  ####For copula parameters####
  copula_score=list()
  dlcopdpar=list()
  for (i in 1:ncol(results$dldth)) {
    dldth=results$dldth[,i]
    dldth2=results$dldth2[,i]
    d2ldth=results$d2ldth[,i] ########THESE ARE NOT LOG LIKELIHOOD FUNCTIONS they are f''(x)
    d2ldth2=results$d2ldth2[,i]
    dthdeta=results$dthdeta[,i]
    dzdeta=results$dzdeta[,i]
    dcdu1=results$dzdeta[,i]
    dcdu2=results$dzdeta[,i]
    
    theta_par_temp=start_par[names(start_par)[grepl(paste(i,i+1,sep=","),names(start_par))]]
    
    dldpar=cbind(dldth,dldth2)
    d2ldpar=cbind(d2ldth,d2ldth2)
    dpardeta=cbind(dthdeta,dzdeta)
    dpardeta=dldpar*0+1 ##########Temporary step
    
    dlcopdpar[[i]]=dcdu1*dcdu2*(results$margin_d[,i]*results$margin_d[,i+1])/results$copula_d[,i]
    copula_score[[i]]=score_function(eta=results$eta[[i]][,c("theta","zeta")],dldpar=dldpar,dpardeta=dpardeta,dlcopdpar=dldpar*0,phi=phi,verbose=verbose[2]) ##Returns updated value
  } 
  
  beta_new_cop=list()
  parameter_names=c("theta","zeta")
  
  #####TEMPORARY UNTIL WE FEED COP MM THROUGH RESULTS LIST
  mm_cop=list()
  for (i in 1:length(copula_score)){
    mm_cop[[i]]=list()
    mm_cop[[i]]$theta=rep(1,nrow(copula_score[[i]]$z_k))
    mm_cop[[i]]$zeta=rep(1,nrow(copula_score[[i]]$z_k))
  }
  
  #-results$eta_nomargin[[1]][,c("mu","sigma","nu","tau")]
  for (i in 1:length(copula_score)){
    e_k=copula_score[[i]]$z_k  
    beta_new_cop[[i]]=list()
    for (j in 1:ncol(copula_score[[i]]$z_k)) {
      mm=mm_cop[[i]][[j]]
      X=as.matrix(mm)
      W=diag(copula_score[[i]]$w_k[,j])
      e_k_par=e_k[,j]
      beta_new_cop[[i]][[parameter_names[j]]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
      edge=paste(i,i+1,sep=",")
      rownames(beta_new_cop[[i]][[parameter_names[j]]])= paste(parameter_names[j],edge,sep=" ")
      #names(beta_new_cop[[j]])= paste(parameter_names[j],margin_number,rownames(beta_new_cop[[j]]),sep=".")
    }
    #beta_new_cop[[i]]
  }
  
  ####For gamlss parameters####
  for (i in 1:ncol(results$margin_d)) {
  
    if (i==1) {dlcopdpar_final=dlcopdpar[[i]]} ###
    else if (i==length(results$derivatives_calculated_all_margins)) {
      dlcopdpar_final=c(dlcopdpar_final,dlcopdpar[[i]]+dlcopdpar[[i-1]])
    } 
    else {
      dlcopdpar_final=c(dlcopdpar_final,dlcopdpar[[i-1]])
    }
  }
  
  dlcopdpar_final=cbind(dlcopdpar_final,dlcopdpar_final,dlcopdpar_final,dlcopdpar_final)
  
  for (i in 1:length(results$derivatives_calculated_all_margins)) {
    
    dldpar=results$derivatives_calculated_all_margins[[i]][,c("dldm","dldd","dldv","dldt")] #First derivative of log likelihood w.r.t. parameters
    d2ldpar=results$derivatives_calculated_all_margins[[i]][,c("d2ldm2","d2ldd2","d2ldv2","d2ldt2")] #Quasi newton second derivative of log function w.r.t parameters
    dpardeta=results$derivatives_calculated_all_margins[[i]][,c("mu.dr","sigma.dr","nu.dr","tau.dr")] #Derivative of inverse link function
    
  }
  
  margin_score=score_function(eta=results$eta_nomargin[[1]][,c("mu","sigma","nu","tau")],dldpar=dldpar,dpardeta=dpardeta,dlcopdpar=dlcopdpar_final,phi=phi,verbose=verbose[1]) ##Returns updated value
  
  
  parameter_names=c("mu","sigma","nu","tau")
  e_k=margin_score$z_k#-results$eta_nomargin[[1]][,c("mu","sigma","nu","tau")]
  beta_r=list()
  for (margin_number in 1:length(results$mm_all)) {
    beta_r[[margin_number]]=list()
    for (parameter_num in 1:length(results$mm_all[[margin_number]])) {
      
      mm=results$mm_all[[margin_number]][[parameter_num]]
      X=as.matrix(mm)
      W=diag(margin_score$w_k[,parameter_num])
      e_k_par=e_k[,parameter_num]
      beta_r[[margin_number]][[parameter_num]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
      rownames(beta_r[[margin_number]][[parameter_num]])= paste(parameter_names[parameter_num],margin_number,rownames(beta_r[[margin_number]][[parameter_num]]),sep=".")
    }
  }
  
  
  ################
  
  
  beta_new=matrix(ncol=1,nrow=0)
  colnames(beta_new)=c("Estimate")

  for (i in 1:length(beta_r)) {
    for (j in 1:length(beta_r[[i]]))
      beta_new=rbind(beta_new,beta_r[[i]][[j]])
  }
  for (i in 1:length(beta_new_cop)) {
    for (j in 1:length(beta_new_cop[[i]]))
      beta_new=rbind(beta_new,beta_new_cop[[i]][[j]])
  }
  
  
  
  end_par=as.vector(beta_new)
  names(end_par)=rownames(beta_new)
  end_par=end_par[names(start_par)]
  names(end_par)==names(start_par)
  
  return_list=list(input_par,end_par)
  names(return_list)=c("input_par","end_par")
  return(return_list)
}

########WORKING NEWTON RAPHSON ITERATIONS

#Starting parameters for iterations
start_par<-extract_parameters_from_fitted(fitted_margins,fitted_copulas,copula_link)
first_input_par=start_par ####

#Newton Raphson iterations
end_par_matrix=matrix(0,ncol=length(start_par),nrow=0)

first_run=TRUE
change=1
run_counter=1
while (abs(change) > .05 ) {
  if(!first_run) {start_log_lik=results$log_lik_results} else {input_par=first_input_par}
  #print(input_par)
  results= calc_joint_likelihood(input_par =input_par  #optim_results$par
                               ,start_fitted_margins = fitted_margins
                               ,margin_dist = ZISICHEL(
                                 mu.link = "log",
                                 sigma.link = "log",
                                 nu.link = "identity",
                                 tau.link = "logit"
                               )
                               , copula_dist="t"
                               , copula_link=copula_link
                               , return_option="list"
                               , dataset=dataset
                               ,verbose=TRUE
                                )

  end_log_lik=results$log_lik_results

  if(first_run) {change=1} else {
    change=end_log_lik["Total"]-start_log_lik["Total"]
    log_lik_output=c(start_log_lik["Total"],end_log_lik["Total"],change)
    names(log_lik_output) = c("Start LogLik","End LogLik","Change")
    print(log_lik_output)
  }
  
  iteration_out=newton_raphson_iteration(results,input_par,phi=.5,verbose=c(TRUE,TRUE))
  cbind(iteration_out[[1]],iteration_out[[2]])
  end_par=iteration_out[[2]]
  input_par=end_par
  first_run=FALSE
  end_par_matrix=rbind(end_par_matrix,end_par)
  run_counter=run_counter+1
  #Showing charts
  colnames(end_par_matrix)=names(end_par)
  pars=round(sqrt(ncol(end_par_matrix)),0)+1
  plot.new();par(mfrow=c(pars,pars-1))
  for (par_name in colnames(end_par_matrix)) {
    plot(1:nrow(end_par_matrix),end_par_matrix[,par_name],main=par_name)    
  }
}

####TESTING dpardeta for copula =1

c(mean(copula_link$theta.linkinv(end_par_matrix[,"theta 1,2"]))
,mean(copula_link$theta.linkinv(end_par_matrix[,"theta 2,3"]))
,mean(copula_link$zeta.linkinv(end_par_matrix[,"zeta 1,2"]))
,mean(copula_link$zeta.linkinv(end_par_matrix[,"zeta 2,3"])))

new_versus_original=cbind(start_par,end_par_matrix[nrow(end_par_matrix),])
colnames(new_versus_original)=c("Gamlss+VineCop","Manual")
print(new_versus_original)

#### Exploratory ####

#### Theta by gender ####

plotDist(dataset[dataset$gender==1,],dist="ZISICHEL")
fitted_margins<-fit_margins(mu.formula=formula(response~age),dataset=dataset[dataset$gender==1,],family=ZISICHEL())
#AIC(fitted_margins[[1]])
#plot(fitted_margins[[1]])
#term.plot(fitted_margins[[2]])

fitted_copulas<-fit_copulas(fitted_margins,copula_dist="t",method="vine")
summary(fitted_copulas)
contour(fitted_copulas)








