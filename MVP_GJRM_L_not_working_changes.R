#The goal of this code is to develop a prototype for a 'GJRM-L' function which allows a user to specify a model
#in the form of a univariate longitudinal model e.g. lme4, gamlss with random, and have it estimated as a 
#copula-based joint regression model.
#GJRM-L(mu.formula, sigma.formula,nu.formula,zeta.forula, tau.formula, other.formula,
#margin.family=list(), copula.family=list()) #All optional except mu.formula

#### 0. Load packages and set environment options ####
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
  
  dataset=dataset[order(dataset$subject,dataset$time),]
  
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

fit_copula_single = function(dataset,mm_mar,input_par,dFUN=dZISICHEL,pFUN=pZISICHEL,copula_dist="t") {
  #### MARGINAL MODEL: Get linear predictors and their inverse -> **Creates eta, eta_link_inv** from **mm_mar** ####
  margin_names=unique(dataset$time)
  num_margins = length(margin_names)
  #observation_count = nrow(dataset[dataset$time==margin_names[1],])
  response=dataset$response
  margin_parameters=names(mm_mar)
  
  #mm_all=list()  
  eta = eta_linkinv = matrix(nrow = nrow(dataset), ncol = length(margin_parameters))
  
  colnames(eta)=colnames(eta_linkinv)= margin_parameters
  
  for (parameter in margin_parameters) {
    par_temp <- input_par[grepl(paste(parameter,sep="."), names(input_par))]
    mm=mm_mar[[parameter]]
    #mm_all[[parameter]]=mm
    eta[, parameter] = rowSums(mm * matrix(rep(par_temp,each=nrow(mm)),ncol=length(par_temp),dimnames=list(NULL,c(names(par_temp))))) #These are the linear predictors
    FUN=margin_dist[[names(margin_dist)[grepl(paste(parameter,"linkinv",sep="."),names(margin_dist))]]]
    eta_linkinv[, parameter]=FUN(eta[,parameter])
  }
  
  #### MARGINAL MODEL: Get density and probability functions for margins -> ** Creates margin_p and margin_d ####
  
  margin_d = dFUN(x =   response
                  ,mu=    eta_linkinv[,"mu"]
                  ,sigma= eta_linkinv[,"sigma"]
                  ,nu=    eta_linkinv[,"nu"]
                  ,tau=   eta_linkinv[,"tau"])
  
  margin_p = pFUN(      response
                        ,mu=    eta_linkinv[,"mu"]
                        ,sigma= eta_linkinv[,"sigma"]
                        ,nu=    eta_linkinv[,"nu"]
                        ,tau=   eta_linkinv[,"tau"])
  
  margin_p_cop_input=cbind(margin_p[dataset$time %in% margin_names[1:(num_margins-1)]],margin_p[dataset$time %in% margin_names[2:(num_margins)]])
  
  cop_fit=BiCopEst(margin_p_cop_input[,1],margin_p_cop_input[,2],family=as.numeric(BiCopName(copula_dist)),se=TRUE)
  
  return(cop_fit)
  ####
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
  return(-1/(x-2))
}

calc_joint_likelihood <- function(input_par,mm_mar,margin_dist,copula_dist,return_option="list",copula_link,dataset,mm_cop,verbose=TRUE,calc_d2=FALSE,dFUN=dZISICHEL,pFUN=pZISICHEL)  {
  
  #return_option="list";start_fitted_margins=fitted_margins;input_par=start_par;copula_dist="t";margin_dist = ZISICHEL(mu.link = "log",sigma.link = "log",nu.link = "identity",tau.link = "logit");verbose=TRUE;dFUN=dZISICHEL;pFUN=pZISICHEL;calc_d2=FALSE
  #input_par =input_par;mm_mar = mm_mar;margin_dist = ZISICHEL(mu.link = "log",sigma.link = "log",    nu.link = "identity",    tau.link = "logit"); copula_dist="t"; copula_link=copula_link; mm_cop = mm_cop; return_option="list"; dataset=dataset; verbose=FALSE;calc_d2 = calc_d2
  
  #tryCatch(
  #  {
  if(return_option=="list" & verbose==TRUE) {
    print("Starting parameters:")
    print(as.matrix(input_par))  
  }
  
  #### MARGINAL MODEL: Get linear predictors and their inverse -> **Creates eta, eta_link_inv** from **mm_mar** ####
  margin_names=unique(dataset$time)
  num_margins = length(margin_names)
  response=dataset$response
  margin_parameters=names(mm_mar)
  
  eta = eta_linkinv = matrix(nrow = nrow(dataset), ncol = length(margin_parameters))
  colnames(eta)=colnames(eta_linkinv)= margin_parameters
  
  for (parameter in margin_parameters) {
    #parameter="mu"
    par_temp <- input_par[grepl(paste(parameter,sep="."), names(input_par))]
    mm=mm_mar[[parameter]]
    eta[, parameter] = rowSums(mm * matrix(rep(par_temp,each=nrow(mm)),ncol=length(par_temp),dimnames=list(NULL,c(names(par_temp)))))
    FUN=margin_dist[[names(margin_dist)[grepl(paste(parameter,"linkinv",sep="."),names(margin_dist))]]]
    eta_linkinv[, parameter]=FUN(eta[,parameter])
  }
  
  #### MARGINAL MODEL: Get density and probability functions for margins -> ** Creates margin_p and margin_d ####
  margin_p = pFUN(    response ,mu=eta_linkinv[,"mu"] ,sigma=eta_linkinv[,"sigma"] ,nu=eta_linkinv[,"nu"],tau=   eta_linkinv[,"tau"])
  margin_d = dFUN(x = response ,mu=eta_linkinv[,"mu"] ,sigma=eta_linkinv[,"sigma"] ,nu=eta_linkinv[,"nu"],tau=   eta_linkinv[,"tau"])
  
  #### Get density of the copula function for each pair of margins ####
  
  #margin_names; num_margins; margin_parameters
  
  #Calculate new eta based on input parameters and model matrix
  copula_parameters=c("theta","zeta")
  
  eta_cop=eta_cop_link_inv=matrix(nrow=nrow(mm_cop$theta),ncol=length(copula_parameters))
  colnames(eta_cop)=colnames(eta_cop_link_inv)=copula_parameters
  for (parameter in copula_parameters) {
    #parameter="theta"
    par_temp = input_par[grepl(parameter, names(input_par))]
    mm=mm_cop[[parameter]]
    eta_cop[,parameter] = rowSums(mm * matrix(rep(par_temp,each=nrow(mm)),ncol=length(par_temp),dimnames=list(NULL,c(names(par_temp))))) #Calc linear predictor
    FUN=copula_link[[names(copula_link)[grepl(paste(parameter,"linkinv",sep="."),names(copula_link))]]]
    eta_cop_link_inv[,parameter]=FUN(eta_cop[,parameter])
  }
  
  copula_d=dcdu1=dcdu2=dldth=d2ldth=dldz=d2ldz=dthdeta=dzdeta=matrix(0,nrow=nrow(eta_cop),ncol=1)
  
  margin_p_cop_input=cbind(margin_p[dataset$time %in% margin_names[1:(num_margins-1)]],margin_p[dataset$time %in% margin_names[2:(num_margins)]])
  
  copula_d=BiCopPDF(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"])
  dldth=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="par",log=TRUE)
  dcdu1=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="u1",log=FALSE)
  dcdu2=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="u2",log=FALSE)
  
  if(calc_d2==TRUE) {
    d2ldth=BiCopDeriv2( margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="par")
  }
  dthdeta=copula_link$dthdeta(eta_cop[,"theta"])
  
  if (length(copula_parameters)>1) {
    dldz=BiCopDeriv(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="par2",log=TRUE)
    if(calc_d2==TRUE) {
      d2ldz=BiCopDeriv2(margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="par2")
    }
    dzdeta=copula_link$dzdeta(eta_cop[,"zeta"])
  }
  
  #### Return Log Likelihood ####
  
  log_lik_results<-c(sum(log(margin_d)),sum(log(copula_d)),sum(log(margin_d))+sum(log(copula_d)))
  names(log_lik_results)=c("Margin","Copula","Total")
  
  #### Extract derivatives of margins ####
  all_derivatives=c(names(margin_dist)[grepl("dl", names(margin_dist))],names(margin_dist)[grepl("dr", names(margin_dist))])
  
  if(calc_d2==TRUE) {all_derivatives=c(all_derivatives,names(margin_dist)[grepl("d2l", names(margin_dist))])}
  derivatives_calculated_all_margins=list()
  
  #Margin derivatives
  
  #for (margin_number in 1:(length_fitted_margins)) {
  i=1 #DON'T DELETE THIS - it's an internal counter
  derivatives_calculated=matrix(nrow=nrow(eta),ncol=length((all_derivatives)))
  for (d_fun in all_derivatives) { 
    if (grepl("dr", d_fun) == TRUE) {
      derivatives_calculated[,i]=margin_dist[[d_fun]](eta[,grepl(sub(".dr","",d_fun), colnames(eta))])
    } else {
      derivatives_calculated[,i]=margin_dist[[d_fun]](y=response,mu=eta_linkinv[,"mu"],sigma=eta_linkinv[,"sigma"],nu=eta_linkinv[,"nu"],tau=eta_linkinv[,"tau"])
    }
    i=i+1
  }
  colnames(derivatives_calculated)=all_derivatives
  derivatives_calculated_all_margins=derivatives_calculated
  #}
  
  if (calc_d2==TRUE) {
    derivatives_copula = list(dldth,dldz,d2ldth,d2ldz,dcdu1,dcdu2,dthdeta,dzdeta)
    names(derivatives_copula)=c("dldth","dldz","d2ldth","d2ldz","dcdu1","dcdu2","dthdeta","dzdeta")
  } else {
    derivatives_copula = list(dldth,dldz,dcdu1,dcdu2,dthdeta,dzdeta)
    names(derivatives_copula)=c("dldth","dldz","dcdu1","dcdu2","dthdeta","dzdeta")
  }
  
  return_list=list(log_lik_results,response,margin_d,margin_p,copula_d,derivatives_copula,margin_dist,copula_dist,eta,eta_linkinv,derivatives_calculated_all_margins,mm_mar,start_par,eta_cop,eta_cop_link_inv,mm_cop)
  names(return_list)=c("log_lik_results","response","margin_d","margin_p","copula_d","derivatives_copula","margin_dist","copula_dist","eta","eta_linkinv","derivatives_calculated_all_margins","mm_mar","start_par","eta_cop","eta_cop_link_inv","mm_cop")
  
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

extract_parameters_from_fitted <- function(fitted_margins,fitted_copulas=NA,copula_link) {
  mu=list()
  sigma=list()
  nu=list()
  tau=list()
  for (i in 1:length(fitted_margins)) {
    mu[[i]]=coefAll(fitted_margins[[i]])[1]
    sigma[[i]]=coefAll(fitted_margins[[i]])[2]
    nu[[i]]=coefAll(fitted_margins[[i]])[3]
    tau[[i]]=coefAll(fitted_margins[[i]])[4]
    #Add back in for multiple margins
    #names(mu[[i]])=paste(names(mu[[i]]),i,sep=".")
    #names(sigma[[i]])=paste(names(sigma[[i]]),i,sep=".")
    #names(nu[[i]])=paste(names(nu[[i]]),i,sep=".")
    #names(tau[[i]])=paste(names(tau[[i]]),i,sep=".")
  }  
  
  if (!(is.na(fitted_copulas))) {
    theta=zeta=vector()
    for (i in 1:length(fitted_copulas)) {
      theta[i]=copula_link$theta.linkfun(fitted_copulas[[i]]$par)
      zeta[i]=copula_link$zeta.linkfun(fitted_copulas[[i]]$par2)
      names(theta)[i]=paste("theta",names(fitted_copulas)[i])
      names(zeta)[i]=paste("zeta",names(fitted_copulas)[i])
    }
    return_list=(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau),theta,zeta))
  } else {
    return_list=(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau)))
  }
  
  return(return_list)
}

score_function <- function(eta,dldpar,dpardeta,dlcopdpar,d2ldpar,response,phi=1,step_size=1,verbose=TRUE) {
  #u_k
  dlcopdpar=dlcopdpar*((dldpar)*response)
  
  dldpar=dldpar+dlcopdpar #commented out for now until we come back and fix it
  
  u_k=dldeta = dldpar * dpardeta
  
  #f_k and w_k
  #f_k_quasi_newton=(-(dldpar*dldpar))
  f_k=(d2ldpar)
  
  w_k=-f_k*(dpardeta*dpardeta)
  
  #w_k=matrix(1,ncol=ncol(f_k),nrow=nrow(f_k)) ####Weights of 1
  
  z_k=(1-phi)*eta+phi*(eta+step_size*(u_k/w_k))
  
  if(verbose==TRUE) {
    steps_mean=round(rbind(colMeans(eta)
                           ,colMeans(dldpar-dlcopdpar)
                           ,colMeans(dlcopdpar)
                           ,colMeans(dpardeta)
                           ,colMeans(dpardeta*dpardeta)
                           ,colMeans(f_k)
                           ,colMeans(w_k)
                           ,colMeans(u_k)
                           ,colMeans((1/w_k)*u_k)
                           ,colMeans(z_k)
    ),8)
    rownames(steps_mean)=c("eta","dldpar","dlcopdpar","dpardeta","dpardeta2","f_k","w_k","u_k","(1/w_k)*u_k","z_k")
    print(steps_mean)  
  }
  return_list=list(colMeans(z_k),u_k,f_k,w_k,z_k)
  names(return_list)=c("par","u_k","f_k","w_k","z_k")
  return(return_list)
}

newton_raphson_iteration=function(results,input_par,phi=1,step_size=1,verbose=c(FALSE,FALSE),calc_d2=FALSE) {
  
  #results=results;input_par=input_par;verbose=c(TRUE,TRUE); calc_d2=FALSE
  
  steps_all=list()
  par_out=list()
  margin_dist=results$margin_dist
  copula_dist=results$copula_dist
  
  
  ####COPULA DERIVATIVES####
  ### Extract COPULA DERIVATIVES from results list (from calc_joint_likelihood)
  
  dlcopdpar=list()
  dldth=results$derivatives_copula$dldth
  dldz=results$derivatives_copula$dldz
  dthdeta=results$derivatives_copula$dthdeta
  dzdeta=results$derivatives_copula$dzdeta
  dcdu1=results$derivatives_copula$dcdu1
  dcdu2=results$derivatives_copula$dcdu2
  
  dldpar=cbind(dldth,dldz)
  dpardeta=cbind(dthdeta,dzdeta)
  
  dlcopdpar=(dcdu1+dcdu2)/(results$copula_d) #+ (dcdu2)/results$copula_d[,i]  #*(results$margin_d[,i]*results$margin_d[,i+1])
  
  eta_input=cbind(results$eta_cop[,"theta"],results$eta_cop[,"zeta"])
  colnames(eta_input)=names(results$eta_cop)
  
  if(calc_d2 == TRUE) {
    
    count=1
    d2names=names(results$derivatives_copula)[grepl("d2",names(results$derivatives_copula))]
    for (names in d2names) {
      if (count==1) {
        d2cdpar=matrix(results$derivatives_copula[[names]])
      } else {
        d2cdpar=cbind(d2cdpar,results$derivatives_copula[[names]])
      }
      count=count+1
    }
    colnames(d2cdpar)=d2names
    
    d2ldpar=(d2cdpar/copula_d)-(dldpar^2)
    
  } else {
    d2ldpar=-(dldpar*dldpar)
  }
  
  #plot.new();par(mfrow=c(1,1));plot(d2ldpar[,1],-(dldpar^2)[,1])
  #abline(a=1,b=1,col="red")
  
  copula_score=score_function(eta=eta_input,dldpar=dldpar,dpardeta=dpardeta,dlcopdpar=dldpar*0,d2ldpar=d2ldpar,response=dldpar*0,phi=phi,step_size=step_size,verbose=verbose[2]) ##Returns updated value
  
  beta_new_cop=list()
  cop_parameter_names=c("theta","zeta")
  
  hess_list=list()
  
  e_k=copula_score$z_k  
  beta_new_cop=list()
  for (j in 1:ncol(copula_score$z_k)) {
    mm=mm_cop[[j]]
    X=as.matrix(mm)
    W=diag(copula_score$w_k[,j])
    e_k_par=e_k[,j]
    beta_new_cop[[cop_parameter_names[j]]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
    rownames(beta_new_cop[[cop_parameter_names[j]]])= paste(paste(cop_parameter_names[j],sep=" "),colnames(mm),sep=".")
    #names(beta_new_cop[[j]])= paste(parameter_names[j],margin_number,rownames(beta_new_cop[[j]]),sep=".")
    if (calc_d2 == TRUE) {
      hess_list[[(cop_parameter_names[j])]]=-t(X)%*%W%*%X
    }
  }
  
  one=cbind(dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],c("time","subject")],dlcopdpar)
  two=cbind(dataset[!(dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)]),c("time","subject")],0)
  colnames(two)=colnames(one)
  one_two=rbind(one,two)
  
  one=cbind(dataset[!(dataset$time %in% unique(dataset$time)[2:(length(unique(dataset$time)))]),c("time","subject")],0)
  two=cbind(dataset[dataset$time %in% unique(dataset$time)[2:(length(unique(dataset$time)))],c("time","subject")],dlcopdpar)
  colnames(one)=colnames(two)
  two_three=rbind(one,two)
  
  one_two=one_two[order(one_two$subject,one_two$time),]
  two_three=two_three[order(two_three$subject,two_three$time),]
  
  ###CHECKS ON ORDERING
  all(one_two$time==two_three$time)
  all(one_two$time==dataset$time)
  all(one_two$subject==two_three$subject)
  all(one_two$subject==dataset$subject)
  
  dlcopdpar_one_col=(one_two$dlcopdpar+two_three$dlcopdpar)*results$margin_d
  
  
  ####MARGIN DERIVATIVES####
  margin_parameters=colnames(results$eta)
  
  dlcopdpar_final=dlcopdpar_one_col
  response_final=results$response
  for (i in 2:length(margin_parameters)) {
    dlcopdpar_final=cbind(dlcopdpar_final,dlcopdpar_one_col)
    response_final=cbind(response_final,results$response)
  }
  
  dldpar=results$derivatives_calculated_all_margins[,c("dldm","dldd","dldv","dldt")] #First derivative of log likelihood w.r.t. parameters
  dpardeta=results$derivatives_calculated_all_margins[,c("mu.dr","sigma.dr","nu.dr","tau.dr")] #Derivative of inverse link function
  
  if(calc_d2 == TRUE) {
    d2ldpar=results$derivatives_calculated_all_margins[,c("d2ldm2","d2ldd2","d2ldv2","d2ldt2")]
  } else {
    d2ldpar=-(dldpar*dldpar)
  }
  
  margin_score=score_function(eta=results$eta[,margin_parameters],dldpar=dldpar,dpardeta=dpardeta,dlcopdpar=dlcopdpar_final,d2ldpar=d2ldpar, response=response_final,phi=phi,step_size=step_size,verbose=verbose[1]) ##Returns updated value
  
  e_k=margin_score$z_k#-results$eta_nomargin[[1]][,c("mu","sigma","nu","tau")]
  beta_r=list()
  #for (margin_number in 1:length(results$mm_all)) {
  for (parameter in margin_parameters) {
    mm=results$mm_mar[[parameter]]
    X=as.matrix(mm)
    W=diag(margin_score$w_k[,which(margin_parameters==parameter)])
    e_k_par=e_k[,which(margin_parameters==parameter)]
    beta_r[[parameter]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
    rownames(beta_r[[parameter]])= paste(parameter,rownames(beta_r[[parameter]]),sep=".")
    if (calc_d2 == TRUE) {
      hess_list[[parameter]]=-t(X)%*%W%*%X
    }
  }
  
  beta_new=matrix(ncol=1,nrow=0)
  colnames(beta_new)=c("Estimate")
  
  for (parameter in names(beta_r)) {
    beta_new=rbind(beta_new,beta_r[[parameter]])
  }
  for (parameter in names(beta_new_cop)) {
    beta_new=rbind(beta_new,beta_new_cop[[parameter]])
  }
  
  end_par=as.vector(beta_new)
  names(end_par)=rownames(beta_new)
  end_par=end_par[names(start_par)]
  
  if (calc_d2==TRUE) {
    
    
    cbind(results$derivatives_copula$d2ldth,results$derivatives_copula$d2ldz)
    
    
    
    se_par=vector()
    for (names in c(margin_parameters)) {
      n=nrow(margin_score$z_k)
      se_var=(diag(solve(-hess_list[[names]])))
      #se_var=diag(hess_list[[names]])/n
      se=sqrt(se_var)#/sqrt(n)
      se_par=c(se_par,se)
    }
    for (names in cop_parameter_names) {
      n=nrow(copula_score$z_k)
      se_var=(diag(solve(-hess_list[[names]])))
      se=sqrt(se_var)#/sqrt(n)
      se_par=c(se_par,se)
    }
    return_list=list(input_par,end_par,se_par)
    names(return_list)=c("input_par","end_par","se_par")
  } else {
    return_list=list(input_par,end_par)
    names(return_list)=c("input_par","end_par")
  }
  
  return(return_list)
}

#### 0. Load data subset and transform to standard longitudinal dataset #######

load("Data/rand_mvt.rds")
head(rand_mvt)

# Basic data setup
response = rand_mvt[,4:(4+2)]#[,4:18] ####Currently limiting to just 5 margins for simplicity
covariates=list()
covariates[[1]] = as.data.frame(rand_mvt[,19]) #Age 19 - changed to age at start to avoid correlation with time
covariates[[2]] = as.data.frame(rand_mvt[,34:(34+4)]) #Time 34:48
covariates[[3]] = as.data.frame(rand_mvt[,3]) #Gender

# Setup data as longitudinal file
dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
head(dataset)

##### 0. OR simualte a dataset

# set up D-vine copula model with mixed pair-copulas
d <- 3
dd <- d*(d-1)/2
order <- 1:d
family <- c(2, 2, 0)
par <- c(0.2, 0.9, .9)
par2 <- c(5,10,7)

# transform to R-vine matrix notation
RVM <- D2RVine(order, family, par, par2)
contour(RVM)

n=100;t=d
copsim=RVineSim(n*t,RVM)

margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
for ( i in 1:ncol(copsim)) {
  margin[,i]=qZISICHEL(copsim[,i])
}

response = as.data.frame(margin)
covariates=list()
covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

##########

#plotDist(dataset,"ZISICHEL")
#####################################################################################
################################## INPUT: Parameters ################################
#####################################################################################

#Formulas
mu.formula = formula("response ~ as.factor(time)+as.factor(gender)")
sigma.formula = formula("~ as.factor(time)+age")
nu.formula = formula("~ as.factor(time)+as.factor(gender)")
tau.formula = formula("~ age")
theta.formula=formula("response~as.factor(time)")
zeta.formula=formula("response~as.factor(time)")

#mu.formula = formula("response ~ 1")
#sigma.formula = formula("~ 1")
#nu.formula = formula("~ 1")
#tau.formula = formula("~ 1")
#theta.formula=formula("response~as.factor(time)")
#zeta.formula=formula("response~1")

#Link functions and distributions
margin_dist = ZISICHEL(
  mu.link = "log",
  sigma.link = "log",
  nu.link = "identity",
  tau.link = "logit"
)
dFUN=dZISICHEL;pFUN=pZISICHEL
copula_dist="t"
copula_link=list(logit,logit_inv,log_2plus,log_2plus_inv,dlogit,dlog_2plus)
names(copula_link)=c("theta.linkfun","theta.linkinv","zeta.linkfun","zeta.linkinv","dthdeta","dzdeta")

##### Setup model matrices and parameter vector for optimisation #####

#Fit gamlss marginal model
gamlss_model <- gamlss(formula = mu.formula
                       , sigma.formula = sigma.formula
                       , nu.formula = nu.formula
                       , tau.formula = tau.formula
                       , family=margin_dist,data=dataset,method=RS(100))

plot(gamlss_model)

fitted_copulas=fitted_margins=list()
fitted_margins[[1]]=gamlss_model

#Extract margin model matrix
mm_mar=list()
for (parameter in c("mu","sigma","nu","tau")) {
  #print(head(model.matrix(fitted_margins[[1]],what=parameter)))
  mm_mar[[parameter]]= model.matrix(fitted_margins[[1]],what=parameter)
  #print(head(mm_mar[[parameter]]))
}

#Fit starting copula with just an intercept for each parameter
start_par<-extract_parameters_from_fitted(fitted_margins,fitted_copulas=NA,copula_link)
copula_model<-fit_copula_single(dataset=dataset,mm_mar=mm_mar,input_par=start_par,dFUN=dFUN,pFUN=pFUN,copula_dist=copula_dist)
fitted_copulas[[1]]=copula_model
start_par<-extract_parameters_from_fitted(fitted_margins,fitted_copulas=fitted_copulas,copula_link)

#Generate the copula model matrix
generate_cop_model_matrix = function(dataset,theta.formula=formula(~1),zeta.formula=formula(~1),time) {
  #time="time"
  margins=unique(dataset[,time])
  num_copulas=length(margins)-1
  first=TRUE
  for (margin_number in 1:num_copulas) {
    time_1=dataset[dataset[,time]==margins[margin_number],]
    time_2=dataset[dataset[,time]==margins[margin_number+1],]
    names(time_2)=paste(names(time_2),2,sep=".")
    add_to_cop_data_matrix=cbind(time_1,time_2)
    if (first==TRUE) {
      cop_data_matrix=add_to_cop_data_matrix; first=FALSE
    } else {
      cop_data_matrix=rbind(cop_data_matrix,add_to_cop_data_matrix)
    } 
  }
  #matrix_of_ones=matrix(1,nrow=nrow(cop_data_matrix),ncol=1)
  #colnames(matrix_of_ones)="(Intercept)"
  cop_model_matrix=list()
  #cop_model_matrix[["theta"]]=cbind(matrix_of_ones,model.frame(formula,data=cop_data_matrix))
  #cop_model_matrix[["zeta"]]=cbind(matrix_of_ones,model.frame(zeta.formula,data=cop_data_matrix))
  cop_model_matrix[["theta"]]=model.matrix(theta.formula,data=cop_data_matrix,contrasts.arg = sapply(cop_data_matrix,is.factor))
  cop_model_matrix[["zeta"]]=model.matrix(zeta.formula,data=cop_data_matrix)
  
  return_list=(cop_model_matrix)
  #names(return_list)=("cop_model_matrix")
  return(return_list)
}
#mm_cop=generate_cop_model_matrix(dataset=dataset,formula=theta.formula,zeta.formula=zeta.formula,time="time")

#hacky model matrix
mm_cop=list()
mm_cop[["theta"]]=model.matrix(gamlss(formula=theta.formula,data=(dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],])))
mm_cop[["zeta"]]=model.matrix(gamlss(formula=zeta.formula,data=(dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],])))

#Create parameter vector from model matrix for copulas
temp_cop_parameter_names=c()
for (parameter in names(mm_cop)) {
  temp_cop_parameter_names=c(temp_cop_parameter_names,paste(parameter,colnames(mm_cop[[parameter]]),sep="."))
}

#Setting up starting parameter vector with correct factors from mm_mar and mm_cop
theta_par_loc=grepl("theta",names(start_par))
zeta_par_loc=grepl("zeta",names(start_par))

temp_cop_start_par=vector(length=length(temp_cop_parameter_names))
names(temp_cop_start_par)=temp_cop_parameter_names

temp_cop_start_par[grepl("Intercept",names(temp_cop_start_par))&grepl("theta",names(temp_cop_start_par))]=start_par[theta_par_loc][1] #Theta starting value
temp_cop_start_par[grepl("Intercept",names(temp_cop_start_par))&grepl("zeta",names(temp_cop_start_par))]=start_par[zeta_par_loc][1] #Theta starting value

start_par=c(start_par[!(theta_par_loc | zeta_par_loc)],temp_cop_start_par)

#Starting joint likelihood
results_start= calc_joint_likelihood(input_par =start_par  #optim_results$par
                                     , mm_mar= mm_mar
                                     , margin_dist = margin_dist
                                     , copula_dist=copula_dist
                                     , copula_link=copula_link
                                     , mm_cop = mm_cop
                                     , return_option="list"
                                     , dataset=dataset
                                     , verbose=FALSE
                                     , calc_d2=TRUE #For ease of debugging - turn off for later versions
)

print(results_start$log_lik_results)

#####################################################################################
################################ INPUT: RUN ITERATIONS ##############################
#####################################################################################
first_input_par=start_par
end_par_matrix=matrix(0,ncol=length(start_par),nrow=0)
end_loglik_matrix=matrix(0,ncol=1,nrow=0)
first_run=TRUE
change=1
run_counter=1
phi=1
step_adjustment=0.5^(1/2)
step_size=1
verbose_option=c(FALSE,FALSE)
stopifnegative=FALSE
calc_d2=FALSE
while ((abs(change) > .1*phi*step_adjustment) & run_counter <= 100) { #
  if(!first_run) {start_log_lik=results$log_lik_results} else {input_par=first_input_par; phi_inner=phi}
  
  #print(input_par)
  results= calc_joint_likelihood(input_par =input_par  #optim_results$par
                                 ,mm_mar = mm_mar
                                 ,margin_dist = ZISICHEL(
                                   mu.link = "log",
                                   sigma.link = "log",
                                   nu.link = "identity",
                                   tau.link = "logit"
                                 )
                                 , copula_dist="t"
                                 , copula_link=copula_link
                                 ,mm_cop = mm_cop
                                 , return_option="list"
                                 , dataset=dataset
                                 , verbose=FALSE
                                 ,calc_d2 = calc_d2
  )
  
  end_log_lik=results$log_lik_results
  
  if(first_run) {change=1} else {
    change=end_log_lik["Total"]-start_log_lik["Total"]
    log_lik_output=c(start_log_lik["Total"],end_log_lik["Total"],change,phi_inner)
    names(log_lik_output) = c("Start LogLik","End LogLik","Change","phi")
    print(log_lik_output)
  }
  
  iteration_out=newton_raphson_iteration(results
                                         ,input_par
                                         ,phi=phi_inner
                                         ,step_size=step_size
                                         ,verbose=verbose_option
                                         ,calc_d2=calc_d2)
  cbind(iteration_out[[1]],iteration_out[[2]])
  end_par=iteration_out[[2]]
  input_par=end_par
  first_run=FALSE
  end_par_matrix=rbind(end_par_matrix,end_par)
  end_loglik_matrix=rbind(end_loglik_matrix,end_log_lik[["Total"]])
  run_counter=run_counter+1
  #Showing charts
  colnames(end_par_matrix)=names(end_par)
  pars=round(sqrt(ncol(end_par_matrix)+1),0)
  plot.new();par(mfrow=c(pars,pars+1))
  for (par_name in colnames(end_par_matrix)) {
    #scale_range=0.5
    #if(abs((end_par_matrix[nrow(end_par_matrix),par_name]-start_par[par_name]))>abs(scale_range*start_par[par_name])){
    #            
    #} else {
    #  lims=c(start_par[par_name]-abs(start_par[par_name])*scale_range,start_par[par_name]+abs(start_par[par_name])*scale_range)
    #  plot(1:nrow(end_par_matrix),end_par_matrix[,par_name],main=par_name,ylim=lims)
    #}
    plot(1:nrow(end_par_matrix),end_par_matrix[,par_name],main=par_name)  
  }
  
  plot(1:nrow(end_par_matrix),end_loglik_matrix[,1],main="loglik")
  
  if(stopifnegative==TRUE) {
    if(length(end_loglik_matrix)>5) {
      if(all(end_loglik_matrix[(length(end_loglik_matrix)-3):(length(end_loglik_matrix)-1)]-end_loglik_matrix[(length(end_loglik_matrix)-2):length(end_loglik_matrix)]>0)){break}
    }
  }
  
  phi_inner=phi*(step_adjustment^min(2,run_counter-1))
}

print("Separate likelihood")
print(results_start$log_lik_results)
print("Joint likelihood")
print(results$log_lik_results)

###############

final_results= calc_joint_likelihood(input_par =end_par  #optim_results$par
                                     ,mm_mar = mm_mar
                                     ,margin_dist = ZISICHEL(
                                       mu.link = "log",
                                       sigma.link = "log",
                                       nu.link = "identity",
                                       tau.link = "logit"
                                     )
                                     , copula_dist="t"
                                     , copula_link=copula_link
                                     ,mm_cop = mm_cop
                                     , return_option="list"
                                     , dataset=dataset
                                     , verbose=FALSE
                                     , calc_d2=calc_d2
)


final_iteration_out=newton_raphson_iteration(final_results,end_par,phi=phi_inner,step_size=step_size,verbose=verbose_option,calc_d2=TRUE)

new_SE=final_iteration_out$se_par
old_SE=sqrt(diag(vcov(fitted_margins[[1]])))
old_SE=c(old_SE,c(fitted_copulas$`1,2`$se,fitted_copulas$`2,3`$se, fitted_copulas$`1,2`$se2,fitted_copulas$`2,3`$se2))

z_score_old=(start_par[0:length(old_SE)]-0)/old_SE
z_score_new=(end_par-0)/new_SE

p_new=pNO(-abs(z_score_new))
p_old=pNO(-abs(z_score_old))

SEs=cbind(start_par,old_SE,p_old,end_par,new_SE,p_new)
print(round(SEs,5))












####ARCHIVE BEYOND HERE####

#### Manual gradient and hessian matrix

results= calc_joint_likelihood(input_par =start_par
                               ,start_fitted_margins = fitted_margins
                               ,margin_dist = ZISICHEL(
                                 mu.link = "log",
                                 sigma.link = "log",
                                 nu.link = "identity",
                                 tau.link = "logit"
                               )
                               , copula_dist="t"
                               , copula_link=copula_link
                               ,mm_cop = mm_cop
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
             , return_option="log_lik"
             ,dataset=dataset
)
grad_l

hessian_l<-hessian(calc_joint_likelihood,
                   , x=input_par
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
                   ,dataset=dataset
)
hessian_l

########LATER
new_SE=sqrt(diag(solve(-hessian_l)))
old_SE=sqrt(diag(vcov(fitted_margins[[1]])))
old_SE=c(old_SE,c(fitted_copulas$`1,2`$se,fitted_copulas$`2,3`$se, fitted_copulas$`1,2`$se2,fitted_copulas$`2,3`$se2))
SEs=cbind(start_par,old_SE,end_par,new_SE)
SEs

#### Exploratory

#### Theta by gender

plotDist(dataset[dataset$gender==1,],dist="ZISICHEL")
fitted_margins<-fit_margins(mu.formula=formula(response~age),dataset=dataset[dataset$gender==1,],family=ZISICHEL())
#AIC(fitted_margins[[1]])
#plot(fitted_margins[[1]])
#term.plot(fitted_margins[[2]])

fitted_copulas<-fit_copulas(fitted_margins,copula_dist="t",method="vine")
summary(fitted_copulas)
contour(fitted_copulas)


#### 3. Separate margin and copula models

#Goal

fitted_margins=list()
fitted_margins[[1]]=gamlss_model

# Starting fit for margins and copula
#best_marginal_fits <- find_best_marginal_fits(dataset,type="counts") #Find best fit margins
#fitted_margins<-fit_margins(mu.formula=formula(response~age+as.factor(gender)),dataset=dataset,margin_dist=ZISICHEL())
#AIC(fitted_margins[[1]])
#plot(fitted_margins[[1]])
#term.plot(fitted_margins[[2]])

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


#### 2. Benchmark models (optional)

###### GLMM ZISICHEL logLik(gamlss_glmm_model) = -11187.34 (df=1042.17) AIC: 24459.02 | BIC: 31322.54
#margin_model_formula_glmm=formula(response~(age)+as.factor(gender)+as.factor(time)+random(as.factor(subject)))
#gamlss_glmm_model <- gamlss(formula = margin_model_formula_glmm,family="ZISICHEL",data=dataset) #So this doesn't run automatically...
#summary(gamlss_glmm_model)
#plot(gamlss_glmm_model)
#term.plot(gamlss_glmm_model)

#### Extract parameters and SE from separate and jointly optimised datasets -> all_out_combined
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



