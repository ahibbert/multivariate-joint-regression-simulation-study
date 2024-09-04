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

SEs_ALL=list()
se_count=1

true_val_matrix=matrix(nrow=0,ncol=4)
for (a in c(3)) {
  for (b in c(1)) {
    for (c in c(1)) {
      for (d in c(-10:10*.1)) { #-5:15*.1
        true_val_matrix=rbind(true_val_matrix,matrix(c(a,b,c,d),ncol=4,byrow=TRUE))
      }}}}

#true_val_matrix=matrix(c(1,1,1,1),ncol=4,byrow=TRUE)


for (run_counter in 1:nrow(true_val_matrix)) {
  true_val=true_val_matrix[run_counter,]
  #### 0. Load RAND data subset and transform to standard longitudinal dataset #######
  #copula.link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv)
  copula.link=list(log,exp,dloginv=exp); two_par_cop=FALSE
  if(two_par_cop) {
    names(copula.link)=c("theta.linkfun","theta.linkinv","dthdeta","zeta.linkfun","zeta.linkinv","dzdeta")} else {names(copula.link)=c("theta.linkfun","theta.linkinv","dthdeta")}
  #true_val=c(.3,.2,.1,.3,.2,.1,-0.8,-2.197,.6,.2)
  #true_val=c(4,1,1,1)
  
  dataset=loadDataset(simOption=3
                      ,plot_dist=FALSE
                      ,n=500
                      ,d=3
                      ,copula.family=BiCopName("C")
                      ,copula.link=copula.link
                      ,qFUN=qEXP
                      ,par.copula=c(rep(true_val[4],2))
                      ,par.margin=true_val[1:3]) # c(0.5,0.5,0.5,0.5,1,1,1,1))
  
  plotDist(dataset,"EXP")
  
  ###########################
  ####### USER INPUT ########
  ###########################
  
  #Formulas
  setup=
  GJRM_L_SETUP(
    mu.formula = ("response ~ time+as.factor(gender)"),#mu.formula = formula("response ~ as.factor(time)+as.factor(gender)+age")
    sigma.formula = ("~ 1"),#sigma.formula = formula("~ as.factor(time)+age")
    nu.formula = ("~ 1"),#nu.formula = formula("~ as.factor(time)+as.factor(gender)")
    tau.formula = ("~ 1"),#tau.formula = formula("~ age")
    theta.formula=("~1"),#theta.formula=formula("response~as.factor(gender)")
    zeta.formula=("~1"),#zeta.formula=formula("response~1")
    margin.family=EXP(mu.link = "log"), #ZISICHEL(mu.link = "log",sigma.link = "log",nu.link = "identity",tau.link = "logit"),
    dFUN=dEXP,pFUN=pEXP,
    copula.family="C",
    copula.link=copula.link,
    start="nofit",
    calc_results=TRUE
  )
  
  #### Newton Raphson optimisation ####
  optim=GJRM_L_OPTIM(
      dataset,
      start_par=setup$start_par,
      mm_mar=setup$mm_mar,
      mm_cop=setup$mm_cop,
      margin_dist=setup$margin.family,
      copula_dist=setup$copula.family,
      copula_link=setup$copula.link,
      dFUN=setup$dFUN,
      pFUN=setup$pFUN,
      verbose=FALSE,
      calc_d2=TRUE,
      phi=.5,
      step_size=1,
      step_adjustment=.5,
      step_reduc_count=3,
      crit_wk=0.0000001,
      true_val=true_val,
      stopifnegative=TRUE,
      stopval=.1,
      max_runs=100,
      plot_optim=TRUE,
      plot_end_only=FALSE
  )
  
  post_optim=GJRM_POST_OPTIM(optim,setup,dataset)
  
  end_par=post_optim$end_par
  new_SE=post_optim$new_SE
  
  margin.family=setup$margin.family
  old_SE=sqrt(diag(vcov(setup$fitted_margins[[1]])))
  old_SE=c(old_SE,c(setup$fitted_copulas$`1,2`$se,setup$fitted_copulas$`2,3`$se, setup$fitted_copulas$`1,2`$se2,setup$fitted_copulas$`2,3`$se2))
  z_old=(setup$start_par[0:length(old_SE)]-0)/old_SE
  z_new=(end_par-0)/new_SE
  SEs=cbind(setup$start_par,old_SE,z_old,end_par,new_SE,z_new,optim$true_val)
  rownames(SEs)=names(optim$true_val)
  colnames(SEs)=c("Start Par","Old SE","Z Old","End Par","New SE","Z New","True Val")
  print(round(SEs,3))
  
  SEs_ALL[[se_count]]=list()
  SEs_ALL[[se_count]][["SE"]]=SEs
  SEs_ALL[[se_count]][["LogLik"]]=optim$end_loglik_matrix
  SEs_ALL[[se_count]][["Converged"]]=(optim$converged)
  se_count=se_count+1
}






#################### POST RUN ANALYSIS #################

rownames(SEs_ALL[[1]]["SE"])

parameter_estimates=list()
for (run in 1:length(SEs_ALL)) {
  for (parameter in rownames(SEs_ALL[[run]]$SE)) {
    parameter_estimates[[parameter]]=rbind(parameter_estimates[[parameter]],c(SEs_ALL[[run]]$SE[parameter,],SEs_ALL[[run]][["Converged"]]))
  }
}

logliks=matrix(ncol=2,nrow=length(SEs_ALL))
for (run in 1:length(SEs_ALL)) {
  logliks[run,]=c(SEs_ALL[[run]]$LogLik[1,3],SEs_ALL[[run]]$LogLik[nrow(SEs_ALL[[run]]$LogLik),3])
}

conv=matrix(ncol=1,nrow=length(SEs_ALL))
for (run in 1:length(SEs_ALL)) {
  conv[run,]=c(SEs_ALL[[run]]$Converged)
}

plot.new(); par(mfrow=c(1,1))
lm_fit=lm(logliks[,2]~logliks[,1])
plot(logliks[,1],logliks[,2],main=paste("Log Likelihoods, intercept: ", lm_fit$coefficients[1], " gradient: ",lm_fit$coefficients[2] ,sep=""))
abline(a=lm_fit[1],b=lm_fit[2],col="red")

plot.new()
par(mfrow=c(4,2))
for (parameter in names(parameter_estimates)) {
  JointReg=parameter_estimates[[parameter]][,"End Par"]
  #SepReg=parameter_estimates[[parameter]][,"Start Par"]
  True=parameter_estimates[[parameter]][,"True Val"]
  JointConvStatus=parameter_estimates[[parameter]][,8]
  #plot(True,SepReg,main=paste(parameter," Start",sep=""))
  #lm_fit=lm(SepReg~as.numeric(True))$coefficients
  #lm_fit[is.na(lm_fit)]=0
  #abline(a=lm_fit[1],b=lm_fit[2],col="red")
  plot(True,JointReg,main=paste(parameter," Joint",sep=""),col=ifelse(conv=="Converged","black","red"))
  lm_fit=lm(JointReg~as.numeric(True))$coefficients
  lm_fit[is.na(lm_fit)]=0
  abline(a=lm_fit[1],b=lm_fit[2],col="red")
}

#plot.new()
#par(mfrow=c(4,2))
for (parameter in names(parameter_estimates)) {
  JointReg=parameter_estimates[[parameter]][,"New SE"]
  #SepReg=parameter_estimates[[parameter]][,"Old SE"]
  True=parameter_estimates[[parameter]][,"True Val"]
  JointConvStatus=parameter_estimates[[parameter]][,8]
  #plot(True,SepReg,main=paste(parameter," Start",sep=""))
  #lm_fit=lm(SepReg~as.numeric(True))$coefficients
  #lm_fit[is.na(lm_fit)]=0
  #abline(a=lm_fit[1],b=lm_fit[2],col="red")
  plot(True,JointReg,main=paste(parameter," Joint",sep=""),col=ifelse(conv=="Converged","black","red"))
  #lm_fit=lm(JointReg~as.numeric(True))$coefficients
  #lm_fit[is.na(lm_fit)]=0
  #abline(a=lm_fit[1],b=lm_fit[2],col="red")
}

plot.new()
par(mfrow=c(4,2))
for (parameter in names(parameter_estimates)) {
  JointReg=logliks[,2]
  SepReg=logliks[,1]
  True=parameter_estimates[[parameter]][,"True Val"]
  JointConvStatus=parameter_estimates[[parameter]][,8]
  plot(True,SepReg,main=paste(parameter," Start",sep=""))
  #lm_fit=lm(SepReg~as.numeric(True))$coefficients
  #lm_fit[is.na(lm_fit)]=0
  #abline(a=lm_fit[1],b=lm_fit[2],col="red")
  plot(True,JointReg-SepReg,main=paste(parameter," Joint",sep=""))
  #lm_fit=lm(JointReg~as.numeric(True))$coefficients
  #lm_fit[is.na(lm_fit)]=0
  #abline(a=lm_fit[1],b=lm_fit[2],col="red")
}


plot(as.numeric(JointReg),as.numeric(SepReg))


#### Post optimisation calculation of SEs ####

#Final SEs


# #######ADD GAMLSS GLMM COMPARISON #######
# gamlss_glmm_model <- gamlss(formula = formula("response ~ time+as.factor(gender)+random(as.factor(subject))")
#                        #, sigma.formula = sigma.formula
#                        #, nu.formula = nu.formula
#                        #, tau.formula = tau.formula
#                        , family=setup$margin.family
#                        , data=dataset,method=mixed(5,100))
# 
# gamlss_coeffs=vector()
# for (name in c("mu","sigma","nu",'tau')) {
#   gamlss_coeffs=c(gamlss_coeffs,coef(gamlss_glmm_model,what=name))
# }
# 
# gamlss_vcov=sqrt(diag(vcov(gamlss_glmm_model)))
# gamlss_coeffs=gamlss_coeffs[!grepl("random",names(gamlss_coeffs))]
# 
# SEs_plus_gamlss=cbind(SEs,cbind(c(gamlss_coeffs,rep(NA,nrow(SEs)-length(gamlss_coeffs))),c(gamlss_vcov,rep(NA,nrow(SEs)-length(gamlss_coeffs)))))
# 
# colnames(SEs_plus_gamlss)=c("Start Par","Old SE","Z Old","End Par","New SE","Z New","True Val","GAMLSS Coeff","GAMLSS SE")
# print(round(SEs_plus_gamlss,3))






