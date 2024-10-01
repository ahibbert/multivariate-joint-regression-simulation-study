########### Starting parameters
#set.seed(1000)
options(scipen=999)
set.seed(1000)
library(VineCopula)
library(gamlss)
library(Rcpp)
sourceCpp("test.cpp")
source("common_functions.R")

#### INNER ITERATION ##################################################################

#1. Choose distribution
#copula_dist="C";margin_dist=EXP(); mu=10; sigma=NA;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5;  #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="C";margin_dist=PO(); mu=10; sigma=NA;nu=NA; tau=NA; theta=2; zeta=NA; simOption=5;  #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="C";margin_dist=GA(); mu=2; sigma=.5;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="C";margin_dist=NO(); mu=3; sigma=2;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="N";margin_dist=NO(); mu=3; sigma=2;nu=NA; tau=NA; theta=0.75; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#simOption=6
#copula_dist="C";margin_dist=NBI(); mu=3; sigma=1;nu=NA; tau=NA; theta=5; zeta=NA#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="N";margin_dist=NBI(); mu=3; sigma=1;nu=NA; tau=NA; theta=5; zeta=NA#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="C"; margin_dist=ST1(); mu=5;sigma=exp(1);nu=5;tau=5;theta=5;zeta=NA
#copula_dist="N"; margin_dist=ST1(); mu=1;sigma=exp(1);nu=10;tau=10;theta=5;zeta=NA

copula_dist="C"; margin_dist=ZISICHEL(); mu=1;sigma=exp(1);nu=-2;tau=.1;theta=5;zeta=NA;
copula_dist="N"; margin_dist=ZISICHEL(); simOption=1

#2. Set simulation parameters and simulate
n=100; d=3;  #simulation parameters

### 2. or choose a dataset #########APPLICATIONS DATASET
# copula_dist="C"; margin_dist=SICHEL();
# dataset=loadDataset(simOption=1, n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
#                    , par.margin=c(mu,sigma,nu,tau), par.copula=rep(theta,d-1))
# 
# plotDist(dataset,margin_dist)

input_par=c(mu,sigma,nu,tau,theta,zeta); names(input_par)=c("mu","sigma","nu","tau","theta","zeta")

no_dl_outer=w_dl_outer=matrix(nrow=0,ncol=4+length(input_par[!is.na(input_par)]))
true_val=matrix(nrow=0,ncol=length(input_par[!is.na(input_par)]))
colnames(no_dl_outer)=colnames(w_dl_outer)=c("iterations",names(input_par[!is.na(input_par)]),"ll margin","ll copula", "ll total")

plot_runs=TRUE
for (run_counter in 1:20) {
  tryCatch({
  start=Sys.time()
  print(paste("################################# RUN:",run_counter))
  
  #Generate dataset
  dataset=loadDataset(simOption=simOption,n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
                      , par.margin=input_par[c("mu","sigma","nu","tau")], par.copula=rep(input_par[c("theta")],d-1))
  
  #plotDist(dataset,margin_dist)
  #3. Fit model
  no_dl=fit_jointreg(dataset, margin_dist,copula_dist,
                     mu.formula = ("response ~ 1"), sigma.formula = ("~ 1"), nu.formula = ("~ 1"), tau.formula = ("~ 1"),
                     theta.formula=("~1"), zeta.formula=("~1"),
                     include_dlcopdpar=FALSE,
                     verbose=3, plot_results=FALSE#,  true_val=par_to_eta(input_par,copula_dist,margin_dist)
                     ,start_step_size = .5, step_adjustment = .5, max_steps=5, inner_stop_crit = .05, outer_stop_crit = 0.1
                     )
  
    #source("common_functions.R")
  w_dl=fit_jointreg(dataset, margin_dist, copula_dist,
                 mu.formula = ("response ~ 1"), sigma.formula = ("~ 1"), nu.formula = ("~ 1"), tau.formula = ("~ 1"),
                 theta.formula=("~1"), zeta.formula=("~1"),
                 include_dlcopdpar=TRUE,
                 verbose=3,plot_results=FALSE#, true_val=par_to_eta(input_par,copula_dist,margin_dist)
                 ,start_step_size = .5, step_adjustment = .5, max_steps=5, inner_stop_crit = .05, outer_stop_crit = 0.1
                 )
    
  end=Sys.time()
  print(paste("Time taken:",end-start))
  
  if(plot_runs==TRUE) {
    plot.new()
    par(mfrow=c(1,3))
    
    w_iter=nrow(w_dl[[2]])
    n_iter=nrow(no_dl[[2]])
    
    
    w_diff=round(w_dl[[2]][nrow(w_dl[[2]]),]-no_dl[[2]][nrow(no_dl[[2]]),],2)
    
    if(w_iter>n_iter) {
      for (type_plot in c("marginal","copula","joint")) {
        plot(1:nrow(w_dl[[2]]),w_dl[[2]][,type_plot],ylim=range(c(w_dl[[2]][,type_plot],no_dl[[2]][,type_plot])),type='l',col="red",main=paste(type_plot,w_diff[type_plot]))
        lines(1:nrow(no_dl[[2]]),no_dl[[2]][,type_plot],col="black")
      }
    } else {
      for (type_plot in c("marginal","copula","joint")) {
        plot(1:nrow(no_dl[[2]]),no_dl[[2]][,type_plot],ylim=range(c(w_dl[[2]][,type_plot],no_dl[[2]][,type_plot])),type='l',col="black",main=paste(type_plot,w_diff[type_plot]))
        lines(1:nrow(w_dl[[2]]),w_dl[[2]][,type_plot],col="red")
      }
    }
  }
  
  no_dl_outer=rbind(no_dl_outer,c(nrow(no_dl$par_history),no_dl$par_history[nrow(no_dl$par_history),],no_dl$log_lik_history[nrow(no_dl$log_lik_history),]))
  w_dl_outer=rbind(w_dl_outer,c(nrow(w_dl$par_history),w_dl$par_history[nrow(w_dl$par_history),],w_dl$log_lik_history[nrow(w_dl$log_lik_history),]))
  true_val=rbind(true_val,input_par[!is.na(input_par)])
  }
  , error=function(e) {
    print(e)
  }
  ) 
  
}

par(mfrow=c(3,3))
boxplot(c(no_dl_outer[,"iterations"], w_dl_outer[,"iterations"])~c(rep("Separate",nrow(no_dl_outer)),rep("Joint",nrow(w_dl_outer))),ylab="Iterations",xlab="Optimisation Method",main="Iterations")

par=names(input_par[!is.na(input_par)])
for (par_name in par) {
  if(par_name %in% c("theta","zeta")) {
    true=get_copula_dist(copula_dist)$copula_link[[paste(par_name,".linkfun",sep="")]](input_par[par_name])
  } else {
    true=margin_dist[[paste(par_name,".linkfun",sep="")]](input_par[par_name])
  }
  boxplot(c(no_dl_outer[,par_name], w_dl_outer[,par_name])~c(rep("Separate",nrow(no_dl_outer)),rep("Joint",nrow(w_dl_outer))),ylab="Parameter value",xlab="Optimisation Method",main=par_name
          ,ylim=range(c(no_dl_outer[,par_name], w_dl_outer[,par_name],true)))
  abline(h=true,col="red")
}
boxplot(c(no_dl_outer[,"ll margin"], w_dl_outer[,"ll margin"])~c(rep("Separate",nrow(no_dl_outer)),rep("Joint",nrow(w_dl_outer))),ylab="Parameter value",xlab="Optimisation Method",main="LogLik - Margin")
boxplot(c(no_dl_outer[,"ll copula"], w_dl_outer[,"ll copula"])~c(rep("Separate",nrow(no_dl_outer)),rep("Joint",nrow(w_dl_outer))),ylab="Parameter value",xlab="Optimisation Method",main="LogLik - Copula")
boxplot(c(no_dl_outer[,"ll total"], w_dl_outer[,"ll total"])~c(rep("Separate",nrow(no_dl_outer)),rep("Joint",nrow(w_dl_outer))),ylab="Parameter value",xlab="Optimisation Method",main="LogLik - Overall")

hist(w_dl_outer[,"ll total"]-no_dl_outer[,"ll total"],xlab="Change in LL",ylab="No. Simualtions",main=paste("Change in overall LogLik",round(mean(w_dl_outer[,"ll total"]-no_dl_outer[,"ll total"]),2)))
