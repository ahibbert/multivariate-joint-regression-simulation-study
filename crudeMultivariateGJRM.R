#######################Required functions####################
source("common_functions.R")
library(copula)
library(gamlss)
options(scipen=999)

############################1. Generate Data#################

dist="NO";mu_vector=c(1,2,3);sigma_vector=c(1,2,3);rho_vector=c(.4,.5,.6);n=1000

#BIVARIATE### dist="NO";mu_vector=c(1,2);sigma_vector=c(1,2);rho_vector=c(.7);n=1000

rho_matrix=matrix(0,ncol=length(mu_vector),nrow=length(mu_vector))
rho_matrix[lower.tri(rho_matrix,diag=FALSE)]<-rho_vector
rho_matrix[upper.tri(rho_matrix,diag=FALSE)]<-rho_vector

dataset<-generateMvtDist(dist,mu_vector,sigma_vector,rho_matrix)

plotDist(dataset,dist)

dataset_matrix<-dataset[dataset[,"time"]==0,"random_variable"]
for(i in 1:length(mu_vector)) {
  dataset_matrix<-cbind(dataset_matrix,dataset[dataset[,"time"]==i,"random_variable"])  
}

orig_par<- par <- c(mu_vector,sigma_vector,rho_vector)

#copula_link<-BI(mu.link="logit")
copula_link<-NO(mu.link="identity")
margin_pdf=NO(mu.link="identity",sigma.link="log");copula_link<-BI(mu.link="logit")
m=ncol(dataset_matrix)
eta_par <- c(
  margin_pdf$mu.linkfun(par[1:m]),
  sigma_vector=margin_pdf$sigma.linkfun(par[(m+1):(2*m)]),
  rho_vector=copula_link$mu.linkfun(par[(2*m+1):length(par)]))

##################################################

calcJointLikelihood <- function(par,dataset_matrix,margin_pfun,margin_dfun,copula_pdf,margin_pdf,copula_link) {
  
  library(copula)
  library(gamlss)
  
  #Test parameters for stepthrough: par=eta_par;margin_pdf=dNO;copula_pdf=normalCopula; i=2; j=1;copula_link
  
  #Extracting vectors for each parameter
  #par=c(mu_vector,sigma_vector,rho_vector)
  m=ncol(dataset_matrix)
  dist<-margin_pdf
  mu_vector=dist$mu.linkinv(par[1:m])
  sigma_vector=dist$sigma.linkinv(par[(m+1):(2*m)])
  rho_vector=copula_link$mu.linkinv(par[(2*m+1):length(par)])
  rho_matrix=matrix(0,ncol=m,nrow=m)
  rho_matrix[lower.tri(rho_matrix,diag=FALSE)]<-rho_vector
  rho_matrix[upper.tri(rho_matrix,diag=FALSE)]<-rho_vector
  
  #Initialising results matrices
  margin_p <- margin_d <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  copula_p <- copula_d <- matrix(0,ncol=(m*(m-1))/2,nrow=nrow(dataset_matrix))
  
  ####Fit each margin with the given parameters using gamlss and passed in distribution cfunction margin_pdf
  for (i in 1:ncol(dataset_matrix)) {
    margin_p[,i] <- margin_pfun(dataset_matrix[,i],mu=mu_vector[i],sigma=((sigma_vector[i])))
    margin_d[,i] <- margin_dfun(dataset_matrix[,i],mu=mu_vector[i],sigma=((sigma_vector[i])))
    #########Update with a generic gamlss function ideally
  }
  
  ####Fit each copula with   
  z=1
  colnames_temp=c()
  for (i in 1:ncol(dataset_matrix)) {
    for (j in 1:ncol(dataset_matrix)) {
      if (i>j) {
        copula_d[,z]<-dCopula(margin_p[,c(i,j)],copula=copula_pdf(rho_matrix[i,j]))
        colnames_temp[z]<-paste(i,",",j,sep="")
        z=z+1
      }
    }
  }
  colnames(copula_d)<-colnames_temp
  log_likelihood=sum(log(margin_d))+sum(log(copula_d))
  
  return(log_likelihood)
  
}

iterateParametersNewtownRaphson <- function(par,dataset_matrix,margin_pdf,margin_pfun,margin_dfun,copula_pdf,step_size,copula_link) {
  
  library(copula)
  library(gamlss)
  library(VineCopula)
  
  #Test parameters for stepthrough: 
  ######################################DELETE AFTER TESTING
  
  #par=eta_par*.5;margin_pdf=NO(mu.link="identity",sigma.link="log");margin_pfun=pNO; margin_dfun=dNO; copula_pdf=normalCopula;step_size=1;i=2;j=1
  
  #Extracting vectors for each parameter
  m=ncol(dataset_matrix)
  mu_vector=par[1:m]
  sigma_vector=par[(m+1):(2*m)]
  rho_vector=par[(2*m+1):length(par)]
  rho_matrix=matrix(0,ncol=m,nrow=m)
  rho_matrix[lower.tri(rho_matrix,diag=FALSE)]<-rho_vector
  rho_matrix[upper.tri(rho_matrix,diag=FALSE)]<-rho_vector
  
  #Initialising results matrices
  margin_p <- margin_d <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  margin_dldm <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  margin_dldd <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  margin_d2ldm2 <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  margin_d2ldmdd <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  margin_d2ldd2 <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  
  margin_dmdnu <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  margin_dddnu <- matrix(0,ncol=m,nrow=nrow(dataset_matrix))
  
  dist<-margin_pdf
  
  ###########LIKELIHOOD 
  for (i in 1:ncol(dataset_matrix)) {
    margin_p[,i] <- margin_pfun(dataset_matrix[,i],mu=mu_vector[i],sigma=dist$sigma.linkinv(sigma_vector[i]))
    margin_d[,i] <- margin_dfun(dataset_matrix[,i],mu=mu_vector[i],sigma=dist$sigma.linkinv(sigma_vector[i]))
  }
  
  ###First and second derivatives
  for (i in 1:ncol(dataset_matrix)) { ###For each margin
    for (j in 1:nrow(dataset_matrix)) { ##For each observation in the margin
      sigma=(sigma_vector[i]);mu=mu_vector[i];y=dataset_matrix[j,i]
      margin_dldm[j,i]<-    dist$dldm(y,dist$mu.linkinv(mu),(dist$sigma.linkinv(sigma)))
      margin_dldd[j,i]<-    dist$dldd(y,dist$mu.linkinv(mu),(dist$sigma.linkinv(sigma)))
      margin_d2ldm2[j,i]<-  dist$d2ldm2((dist$sigma.linkinv(sigma)))
      margin_d2ldmdd[j,i]<- dist$d2ldmdd(y)
      margin_d2ldd2[j,i]<-  dist$d2ldd2((dist$sigma.linkinv(sigma)))
      margin_dmdnu[j,i]<- dist$mu.dr((mu))
      margin_dddnu[j,i]<- dist$sigma.dr((sigma))
    }
  }
  
  
  dCopulaWrapper <- function(x,u,copula_pdf) {
    return(dCopula(u,copula=copula_pdf(copula_link$mu.linkinv(x)),log=TRUE))
    #dCopulaWrapper(x=rho_matrix[i,j],u=margin_p[,c(i,j)],copula_pdf)
    ###TEST #### dCopula(margin_p[,c(i,j)],copula=copula_pdf(copula_link$mu.linkinv(rho_matrix[i,j])),log=TRUE) == dCopulaWrapper(x=rho_matrix[i,j],u=margin_p[,c(i,j)],copula_pdf)
  }
  
  z=1
  colnames_temp=c()
  copula_p <- copula_d <- copula_dtdnu<- copula_dfdt <- copula_dfdt2 <-copula_dldt <- copula_dldt2 <- nd_copula_dldt <- nd_copula_dldt2 <- matrix(NA,ncol=m*(m-1)/2,nrow=nrow(dataset_matrix))
  for (i in 1:ncol(dataset_matrix)) {
    for (j in 1:ncol(dataset_matrix)) {
      if (i>j) {
        copula_d[,z]<-dCopula(margin_p[,c(i,j)],copula=copula_pdf(copula_link$mu.linkinv(rho_matrix[i,j])))
        
        for (l in 1:nrow(margin_p[,c(i,j)])) {
          nd_copula_dldt[l,z]<- grad(func=dCopulaWrapper,x=rho_matrix[i,j],u=margin_p[l,c(i,j)],copula_pdf=copula_pdf)
          nd_copula_dldt2[l,z]<- hessian(func=dCopulaWrapper,x=rho_matrix[i,j],u=margin_p[l,c(i,j)],copula_pdf=copula_pdf)
        }

        copula_dfdt[,z]<-BiCopDeriv(margin_p[,c(i)],margin_p[,c(j)],family=1,par=copula_link$mu.linkinv(rho_matrix[i,j]),deriv="par")
        copula_dfdt2[,z]<-BiCopDeriv2(margin_p[,c(i)],margin_p[,c(j)],family=1,par=copula_link$mu.linkinv(rho_matrix[i,j]),deriv="par")
        copula_dtdnu[,z]<-copula_link$mu.dr(rho_matrix[i,j])
        colnames_temp[z]<-paste(i,",",j,sep="")
        z=z+1
      }
    }
  }
  colnames(copula_d)<-colnames_temp
  colnames(copula_dfdt)<-colnames_temp
  colnames(copula_dfdt2)<-colnames_temp
  
  copula_dldt=(1/copula_d)*copula_dfdt
  copula_dldt2=(1/(copula_d^2))*((copula_d*copula_dfdt2)-(copula_dfdt^2))
  
  ##############WHILE USING NUMERICAL DERIVATIVES
  copula_dldt=nd_copula_dldt
  copula_dldt2=nd_copula_dldt2

  f=cbind(margin_d2ldm2,margin_d2ldd2,copula_dfdt2)
  dalldnu=cbind(margin_dmdnu,margin_dddnu,copula_dtdnu)
  w=-1*f*dalldnu*dalldnu
  u=cbind(margin_dldm*margin_dmdnu,margin_dldd*margin_dddnu,copula_dldt*(copula_dtdnu))
  colMeans((1/w)*u)
  #eta_par
  
  z=par + colMeans((1/w)*u)
  
  new_par<- step_size*z+(1-step_size)*par
  
  results <- list()

  results[[1]] <- (par)
  results[[2]] <- new_par
  
  log_likelihood_orig=calcJointLikelihood(par=par,dataset_matrix,margin_pfun,margin_dfun,copula_pdf,margin_pdf,copula_link)
  results[[3]] <- log_likelihood_orig
  
  log_likelihood_new=calcJointLikelihood(par=new_par,dataset_matrix,margin_pfun,margin_dfun,copula_pdf,margin_pdf,copula_link)
  results[[4]] <- log_likelihood_new
  
  results[[5]] <- colMeans(u)
  results[[6]] <- colMeans(f)
  
  par_trans_end=c(
    dist$mu.linkinv(z[1:m]) ,
    dist$sigma.linkinv(z[(m+1):(2*m)]),
    copula_link$mu.linkinv(z[(2*m+1):length(par)]) 
  )
  
  results[[7]] <- par_trans_end
  
  #########par=new_par
  
  names(results) <- c("par_start","par_end","ll_orig","ll_new","grad","hessian","par_trans_end")
  
  return(results)
  
  }

gamlssStartingFit <- function(dataset_matrix,margin_pdf) {
  
  fits<-list()
  m=ncol(dataset_matrix)
  pars<-matrix(ncol=2,nrow=m)
  for (i in 1:m) {
    fits[[i]]<-coefAll(gamlss(dataset_matrix[,i]~1,family=margin_pdf))
    pars[i,]<-c(fits[[i]]$mu,fits[[i]]$sigma)
  }
  
  return(c(pars[,1],pars[,2]))
  
}

tauStartingFit <- function(dataset_matrix,cor_method="kendall") {
  
  start_cor <- vector(length=ncol(dataset_matrix))
  z=1
  for (i in 1:ncol(dataset_matrix)) {
    for (j in 1:ncol(dataset_matrix)) {
      if (i>j) {
        start_cor[z]<-cor(dataset_matrix[,i],dataset_matrix[,j],method=cor_method)
        z=z+1
        }
    }
  }
  return(start_cor)
}

#################2. Run Functions#################################

calcJointLikelihood(eta_par,dataset_matrix,margin_pfun=pNO,margin_dfun=dNO,copula_pdf=normalCopula,margin_pdf,copula_link)

test_spread=.25
start_par=eta_par*(1+runif(length(par),-1*test_spread,1*test_spread))
print(start_par)
step_size=.5; change_in_ll <- 10; results<-list();results$par_end=start_par
z=0
while (change_in_ll>.1) {
  results<-iterateParametersNewtownRaphson(results$par_end,dataset_matrix,margin_pdf,margin_pfun=pNO,margin_dfun=dNO,copula_pdf=normalCopula,step_size,copula_link)
  change_in_ll <- results$ll_new-results$ll_orig
  print(paste("Log Likelihood: ", results$ll_new," | Par Est: ",sep=""))
  print(c(results$par_trans_end))
  z=z+1
  
  if(z<5) {step_size=step_size*.5}
}


test_par<-c(1,2,3,log(1),log(2),log(3),copula_link$mu.linkfun(c(.3,.5,.7)))


rho_test=1:9*.1
error=-10:10*.05
rhoLikelihood=matrix(NA,nrow=length(rho_test)*length(error),ncol=3)
z=1
for (j in 1:length(error)) {
  for (i in 1:length(rho_test)) {
    rhoLikelihood[z,1]=rho_test[i]
    rhoLikelihood[z,2]=error[j]
    rhoLikelihood[z,3]=calcJointLikelihood(c(eta_par[1:(m*2)]*(1+error[j]),copula_link$mu.linkfun(rho_test[i])),dataset_matrix,margin_pfun=pNO,margin_dfun=dNO,copula_pdf=normalCopula,margin_pdf,copula_link)
    z=z+1
  }
}
par(mfrow=c(sqrt(length(error))+1,sqrt(length(error))+1))
for (j in 1:length(error)) {
  plot(rhoLikelihood[rhoLikelihood[,2]==error[j],c(1,3)],main=paste("Error:",error[j]," | rho: ", copula_link$mu.linkinv(eta_par[(2*m+1):length(par)])),xlab="Estimated Rho",ylab="Likelihood")
}












######################

library(numDeriv)
hessian(calcJointLikelihood,rep(.5,length(par)),dataset_matrix=dataset_matrix,margin_pdf=dNO,copula_pdf=normalCopula)








###Returns a list of fitted gamlss objects to the margins
fitMargins <- function(dataset_matrix,dist) {
  
  marginFits <- list()
  
  for (i in 1:ncol(dataset_matrix)) {
    marginFits[[i]] <- gamlss(dataset_matrix[,i]~1,family=dist)
  }
  
  return(marginFits)
  
}

###Actually get the hessians
marginFits<- fitMargins(dataset_matrix,dist)
get.K()

marginFits[[1]]


library(numDeriv)

help("numDeriv-package")