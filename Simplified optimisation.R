###########MVP OPTIMISATION

library(VineCopula)
library(gamlss)
source("common_functions.R")

#Defining link functions for copula manually
copula_link=copula.link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv)
names(copula_link)=names(copula.link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")
copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE
if(two_par_cop) {
  names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}

dataset=loadDataset(simOption=4
                    ,plot_dist=FALSE
                    ,n=1000
                    ,d=3
                    ,copula.family=BiCopName("C")
                    ,copula.link=copula_link
                    ,qFUN=qEXP
                    ,par.copula=c(2,2)
                    ,par.margin=c(1,1,-1,-3))

plotDist(dataset,dist="EXP")

#Start parameters
input_par=c(0.5,1); names(input_par)=c("mu","theta");margin_dist=EXP(); copula_dist="C"

###OUTER ITERATION - currently in non-link form as it seems to work

change=1
while (abs(change)>0.0001) {
  outer_optim_output=optim_outer(par=input_par,dataset,margin_dist,copula_dist,copula_link=copula_link)
  input_par=outer_optim_output$par_end
  change=sum(outer_optim_output$par_change)
}
  
###INNER ITERATION

#Fit one parameter at a time like RS maybe?

for(par_name in names(input_par)) {
  par_value=input_par[par_name]
  
  
}


score_function
inner_iterations
backfitting
###pseudo code
for (parameter in parameters) {
  mm=mm_cop[[parameter]]
  X=as.matrix(mm)
  W=diag(copula_score$w_k[,j])
  e_k_par=e_k[,j]
  beta_new_cop[[copula_parameters[j]]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
  rownames(beta_new_cop[[parameter]])= paste(paste(parameter,sep=" "),colnames(mm),sep=".")
  #names(beta_new_cop[[j]])= paste(parameter_names[j],margin_number,rownames(beta_new_cop[[j]]),sep=".")
  if (calc_d2 == TRUE) {
    hess_list[[(copula_parameters[j])]]=-t(X)%*%W%*%X
  }
}





