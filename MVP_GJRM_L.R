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
simOption=2

#### 0. Load RAND data subset and transform to standard longitudinal dataset #######
loadDataset(simOption)

#### INPUT: Parameters ####

#Formulas

mu.formula = formula("response ~ 1")#mu.formula = formula("response ~ as.factor(time)+as.factor(gender)+age")
sigma.formula = formula("~ 1")#sigma.formula = formula("~ as.factor(time)+age")
nu.formula = formula("~ 1")#nu.formula = formula("~ as.factor(time)+as.factor(gender)")
tau.formula = formula("~ 1")#tau.formula = formula("~ age")
theta.formula=formula("response~1")#theta.formula=formula("response~as.factor(gender)")
zeta.formula=formula("response~1")#zeta.formula=formula("response~1")

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
                       , family=margin_dist,data=dataset,method=RS(50))

plot(gamlss_model)

fitted_copulas=fitted_margins=list()
fitted_margins[[1]]=gamlss_model

#Extract margin model matrix
mm_mar=list()
for (parameter in c("mu","sigma","nu","tau")) {
  print(head(model.matrix(fitted_margins[[1]],what=parameter)))
  mm_mar[[parameter]]= model.matrix(fitted_margins[[1]],what=parameter)
  print(head(mm_mar[[parameter]]))
}

#Fit starting copula with just an intercept for each parameter
start_par<-extract_parameters_from_fitted(fitted_margins,fitted_copulas=NA,copula_link)
copula_model<-fit_copula_single(dataset=dataset,mm_mar=mm_mar,input_par=start_par,dFUN=dFUN,pFUN=pFUN,copula_dist=copula_dist)
fitted_copulas[[1]]=copula_model
start_par<-extract_parameters_from_fitted(fitted_margins,fitted_copulas=fitted_copulas,copula_link)


#hacky model matrix
#mm_cop=generate_cop_model_matrix(dataset=dataset,formula=theta.formula,zeta.formula=zeta.formula,time="time")
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
)

print(results_start$log_lik_results)

#### Newton Raphson optimisation ####

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
stopifnegative=TRUE
while ((abs(change) > .1*phi) & run_counter <= 100) { #
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
                                )

  end_log_lik=results$log_lik_results

  if(first_run) {change=1} else {
    change=end_log_lik["Total"]-start_log_lik["Total"]
    log_lik_output=c(start_log_lik["Total"],end_log_lik["Total"],change,phi_inner)
    names(log_lik_output) = c("Start LogLik","End LogLik","Change","phi")
    print(log_lik_output)
  }
  
  iteration_out=newton_raphson_iteration(results,input_par,phi=phi_inner,step_size=step_size,verbose=verbose_option)
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

final_results= calc_joint_likelihood(input_par =end_par_matrix[which.max(end_loglik_matrix),] 
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
                               , calc_d2=TRUE
)


final_iteration_out=newton_raphson_iteration(final_results,end_par,phi=phi_inner,step_size=step_size,verbose=verbose_option,calc_d2=TRUE)

new_SE=final_iteration_out$se_par
old_SE=sqrt(diag(vcov(fitted_margins[[1]])))
old_SE=c(old_SE,c(fitted_copulas$`1,2`$se,fitted_copulas$`2,3`$se, fitted_copulas$`1,2`$se2,fitted_copulas$`2,3`$se2))

z_score_old=(start_par[0:length(old_SE)]-0)/old_SE
z_score_new=(end_par-0)/new_SE
SEs=cbind(start_par,old_SE,z_score_old,end_par,new_SE,z_score_new)
print(round(SEs,3))













