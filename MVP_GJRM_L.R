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

#### 0. Load RAND data subset and transform to standard longitudinal dataset #######
#copula.link=list(logit,logit_inv,dlogit,log_2plus,log_2plus_inv,dlog_2plus)

copula.link=list(log,exp,dlog_inv)
if(length(copula.link)==3) {
  names(copula.link)=c("theta.linkfun","theta.linkinv","dthdeta")
} else {names(copula.link)=c("theta.linkfun","theta.linkinv","dthdeta","zeta.linkfun","zeta.linkinv","dzdeta")
}

dataset=loadDataset(simOption=3
                    ,plot_dist=FALSE
                    ,n=100
                    ,d=5
                    ,copula.family=3
                    ,copula.link=copula.link
                    ,qFUN=qZISICHEL
                    ,par.copula=c(0.5,0.5,0.5,0.5))

###########################
####### USER INPUT ########
###########################

#Formulas

setup=
GJRM_L_SETUP(
  mu.formula = ("response ~ time+as.factor(gender)"),#mu.formula = formula("response ~ as.factor(time)+as.factor(gender)+age")
  sigma.formula = ("~ time+as.factor(gender)"),#sigma.formula = formula("~ as.factor(time)+age")
  nu.formula = ("~ 1"),#nu.formula = formula("~ as.factor(time)+as.factor(gender)")
  tau.formula = ("~ 1"),#tau.formula = formula("~ age")
  theta.formula=("~1"),#theta.formula=formula("response~as.factor(gender)")
  zeta.formula=("~1"),#zeta.formula=formula("response~1")
  margin.family=ZISICHEL(mu.link = "log",sigma.link = "log",nu.link = "identity",tau.link = "logit"),
  dFUN=dZISICHEL,pFUN=pZISICHEL,
  copula.family="C",
  copula.link=copula.link,
  start="fit"
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
    verbose=TRUE,
    calc_d2=FALSE,
    phi=0.5,
    step_size=1,
    step_adjustment=0.5,
    step_reduc_count=3,
    crit_wk=0.0000001,
    true_val=NA,
    stopifnegative=FALSE,
    stopval=.1
)

post_optim=GJRM_POST_OPTIM(optim,setup,dataset)

margin.family=setup$margin.family
old_SE=sqrt(diag(vcov(setup$fitted_margins[[1]])))
old_SE=c(old_SE,c(setup$fitted_copulas$`1,2`$se,setup$fitted_copulas$`2,3`$se, setup$fitted_copulas$`1,2`$se2,setup$fitted_copulas$`2,3`$se2))
z_old=(start_par[0:length(old_SE)]-0)/old_SE
z_new=(end_par-0)/new_SE
SEs=cbind(start_par,old_SE,z_old,end_par,new_SE,z_new,optim$true_val)
rownames(SEs)=names(optim$true_val)
colnames(SEs)=c("Start Par","Old SE","Z Old","End Par","New SE","Z New","True Val")
print(round(SEs,3))



#### Post optimisation calculation of SEs ####

#Final SEs


#######ADD GAMLSS GLMM COMPARISON #######
gamlss_glmm_model <- gamlss(formula = formula("response ~ time+as.factor(gender)+random(as.factor(subject))")
                       , sigma.formula = sigma.formula
                       , nu.formula = nu.formula
                       , tau.formula = tau.formula
                       , family=margin_dist
                       , data=dataset,method=mixed(5,100))

gamlss_coeffs=vector()
for (name in c("mu","sigma","nu",'tau')) {
  gamlss_coeffs=c(gamlss_coeffs,coef(gamlss_glmm_model,what=name))
}

gamlss_vcov=sqrt(diag(vcov(gamlss_glmm_model)))
gamlss_coeffs=gamlss_coeffs[!grepl("random",names(gamlss_coeffs))]

SEs_plus_gamlss=cbind(SEs,cbind(c(gamlss_coeffs,rep(NA,nrow(SEs)-length(gamlss_coeffs))),c(gamlss_vcov,rep(NA,nrow(SEs)-length(gamlss_coeffs)))))

colnames(SEs_plus_gamlss)=c("Start Par","Old SE","Z Old","End Par","New SE","Z New","True Val","GAMLSS Coeff","GAMLSS SE")
print(round(SEs_plus_gamlss,3))






