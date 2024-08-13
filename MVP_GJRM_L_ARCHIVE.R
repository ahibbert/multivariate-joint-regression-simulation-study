

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



