########### Starting parameters
set.seed(1000)
options(scipen=999)
library(VineCopula)
library(gamlss)
source("common_functions.R")

### OUTER ITERATIONS ####

######## Select distributions and range for parameters
copula_dist="C"
#margin_dist=EXP();      min_par=c(1,1);         max_par=c(10,10) #EXP #Note these can be times by 2 or half as start parameters so be careful
#margin_dist=GA();       min_par=c(1,1,2);       max_par=c(10,10,10) #GA
margin_dist=ZISICHEL(); min_par=c(1,1,-2,.1,2); max_par=c(10,10,2,.5,10) # ZISICHEL

num_par=length(margin_dist$parameters)+if(copula_dist=="t") {2} else {1}
copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE #copula_link=copula.link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv)
if(two_par_cop) {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}

final_lik_separate=final_lik_joint=matrix(ncol=5,nrow=0)
final_par_separate=final_par_joint=matrix(ncol=2+num_par,nrow=0)
final_true_val=matrix(ncol=num_par,nrow=0)

n=1000; d=3; crit_lik_change=0.05
for (random_simulation in 1:100) {
  
  print(random_simulation)
  plot_results=FALSE
  
  sim_parameters=runif(num_par,min_par,max_par)
  
  a=sim_parameters[length(sim_parameters)];b=sim_parameters[(1:(length(sim_parameters)-1))];d=3
  true_val=c(b,a)
  
  start_par=true_val*runif(length(true_val),0.5,2)
  names(start_par)=c(names(margin_dist$parameters),"theta"); input_par=start_par
  dataset=loadDataset(simOption=5,n=n,d=d, copula.family=BiCopName(copula_dist), copula.link=copula_link, qFUN=paste("q",margin_dist$family[1],sep="")
                      , par.copula=rep(a,d-1), par.margin=b,plot_dist = FALSE)
  #plotDist(dataset,dist=margin_dist)
  
  
  #####Find close starting parameters (just outer iterations with no covariates)
  log_lik_history=matrix(ncol=3+2,nrow=0)
  par_history=matrix(ncol=length(input_par)+2,nrow=0)
  
  ### Run fit for separate and joint optimisation
  for (use_dlcopdpar in c(FALSE,TRUE)) {
    
    ### CORE ITERATION
    input_par=start_par
    change=1;log_lik_start=0;log_lik_change=1000;run_counter=1;step_size=1;step_adjustment=0.5
    while (abs(log_lik_change)>crit_lik_change) {
      
      step_size=step_size*(step_adjustment^min(5,run_counter))
      
      #Run optimisation
      outer_optim_output=optim_outer(par=input_par,dataset,margin_dist,copula_dist,copula_link=copula_link,use_dlcopdpar=use_dlcopdpar,verbose=FALSE,step_size=step_size)
      
      #Capture outputs
      input_par=outer_optim_output$par_end
      change=sum(outer_optim_output$par_change)
      #print(outer_optim_output$log_lik)
      run_counter=run_counter+1
      log_lik=outer_optim_output$log_lik["joint"]
      log_lik_change=log_lik-log_lik_start
      log_lik_start=log_lik
      
      #Capture changes in parameters
      copula_deriv=if(use_dlcopdpar==TRUE){1}else{0}
      
      log_lik_history=rbind(log_lik_history,c(copula_deriv,run_counter,outer_optim_output$log_lik))
      par_history=rbind(par_history,c(copula_deriv,run_counter,input_par))
    }
    colnames(log_lik_history)[1:2]=colnames(par_history)[1:2]=c("use_dlcopdpar","run_counter")
  }
  
  #Plot likelihood and parameters
  if(plot_results==TRUE) {
    plot.new()
    par(mfrow=c(1,3))
    for (i in colnames(log_lik_history)[3:5]) {
      plot( log_lik_history[log_lik_history[,"use_dlcopdpar"]==1,i],xlab="LogLik",ylab="Iteration",main=i,type = "l",col="blue",xlim=c(1,max(log_lik_history[,"run_counter"])),ylim=range(log_lik_history[,i]))
      lines(log_lik_history[log_lik_history[,"use_dlcopdpar"]==0,i],xlab="LogLik",ylab="Iteration",main=i,type = "l",col="red",xlim=c(1,max(log_lik_history[,"run_counter"])),ylim=range(log_lik_history[,i]))
      legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))
      
    }
    
    par(mfrow=c(1,(ncol(par_history)-2)))
    for (i in 1:(ncol(par_history)-2)) {
      par_dlcop=par_history[par_history[,"use_dlcopdpar"]==1,]
      par_nodlcop=par_history[par_history[,"use_dlcopdpar"]==0,]
      plot(par_dlcop[,i+2],col="blue",type="l",main=colnames(par_history)[i+2],ylim=range(c(par_dlcop[,i+2],par_nodlcop[,i+2],true_val[i])),ylab="Parameter estimate")
      lines(par_nodlcop[,i+2],col="red",type="l")
      abline(h=true_val[i])
      legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))
    }
  }

  lik0=log_lik_history[log_lik_history[,"use_dlcopdpar"]==0,]
  par0=par_history[par_history[,"use_dlcopdpar"]==0,]
  lik1=log_lik_history[log_lik_history[,"use_dlcopdpar"]==1,]
  par1=par_history[par_history[,"use_dlcopdpar"]==1,]
  
  final_lik_separate=rbind(final_lik_separate,lik0[nrow(lik0),])
  final_lik_joint=rbind(final_lik_joint,lik1[nrow(lik1),])
  final_par_separate=rbind(final_par_separate,par0[nrow(par0),])
  final_par_joint=rbind(final_par_joint,par1[nrow(par1),])
  final_true_val=rbind(final_true_val,true_val)
}

plot.new()
par(mfrow=c(length(start_par),2))
for (i in 1:length(start_par)) {
  plot(final_par_separate[,2+i]~final_true_val[,i],xlim=range(c(final_true_val[,i])),ylim=range(final_true_val[,i]),main=paste("Separate Optimisation:",names(start_par)[i]),xlab="True Value",ylab="Estimate")
  abline(0,1,col="red")
  plot(final_par_joint[,2+i]~final_true_val[,i],xlim=range(c(final_true_val[,i])),ylim=range(final_true_val[,i]),main=paste("Joint Optimisation:",names(start_par)[i]),xlab="True Value",ylab="Estimate")
  abline(0,1,col="red")
}

#Remove extreme outliers (optimisation problems)
outliers=which(final_lik_joint[,5]-final_lik_separate[,5]>100)
final_lik_joint=final_lik_joint[-outliers,]
final_lik_separate=final_lik_separate[-outliers,]

plot.new()
par(mfrow=c(1,2))
lm_coef=lm(final_lik_joint[,5]~final_lik_separate[,5])$coefficients
plot(final_lik_joint[,5]~final_lik_separate[,5],xlab="Separate",ylab="Joint",main=paste(c("Joint v Separate Likelihood: ","Intercept:","Slope:"), c("",lm_coef)),xlim=c(-12000,0))
abline(a=0,b=1,col="red")
hist(final_lik_joint[,5]-final_lik_separate[,5],main="Likelihood Improvement",xlim=c(-100,100),breaks=100)

plot(final_lik_joint[,5]/final_lik_separate[,5]~final_true_val[,2])

#### INNER ITERATION ####

margin_dist=EXP(); copula_dist="C"
copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE #copula_link=copula.link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv)
if(two_par_cop) {
  names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}
true_val=c(1,1,1,1); n=1000; d=3
dataset=loadDataset(simOption=3,n=n,d=d, copula.family=BiCopName("C"), copula.link=copula_link, qFUN=qEXP
                    , par.copula=rep(1,d-1), par.margin=c(1,1,1,1),plot_dist = FALSE)
plotDist(dataset,"EXP")

#Fit one parameter at a time like RS maybe?

mm=create_model_matrices(
  mu.formula = ("response ~ time+as.factor(gender)"),
  sigma.formula = ("~ 1"),
  nu.formula = ("~ 1"),
  tau.formula = ("~ 1"),
  theta.formula=("~1"),
  zeta.formula=("~1"),
  margin.family=margin_dist,copula.family=copula_dist,copula.link=copula_link
)

###CALCULATE eta,dldpar,d2ldpar,dpardeta

#First run set all covariates except intercept to zero and setup starting vector
input_par=c(2,3); names(input_par)=c("mu","theta")
par_eta=input_par
par_cov=as.numeric(vector())
for (par_name in names(mm)) {
    par_cov_single=as.numeric(vector(length=length(colnames(mm[[par_name]]))))
    names(par_cov_single)=paste(par_name,colnames(mm[[par_name]]),sep=".")
    par_cov_single[1]=par_eta[par_name]
    if(length(par_cov_single)>1) {
      par_cov_single[2:length(par_cov_single)]=0
    }
    par_cov=c(par_cov,par_cov_single)
}

include_dlcopdpar=TRUE
first_outer_run=TRUE
outer_log_lik_change=outer_start_log_lik=outer_end_log_lik=0
inner_stop_crit=outer_stop_crit=0.1
start_step_size=0.5; step_adjustment=1; max_steps=5
while (first_outer_run==TRUE | abs(outer_log_lik_change)>outer_stop_crit) {
  print("STARTING NEW OUTER ITERATION")
  if(first_outer_run==FALSE) {
    print(c(outer_start_log_lik,outer_end_log_lik,outer_log_lik_change))
  }
  first_outer_run=FALSE
  
  for (par_name in names(mm)) {
    print(paste("INNTER ITERATION: Parameter:",par_name))
    first_inner_run=TRUE; change_log_lik=0
    run_counter=1
    step_size=start_step_size
    while (first_inner_run==TRUE | abs(change_log_lik)>inner_stop_crit) {
    
      first_inner_run=FALSE
      eta_out=calc_eta(par_cov,mm,margin_dist,copula_link)
      eta=eta_out$eta; eta_dr=eta_out$eta_dr; eta_inv=eta_out$eta_inv
      
      #Given par_cov and mm, calculate eta, dldpar, d2ldpar, dpardeta
      
      calc_lik_out=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist)
      log_lik=calc_lik_out$log_lik; margin_d=calc_lik_out$margin_d; margin_p=calc_lik_out$margin_p; margin_deriv=calc_lik_out$margin_deriv; copula_d=calc_lik_out$copula_d; copula_p=calc_lik_out$copula_p; Fx_1_2=calc_lik_out$Fx_1_2;order_copula=calc_lik_out$order_copula
      
      par1=eta_inv[["theta"]]
      if(!"zeta" %in% names(eta_inv)) {par2=eta_inv[["theta"]]*0} else {par2=eta_inv[["zeta"]]}
      dldth=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par",log=TRUE)
      dcdth=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par",log=FALSE)
      d2cdth=BiCopDeriv2( Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par")
      d2ldth2=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
      if("zeta" %in% names(eta_inv)) {
        dldz=BiCopDeriv(    Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par2",log=TRUE)
        dcdz=BiCopDeriv(    Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par2",log=FALSE)
        d2cdz=BiCopDeriv2(  Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par2")  
        d2ldz2=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)
        
        d2cdthdz=BiCopDeriv2(  Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="par1par2")  
        d2ldthdz=(d2cdthdz*copula_d-dcdth*dcdz)/(copula_d^2)
      }
      dcdu1=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="u1",log=FALSE)
      dcdu2=BiCopDeriv(   Fx_1_2[,1],Fx_1_2[,2],family = as.numeric(BiCopName(copula_dist)),par=par1,par2=par2,deriv="u2",log=FALSE)
      
      dldth[!is.finite(dldth)]=0
      d2ldth2[!is.finite(d2ldth2)]=0
      if("zeta" %in% names(eta_inv)) {
        dldz[!is.finite(dldz)]=0
        d2ldz2[!is.finite(d2ldz2)]=0
        d2ldthdz[!is.finite(d2ldthdz)]=0
      }
      
      ### Calculate copula derivatives w.r.t margin parameters
      if(!par_name %in% c("mu","sigma","nu","tau")) {
        if(!"zeta" %in% names(eta_inv)) {
          d1=dldth
          d2=d2ldth2
        } else {
          d1=cbind(dldth,dldz)
          d2=cbind(d2ldth2,d2ldz2)
        }
      } else {
        margin_par=names(mm)[names(mm) %in% c("mu","sigma","nu","tau")]
        response=dataset$response
        #Extract margin calculations for F(x), f(x), response and derivatives at time 1 and time 2, join to copula values for time 1 and time 2
        margin_deriv_1=margin_deriv_2=margin_deriv_2cross=matrix(ncol=length(margin_par),nrow=length(response))
        for (i in 1:length(margin_par)) {
          margin_deriv_1[,i]=margin_deriv[grepl("dld",names(margin_deriv))][[i]]
          #margin_deriv_2[,i]=margin_deriv[grepl("d2ld",names(margin_deriv))&endsWith(names(margin_deriv),"2")][[i]]
        }
        colnames(margin_deriv_1)=paste("dld",margin_par,sep="")
        #colnames(margin_deriv_2)=paste("d2ld",names(margin_par),sep="")
        
        order_margin=dataset[,c("time","subject")]
        margin_components=cbind(order_margin,response,margin_p,margin_d,margin_deriv_1)
        margin_components_Ft_plus=margin_components
        margin_components_Ft_plus$time=margin_components_Ft_plus$time-1
        margin_plus=merge(margin_components,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)
        
        copula_components=cbind(order_copula,dcdu1,dcdu2,copula_d)
        copula_merged=merge(copula_components,margin_plus,by.x=c("time1","subject1"),by.y=c("time","subject"),all.x=TRUE)
        
        #Calculate copula derivative with respect to marginal parameters
        input=copula_merged
        dlcopdpar=matrix(0,nrow=nrow(input),ncol=length(margin_par))
        i=1
        for (par_name in margin_par) {
          
          #Take parameters from input for clarity
          dc_tplus_du_t=input[,"dcdu1"]
          dc_tplus_du_tplus=input[,"dcdu2"]
          l_t=input[,paste(paste("dld",par_name,sep=""),".x",sep="")]
          l_t_plus=input[,paste(paste("dld",par_name,sep=""),".y",sep="")]
          x_t=input[,"response.x"]
          x_t_plus=input[,"response.y"]
          f_t=input[,"margin_d.x"]
          f_t_plus=input[,"margin_d.y"]
          du_t_dmu=x_t*f_t*l_t
          du_t_plus_dmu=x_t_plus*f_t_plus*l_t_plus
          c_tplus=input[,"copula_d"]
          
          du_t_dmu=x_t*f_t*l_t
          du_t_plus_dmu=x_t_plus*f_t_plus*l_t_plus
          
          dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
          dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
          dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
          dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
          dcdmu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
          dcdmu_tplus[is.nan(dcdmu_tplus)|is.na(dcdmu_tplus)]=0
          
          dlcopdpar[,i]=dcdmu_tplus
          
        }
        colnames(dlcopdpar)=paste("dlcopd",margin_par,sep="")
        
        par_dlcopdpar=dlcopdpar[,paste("dlcopd",margin_par,sep="")]
        merged_dlcopdpar=merge(cbind(order_copula,par_dlcopdpar),cbind(order_copula,par_dlcopdpar),by.x=c("time2","subject2"),by.y=c("time1","subject1"),all=TRUE)
        merged_dlcopdpar[is.na(merged_dlcopdpar)]=0
        d1_cop=0.5*(merged_dlcopdpar[,"par_dlcopdpar.x"]+merged_dlcopdpar[,"par_dlcopdpar.y"])
        
        margin_deriv_subnames=c("m","d","v","t")
        names(margin_deriv_subnames)=c("mu","sigma","nu","tau")
        d1=margin_deriv[grepl(paste("dld",margin_deriv_subnames[par_name],sep=""),names(margin_deriv))][[1]]
        d2=margin_deriv[grepl(paste("d2ld",margin_deriv_subnames[par_name],sep=""),names(margin_deriv))][[1]]
        d1_plus_cop=d1 + if(include_dlcopdpar==TRUE){d1_cop} else {0*d1_cop}
        d1=d1_plus_cop
      }
      
      score=score_function_v2(eta=eta[[par_name]],dldpar=d1,d2ldpar=-(d1*d1),dpardeta=eta_dr[[par_name]])
      
      X=as.matrix(mm[[par_name]])
      W=diag(as.vector(score$w_k))
      z_k=score$z_k
      beta_start=par_cov[paste(paste(par_name,sep=" "),colnames(mm[[par_name]]),sep=".")]
      beta=beta_start*(1-step_size) + (step_size)*(as.vector(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z_k))
      names(beta)=paste(paste(par_name,sep=" "),colnames(mm[[par_name]]),sep=".")
      
      par_cov_new=par_cov
      par_cov_new[names(beta)]=beta
      
      eta_out=calc_eta(par_cov_new,mm,margin_dist,copula_link)
      
      #Update parameters
      eta=eta_out$eta; eta_dr=eta_out$eta_dr; eta_inv=eta_out$eta_inv
      par_cov=par_cov_new
      
      calc_lik_out_end=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist)
      
      print("Start Log Lik:")
      print(calc_lik_out$log_lik)
      print("End Log Lik")
      print(calc_lik_out_end$log_lik)
      
      change_log_lik=calc_lik_out_end$log_lik["joint"]-calc_lik_out$log_lik["joint"]
      step_size = (step_adjustment^min(run_counter,max_steps))*start_step_size
      run_counter=run_counter+1
    }
  }
  outer_start_log_lik=outer_end_log_lik
  outer_end_log_lik=calc_lik_out_end$log_lik["joint"]
  outer_log_lik_change=outer_start_log_lik-outer_end_log_lik
  
  if(abs(outer_log_lik_change)<=outer_stop_crit) {
    print("OUTER CONVERGED")
  }
}

print(c((par_cov),true_val))