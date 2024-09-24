########### Starting parameters
#set.seed(1000)
options(scipen=999)
library(VineCopula)
library(gamlss)
source("common_functions.R")

# Choose distributions

#copula_dist="C";margin_dist=EXP();   #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="C";margin_dist=GA();   #    min_par=c(1,1,2);       max_par=c(10,10,10)
copula_dist="C"; margin_dist=ST1();
#margin_dist=ZISICHEL(); min_par=c(1,1,-2,.1,2); max_par=c(10,10,2,.5,10) 

mu_range    =c(0)
sigma_range =c(10)
nu_range    =c(5)
tau_range   =c(10)
theta_range =c(10)
zeta_range  =c(NA)

runs=100; n=100; d=10
start_step_size=.25;step_adjustment=.5;max_steps=5; crit_lik_change=.05*start_step_size

output=list(); counter=1

for (mu in mu_range) {
  for (sigma in sigma_range) {
    for (nu in nu_range) {
      for (tau in tau_range) {
        for (theta in theta_range) {
          for (zeta in zeta_range) {
            for (run_counter in 1:runs) {
              
              tryCatch({
                print(c(mu,sigma,nu,tau,theta,zeta))
              parameters=c(mu,sigma,nu,tau,theta,zeta); names(parameters)=c("mu","sigma","nu","tau","theta","zeta")
              parameters=parameters[c(names(margin_dist$parameters), get_copula_dist(copula_dist)$parameters)]
              
              print(c(counter,run_counter,parameters))  
            
              dataset=loadDataset(n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
                                  , par.margin=c(mu,sigma,nu,tau), par.copula=rep(theta,d-1))
              
              start_par=get_starting_values(copula_dist,margin_dist,dataset)
              
              fit_withdl=fit_jointreg_nocov(input_par=start_par,margin_dist,copula_dist,data=dataset,use_dlcopdpar = TRUE,plot_results = TRUE,
                                            start_step_size=start_step_size,step_adjustment=step_adjustment,max_steps=max_steps, crit_lik_change=crit_lik_change)
              fit_nodl=fit_jointreg_nocov(input_par=start_par,margin_dist,copula_dist,data=dataset,use_dlcopdpar = FALSE,plot_results = FALSE,
                                          start_step_size=start_step_size,step_adjustment=step_adjustment,max_steps=max_steps,crit_lik_change=crit_lik_change)
              
              return_list=list(fit_nodl,fit_withdl,parameters)
              names(return_list)=c("fit_nodl","fit_withdl","parameters")
              
              output[[counter]]=return_list
              counter=counter+1
              }, error=function(e) {
                print("Error in optimisation")
                output[[counter]]=list("error")
                counter=counter+1
              })
              
            }
          }
        }
      }
    }
  }
}

#### Plotting ####
out_w=out_l=matrix(ncol=3+length(start_par),nrow=length(output))
out_true=matrix(ncol=length(start_par),nrow=length(output))
for (i in (1:length(output))) {
  max_runs_nodl=max(output[[i]]$fit_nodl$log_lik_history[,"run_counter"])
  max_runs_withdl=max(output[[i]]$fit_withdl$log_lik_history[,"run_counter"])
  ll_n=output[[i]]$fit_nodl$log_lik_history[max_runs_nodl,]
  ll_w=output[[i]]$fit_withdl$log_lik_history[max_runs_withdl,]
  par_n=output[[i]]$fit_nodl$par_history[max_runs_nodl,]
  par_w=output[[i]]$fit_withdl$par_history[max_runs_withdl,]
  true=output[[i]]$parameters
  out_l[i,]=c(par_n[3:length(par_n)],ll_n[3:5])
  out_w[i,]=c(par_w[3:length(par_w)],ll_w[3:5])
  out_true[i,]=true
}
colnames(out_l)=colnames(out_w)=c(names(par_n)[3:length(par_n)],names(ll_n)[3:length(ll_n)]); colnames(out_true)=names(true)

library(ggplot2); library(ggpubr)
plot_list=list(); i=1
#plot_par="mu"; against_par="theta"; val_range=mu_range
#plot_par="sigma"; against_par="theta"; val_range=sigma_range
#plot_par="theta"; against_par="sigma"; val_range=theta_range
#plot_par="nu"; against_par="sigma"; val_range=nu_range
plot_par="tau"; against_par="sigma"; val_range=tau_range

for (par_val in val_range) {
  
  # create a data frame
  variety=as.factor(rep(out_true[out_true[,plot_par]==par_val,against_par],2))
  note=c(out_l[out_true[,plot_par]==par_val,plot_par],out_w[out_true[,plot_par]==par_val,plot_par])
  treatment=rep(c("separate","joint"),each=length(note)/2)
  data=data.frame(variety, treatment ,  note)
  
  # grouped boxplot
  plot_list[[i]] = ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
    geom_boxplot() +
    xlab(against_par) + ylab(plot_par) + ggtitle(paste(paste(plot_par,"=",sep=""),par_val)) +
    geom_hline(yintercept=par_val, linetype="dashed", color = "red")

  i=i+1
}
ggarrange(plotlist=plot_list,common.legend=TRUE,ncol=3,nrow=3)

plot_par="joint"; against_par="mu"; val_range=theta_range
# create a data frame
variety=rep(paste(out_true[,1],out_true[,2],out_true[,3]),2)
note=c(out_l[,plot_par],out_w[,plot_par])
treatment=rep(c("separate","joint"),each=length(note)/2)
data=data.frame(variety, treatment ,  note)
# grouped boxplot
ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
  geom_boxplot() +
  xlab(against_par) + ylab(plot_par) + ggtitle(paste(paste(plot_par,"=",sep=""),par_val)) +
  facet_wrap(~variety, scale="free")


####################################################

### OUTER ITERATIONS ####

######## Select distributions and range for parameters


num_par=length(margin_dist$parameters)+if(copula_dist=="t") {2} else {1}
 #copula_link=copula.link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv)

final_lik_separate=final_lik_joint=matrix(ncol=5,nrow=0)
final_par_separate=final_par_joint=matrix(ncol=2+num_par,nrow=0)
final_true_val=matrix(ncol=num_par,nrow=0)

n=1000; d=3; crit_lik_change=0.05; start_step_size=.5; plot_results=TRUE; step_adjustment=.9; max_steps=5
for (random_simulation in 1:100) {
  
  print(random_simulation)

  sim_parameters=runif(num_par,min_par,max_par)
  
  a=sim_parameters[length(sim_parameters)];b=sim_parameters[(1:(length(sim_parameters)-1))];d=3
  true_val=c(b,a)
  
  start_par=true_val*runif(length(true_val),0.5,2)
  names(start_par)=c(names(margin_dist$parameters),"theta"); input_par=start_par
  

  log_lik_history=matrix(ncol=3+2,nrow=0)
  par_history=matrix(ncol=length(input_par)+2,nrow=0)
  
  ### Run fit for separate and joint optimisation
  for (use_dlcopdpar in c(FALSE,TRUE)) {
    copula_deriv=if(use_dlcopdpar==TRUE){1}else{0}
    ### CORE ITERATION
    input_par=start_par
    change=1;log_lik_start=0;log_lik_change=1000;run_counter=1;step_size=start_step_size;
    while (abs(log_lik_change)>crit_lik_change) {
      step_size=step_size*(step_adjustment^min(max_steps,run_counter))
      par_history=rbind(par_history,c(copula_deriv,run_counter,input_par))
      
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
      
      
      log_lik_history=rbind(log_lik_history,c(copula_deriv,run_counter,outer_optim_output$log_lik))
      
    }
    par_history=rbind(par_history,c(copula_deriv,run_counter,input_par))
    colnames(log_lik_history)[1:2]=colnames(par_history)[1:2]=c("use_dlcopdpar","run_counter")
  }

  #Plot likelihood and parameters
  if(plot_results==TRUE) {
    #plot.new()
    #par(mfrow=c(1,3))
    #for (i in colnames(log_lik_history)[3:5]) {
    #  plot( log_lik_history[log_lik_history[,"use_dlcopdpar"]==1,i],xlab="LogLik",ylab="Iteration",main=i,type = "l",col="blue",xlim=c(1,max(log_lik_history[,"run_counter"])),ylim=range(log_lik_history[,i]))
    #  lines(log_lik_history[log_lik_history[,"use_dlcopdpar"]==0,i],xlab="LogLik",ylab="Iteration",main=i,type = "l",col="red",xlim=c(1,max(log_lik_history[,"run_counter"])),ylim=range(log_lik_history[,i]))
    #  legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))
    #  
    #}
    
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
outliers=which(abs(final_lik_joint[,5]-final_lik_separate[,5])>100)
final_lik_joint=final_lik_joint[-outliers,]
final_lik_separate=final_lik_separate[-outliers,]

plot.new()
par(mfrow=c(1,2))
lm_coef=lm(final_lik_joint[,5]~final_lik_separate[,5])$coefficients
plot(final_lik_joint[,5]~final_lik_separate[,5],xlab="Separate",ylab="Joint",main=paste(c("Joint v Separate Likelihood: ","Intercept:","Slope:"), c("",lm_coef)),xlim=c(-12000,0))
abline(a=0,b=1,col="red")
hist(final_lik_joint[,5]-final_lik_separate[,5],main="Likelihood Improvement",xlab="Likelihood Change",breaks=10)


mean(final_lik_joint[,"run_counter"])
mean(final_lik_separate[,"run_counter"])
par(mfrow=c(1,2))
hist(final_lik_joint[,"run_counter"])
hist(final_lik_separate[,"run_counter"])
#### INNER ITERATION ####

margin_dist=EXP(); copula_dist="C"; true_val=c(1,1,1,3); 
margin_dist=GA(); copula_dist="C"; true_val=c(1,1,1,1,3); 
margin_dist=ZISICHEL(); copula_dist="C"; true_val=c(1,1,1,1,-1,-2.20,3); 

copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE #copula_link=copula.link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv)
if(two_par_cop) {
  names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}
n=500; d=3
dataset=loadDataset(simOption=3,n=n,d=d, copula.family=BiCopName(copula_dist), copula.link=copula_link, qFUN=eval(parse(text=paste("q",margin_dist$family[1],sep="")))
                    , par.copula=rep(3,d-1), par.margin=c(1,1,1,1),plot_dist = FALSE)
plotDist(dataset,margin_dist$family[1])

#Fit one parameter at a time like RS maybe?
copula_link=get_copula_dist(copula_dist)$copula_link
mm=create_model_matrices(
  mu.formula = ("response ~ 1"),
  sigma.formula = ("~ 1"),
  nu.formula = ("~ 1"),
  tau.formula = ("~ 1"),
  theta.formula=("~1"),
  zeta.formula=("~1"),
  margin.family=margin_dist,copula.family=copula_dist,copula.link=copula_link
)

##Maybe do one backfit iteration first?

###CALCULATE eta,dldpar,d2ldpar,dpardeta

#First run set all covariates except intercept to zero and setup starting vector
input_par=get_starting_values(copula_dist,margin_dist,dataset); names(input_par)=c(names(margin_dist$parameters),"theta")
#input_par=log(input_par)

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
inner_stop_crit=outer_stop_crit=.1
start_step_size=1; step_adjustment=.9; max_steps=10
log_lik_history=matrix(ncol=3,nrow=0)
par_history=matrix(ncol=length(par_cov),nrow=0); colnames(par_history)=names(par_cov)
outer_run_counter=1
while (first_outer_run==TRUE | (abs(outer_log_lik_change)>outer_stop_crit) & (outer_log_lik_change>0)) {
  print("STARTING NEW OUTER ITERATION")
  if(first_outer_run==FALSE) {
    print(c(outer_start_log_lik,outer_end_log_lik,outer_log_lik_change))
  }
  first_outer_run=TRUE
  
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
      
      if(first_outer_run==TRUE) {
        outer_start_log_lik=log_lik["joint"]; first_outer_run=FALSE
      }
      
      
      log_lik_history=rbind(log_lik_history,calc_lik_out$log_lik)
      par_history=rbind(par_history,par_cov)
      
      Fx_1_2[Fx_1_2>1]=1;Fx_1_2[Fx_1_2<0]=0
      
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
        for (inner_par_name in margin_par) {
          
          #Take parameters from input for clarity
          dc_tplus_du_t=input[,"dcdu1"]
          dc_tplus_du_tplus=input[,"dcdu2"]
          l_t=input[,paste(paste("dld",inner_par_name,sep=""),".x",sep="")]
          l_t_plus=input[,paste(paste("dld",inner_par_name,sep=""),".y",sep="")]
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
          i=i+1
        }
        colnames(dlcopdpar)=paste("dlcopd",margin_par,sep="")
        
        par_dlcopdpar=dlcopdpar[,paste("dlcopd",margin_par,sep="")]
        merged_dlcopdpar=merge(cbind(order_copula,par_dlcopdpar),cbind(order_copula,par_dlcopdpar),by.x=c("time2","subject2"),by.y=c("time1","subject1"),all=TRUE)
        merged_dlcopdpar[is.na(merged_dlcopdpar)]=0
        
        
        x_comp=grepl("dlcopd",colnames(merged_dlcopdpar))&grepl(".x",colnames(merged_dlcopdpar))
        y_comp=grepl("dlcopd",colnames(merged_dlcopdpar))&grepl(".y",colnames(merged_dlcopdpar))
        
        d1_cop=0.5*(merged_dlcopdpar[,x_comp]+merged_dlcopdpar[,y_comp])
        
        margin_deriv_subnames=c("m","d","v","t")
        names(margin_deriv_subnames)=c("mu","sigma","nu","tau")
        d1=as.matrix(margin_deriv[grepl(paste("dld",margin_deriv_subnames[par_name],sep=""),names(margin_deriv))][[1]])
        colnames(d1)=paste("dld",par_name,sep="")
        d2=as.matrix(margin_deriv[grepl(paste("d2ld",margin_deriv_subnames[par_name],sep=""),names(margin_deriv))][[1]])
        colnames(d2)=paste("d2ld",par_name,sep="")
        d1_plus_cop=d1 + if(include_dlcopdpar==TRUE){d1_cop} else {0*d1_cop}
        d1=d1_plus_cop
        
        d1=d1[,grepl(par_name,colnames(d1))]
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
      
      #plot_count=3+length(par_cov)
      #sides=round(sqrt(plot_count))
      
      #par(mfrow=c(sides+1,sides))
      #plot(log_lik_history[,3],type="l",main="LogLik - Overall")
      #plot(log_lik_history[,1],type="l",main="LogLik - Margin")
      #plot(log_lik_history[,2],type="l",main="LogLik - Copula")
      
      #for(i in 1:length(colnames(par_history))) {
      #  plot(par_history[,i],type="l",ylim=range(c(par_history[,i],true_val[i])),main=colnames(par_history)[i],xlab="Iteration",ylab="Parameter estimate")
      #  abline(h=true_val[i],col="red")
      #}
      
      print("End Log Lik")
      print(calc_lik_out_end$log_lik)
      
      change_log_lik=calc_lik_out_end$log_lik["joint"]-calc_lik_out$log_lik["joint"]
      
      if(change_log_lik<0) {change_log_lik=0}
      
      step_size = (step_adjustment^min(run_counter,max_steps))*start_step_size
      run_counter=run_counter+1
      
      outer_run_counter=outer_run_counter+1
    }
  }
  outer_end_log_lik=calc_lik_out_end$log_lik["joint"]
  outer_log_lik_change=outer_end_log_lik-outer_start_log_lik
  
  if(abs(outer_log_lik_change)<=outer_stop_crit) {
    print("OUTER CONVERGED")
  }
}


