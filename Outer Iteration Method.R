

#### OUTER ITERATION ####
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
start_step_size=1;step_adjustment=.5;max_steps=5; crit_lik_change=.05*start_step_size

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

#### Plotting
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


#### OUTER ITERATIONS #### OLD ####

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

####
####



