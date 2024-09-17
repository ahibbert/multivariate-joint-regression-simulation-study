###########MVP OPTIMISATION

options(scipen=999)
library(VineCopula)
library(gamlss)
source("common_functions.R")

#Defining link functions for copula manually

copula_link=copula.link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv)
names(copula_link)=names(copula.link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")
copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE
if(two_par_cop) {
  names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}

dataset=loadDataset(simOption=3
                    ,plot_dist=FALSE
                    ,n=1000
                    ,d=3
                    ,copula.family=BiCopName("C")
                    ,copula.link=copula_link
                    ,qFUN=qEXP
                    ,par.copula=c(3,3)
                    ,par.margin=c(1,1,1))

plotDist(dataset,dist="EXP")

###FIND BEST STARTING PARAMETERS - currently in non-link form as it seems to work

#START PARAMETERS
input_par=c(1,3); names(input_par)=c("mu","theta");
margin_dist=EXP(); copula_dist="C"
copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE #copula_link=copula.link=list(logit,logit_inv,dlogit_inv,log_2plus,log_2plus_inv,dlog_2plus_inv)
if(two_par_cop) {
  names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}

#####Find close starting parameters (just outer iterations with no covariates)
#change=1;log_lik_start=0;log_lik_change=1000
#while (abs(log_lik_change)>0.1) {
#  outer_optim_output=optim_outer(par=input_par,dataset,margin_dist,copula_dist,copula_link=copula_link)
#  input_par=outer_optim_output$par_end
#  change=sum(outer_optim_output$par_change)
#  log_lik=outer_optim_output$log_lik["joint"]
#  log_lik_change=log_lik-log_lik_start
#  log_lik_start=log_lik
#}

###INNER ITERATION

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
input_par=c(1,3); names(input_par)=c("mu","theta")
par_eta=outer_optim_output$par_eta_end
par_cov=as.numeric(vector())
for (par_name in names(mm)) {
    par_cov_single=vector(length=length(colnames(mm[[par_name]])))
    names(par_cov_single)=paste(par_name,colnames(mm[[par_name]]),sep=".")
    par_cov_single[1]=par_eta[par_name]
    par_cov_single[2:length(par_cov_single)]=0
    par_cov=c(par_cov,par_cov_single)
}

print("THIS FUNCTION ASSUMES RESPONSE IS ORDERED AS TIME, SUBJECT | PAR INPUT MUST BE NAMED")

include_dlcopdpar=FALSE
first_outer_run=TRUE
outer_log_lik_change=outer_start_log_lik=outer_end_log_lik=0
inner_stop_crit=outer_stop_crit=0.1
while (first_outer_run==TRUE | abs(outer_log_lik_change)>outer_stop_crit) {
  print("STARTING NEW OUTER ITERATION")
  if(first_outer_run==FALSE) {
    print(c(outer_start_log_lik,outer_end_log_lik,outer_log_lik_change))
  }
  first_outer_run=FALSE
  

  for (par_name in names(mm)) {
    print(paste("INNTER ITERATION: Parameter:",par_name))
    first_inner_run=TRUE; change_log_lik=0
    
    step_size=0.25
    step_adjustment=0.5
    run_counter=1
    max_steps=5
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
      
      
      ### Calculate copula derivatives w.r.t margin parameters
      if(!par_name %in% c("mu","sigma","nu","tau")) {
        if(!"zeta" %in% names(eta_inv)) {
          d1=sum(dldth)
          d2=sum(d2ldth2)
        } else {
          d1=colSums(cbind(dldth,dldz))
          d2=colSums(cbind(d2ldth2,d2ldz2))
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
          
          dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
          dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
          dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
          dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
          dcdu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
          dcdu_tplus[is.nan(dcdu_tplus)|is.na(dcdu_tplus)]=0
          
          dlcopdpar[,i]=dcdu_tplus; i=i+1
          
          #num_deriv=margin_copula_merged_2[,"num_dlcopdpar_ordered.Ft"]
          #num_deriv_nolog=margin_copula_merged_2[,"num_dlcopdpar_nolog_ordered.Ft"]
          
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
      
      score=score_function_v2(eta=eta[[par_name]],dldpar=d1,d2ldpar=d2,dpardeta=eta_dr[[par_name]])
      
      X=as.matrix(mm[[par_name]])
      W=diag(as.vector(score$w_k))
      z_k=score$z_k
      beta_start=par_cov[paste(paste(par_name,sep=" "),colnames(mm[[par_name]]),sep=".")]
      beta=beta_start*step_size + (1-step_size)*(as.vector(solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%z_k))
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
      step_size = (step_adjustment^min(run_counter,max_steps))*step_size
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

print(par_cov)