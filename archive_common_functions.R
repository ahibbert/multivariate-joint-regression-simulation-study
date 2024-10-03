generate_cop_model_matrix = function(dataset,theta.formula=formula(~1),zeta.formula=formula(~1),time) {
  #time="time"
  margins=unique(dataset[,time])
  num_copulas=length(margins)-1
  first=TRUE
  for (margin_number in 1:num_copulas) {
    time_1=dataset[dataset[,time]==margins[margin_number],]
    time_2=dataset[dataset[,time]==margins[margin_number+1],]
    names(time_2)=paste(names(time_2),2,sep=".")
    add_to_cop_data_matrix=cbind(time_1,time_2)
    if (first==TRUE) {
      cop_data_matrix=add_to_cop_data_matrix; first=FALSE
    } else {
      cop_data_matrix=rbind(cop_data_matrix,add_to_cop_data_matrix)
    } 
  }
  #matrix_of_ones=matrix(1,nrow=nrow(cop_data_matrix),ncol=1)
  #colnames(matrix_of_ones)="(Intercept)"
  cop_model_matrix=list()
  #cop_model_matrix[["theta"]]=cbind(matrix_of_ones,model.frame(formula,data=cop_data_matrix))
  #cop_model_matrix[["zeta"]]=cbind(matrix_of_ones,model.frame(zeta.formula,data=cop_data_matrix))
  cop_model_matrix[["theta"]]=model.matrix(theta.formula,data=cop_data_matrix,contrasts.arg = sapply(cop_data_matrix,is.factor))
  cop_model_matrix[["zeta"]]=model.matrix(zeta.formula,data=cop_data_matrix)
  
  return_list=(cop_model_matrix)
  #names(return_list)=("cop_model_matrix")
  return(return_list)
}

calc_joint_likelihood <- function(input_par,mm_mar,margin_dist,copula_dist,return_option="list",copula_link,dataset,mm_cop,verbose=TRUE,calc_d2=FALSE,dFUN=dZISICHEL,pFUN=pZISICHEL)  {
  
  #return_option="list";input_par=setup$start_par;copula_dist="C";margin_dist = setup$margin.family;
  #verbose=TRUE;dFUN=setup$dFUN;pFUN=setup$pFUN;calc_d2=TRUE;copula_link=setup$copula.link; mm_mar=setup$mm_mar;mm_cop=setup$mm_cop
  
  #tryCatch(
  #  {
  if(return_option=="list" & verbose==TRUE) {
    print("Starting parameters:")
    print(as.matrix(input_par))  
  }
  
  #### MARGINAL MODEL: Get linear predictors and their inverse -> **Creates eta, eta_link_inv** from **mm_mar** ####
  margin_names=unique(dataset$time)
  num_margins = length(margin_names)
  response=dataset$response
  margin_parameters=names(mm_mar)
  
  eta = eta_linkinv = matrix(nrow = nrow(dataset), ncol = length(margin_parameters))
  colnames(eta)=colnames(eta_linkinv)= margin_parameters
  
  for (parameter in margin_parameters) {
    par_temp <- input_par[grepl(paste(parameter,sep="."), names(input_par))]
    mm=mm_mar[[parameter]]
    eta[, parameter] = rowSums(mm * matrix(rep(par_temp,each=nrow(mm)),ncol=length(par_temp),dimnames=list(NULL,c(names(par_temp))))) #These are the linear predictors
    FUN=margin_dist[[names(margin_dist)[grepl(paste(parameter,"linkinv",sep="."),names(margin_dist))]]]
    eta_linkinv[, parameter]=FUN(eta[,parameter])
  }
  
  #### MARGINAL MODEL: Get density and probability functions for margins -> ** Creates margin_p and margin_d ####
  
  input_list=list(q=response,x=response) 
  
  for (parameter in margin_parameters) {
    input_list[[parameter]]=eta_linkinv[,parameter] 
  }
  
  #### Get density of the copula function for each pair of margins ####
  
  #margin_names; num_margins; margin_parameters
  
  #Calculate new eta based on input parameters and model matrix
  copula_parameters=names(mm_cop)
  if(length(mm_cop)>1){two_par_cop=TRUE} else {two_par_cop=FALSE}
  
  eta_cop=eta_cop_link_inv=matrix(nrow=nrow(mm_cop$theta),ncol=length(copula_parameters))
  colnames(eta_cop)=colnames(eta_cop_link_inv)=copula_parameters
  for (parameter in copula_parameters) {
    #parameter="theta"
    par_temp = input_par[grepl(parameter, names(input_par))]
    mm=mm_cop[[parameter]]
    eta_cop[,parameter] = rowSums(mm * matrix(rep(par_temp,each=nrow(mm)),ncol=length(par_temp),dimnames=list(NULL,c(names(par_temp))))) #Calc linear predictor
    FUN=copula_link[[names(copula_link)[grepl(paste(parameter,"linkinv",sep="."),names(copula_link))]]]
    eta_cop_link_inv[,parameter]=FUN(eta_cop[,parameter])
  }
  
  ####NUMERICAL DERIVATIVES FOR DLCOPDPAR IN linear predictor form - TO DELETE
  
  fx_eps=list(); count=1; eps_ch=.000001;
  zero_mat_n=matrix(0,nrow=nrow(dataset)/num_margins)
  for (eps in -eps_ch + 0:1*eps_ch*2) {
    
    args=names(input_list)[names(input_list)%in%formalArgs(pFUN)]
    input_list$mu = input_list$mu+eps
    margin_p = do.call(pFUN,args=input_list[args])
    
    margin_p_cop_input_u1=margin_p_cop_input_u2=matrix(ncol=2,nrow=0)
    order_copula=order_copula_u2=matrix(ncol=4,nrow=0)
    for (i in 1:(num_margins)) {
      if(i!=1) {
        #if(i==num_margins){margin_p_cop_input_u1=rbind(margin_p_cop_input_u1,cbind(zero_mat_n,zero_mat_n))}
        margin_p_cop_input_u2=rbind(margin_p_cop_input_u2,cbind(margin_p[dataset$time == margin_names[i-1]],margin_p[dataset$time == margin_names[i]]))
        margin_p_cop_input_u2[margin_p_cop_input_u2>1]=1
        margin_p_cop_input_u2[margin_p_cop_input_u2<0]=0
        order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i-1],c("time","subject")],dataset[dataset$time == margin_names[i],c("time","subject")]))
      }
      if(i!=num_margins) {
        #if(i==1) { margin_p_cop_input_u2=rbind(margin_p_cop_input_u2,cbind(zero_mat_n,zero_mat_n))}
        margin_p_cop_input_u1=rbind(margin_p_cop_input_u1,cbind(margin_p[dataset$time == margin_names[i]],margin_p[dataset$time == margin_names[i+1]]))
        margin_p_cop_input_u1[margin_p_cop_input_u1>1]=1
        margin_p_cop_input_u1[margin_p_cop_input_u1<0]=0
        order_copula_u2=rbind(order_copula_u2,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
      }
      
    }
    
    names(order_copula)=c("time1","subject1","time2","subject2")
    fx_eps[[count]]= 
      c(log(BiCopPDF(  margin_p_cop_input_u1[,1],margin_p_cop_input_u1[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0})),zero_mat_n) + c(zero_mat_n,log(BiCopPDF(  margin_p_cop_input_u2[,1],margin_p_cop_input_u2[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0})))
    count=count+1
  }
  num_dlcopdpar=(fx_eps[[2]]-fx_eps[[1]])/(2*eps_ch)
  
  fx_eps=list(); count=1; eps_ch=.000001;
  zero_mat_n=matrix(0,nrow=nrow(dataset)/num_margins)
  for (eps in -eps_ch + 0:1*eps_ch*2) {
    
    args=names(input_list)[names(input_list)%in%formalArgs(pFUN)]
    input_list$mu = input_list$mu+eps
    margin_p = do.call(pFUN,args=input_list[args])
    
    margin_p_cop_input_u1=margin_p_cop_input_u2=matrix(ncol=2,nrow=0)
    order_copula=order_copula_u2=matrix(ncol=4,nrow=0)
    for (i in 1:(num_margins)) {
      if(i!=1) {
        #if(i==num_margins){margin_p_cop_input_u1=rbind(margin_p_cop_input_u1,cbind(zero_mat_n,zero_mat_n))}
        margin_p_cop_input_u2=rbind(margin_p_cop_input_u2,cbind(margin_p[dataset$time == margin_names[i-1]],margin_p[dataset$time == margin_names[i]]))
        
        margin_p_cop_input_u2[margin_p_cop_input_u2>1]=1
        margin_p_cop_input_u2[margin_p_cop_input_u2<0]=0
        
        order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i-1],c("time","subject")],dataset[dataset$time == margin_names[i],c("time","subject")]))
      }
      if(i!=num_margins) {
        #if(i==1) { margin_p_cop_input_u2=rbind(margin_p_cop_input_u2,cbind(zero_mat_n,zero_mat_n))}
        margin_p_cop_input_u1=rbind(margin_p_cop_input_u1,cbind(margin_p[dataset$time == margin_names[i]],margin_p[dataset$time == margin_names[i+1]]))
        margin_p_cop_input_u1[margin_p_cop_input_u1>1]=1
        margin_p_cop_input_u1[margin_p_cop_input_u1<0]=0
        order_copula_u2=rbind(order_copula_u2,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
      }
      
    }
    
    names(order_copula)=c("time1","subject1","time2","subject2")
    fx_eps[[count]]=
      (c((BiCopPDF(  margin_p_cop_input_u1[,1],margin_p_cop_input_u1[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0})),zero_mat_n) +
         c(zero_mat_n,(BiCopPDF(  margin_p_cop_input_u2[,1],margin_p_cop_input_u2[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0}))))
    count=count+1
  }
  num_dlcopdpar_nolog=(fx_eps[[2]]-fx_eps[[1]])/(2*eps_ch)
  
  # ###OK for dcdu1 and dcdu2
  # fx_eps=list(); count=1; eps_ch=.00001
  # for (eps in -eps_ch + 0:1*eps_ch*2) {
  #   
  #   args=names(input_list)[names(input_list)%in%formalArgs(pFUN)]
  #   #input_list$mu = input_list$mu+eps
  #   margin_p = do.call(pFUN,args=input_list[args])
  #   
  #   margin_p_cop_input=matrix(ncol=2,nrow=0)
  #   order_copula=matrix(ncol=4,nrow=0)
  #   for (i in 1:(num_margins-1)) {
  #     margin_p_cop_input=rbind(margin_p_cop_input,cbind(margin_p[dataset$time == margin_names[i]],margin_p[dataset$time == margin_names[i+1]]))
  #     order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
  #   }
  #   names(order_copula)=c("time1","subject1","time2","subject2")
  #   fx_eps[[count]]=(BiCopPDF(  margin_p_cop_input[,1]+eps,margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0}))
  #   count=count+1
  # }
  # num_dcdu1=(fx_eps[[2]]-fx_eps[[1]])/(2*eps_ch)
  # 
  # ###OK for dcdu1 and dcdu2
  # fx_eps=list(); count=1; eps_ch=.00001
  # for (eps in -eps_ch + 0:1*eps_ch*2) {
  #   
  #   args=names(input_list)[names(input_list)%in%formalArgs(pFUN)]
  #   #input_list$mu = input_list$mu+eps
  #   margin_p = do.call(pFUN,args=input_list[args])
  #   
  #   margin_p_cop_input=matrix(ncol=2,nrow=0)
  #   order_copula=matrix(ncol=4,nrow=0)
  #   for (i in 1:(num_margins-1)) {
  #     margin_p_cop_input=rbind(margin_p_cop_input,cbind(margin_p[dataset$time == margin_names[i]],margin_p[dataset$time == margin_names[i+1]]))
  #     order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
  #   }
  #   names(order_copula)=c("time1","subject1","time2","subject2")
  #   fx_eps[[count]]=(BiCopPDF(  margin_p_cop_input[,1],margin_p_cop_input[,2]+eps,family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0}))
  #   count=count+1
  # }
  # num_dcdu2=(fx_eps[[2]]-fx_eps[[1]])/(2*eps_ch)
  
  #########################
  
  args=names(input_list)[names(input_list)%in%formalArgs(pFUN)]
  margin_p = do.call(pFUN,args=input_list[args])
  args=names(input_list)[names(input_list)%in%formalArgs(dFUN)]
  margin_d = do.call(dFUN,args=input_list[args])
  
  ###order of mm_cop and other 
  copula_d=dcdu1=dcdu2=dldth=d2ldth=dldz=d2ldz=dthdeta=dzdeta=matrix(0,nrow=nrow(eta_cop),ncol=1)
  
  margin_p_cop_input=matrix(ncol=2,nrow=0)
  order_copula=matrix(ncol=4,nrow=0)
  for (i in 1:(num_margins-1)) {
    margin_p_cop_input=rbind(margin_p_cop_input,cbind(margin_p[dataset$time == margin_names[i]],margin_p[dataset$time == margin_names[i+1]]))
    order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
  }
  names(order_copula)=c("time1","subject1","time2","subject2")
  
  copula_d=BiCopPDF(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0})
  
  #Check that optimal value is actually true value, yes
  #sum(log(BiCopPDF(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"]*0+exp(1.1),par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0})))
  
  dldth=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{0},deriv="par",log=TRUE)
  
  ###############
  ###############
  
  if(calc_d2==TRUE) {
    dcdth=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{0},deriv="par",log=FALSE)
    d2cdth=BiCopDeriv2( margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{0},deriv="par")
    d2ldth=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
  }
  
  if(calc_d2==TRUE & two_par_cop) {
    dcdz=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{0},deriv="par2",log=FALSE)
    d2cdz=BiCopDeriv2( margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{0},deriv="par2")
    d2ldz=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)
  }
  
  dcdu1=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0},deriv="u1",log=FALSE)
  dcdu2=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0},deriv="u2",log=FALSE)
  dthdeta=copula_link$dthdeta(eta_cop[,"theta"])
  if (two_par_cop) {dzdeta=copula_link$dzdeta(eta_cop[,"zeta"])}
  
  ##################NUMDERIV CHECK  - DELETE AFTER
  ##################
  
  eps=.0001
  par_change=eta_cop_link_inv[,"theta"]+eps
  eps_plus=log(BiCopPDF(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=par_change,par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0}))
  par_change=eta_cop_link_inv[,"theta"]-eps
  eps_minus=log(BiCopPDF(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=par_change,par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0}))
  num_dldth=(eps_plus-eps_minus)/(2*eps)
  
  
  ###NOT UPDATED YET
  eps=.0001
  par_change=eta_cop_link_inv[,"theta"]+eps
  eps_plus=(BiCopDeriv(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=par_change,par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0},log=TRUE))
  par_change=eta_cop_link_inv[,"theta"]-eps
  eps_minus=(BiCopDeriv(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.numeric(BiCopName(copula_dist)),par=par_change,par2=if(two_par_cop){eta_cop_link_inv[,"zeta"]}else{eta_cop_link_inv[,"theta"]*0},log=TRUE))
  num_d2ldth=(eps_plus-eps_minus)/(2*eps)
  
  #plot.new()
  #par(mfrow=c(2,2))
  #plot(dldth,num_dldth)
  #plot(d2ldth,num_d2ldth)
  
  #### Return Log Likelihood ####
  
  log_lik_results<-c(sum(log(margin_d)),sum(log(copula_d)),sum(log(margin_d))+sum(log(copula_d)))
  names(log_lik_results)=c("Margin","Copula","Total")
  
  #### Extract derivatives of margins ####
  all_derivatives=c(names(margin_dist)[grepl("dl", names(margin_dist))],names(margin_dist)[grepl("dr", names(margin_dist))])
  
  if(calc_d2==TRUE) {all_derivatives=c(all_derivatives,names(margin_dist)[grepl("d2l", names(margin_dist))])}
  derivatives_calculated_all_margins=list()
  
  #Margin derivatives
  
  #for (margin_number in 1:(length_fitted_margins)) {
  i=1 #DON'T DELETE THIS - it's an internal counter
  derivatives_calculated=matrix(nrow=nrow(eta),ncol=length((all_derivatives)))
  for (d_fun in all_derivatives) { 
    if (grepl("dr", d_fun) == TRUE) {
      derivatives_calculated[,i]=margin_dist[[d_fun]](eta[,grepl(sub(".dr","",d_fun), colnames(eta))])
    } else {
      
      input_list=list(y=response)
      for (parameter in margin_parameters) {
        input_list[[parameter]]=eta_linkinv[,parameter]
      }
      args=names(input_list)[names(input_list)%in%formalArgs(margin_dist[[d_fun]])]
      derivatives_calculated[,i]=do.call(margin_dist[[d_fun]],args=input_list[args])
      
    }
    i=i+1
  }
  colnames(derivatives_calculated)=all_derivatives
  derivatives_calculated_all_margins=derivatives_calculated
  #}
  
  if (calc_d2==TRUE & two_par_cop) {
    derivatives_copula = list(dldth,dldz,d2ldth,d2ldz,dcdu1,dcdu2,dthdeta,dzdeta)
    names(derivatives_copula)=c("dldth","dldz","d2ldth","d2ldz","dcdu1","dcdu2","dthdeta","dzdeta")
  } else if (calc_d2==TRUE & !two_par_cop) {
    derivatives_copula = list(dldth,d2ldth,dcdu1,dcdu2,dthdeta)
    names(derivatives_copula)=c("dldth","d2ldth","dcdu1","dcdu2","dthdeta")
  } else if (calc_d2==FALSE & two_par_cop) {
    derivatives_copula = list(dldth,dldz,dcdu1,dcdu2,dthdeta,dzdeta)
    names(derivatives_copula)=c("dldth","dldz","dcdu1","dcdu2","dthdeta","dzdeta")
  } else if (calc_d2==FALSE & !two_par_cop) {
    derivatives_copula = list(dldth,dcdu1,dcdu2,dthdeta)
    names(derivatives_copula)=c("dldth","dcdu1","dcdu2","dthdeta")
  }
  
  return_list=list(log_lik_results,response,margin_d,margin_p,copula_d, order_copula,derivatives_copula,margin_dist,copula_dist,eta,eta_linkinv,derivatives_calculated_all_margins,mm_mar,input_par,eta_cop,eta_cop_link_inv,mm_cop,num_dlcopdpar,num_dlcopdpar_nolog)
  names(return_list)=c("log_lik_results","response","margin_d","margin_p","copula_d","order_copula","derivatives_copula","margin_dist","copula_dist","eta","eta_linkinv","derivatives_calculated_all_margins","mm_mar","input_par","eta_cop","eta_cop_link_inv","mm_cop","num_dlcopdpar","num_dlcopdpar_nolog")
  
  if (return_option == "list") {
    if (verbose==TRUE) {
      print("Log Likelihoods:")
      print(log_lik_results)
    }
    return(return_list)}
  if (return_option == "log_lik") {return(as.numeric(log_lik_results["Total"]))}
  if (return_option == "log_lik_cop") {return(sum(log_lik_results[grepl(",",names(log_lik_results))]))}
  #}, 
  #error = function(err) {
  #  print("Function error")
  #  return(NA)
  #}
  #)
}

extract_parameters_from_fitted <- function(fitted_margins,fitted_copulas=NA,copula_link) {
  mu=list()
  sigma=list()
  nu=list()
  tau=list()
  for (i in 1:length(fitted_margins)) {
    mu[[i]]=coefAll(fitted_margins[[i]])[1]
    sigma[[i]]=coefAll(fitted_margins[[i]])[2]
    nu[[i]]=coefAll(fitted_margins[[i]])[3]
    tau[[i]]=coefAll(fitted_margins[[i]])[4]
  }  
  
  if (!(is.na(fitted_copulas))) {
    theta=zeta=vector()
    for (i in 1:length(fitted_copulas)) {
      theta[i]=copula_link$theta.linkfun(fitted_copulas[[i]]$par)
      names(theta)[i]=paste("theta",names(fitted_copulas)[i])
      cop_par=c(theta)
      if ((fitted_copulas[[i]]$par2)!=0) {
        zeta[i]=copula_link$zeta.linkfun(fitted_copulas[[i]]$par2)
        names(zeta)[i]=paste("zeta",names(fitted_copulas)[i])
        cop_par=c(theta,zeta)
      }
      
    }
    return_list=(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau),cop_par))
  } else {
    return_list=(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau)))
  }
  
  return(return_list)
}

score_function <- function(eta,dldpar,d2ldpar,dpardeta,dlcopdpar,response=NA,phi=1,step_size=1,verbose=TRUE,crit_wk=0.0000001) {
  
  ###TRY USING SECOND DERIVATIVE FOR COPULA
  dldpar=dldpar+dlcopdpar
  
  u_k=dldeta = dldpar * dpardeta
  f_k=d2ldpar
  w_k=-f_k*(dpardeta*dpardeta)
  
  #Stop if weights are too small
  w_k[abs(w_k)<crit_wk]=1
  u_k[abs(w_k)<crit_wk]=0
  
  w_k[abs(u_k)<crit_wk]=1
  u_k[abs(u_k)<crit_wk]=0
  
  z_k=(1-phi)*eta+phi*(eta+step_size*(u_k/w_k))
  
  if(verbose==TRUE) {
    steps_mean=round(rbind(colMeans(as.matrix(eta))
                           ,colMeans(as.matrix(dldpar-dlcopdpar))
                           ,colMeans(as.matrix(dlcopdpar))
                           ,colMeans(as.matrix(dpardeta))
                           ,colMeans(as.matrix(dpardeta*dpardeta))
                           ,colMeans(as.matrix(f_k))
                           ,colMeans(as.matrix(w_k))
                           ,colMeans(as.matrix(u_k))
                           ,colMeans(as.matrix(u_k/w_k))
                           ,colMeans(as.matrix(z_k))
    ),8)
    rownames(steps_mean)=c("eta","dldpar","dlcopdpar","dpardeta","dpardeta2","f_k","w_k","u_k","(1/w_k)*u_k","z_k")
    print(steps_mean)  
  }
  return_list=list(colMeans(as.matrix(z_k)),as.matrix(u_k),as.matrix(f_k),as.matrix(w_k),as.matrix(z_k))
  names(return_list)=c("par","u_k","f_k","w_k","z_k")
  return(return_list)
}

newton_raphson_iteration=function(results,input_par,phi=1,step_size=1,verbose=c(FALSE,FALSE),calc_d2=FALSE,crit_wk=.00001) {
  
  #Start from starting parameters
  #results=calc_joint_likelihood(input_par=setup$start_par,mm_mar=setup$results$mm_mar,margin_dist=setup$margin.family,copula_dist=setup$copula.family,return_option="list",copula_link=setup$copula.link,dataset=dataset,mm_cop=setup$results$mm_cop,verbose=TRUE,calc_d2=TRUE,dFUN=setup$dFUN,pFUN=setup$pFUN)
  #input_par=setup$start_par;verbose=c(TRUE,TRUE); calc_d2=TRUE;phi=1;step_size=1;crit_wk=.00001
  
  #Start from last optim
  #results=optim$results;input_par=optim$input_par;verbose=c(TRUE,TRUE); calc_d2=TRUE;phi=1;step_size=1;crit_wk=.00001
  
  steps_all=list()
  par_out=list()
  margin_dist=results$margin_dist
  copula_dist=results$copula_dist
  mm_cop=results$mm_cop
  mm_mar=results$mm_mar
  num_margins=length(unique(dataset$time))
  
  if(length(names(mm_cop))==1) {
    two_par_cop=FALSE
  } else {
    two_par_cop=TRUE
  }
  
  ##### COPULA ITERATION ####
  
  # Extract COPULA DERIVATIVES from results list (from calc_joint_likelihood)
  
  dlcopdpar=list()
  dldth=results$derivatives_copula$dldth
  
  dthdeta=results$derivatives_copula$dthdeta
  
  dcdu1=results$derivatives_copula$dcdu1
  dcdu2=results$derivatives_copula$dcdu2
  
  eta=results$eta
  eta_linkinv=results$eta_linkinv
  
  if(calc_d2==TRUE) {
    if(two_par_cop) {
      d2ldpar=cbind(results$derivatives_copula$d2ldth,results$derivatives_copula$d2ldz)  
    } else {
      d2ldpar=matrix(results$derivatives_copula$d2ldth)
      colnames(d2ldpar)="d2ldth"
    }
  } else {
    if(two_par_cop) {
      d2ldpar=cbind(-(dldth^2),-(dldz^2))
    } else {
      d2ldpar=-(dldth^2)  
    }
    
  }
  
  if(two_par_cop) {
    dldz=results$derivatives_copula$dldz
    dzdeta=results$derivatives_copula$dzdeta
    dldpar=cbind(dldth,dldz)
    dpardeta=cbind(dthdeta,dzdeta)
  } else {
    #dldz=results$derivatives_copula$dldz
    #dzdeta=results$derivatives_copula$dzdeta
    dldpar=dldth
    dpardeta=dthdeta
  }
  
  eta_input=results$eta_cop[,names(mm_cop)]
  colnames(eta_input)=names(results$eta_cop)
  copula_score=score_function(eta=eta_input,dldpar=dldpar,d2ldpar=d2ldpar,dpardeta=dpardeta,dlcopdpar=dldpar*0,response=dldpar*0,phi=phi,step_size=step_size,verbose=verbose[2],crit_wk = crit_wk) ##Returns updated value
  
  beta_new_cop=list()
  #parameter_names=c("theta","zeta")
  copula_parameters=names(mm_cop)
  
  hess_list=list()
  e_k=copula_score$z_k  
  beta_new_cop=list()
  for (parameter in copula_parameters) {
    mm=mm_cop[[parameter]]
    mm_merged=merge(results$order_copula,mm,by="row.names")[order(results$order_copula$time1,results$order_copula$subject1),]
    mm_merged=mm_merged[order(mm_merged$time1,mm_merged$subject1),]
    mm_ordered=as.matrix(mm_merged[,colnames(mm)])
    colnames(mm_ordered)=colnames(mm)
    mm=mm_ordered
    
    X=as.matrix(mm)
    j=which(copula_parameters==parameter)
    W=diag(copula_score$w_k[,j])
    e_k_par=e_k[,j]
    beta_new_cop[[copula_parameters[j]]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
    rownames(beta_new_cop[[parameter]])= paste(paste(parameter,sep=" "),colnames(mm),sep=".")
    #names(beta_new_cop[[j]])= paste(parameter_names[j],margin_number,rownames(beta_new_cop[[j]]),sep=".")
    if (calc_d2 == TRUE) {
      hess_list[[(copula_parameters[j])]]=-t(X)%*%W%*%X
    }
  }
  
  ##### MARGIN ITERATION ####
  
  margin_parameters=colnames(results$eta)
  
  dldpar=as.matrix(results$derivatives_calculated_all_margins[,grepl("dld",colnames(results$derivatives_calculated_all_margins))]) #First derivative of log likelihood w.r.t. parameters
  
  colnames(dldpar)=colnames(results$derivatives_calculated_all_margins)[grepl("dld",colnames(results$derivatives_calculated_all_margins))]
  
  if (calc_d2==TRUE) {
    d2ldpar=results$derivatives_calculated_all_margins[,grepl("d2ld",colnames(results$derivatives_calculated_all_margins))]
  } else {
    d2ldpar=-(dldpar^2) #Quasi newton second derivative of log function w.r.t parameters
  }
  
  dpardeta=results$derivatives_calculated_all_margins[,grepl(".dr",colnames(results$derivatives_calculated_all_margins))] #Derivative of inverse link function
  
  #Calculating derivatives of the copula function with regard to marginal parameters
  #OK so what we want for dldcop is x2*dldpar(x2)*margin_d(x2)*dcdu2 + x1*dldpar(x1)*margin_d(x1)*dcdu1
  
  order_margin=dataset[,c("time","subject")]
  
  ### Let's change this to a dplyr join to be confident..
  
  order_num_dlcopdpar=matrix(0,ncol=2,nrow=0)
  for(i in 1:num_margins) {
    order_num_dlcopdpar=rbind(order_num_dlcopdpar,cbind(rep(i,(nrow(dataset)/num_margins)),1:(nrow(dataset)/num_margins)))
  }
  
  num_dlcopdpar_ordered=results$num_dlcopdpar[order(order_num_dlcopdpar[,2],order_num_dlcopdpar[,1])]
  num_dlcopdpar_nolog_ordered=results$num_dlcopdpar_nolog[order(order_num_dlcopdpar[,2],order_num_dlcopdpar[,1])]
  margin_components=cbind(order_margin,results$response,results$margin_d,results$margin_p,results$eta_linkinv,dldpar,num_dlcopdpar_ordered,num_dlcopdpar_nolog_ordered)
  
  #So from the margin we need dF(xit)/dMu
  
  margin_components_Ft=margin_components
  margin_components_Ft_minus=margin_components
  margin_components_Ft_minus$time=margin_components_Ft_minus$time+1
  margin_components_Ft_plus=margin_components
  margin_components_Ft_plus$time=margin_components_Ft_plus$time-1
  
  margin_minus=merge(margin_components_Ft,margin_components_Ft_minus,by=c("time","subject"),all.x=TRUE)
  margin_plus_minus=merge(margin_minus,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)
  
  colnames(margin_plus_minus)=gsub(".x",".Ft",gsub(".y",".Ft_minus",colnames(margin_plus_minus)))
  
  #### OK so .x is F_t, .y is F_t-1, no name is 
  #head(margin_plus_minus[margin_plus_minus$subject==1,])
  
  #copula_components=cbind(results$order_copula,dcdu1,dcdu2,results$copula_d)
  copula_components=cbind(results$order_copula,dcdu1,dcdu2,results$copula_d)
  
  margin_copula_merged=merge(margin_plus_minus,copula_components,by.x=c("time","subject"),by.y=c("time1","subject1"),all.x=TRUE)
  margin_copula_merged_2=merge(margin_copula_merged,copula_components,by.x=c("time","subject"),by.y=c("time2","subject2"),all.x=TRUE)
  
  colnames(margin_copula_merged_2)=gsub(".x",".c_plus",gsub(".y",".c_minus",colnames(margin_copula_merged_2)))
  
  margin_deriv_names=colnames(results$derivatives_calculated_all_margins)[grepl("dld",colnames(results$derivatives_calculated_all_margins))]
  
  dlcopdpar=dldpar*0
  
  margin_copula_merged_2[is.na(margin_copula_merged_2)]=0
  
  input=margin_copula_merged_2
  
  for (i in 1:length(margin_deriv_names)) {
    
    
    input[input$subject==1,]
    
    dc_tplus_du_t=input[,"dcdu1.c_plus"]
    dc_tplus_du_tplus=input[,"dcdu2.c_plus"]
    
    dc_tminus_du_tminus=input[,"dcdu1.c_minus"]
    dc_tminus_du_t=input[,"dcdu2.c_minus"]
    
    l_t=input[,"mu.Ft"]
    l_t_minus=input[,"mu.Ft_minus"]
    l_t_plus=input[,"mu"]
    
    x_t=input[,"results$response.Ft"]
    x_t_minus=input[,"results$response.Ft_minus"]
    x_t_plus=input[,"results$response"]
    
    du_t_dmu=-x_t*exp(-x_t*(1/l_t)) * (l_t^-2)
    du_t_plus_dmu=-x_t_plus*exp(-x_t_plus*(1/l_t_plus)) * (l_t_plus^-2)
    du_t_minus_dmu=-x_t_minus*exp(-x_t_minus*(1/l_t_minus)) * (l_t_minus^-2)
    
    #du_t_dmu=input[,paste(margin_deriv_names[i],".Ft",sep="")]
    #du_tminus_dmu=input[,paste(margin_deriv_names[i],".Ft_minus",sep="")]
    #du_tplus_dmu=input[,paste(margin_deriv_names[i],"",sep="")]
    
    c_tplus=input[,"results$copula_d.c_plus"]
    c_tminus=input[,"results$copula_d.c_minus"]
    
    dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
    dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
    dc_minus_dt_minus_dmu=dc_tminus_du_tminus * du_t_minus_dmu
    dc_minus_dt_dmu=dc_tminus_du_t * du_t_dmu
    
    dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
    dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
    dc_minus_dt_minus_dmu[is.nan(dc_minus_dt_minus_dmu)]=0
    dc_minus_dt_dmu[is.nan(dc_minus_dt_dmu)]=0
    
    dcdu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
    dcdu_tminus=((dc_minus_dt_minus_dmu + dc_minus_dt_dmu) / c_tminus)
    
    dcdu_tplus[is.nan(dcdu_tplus)|is.na(dcdu_tplus)]=0
    dcdu_tminus[is.nan(dcdu_tminus)|is.na(dcdu_tminus)]=0
    
    dlcopdpar[,i]=0.5*(dcdu_tplus+dcdu_tminus)
    
    num_deriv=margin_copula_merged_2[,"num_dlcopdpar_ordered.Ft"]
    num_deriv_nolog=margin_copula_merged_2[,"num_dlcopdpar_nolog_ordered.Ft"]
    
    #dcdu1_f1[is.nan(dcdu1_f1)]=0
    #dcdu2_f2[is.nan(dcdu2_f2)]=0
    
    #dlcopdpar[,i]= 0.5*(dcdu1_f1+dcdu2_f2)
  }
  dlcopdpar[is.na(dlcopdpar)]=0
  
  #plot.new()
  #par(mfrow=c(3,3))
  #plot(num_deriv,dlcopdpar,main="dlcopdpar",xlim=range(num_deriv),ylim=range(num_deriv))
  #plot(num_deriv[dcdu_tplus==0],dlcopdpar[dcdu_tplus==0],main="dlcopdpar - when dcdu_tplus is zero",xlim=range(num_deriv),ylim=range(num_deriv))
  #plot(num_deriv[dcdu_tminus==0],dlcopdpar[dcdu_tminus==0],main="dlcopdpar - when dcdu_tminus is zero",xlim=range(num_deriv),ylim=range(num_deriv))
  #plot(num_deriv_nolog[(dcdu_tminus)==0],(dcdu_tplus*c_tplus)[(dcdu_tminus)==0],main="dcopdpar - when dcdu2 is zero",xlim=range(num_deriv_nolog),ylim=range(num_deriv_nolog))
  #plot(num_deriv_nolog[(dcdu_tplus)==0],(dcdu_tminus*c_tminus)[(dcdu_tplus)==0],main="dcopdpar - when dcdu1 is zero",xlim=range(num_deriv_nolog),ylim=range(num_deriv_nolog))
  #plot(results$num_dcdu1,dcdu_tplus[dcdu_tplus!=0],main="dcdu1")
  #plot(results$num_dcdu2,dcdu_tminus[dcdu_tminus!=0],main="dcdu2")
  
  margin_score=score_function(eta=results$eta[,margin_parameters],dldpar=dldpar,d2ldpar=d2ldpar,dpardeta=dpardeta,dlcopdpar=dlcopdpar,response=response_final,phi=phi,step_size=step_size,verbose=verbose[2],crit_wk = crit_wk) ##Returns updated value
  
  e_k=margin_score$z_k#-results$eta_nomargin[[1]][,c("mu","sigma","nu","tau")]
  beta_r=list()
  #for (margin_number in 1:length(results$mm_all)) {
  for (parameter in margin_parameters) {
    mm=results$mm_mar[[parameter]]
    X=as.matrix(mm)
    j=which(margin_parameters==parameter)
    W=diag(margin_score$w_k[,j])
    e_k_par=e_k[,j]
    beta_r[[parameter]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
    rownames(beta_r[[parameter]])= paste(parameter,rownames(beta_r[[parameter]]),sep=".")
    if (calc_d2 == TRUE) {
      hess_list[[parameter]]=-t(X)%*%W%*%X
    }
  }
  
  #### COMBINE MARGIN AND COPULA ####
  
  beta_new=matrix(ncol=1,nrow=0)
  colnames(beta_new)=c("Estimate")
  
  for (parameter in names(beta_r)) {
    beta_new=rbind(beta_new,beta_r[[parameter]])
  }
  for (parameter in names(beta_new_cop)) {
    beta_new=rbind(beta_new,beta_new_cop[[parameter]])
  }
  
  end_par=as.vector(beta_new)
  names(end_par)=rownames(beta_new)
  end_par=end_par[names(input_par)]
  
  hess_list_solved=list()
  
  if (calc_d2==TRUE) {
    for (parameter in names(hess_list)) {
      hess_list_solved[[parameter]]=diag(sqrt(solve(-hess_list[[parameter]])))
    }
    se_par=unlist(hess_list_solved)[names(end_par)]
    return_list=list(input_par,end_par,se_par)
    names(return_list)=c("input_par","end_par","se_par")
  } else {
    return_list=list(input_par,end_par)
    names(return_list)=c("input_par","end_par")
  }
  
  return(return_list)
  
  ########
}

GJRM_L_SETUP<-function(
    mu.formula = ("response ~ 1"),#mu.formula = formula("response ~ as.factor(time)+as.factor(gender)+age")
    sigma.formula = ("~ 1"),#sigma.formula = formula("~ as.factor(time)+age")
    nu.formula = ("~ 1"),#nu.formula = formula("~ as.factor(time)+as.factor(gender)")
    tau.formula = ("~ 1"),#tau.formula = formula("~ age")
    theta.formula=("~1"),#theta.formula=formula("response~as.factor(gender)")
    zeta.formula=("~1"),#zeta.formula=formula("response~1") 
    margin.family=ZISICHEL(mu.link = "log",sigma.link = "log",nu.link = "identity",tau.link = "logit"),
    dFUN=dZISICHEL,pFUN=pZISICHEL,
    copula.family="t",
    copula.link=NA,
    start="fit",
    calc_results=FALSE
) {
  
  if(copula.family %in% c("t")){two_par_cop=TRUE} else {two_par_cop=FALSE}
  
  #Turn text formula inputs into formulas
  mu.formula = formula(mu.formula)#mu.formula = formula("response ~ as.factor(time)+as.factor(gender)+age")
  sigma.formula = formula(sigma.formula)#sigma.formula = formula("~ as.factor(time)+age")
  nu.formula = formula(nu.formula)#nu.formula = formula("~ as.factor(time)+as.factor(gender)")
  tau.formula = formula(tau.formula)#tau.formula = formula("~ age")
  theta.formula=formula(paste("response",theta.formula,sep=""))#theta.formula=formula("response~as.factor(gender)")
  zeta.formula=formula(paste("response",zeta.formula,sep=""))#zeta.formula=formula("response~1") 
  
  #Fit gamlss marginal models
  print("Fitting starting GAMLSS margin")
  gamlss_model <- gamlss(formula = mu.formula
                         , sigma.formula = sigma.formula
                         , nu.formula = nu.formula
                         , tau.formula = tau.formula
                         , family=margin.family,data=dataset)
  
  fitted_copulas=fitted_margins=list()
  fitted_margins[[1]]=gamlss_model
  
  #Extract margin model matrix
  mm_mar=list()
  for (parameter in gamlss_model$parameters) {
    mm_mar[[parameter]]= model.matrix(gamlss_model,what=parameter)
  }
  
  print("Fitting starting COPULA")
  #Fit starting copula with just an intercept for each parameter
  
  if(start=="fit") {
    gamlss_unif_resid=pnorm(gamlss_model$residuals)
    timepoints=unique(dataset$time)
    copula_response=matrix(ncol=2,nrow=0)
    for (i in timepoints[1:(length(timepoints)-1)]) {
      copula_response=rbind(copula_response,cbind(gamlss_unif_resid[dataset$time==i],gamlss_unif_resid[dataset$time==i+1]))
    }
    fitted_copula=BiCopEst(copula_response[,1],copula_response[,2],family=BiCopName(copula.family))
    start_par<-c(extract_parameters_from_fitted(fitted_margins,fitted_copulas=list(fitted_copula),copula_link=copula.link))
  } else {
    cop_par=log(c(BiCopTau2Par(family=BiCopName(copula.family),tau=
                                 cor( dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],"response"],
                                      dataset[dataset$time %in% unique(dataset$time)[2:(length(unique(dataset$time)))],"response"],
                                      method="kendall"))))
    if(two_par_cop) {cop_par=c(cop_par,2);names(cop_par)=c("theta","zeta")}else{names(cop_par)="theta"}
    start_par<-c(extract_parameters_from_fitted(fitted_margins,fitted_copulas=NA,copula_link=copula.link),cop_par)  
    fitted_copula=NA
  }
  
  #Get copula model matrix
  #mm_cop=generate_cop_model_matrix(dataset=dataset,formula=theta.formula,zeta.formula=zeta.formula,time="time")
  mm_cop=list()
  mm_cop[["theta"]]=model.matrix(gamlss(formula=theta.formula,data=(dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],])))
  if(two_par_cop) {
    mm_cop[["zeta"]]=model.matrix(gamlss(formula=zeta.formula,data=(dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],])))
  }
  
  print("Extracting start parameters and model matrices")
  #Create parameter vector from model matrix for copulas
  temp_cop_parameter_names=c()
  for (parameter in names(mm_cop)) {
    temp_cop_parameter_names=c(temp_cop_parameter_names,paste(parameter,colnames(mm_cop[[parameter]]),sep="."))
  }
  
  #Setting up starting parameter vector with correct factors from mm_mar and mm_cop
  theta_par_loc=grepl("theta",names(start_par))
  if(two_par_cop) {
    zeta_par_loc=grepl("zeta",names(start_par))
  }
  
  temp_cop_start_par=vector(length=length(temp_cop_parameter_names))
  names(temp_cop_start_par)=temp_cop_parameter_names
  
  temp_cop_start_par[grepl("Intercept",names(temp_cop_start_par))&grepl("theta",names(temp_cop_start_par))]=start_par[theta_par_loc][1] #Theta starting value
  if(two_par_cop) {
    temp_cop_start_par[grepl("Intercept",names(temp_cop_start_par))&grepl("zeta",names(temp_cop_start_par))]=start_par[zeta_par_loc][1] #Theta starting value
    start_par=c(start_par[!(theta_par_loc | zeta_par_loc)],temp_cop_start_par)
  } else {
    start_par=c(start_par[!(theta_par_loc)],temp_cop_start_par)
  }
  
  if(calc_results==TRUE) {
    results=calc_joint_likelihood(input_par = start_par  #optim_results$par
                                  , mm_mar = mm_mar
                                  , margin_dist = margin.family
                                  , copula_dist=copula.family
                                  , copula_link=copula.link
                                  , mm_cop = mm_cop
                                  , return_option="list"
                                  , dataset=dataset
                                  , verbose=TRUE
                                  , calc_d2 = TRUE
    )
    return_list=list(start_par,mm_mar,mm_cop,margin.family,copula.family,copula.link,dFUN,pFUN,fitted_margins,fitted_copulas,results)
    names(return_list)=c("start_par","mm_mar","mm_cop","margin.family","copula.family","copula.link","dFUN","pFUN","fitted_margins","fitted_copulas","results")  
    
  } else {
    return_list=list(start_par,mm_mar,mm_cop,margin.family,copula.family,copula.link,dFUN,pFUN,fitted_margins,fitted_copulas)
    names(return_list)=c("start_par","mm_mar","mm_cop","margin.family","copula.family","copula.link","dFUN","pFUN","fitted_margins","fitted_copulas")  
  }
  
  
  return(return_list)
}

GJRM_L_OPTIM <- function(
    dataset,
    start_par,
    mm_mar,
    mm_cop,
    margin_dist,
    copula_dist,
    copula_link,
    dFUN=dZISICHEL,
    pFUN=pZISICHEL,
    verbose=TRUE,
    calc_d2=FALSE,
    phi=0.5,
    step_size=1,
    step_adjustment=0.5,
    step_reduc_count=3,
    crit_wk=0.0000001  
    ,true_val=NA
    ,stopifnegative=TRUE
    ,stopval=.1
    ,max_runs=100
    ,plot_optim=FALSE
    ,plot_end_only=TRUE
) {
  
  if(all(is.na(true_val))) {true_val=start_par}
  
  first_input_par=start_par
  end_par_matrix=matrix(0,ncol=length(setup$start_par),nrow=0)
  end_loglik_matrix=matrix(0,ncol=3,nrow=0)
  first_run=TRUE
  change=1
  run_counter=1
  converged_status="Undefined"
  error="No Error"
  #Note doesn't print first set of starting parameters in plot - to fix
  
  stopifnegative=stopifnegative
  while ((abs(change) > stopval) & run_counter <= max_runs) { #
    if(!first_run) {start_log_lik=results$log_lik_results} else {input_par=first_input_par; phi_inner=phi}
    
    try_result=tryCatch({
      results=calc_joint_likelihood(input_par = input_par  #optim_results$par
                                    , mm_mar = setup$mm_mar
                                    , margin_dist = setup$margin.family
                                    , copula_dist=setup$copula.family
                                    , copula_link=setup$copula.link
                                    , mm_cop = setup$mm_cop
                                    , return_option="list"
                                    , dataset=dataset
                                    , verbose=verbose
                                    ,calc_d2 = calc_d2
      )
      "Success"
    }, error=function(e) {
      print(e)
      print("Convergence Error")
      return("Failure")
    }
    )
    
    if(try_result!="Success") {converged_status="Extreme Value";break}
    
    end_log_lik=results$log_lik_results
    end_loglik_matrix=rbind(end_loglik_matrix,end_log_lik)
    
    if(stopifnegative==TRUE) {
      if(nrow(end_loglik_matrix)>=3) {
        if(all(end_loglik_matrix[(length(end_loglik_matrix)-3):(length(end_loglik_matrix)-1)]-end_loglik_matrix[(length(end_loglik_matrix)-2):length(end_loglik_matrix)]>0))
        {converged_status="Negative";break}
      }
      if((change/end_log_lik["Total"])< -0.05) {converged_status="Negative";break}
    }
    
    
    if(first_run) {change=1} else {
      change=end_log_lik["Total"]-start_log_lik["Total"]
      log_lik_output=c(start_log_lik["Total"],end_log_lik["Total"],change,phi_inner)
      names(log_lik_output) = c("Start LogLik","End LogLik","Change","phi")
      print(log_lik_output)
    }
    
    try_result=tryCatch({
      iteration_out=newton_raphson_iteration(results
                                             ,input_par
                                             ,phi=phi_inner
                                             ,step_size=step_size
                                             ,verbose=c(verbose,verbose)
                                             ,crit_wk=crit_wk
                                             ,calc_d2=calc_d2)
      "Success"
    }, error=function(e) {
      print(e)
      print("Convergence Error")
      return("Extreme Value")
    }
    )
    
    if(try_result!="Success") {converged_status="Extreme Value";break}
    
    cbind(iteration_out[[1]],iteration_out[[2]])
    end_par=iteration_out[[2]]
    input_par=end_par
    first_run=FALSE
    end_par_matrix=rbind(end_par_matrix,end_par)
    
    run_counter=run_counter+1
    #Showing charts
    colnames(end_par_matrix)=names(end_par)
    pars=round(sqrt(ncol(end_par_matrix)+3),0)
    
    
    #true_val=c(0.3,0.2,0.1,0.3,0.2,0.1,-0.8,-2.19,2.5,2) #Temporary
    names(true_val)=colnames(end_par_matrix)
    
    if(!plot_end_only) {
      if(plot_optim==TRUE) {
        plot.new();par(mfrow=c(pars,pars+1))
        for (par_name in colnames(end_par_matrix)) {
          true=true_val[par_name]
          actual=c(setup$start_par[par_name],end_par_matrix[,par_name])
          plot(1:(1+nrow(end_par_matrix)),actual,main=par_name,ylim=c(min(true-.1,min(actual)-.1),max(true+.1,max(actual)+.1)))   
          abline(h=true,col="red")
        }
        plot(1:nrow(end_par_matrix),end_loglik_matrix[,1],main="loglik - margin")
        plot(1:nrow(end_par_matrix),end_loglik_matrix[,2],main="loglik - copula")
        plot(1:nrow(end_par_matrix),end_loglik_matrix[,3],main="loglik - total")
        
      }  
    }
    
    if(abs(change) <= stopval) {converged_status="Converged"} else {converged_status="Not Converged"}
    
    phi_inner=phi*(step_adjustment^min(step_reduc_count,run_counter-1))
  }
  
  if(plot_end_only==TRUE & plot_optim==TRUE&converged_status!="Extreme Value") {
    plot.new();par(mfrow=c(pars,pars+1))
    for (par_name in colnames(end_par_matrix)) {
      true=true_val[par_name]
      actual=c(setup$start_par[par_name],end_par_matrix[,par_name])
      plot(1:(1+nrow(end_par_matrix)),actual,main=par_name,ylim=c(min(true-.1,min(actual)-.1),max(true+.1,max(actual)+.1)))   
      abline(h=true,col="red")
    }
    plot(1:nrow(end_par_matrix),end_loglik_matrix[,1],main="loglik - margin")
    plot(1:nrow(end_par_matrix),end_loglik_matrix[,2],main="loglik - copula")
    plot(1:nrow(end_par_matrix),end_loglik_matrix[,3],main="loglik - total")
  }
  
  print("Separate likelihoods")
  print(end_loglik_matrix[1,] )
  print("Joint likelihood")
  print(end_loglik_matrix[which.max(end_loglik_matrix[,"Total"]),] )
  
  return_list=list(end_par_matrix,end_loglik_matrix,results,true_val,converged_status)
  names(return_list)=c("end_par_matrix","end_loglik_matrix","results","true_val","converged_status")
  return(return_list)
}

GJRM_POST_OPTIM=function(optim,setup,dataset) {
  
  #Extract final parameters
  end_par=optim$end_par_matrix[nrow(optim$end_par_matrix),]
  
  if(optim$converged_status=="Extreme Value") {
    return_list=list(end_par,end_par*0,"Extreme Value Error")
  } else {
    
    results=optim$results
    start_par=setup$start_par
    mm_mar=setup$mm_mar
    mm_cop=setup$mm_cop
    margin_dist=setup$margin.family
    copula_dist=setup$copula.family
    copula_link=setup$copula.link
    dFUN=setup$dFUN
    pFUN=setup$pFUN
    dataset=dataset
    
    #Calculate final likelihood
    final_results= calc_joint_likelihood(input_par =end_par  #optim_results$par
                                         ,mm_mar = mm_mar
                                         ,margin_dist = margin_dist
                                         , copula_dist=copula_dist
                                         , copula_link=copula_link
                                         , mm_cop = mm_cop
                                         , return_option="list"
                                         , dataset=dataset
                                         , verbose=FALSE
                                         , calc_d2 = TRUE
    )
    final_iteration_out=newton_raphson_iteration(final_results,end_par,phi=0.5,step_size=0.5,verbose=c(FALSE,FALSE),calc_d2=TRUE)
    new_SE=final_iteration_out$se_par
    return_list=list(end_par,new_SE,final_results)  
  }
  
  names(return_list)=c("end_par","new_SE","final_results")
  return(return_list)
}




