generateBivDist <- function(n,a,b,c,mu1,mu2,dist) {
  
  if(dist=="GA") {
    #Simulating bivariate random variable according to functional input
    w<-rbeta(n,a,b)
    margin_1<-w*rgamma(n,shape=a+b,scale=mu1)
    margin_2<-w*rgamma(n,shape=a+b,scale=mu2)
    
  }
  if(dist=="NO") {
    require(MASS)
    normData<-mvrnorm(n,mu=c(mu1,mu2),Sigma = matrix(c(a^2,c*a*b,c*a*b,b^2),nrow=2))
    margin_1<-normData[,1]
    margin_2<-normData[,2]
  }
  
  if(dist=="PO") {
    
    #Compound multiple poisson of Stein & Juritz, 1987
    c<-rgamma(n,c)
  
    margin_1=vector(length = n) 
    margin_2=vector(length = n) 
    for (i in 1:n) {
      margin_1[i]=rpois(1,mu1*c[i])
    }
    for (i in 1:n) {
      margin_2[i]=rpois(1,mu2*c[i])
    }
  }

  #Transforming data to format required for random effect models
  patient<-as.factor(seq(1:n))
  dataset<-as.data.frame(rbind(cbind(patient,margin_1,0)
                               ,cbind(patient,margin_2,1)))
  colnames(dataset)<-c("patient","random_variable","time")
  
  dataset<-dataset[order(dataset$patient),]
  
  return(dataset)
}

fitBivModels <-function(dataset,dist,include="ALL",a,b,c,mu1,mu2) {
  
  n=nrow(dataset[dataset$time==0,])
  
  #Data Setup
  gamma_c_mu1<-dataset[dataset$time==0,]
  gamma_c_mu2<-dataset[dataset$time==1,]
  
  if(dist=="GA"){
    actuals<-c( log(a/mu1)
                , log(a/mu2)
                , NA
                , NA
                , NA
                ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="kendall")[1,2]*100
                ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="pearson")[1,2]*100
    )
  }
  if(dist=="NO"){
    actuals<-c( 
      mu1
      , mu2
      , (a*sqrt(1-c^2))/sqrt(n)
      , (b*sqrt(1-c^2))/sqrt(n)
      , sqrt(a^2+b^2-2*a*b*c)/sqrt(n)
      ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="kendall")[1,2]*100
      ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="pearson")[1,2]*100
    )
  }
  if(dist=="PO"){
    
    e_x1 = mu1*c
    e_x2 = mu2*c
    v_x1 = (((mu1^2)*(c)+(mu1*c))/((mu1*c)^2))
    v_x2 = (((mu2^2)*(c)+(mu2*c))/((mu2*c)^2))
    
    actuals<-c( 
      e_x1
      , e_x2
      , sqrt(v_x1)     /sqrt(n)
      , sqrt(v_x2)     /sqrt(n)
      , sqrt(
          (v_x2 + v_x1)
          - log((mu1*mu2*c)/(e_x1*e_x2))
      ) /sqrt(n)
      ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="kendall")[1,2]*100
      ,cor(cbind(gamma_c_mu1$random_variable,gamma_c_mu2$random_variable),method="pearson")[1,2]*100
    )
  }
  
  if(include=="ALL" || include=="non-GJRM" ) {
  
    require(gamlss)
    require(gee)
    require(lme4)
    require(MASS)
    require(gamlss.mx)
    
    ###Non-GJRM models first as GJRM breaks base gamlss
    
    if(dist=="GA") {
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1), data=dataset, family=Gamma(link = "log"), maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=Gamma(link = "log"), maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=GA()) ))
      #model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=GA(), method=CG(1000))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=GA()
                              , g.control = gamlss.control(trace = FALSE,method=CG(1000)), mixture="gq",K=2)))
      
      invisible(capture.output(model_lme4 <- glmer(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset, family=Gamma(link="log"))))
    }
    if(dist=="NO") {
      invisible(capture.output(model_glm <- glm(random_variable~as.factor(time==1), data=dataset, family=gaussian, maxit=1000)))
      invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=gaussian, maxiter=25, corstr = "exchangeable")))
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NO())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=NO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family= NO()
                              , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- lmer(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset)
    }
    
    if(dist=="PO"||dist=="NB") {
      invisible(capture.output(model_glm <- glm.nb(random_variable~as.factor(time==1), data=dataset, maxit=1000)))
      #invisible(capture.output(model_gee<-gee(random_variable~as.factor(time==1), id=patient, data=dataset, family=negative.binomial, maxiter=25, corstr = "exchangeable")))
      
      invisible(capture.output(model_re_nosig <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), data=dataset, family=NBI())))
      #invisible(capture.output(model_re <- gamlss(formula=random_variable~as.factor(time==1)+random(as.factor(patient)), sigma.formula=~as.factor(time==1), data=dataset, family=PO(), method=CG(1000))))
      invisible(capture.output(model_re_np <- gamlssNP(formula=random_variable~as.factor(time==1), sigma.formula=~as.factor(time==1), random=as.factor(dataset$patient), data=dataset, family=NBI()
                              , g.control = gamlss.control(trace = FALSE), mixture="gq",K=2)))
      
      model_lme4 <- glmer.nb(formula=random_variable~as.factor(time==1) + (1|patient), data=dataset)
    }
    
    ###Capturing coefficient values and errors from each model
    summary_glm<-c( summary(model_glm)$coeff[1]
                    ,summary(model_glm)$coeff[2]
                    ,summary(model_glm)$coeff[3]
                    ,summary(model_glm)$coeff[4]
                    , logLik(model_glm)
                    , AIC(model_glm)
                    , BIC(model_glm)
                    , 3
    )
    if(dist=="PO"||dist=="NB"){summary_gee<-c(NA,NA,NA,NA,NA,NA,NA)} else{
    summary_gee<-c( summary(model_gee)$coeff[1]
                    , summary(model_gee)$coeff[2]
                    , summary(model_gee)$coeff[3] 
                    , summary(model_gee)$coeff[4] 
                    , NA
                    , NA
                    , NA
                    , 4
    )}
    
    invisible(capture.output(
      summary_re_nosig<-c( summary(model_re_nosig)[1]
                           ,summary(model_re_nosig)[2]
                           ,if(dist=="PO"){summary(model_re_nosig)[4]}else{summary(model_re_nosig)[4]}
                           ,if(dist=="PO"){summary(model_re_nosig)[5]}else{summary(model_re_nosig)[5]}
                           ,logLik(model_re_nosig)
                           ,AIC(model_re_nosig)
                           ,BIC(model_re_nosig)
                           , model_re_nosig$df.fit
      )
    ))
    
    #invisible(capture.output(
      #summary_re<-c( summary(model_re)[1]
      #               ,summary(model_re)[2]
      #               ,if(dist=="PO"){summary(model_re)[3]}else{summary(model_re)[5]}
      #               ,if(dist=="PO"){summary(model_re)[4]}else{summary(model_re)[6]}
      #               , AIC(model_re)
      #               , BIC(model_re)
      #               , model_re$df.fit
      #)
    #))
    
    invisible(capture.output(
      summary_re_np<-c( summary(model_re_np)[1]
                     ,summary(model_re_np)[2]
                     ,if(dist=="PO"){summary(model_re_np)[6]}else{summary(model_re_np)[6]}
                     ,if(dist=="PO"){summary(model_re_np)[7]}else{summary(model_re_np)[7]}
                     ,logLik(model_re_np)
                     , AIC(model_re_np)
                     , BIC(model_re_np)
                     , model_re_np$df.fit
      )
    ))
    
    
    summary_lme4 <- c(summary(model_lme4)$coefficients[1]
                      ,summary(model_lme4)$coefficients[2]
                      ,summary(model_lme4)$coefficients[3]
                      ,summary(model_lme4)$coefficients[4]
                      ,logLik(model_lme4)
                      , AIC(model_lme4)
                      , BIC(model_lme4)
                      , 4)
    
    ###Calculating effective degrees of freedom from Donohue
    X<-getME(model_lme4,name="X")
    Z<-getME(model_lme4,name="Z")
    U<-cbind(X,Z)
    W<-model_lme4@resp$sqrtrwt #weights(model_lme4,type = "working")
    UWU=(t(as.matrix(U))%*%(diag(as.vector(W)))%*%as.matrix(U))
    dim(UWU)
    D<-getME(model_lme4,name="Lambda")
    
    if(sum(D)==0) {lme_EDF=summary_lme4[length(summary_lme4)]} else {
      D_inv<-solve(D)
      dinv_plus_00<-c(0,0,diag(D_inv))
      lme_EDF=sum(diag(UWU%*%solve(UWU+diag(dinv_plus_00))))
      
    }
    
    summary_lme4 <- c(summary(model_lme4)$coefficients[1]
                      ,summary(model_lme4)$coefficients[2]
                      ,summary(model_lme4)$coefficients[3]
                      ,summary(model_lme4)$coefficients[4]
                      ,logLik(model_lme4)
                      ,-2*logLik(model_lme4)+2*lme_EDF
                      , BIC(model_lme4)
                      ,lme_EDF)
    
  }
  
  if(include=="ALL" || include=="GJRM" ) {
  
    require(GJRM)
    
    #Setting up GJRM equations
    eq.mu.1 <- formula(random_variable~1)
    eq.mu.2 <- formula(random_variable.1~1)
    fl <- list(eq.mu.1, eq.mu.2)
    
    if(dist=="NO"){margin_dist="N"}
    if(dist=="GA"){margin_dist="GA"}
    if(dist=="PO"){margin_dist="NBI"}
    
    model_copula<-    gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "C0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
    model_copula_n<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "N",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_j<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "J0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_g<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "G0",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_f<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "F",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_amh<-gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "AMH",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_fgm<-gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "FGM",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_pl<- gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "PL",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_h<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "HO",data=data.frame(gamma_c_mu1,gamma_c_mu2),model="B")
    model_copula_t<-  gjrm(fl, margins = c(margin_dist,margin_dist) , copula = "T",data=data.frame(gamma_c_mu1,gamma_c_mu2),model ="B")
    
    summary_cop<-c( model_copula$coefficients[1]
                    , model_copula$coefficients[2]
                    , summary(model_copula)$tableP1[2] #SE for time 0
                    , summary(model_copula)$tableP2[2] #SE for time 1
                    ,logLik(model_copula)
                    , 2*5-2*logLik(model_copula)
                    ,BIC(model_copula)
                    , 5
    )
    
    summary_cop_n<-c( model_copula_n$coefficients[1]
                      , model_copula_n$coefficients[2] 
                      , summary(model_copula_n)$tableP1[2] #SE for time 0
                      , summary(model_copula_n)$tableP2[2] #SE for time 1
                      , logLik(model_copula_n)
                      , 2*5-2*logLik(model_copula_n)
                      ,BIC(model_copula_n)
                      , 5
                      
    )
    summary_cop_j<-c( model_copula_j$coefficients[1]
                      , model_copula_j$coefficients[2]
                      , summary(model_copula_j)$tableP1[2] #SE for time 0
                      , summary(model_copula_j)$tableP2[2] #SE for time 1
                      , logLik(model_copula_j)
                      , 2*5-2*logLik(model_copula_j)
                      ,BIC(model_copula_j)
                      , 5
                      
    )
    summary_cop_g<-c( model_copula_g$coefficients[1]
                      , model_copula_g$coefficients[2] 
                      , summary(model_copula_g)$tableP1[2] #SE for time 0
                      , summary(model_copula_g)$tableP2[2] #SE for time 1
                      , logLik(model_copula_g)
                      , 2*5-2*logLik(model_copula_g)
                      ,BIC(model_copula_g)
                      , 5
                      
    )
    summary_cop_f<-c( model_copula_f$coefficients[1]
                      , model_copula_f$coefficients[2]
                      , summary(model_copula_f)$tableP1[2] #SE for time 0
                      , summary(model_copula_f)$tableP2[2] #SE for time 1
                      , logLik(model_copula_f)
                      , 2*5-2*logLik(model_copula_f)
                      ,BIC(model_copula_f)
                      , 5
                      
    )
    summary_cop_amh<-c( model_copula_amh$coefficients[1]
                        , model_copula_amh$coefficients[2]
                        , summary(model_copula_amh)$tableP1[2] #SE for time 0
                        , summary(model_copula_amh)$tableP2[2] #SE for time 1
                        , logLik(model_copula_amh)
                        , 2*5-2*logLik(model_copula_amh)
                        ,BIC(model_copula_amh)
                        , 5
                        
    )
    summary_cop_fgm<-c( model_copula_fgm$coefficients[1]
                        , model_copula_fgm$coefficients[2]
                        , summary(model_copula_fgm)$tableP1[2] #SE for time 0
                        , summary(model_copula_fgm)$tableP2[2] #SE for time 1
                        , logLik(model_copula_fgm)
                        , 2*5-2*logLik(model_copula_fgm)
                        ,BIC(model_copula_fgm)
                        , 5
                        
    )
    summary_cop_pl<-c( model_copula_pl$coefficients[1]
                       , model_copula_pl$coefficients[2]
                       , summary(model_copula_pl)$tableP1[2] #SE for time 0
                       , summary(model_copula_pl)$tableP2[2] #SE for time 1
                       , logLik(model_copula_pl)
                       , 2*5-2*logLik(model_copula_pl)
                       ,BIC(model_copula_pl)
                       , 5
                       
    )
    summary_cop_h<-c( model_copula_h$coefficients[1]
                      , model_copula_h$coefficients[2]
                      , summary(model_copula_h)$tableP1[2] #SE for time 0
                      , summary(model_copula_h)$tableP2[2] #SE for time 1
                      , logLik(model_copula_h)
                      , 2*5-2*logLik(model_copula_h)
                      ,BIC(model_copula_h)
                      , 5
                      
    )
    summary_cop_t<-c( model_copula_t$coefficients[1]
                      , model_copula_t$coefficients[2]
                      , summary(model_copula_t)$tableP1[2] #SE for time 0
                      , summary(model_copula_t)$tableP2[2] #SE for time 1
                      , logLik(model_copula_t)
                      , AIC(model_copula_t)
                      ,BIC(model_copula_t)
                      , 6
                      
    )
  }

  ########### 4. Combining results #########
  
  if(include=="ALL") {
    results_exbias<- rbind(summary_glm,summary_gee,summary_re_nosig,summary_re_np,summary_lme4,summary_cop,summary_cop_n,summary_cop_j,summary_cop_g,summary_cop_f,summary_cop_amh,summary_cop_fgm,summary_cop_pl,summary_cop_h,summary_cop_t,actuals)  
  }
  if(include=="GJRM") {
    results_exbias<- rbind(summary_cop,summary_cop_n,summary_cop_j,summary_cop_g,summary_cop_f,summary_cop_amh,summary_cop_fgm,summary_cop_pl,summary_cop_h,summary_cop_t,actuals)
  }
  if(include=="non-GJRM") {
    results_exbias<- rbind(summary_glm,summary_gee,summary_re_nosig,summary_re_np,summary_lme4,actuals)
  }
  
  results_exbias[,1:4]<- round(results_exbias[,1:4],4)
  results_exbias[,5:8]<- round(results_exbias[,5:8],0)
  
  colnames(results_exbias)<-c("b_1","b_2","se_b1","se_b2","LogLik","AIC","BIC","EDF")
  
  return(results_exbias)
  
}

plotDist <- function (dataset,dist) {
  
  require(gamlss)  
  require(latex2exp)
  require(ggplot2)
  require(ggpubr)
  
  num_margins=length(unique(dataset[,"time"]))
  
  margin_data=list()
  margin_unif=list()
  margin_fit=list()
  
  for (i in 1:num_margins) {
    margin_data[[i]]<-dataset[dataset[,"time"]==i,"response"]
    margin_fit[[i]]<-gamlss(margin_data[[i]]~1,family=dist)
    margin_unif[[i]]<-pnorm(margin_fit[[i]]$residuals)
  }
  
  ##plot.new()
  #par(mfrow=c(1,num_margins))
  
  #for (i in 1:num_margins) {histDist(margin_data[[i]],family=dist,xlab=TeX(paste("$Y_",i,"$")),main=paste("Histogram of margin",i,"and fitted",dist))}
  #invisible(readline(prompt="Press [enter] to continue"))
  
  plots=list()
  
  z=1
  for (i in 1:(num_margins)) {
      for (j in 1:(num_margins)) {
        if(i==j) {
          input_data=data.frame(margin_data[[i]])
          colnames(input_data)<-"X1"
          p <- ggplot(input_data, aes(x=X1)) + 
            geom_histogram() + 
            labs(x = TeX(paste("$Y_",i,"$")))
        }
        if(i!=j) {
          input_data=data.frame(cbind(margin_unif[[i]],margin_unif[[j]]))
          p=ggplot(data=input_data,aes(x=X1,y=X2)) +
            #geom_point(size=0.25,color="black") + 
            geom_density_2d(contour_var="density",bins=20,color="black") + 
            scale_fill_brewer() +
            labs(x = TeX(paste("$Y_",i,"$")), y=TeX(paste("$Y_",j,"$")),fill="density")
        }
        
        plots[[z]]=p
        z=z+1
      }
  }
  ggarrange(plotlist=plots,ncol=num_margins,nrow=num_margins)
  
}

generateMvtDist<-function(n,dist,mu_vector,sigma_vector,rho_vector) {
  
  if(dist=="NO") {
    require(MASS)
    cor_matrix<-diag(rep(1,length(sigma_vector)))
    if (length(rho_vector)==1) {
      cor_matrix[lower.tri(cor_matrix,diag=FALSE)]<-rho_vector[1]
      cor_matrix[upper.tri(cor_matrix,diag=FALSE)]<-rho_vector[1]
      
      cov_matrix<-diag(sigma_vector) %*% cor_matrix %*% diag(sigma_vector)
    } else {
      cor_matrix=cor_matrix+rho_vector
      cov_matrix<-diag(sigma_vector) %*% cor_matrix %*% diag(sigma_vector)
    }
    
    data<-mvrnorm(n,mu=mu_vector,Sigma = cov_matrix)
  }
  
  data_output<-cbind(data[,1],0)
  for (i in 2:ncol(data)) {
    data_output<-rbind(data_output,cbind(data[,i],i-1)) 
  }
  colnames(data_output) <- cbind("random_variable","time")
  return(data_output)
}


create_longitudinal_dataset <- function(response,covariates,labels=NA) {
  num_time_points=ncol(response)
  if(num_time_points <=1) {print('Not enough time points')}
  
  dataset<-matrix(data=NA,ncol=2+length(covariates),nrow=0) 
  subject<-as.factor(seq(1:nrow(response)))
  
  for (t in 1:ncol(response)) {
    
    dataset_temp<-cbind(subject,t,response[,t])
    
    for (i in 1:length(covariates)) {
      if (ncol(covariates[[i]]) == 1 ) {
        covariate_for_time=covariates[[i]]
      } else {
        covariate_for_time=covariates[[i]][,t]
      }
      dataset_temp<-cbind(dataset_temp,covariate_for_time)  
    }
    
    ###Add dataset temp to full table
    dataset <- rbind(dataset,dataset_temp)
  }
  
  if(!all(is.na(labels))) {
    colnames(dataset) <- labels  
  }
  
  dataset=dataset[order(dataset$subject,dataset$time),]
  
  return(dataset)
}

fit_margins <- function(mu.formula = formula(response~1), sigma.formula = ~1, 
                        nu.formula = ~1, tau.formula = ~1,dataset,margin_dist) {
  #i=1
  require(gamlss)
  fits<-list()
  for (i in 1:max(dataset$time)) {
    fits[[i]]<-gamlss(formula=mu.formula,family=margin_dist,data=dataset[dataset$time==i,])
  }
  return(fits)
}

find_best_marginal_fits <- function(dataset,type) {
  require(gamlss)
  fits<-list()
  for (i in 1:max(dataset$time)) {
    fits[[i]]<-fitDist(dataset[dataset$time==i,"response"],type = type)$fits[1:6]
  }
  return(fits)
}

fit_copulas = function(margin_fits,copula_dist="N",method="vine",dataset) {
  require(gamlss)
  require(VineCopula)
  
  fitted_margins_length=length(margin_fits)
  
  num_margins=length(unique(dataset$time))
  observation_count=nrow(dataset[dataset$time==unique(dataset$time)[1],])
  
  fit_unif=matrix(0,nrow=observation_count,ncol=num_margins)
  
  if(fitted_margins_length==1) {
    #margin_fits=fitted_margins #for testing
    for (i in unique(dataset$time)) {
      print(i)
      fit_unif[,i]=pZISICHEL(fitted_margins[[1]]$y[dataset$time==i]
                             ,mu=    margin_fits[[1]]$mu.fv[dataset$time==i]
                             ,sigma= margin_fits[[1]]$sigma.fv[dataset$time==i]
                             ,nu=    margin_fits[[1]]$nu.fv[dataset$time==i]
                             ,tau=   margin_fits[[1]]$tau.fv[dataset$time==i])
    }
  } else {
    
    for (i in 1:num_margins) {
      fit_unif[,i]=pNO(margin_fits[[i]]$residuals)    
    }
    
  }
  
  if (method=="novine") {
    cop_fit=list()  
    for (i in 1:(num_margins-1)){
      cop_fit[[i]]=BiCopEst(fit_unif[,i],fit_unif[,i+1],family=as.numeric(BiCopName(copula_dist)),se=TRUE)
    }
    names(cop_fit)=paste(1:(num_margins-1),(1:(num_margins-1))+1,sep=",")
    return(cop_fit)
  } else {
    vinefit<-RVineStructureSelect(fit_unif)
    return(vinefit)
  }
}

fit_copula_single = function(dataset,mm_mar,input_par,dFUN=dZISICHEL,pFUN=pZISICHEL,copula_dist="t") {
  #### MARGINAL MODEL: Get linear predictors and their inverse -> **Creates eta, eta_link_inv** from **mm_mar** ####
  margin_names=unique(dataset$time)
  num_margins = length(margin_names)
  #observation_count = nrow(dataset[dataset$time==margin_names[1],])
  response=dataset$response
  margin_parameters=names(mm_mar)
  
  #mm_all=list()  
  eta = eta_linkinv = matrix(nrow = nrow(dataset), ncol = length(margin_parameters))
  
  colnames(eta)=colnames(eta_linkinv)= margin_parameters
  
  for (parameter in margin_parameters) {
    par_temp <- input_par[grepl(paste(parameter,sep="."), names(input_par))]
    mm=mm_mar[[parameter]]
    #mm_all[[parameter]]=mm
    eta[, parameter] = rowSums(mm * matrix(rep(par_temp,each=nrow(mm)),ncol=length(par_temp),dimnames=list(NULL,c(names(par_temp))))) #These are the linear predictors
    FUN=margin_dist[[names(margin_dist)[grepl(paste(parameter,"linkinv",sep="."),names(margin_dist))]]]
    eta_linkinv[, parameter]=FUN(eta[,parameter])
  }
  
  #### MARGINAL MODEL: Get density and probability functions for margins -> ** Creates margin_p and margin_d ####
  
  margin_d = dFUN(x =   response
                  ,mu=    eta_linkinv[,"mu"]
                  ,sigma= eta_linkinv[,"sigma"]
                  ,nu=    eta_linkinv[,"nu"]
                  ,tau=   eta_linkinv[,"tau"])
  
  margin_p = pFUN(      response
                        ,mu=    eta_linkinv[,"mu"]
                        ,sigma= eta_linkinv[,"sigma"]
                        ,nu=    eta_linkinv[,"nu"]
                        ,tau=   eta_linkinv[,"tau"])
  
  margin_p_cop_input=cbind(margin_p[dataset$time %in% margin_names[1:(num_margins-1)]],margin_p[dataset$time %in% margin_names[2:(num_margins)]])
  
  cop_fit=BiCopEst(margin_p_cop_input[,1],margin_p_cop_input[,2],family=as.numeric(BiCopName(copula_dist)),se=TRUE)
  
  return(cop_fit)
  ####
}

likelihood_fitted_margin_copula <- function(fitted_margins,fitted_copulas) { 
  
  likelihood=fitted_copulas$logLik
  for (i in 1:length(fitted_margins)) {
    likelihood=likelihood+logLik(fitted_margins[[i]])
  }
  
  df=0
  for ( i in 1:length(fitted_margins)) {
    df=df+fitted_margins[[i]]$df.fit
  }
  df = df+(length(fitted_margins)*(length(fitted_margins)-1))
  
  return(
    paste("Overall LogLik:", as.character(round(likelihood,2)),
          "\nGlobal Deviance:", as.character(round(likelihood,2)*-2),
          "\nDegrees of Freedom:", as.character(round(df,2)))
  )
  
} ###List of gamlss objects, VineCopula object as input

logit <- function(x) {
  return(log(x/(1-x)))
}

logit_inv <- function(x) {
  return(
    exp(x)/(1+exp(x))
  )
}

log_2plus <- function(x) {
  return(
    log(x-2)
  )
}

log_2plus_inv <- function(x) {
  y=exp(x)+2
  #Adjust for error when close to two.
  y[y==2]=y[y==2]+0.00001
  return(y)
}

dlogit <- function(x) {
  return(1/(x-(x^2)))
}

dlog_2plus <- function(x) {
  return(1/(x-2))
}

#Generate the copula model matrix
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
  
  #return_option="list";input_par=input_par;copula_dist="t";margin_dist = ZISICHEL(mu.link = "log",sigma.link = "log",nu.link = "identity",tau.link = "logit");verbose=TRUE;dFUN=dZISICHEL;pFUN=pZISICHEL;calc_d2=TRUE
  
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
  
  margin_p = pFUN(    response ,mu=eta_linkinv[,"mu"] ,sigma=eta_linkinv[,"sigma"] ,nu=eta_linkinv[,"nu"],tau=   eta_linkinv[,"tau"])
  margin_d = dFUN(x = response ,mu=eta_linkinv[,"mu"] ,sigma=eta_linkinv[,"sigma"] ,nu=eta_linkinv[,"nu"],tau=   eta_linkinv[,"tau"])
  
  #### Get density of the copula function for each pair of margins ####
  
  #margin_names; num_margins; margin_parameters
  
  #Calculate new eta based on input parameters and model matrix
  copula_parameters=c("theta","zeta")
  
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
  
  copula_d=dcdu1=dcdu2=dldth=d2ldth=dldz=d2ldz=dthdeta=dzdeta=matrix(0,nrow=nrow(eta_cop),ncol=1)
  
  margin_p_cop_input=cbind(margin_p[dataset$time %in% margin_names[1:(num_margins-1)]],margin_p[dataset$time %in% margin_names[2:(num_margins)]])
  
  copula_d=BiCopPDF(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"])
  dldth=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="par",log=TRUE)
  dcdu1=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="u1",log=FALSE)
  dcdu2=BiCopDeriv(   margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="u2",log=FALSE)
  
  if(calc_d2==TRUE) {
    d2ldth=BiCopDeriv2( margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="par")
  }
  dthdeta=copula_link$dthdeta(eta_cop[,"theta"])
  
  if (length(copula_parameters)>1) {
    dldz=BiCopDeriv(  margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="par2",log=TRUE)
    if(calc_d2==TRUE) {
      d2ldz=BiCopDeriv2(margin_p_cop_input[,1],margin_p_cop_input[,2],family = as.vector(BiCopName(copula_dist)),par=eta_cop_link_inv[,"theta"],par2=eta_cop_link_inv[,"zeta"],deriv="par2")
    }
    dzdeta=copula_link$dzdeta(eta_cop[,"zeta"])
  }
  
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
      derivatives_calculated[,i]=margin_dist[[d_fun]](y=response,mu=eta_linkinv[,"mu"],sigma=eta_linkinv[,"sigma"],nu=eta_linkinv[,"nu"],tau=eta_linkinv[,"tau"])
    }
    i=i+1
  }
  colnames(derivatives_calculated)=all_derivatives
  derivatives_calculated_all_margins=derivatives_calculated
  #}
  
  if (calc_d2==TRUE) {
    derivatives_copula = list(dldth,dldz,d2ldth,d2ldz,dcdu1,dcdu2,dthdeta,dzdeta)
    names(derivatives_copula)=c("dldth","dldz","d2ldth","d2ldz","dcdu1","dcdu2","dthdeta","dzdeta")
  } else {
    derivatives_copula = list(dldth,dldz,dcdu1,dcdu2,dthdeta,dzdeta)
    names(derivatives_copula)=c("dldth","dldz","dcdu1","dcdu2","dthdeta","dzdeta")
  }
  
  return_list=list(log_lik_results,response,margin_d,margin_p,copula_d,derivatives_copula,margin_dist,copula_dist,eta,eta_linkinv,derivatives_calculated_all_margins,mm_mar,start_par,eta_cop,eta_cop_link_inv,mm_cop)
  names(return_list)=c("log_lik_results","response","margin_d","margin_p","copula_d","derivatives_copula","margin_dist","copula_dist","eta","eta_linkinv","derivatives_calculated_all_margins","mm_mar","start_par","eta_cop","eta_cop_link_inv","mm_cop")
  
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
    #Add back in for multiple margins
    #names(mu[[i]])=paste(names(mu[[i]]),i,sep=".")
    #names(sigma[[i]])=paste(names(sigma[[i]]),i,sep=".")
    #names(nu[[i]])=paste(names(nu[[i]]),i,sep=".")
    #names(tau[[i]])=paste(names(tau[[i]]),i,sep=".")
  }  
  
  if (!(is.na(fitted_copulas))) {
    theta=zeta=vector()
    for (i in 1:length(fitted_copulas)) {
      theta[i]=copula_link$theta.linkfun(fitted_copulas[[i]]$par)
      zeta[i]=copula_link$zeta.linkfun(fitted_copulas[[i]]$par2)
      names(theta)[i]=paste("theta",names(fitted_copulas)[i])
      names(zeta)[i]=paste("zeta",names(fitted_copulas)[i])
    }
    return_list=(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau),theta,zeta))
  } else {
    return_list=(c(unlist(mu),unlist(sigma),unlist(nu),unlist(tau)))
  }
  
  return(return_list)
}

score_function <- function(eta,dldpar,dpardeta,dlcopdpar,response,phi=1,step_size=1,verbose=TRUE) {
  #u_k
  dlcopdpar=dlcopdpar*((dldpar)*response)
  
  dldpar=dldpar+dlcopdpar #commented out for now until we come back and fix it
  
  u_k=dldeta = dldpar * dpardeta
  
  #f_k and w_k
  #f_k_quasi_newton=(-(dldpar*dldpar))
  f_k=-(dldpar^2)
  
  w_k=-f_k*(dpardeta*dpardeta)
  
  #w_k=matrix(1,ncol=ncol(f_k),nrow=nrow(f_k)) ####Weights of 1
  
  z_k=(1-phi)*eta+phi*(eta+step_size*(u_k/w_k))
  
  if(verbose==TRUE) {
    steps_mean=round(rbind(colMeans(eta)
                           ,colMeans(dldpar-dlcopdpar)
                           ,colMeans(dlcopdpar)
                           ,colMeans(dpardeta)
                           ,colMeans(dpardeta*dpardeta)
                           ,colMeans(f_k)
                           ,colMeans(w_k)
                           ,colMeans(u_k)
                           ,colMeans((1/w_k)*u_k)
                           ,colMeans(z_k)
    ),8)
    rownames(steps_mean)=c("eta","dldpar","dlcopdpar","dpardeta","dpardeta2","f_k","w_k","u_k","(1/w_k)*u_k","z_k")
    print(steps_mean)  
  }
  return_list=list(colMeans(z_k),u_k,f_k,w_k,z_k)
  names(return_list)=c("par","u_k","f_k","w_k","z_k")
  return(return_list)
}

newton_raphson_iteration=function(results,input_par,phi=1,step_size=1,verbose=c(FALSE,FALSE),calc_d2=FALSE) {
  
  #results=final_results;input_par=end_par_matrix[which.max(end_loglik_matrix),] ;verbose=c(TRUE,TRUE); calc_d2=TRUE
  
  steps_all=list()
  par_out=list()
  margin_dist=results$margin_dist
  copula_dist=results$copula_dist
  
  ### Extract COPULA DERIVATIVES from results list (from calc_joint_likelihood)
  
  dlcopdpar=list()
  dldth=results$derivatives_copula$dldth
  dldz=results$derivatives_copula$dldz
  dthdeta=results$derivatives_copula$dthdeta
  dzdeta=results$derivatives_copula$dzdeta
  dcdu1=results$derivatives_copula$dcdu1
  dcdu2=results$derivatives_copula$dcdu2
  
  dldpar=cbind(dldth,dldz)
  #d2ldpar=cbind(d2ldth,d2ldz)
  dpardeta=cbind(dthdeta,dzdeta)
  #dpardeta=dldpar*0+1 ##########Temporary step
  
  dlcopdpar=(dcdu1+dcdu2)/(results$copula_d) #+ (dcdu2)/results$copula_d[,i]  #*(results$margin_d[,i]*results$margin_d[,i+1])
  
  eta_input=cbind(results$eta_cop[,"theta"],results$eta_cop[,"zeta"])
  colnames(eta_input)=names(results$eta_cop)
  copula_score=score_function(eta=eta_input,dldpar=dldpar,dpardeta=dpardeta,dlcopdpar=dldpar*0,response=dldpar*0,phi=phi,step_size=step_size,verbose=verbose[2]) ##Returns updated value
  
  beta_new_cop=list()
  #parameter_names=c("theta","zeta")
  cop_parameter_names=c("theta","zeta")
  
  hess_list=list()
  
  e_k=copula_score$z_k  
  beta_new_cop=list()
  for (j in 1:ncol(copula_score$z_k)) {
    mm=mm_cop[[j]]
    X=as.matrix(mm)
    W=diag(copula_score$w_k[,j])
    e_k_par=e_k[,j]
    beta_new_cop[[cop_parameter_names[j]]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
    rownames(beta_new_cop[[cop_parameter_names[j]]])= paste(paste(cop_parameter_names[j],sep=" "),colnames(mm),sep=".")
    #names(beta_new_cop[[j]])= paste(parameter_names[j],margin_number,rownames(beta_new_cop[[j]]),sep=".")
    if (calc_d2 == TRUE) {
      hess_list[[(cop_parameter_names[j])]]=-t(X)%*%W%*%X
    }
  }
  
  one=cbind(dataset[dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)],c("time","subject")],dlcopdpar)
  two=cbind(dataset[!(dataset$time %in% unique(dataset$time)[1:(length(unique(dataset$time))-1)]),c("time","subject")],0)
  colnames(two)=colnames(one)
  one_two=rbind(one,two)
  
  one=cbind(dataset[!(dataset$time %in% unique(dataset$time)[2:(length(unique(dataset$time)))]),c("time","subject")],0)
  two=cbind(dataset[dataset$time %in% unique(dataset$time)[2:(length(unique(dataset$time)))],c("time","subject")],dlcopdpar)
  colnames(one)=colnames(two)
  two_three=rbind(one,two)
  
  one_two=one_two[order(one_two$subject,one_two$time),]
  two_three=two_three[order(two_three$subject,two_three$time),]
  
  ###CHECKS ON ORDERING
  all(one_two$time==two_three$time)
  all(one_two$time==dataset$time)
  all(one_two$subject==two_three$subject)
  all(one_two$subject==dataset$subject)
  
  dlcopdpar_one_col=(one_two$dlcopdpar+two_three$dlcopdpar)*results$margin_d
  
  margin_parameters=colnames(results$eta)
  
  dlcopdpar_final=dlcopdpar_one_col
  response_final=results$response
  for (i in 2:length(margin_parameters)) {
    dlcopdpar_final=cbind(dlcopdpar_final,dlcopdpar_one_col)
    response_final=cbind(response_final,results$response)
  }
  
  dldpar=results$derivatives_calculated_all_margins[,c("dldm","dldd","dldv","dldt")] #First derivative of log likelihood w.r.t. parameters
  
  if (calc_d2==TRUE) {
    d2ldpar=results$derivatives_calculated_all_margins[,c("d2ldm2","d2ldd2","d2ldv2","d2ldt2")] #Quasi newton second derivative of log function w.r.t parameters
  }
  
  dpardeta=results$derivatives_calculated_all_margins[,c("mu.dr","sigma.dr","nu.dr","tau.dr")] #Derivative of inverse link function
  
  margin_score=score_function(eta=results$eta[,margin_parameters],dldpar=dldpar,dpardeta=dpardeta,dlcopdpar=dlcopdpar_final*0,response=response_final,phi=phi,step_size=step_size,verbose=verbose[1]) ##Returns updated value
  
  e_k=margin_score$z_k#-results$eta_nomargin[[1]][,c("mu","sigma","nu","tau")]
  beta_r=list()
  #for (margin_number in 1:length(results$mm_all)) {
  for (parameter in margin_parameters) {
    mm=results$mm_mar[[parameter]]
    X=as.matrix(mm)
    W=diag(margin_score$w_k[,which(margin_parameters==parameter)])
    e_k_par=e_k[,which(margin_parameters==parameter)]
    beta_r[[parameter]]=solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%e_k_par
    rownames(beta_r[[parameter]])= paste(parameter,rownames(beta_r[[parameter]]),sep=".")
    if (calc_d2 == TRUE) {
      hess_list[[parameter]]=-t(X)%*%W%*%X
    }
  }
  
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
  end_par=end_par[names(start_par)]
  
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
}

loadDataset <- function(simOption,plot=FALSE) {

  if (simOption==1) {
    load("Data/rand_mvt.rds")
    head(rand_mvt)
    
    # Basic data setup
    response = rand_mvt[,4:(4+2)]#[,4:18] ####Currently limiting to just 5 margins for simplicity
    covariates=list()
    covariates[[1]] = as.data.frame(rand_mvt[,19]) #Age 19:33 - changed to age at start to avoid correlation with time
    covariates[[2]] = as.data.frame(rand_mvt[,34:(34+4)]) #Time 34:48
    covariates[[3]] = as.data.frame(rand_mvt[,3]) #Gender
    
    # Setup data as longitudinal file
    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
  } else if (simOption==2) {
    # set up D-vine copula model with mixed pair-copulas
    d <- 3
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(2, 2, 0)
    par <- c(logit_inv(.4), logit_inv(.4), logit_inv(.4))
    par2 <- c(log_2plus_inv(1.1),log_2plus_inv(1.1),log_2plus_inv(1.1))
    
    
    # transform to R-vine matrix notation
    RVM <- D2RVine(order, family, par, par2)
    contour(RVM)
    
    n=1000;t=d
    copsim=RVineSim(n*t,RVM)
    
    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {
      margin[,i]=qZISICHEL(copsim[,i],mu=exp(0.7*i),sigma=exp(0.3*i),nu=-1.6,tau=0.05)#Update to i*mu/sigma as needed
    }
    
    response = as.data.frame(margin)
    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender
    
    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
    
  }
  if(plot==TRUE) {plotDist(dataset,"ZISICHEL")}
  return(dataset)
}
