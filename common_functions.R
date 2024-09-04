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
    margin_unif[[i]]<-(margin_fit[[i]]$residuals)
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
            geom_density_2d(contour_var="density",bins=10,color="black") + 
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
  rownames(dataset)=1:nrow(dataset)
  
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

fit_copula_single = function(dataset,margin.family,mm_mar,input_par,dFUN=dZISICHEL,pFUN=pZISICHEL,copula_dist="t") {
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
    FUN=margin.family[[names(margin.family)[grepl(paste(parameter,"linkinv",sep="."),names(margin.family))]]]
    eta_linkinv[, parameter]=FUN(eta[,parameter])
  }
  
  #### MARGINAL MODEL: Get density and probability functions for margins -> ** Creates margin_p and margin_d ####
  
  input_list=list(x=response,q=response)
  for (parameter in margin_parameters) {
    input_list[[parameter]]=eta_linkinv[,parameter]
  }
  
  args=names(input_list)[names(input_list)%in%formalArgs(dFUN)]
  margin_d = do.call(dFUN,args=input_list[args])
  args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
  margin_p = do.call(pFUN,args=input_list[args])

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
  return(-1/(x-2))
}

dlog <-function(x) {
  return(1/x)
}

dlog_inv <-function(x) {
  return(exp(x))
}


dlogit_inv <- function(x) {
  return(exp(x)/((1+exp(x))^2))
}

dlog_2plus_inv <- function(x) {
  return(exp(x))
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

score_function <- function(eta,dldpar,d2ldpar,dpardeta,dlcopdpar,response,phi=1,step_size=1,verbose=TRUE,crit_wk=0.0000001) {
  
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
    
    #WHOLE THING
    # dlcopdpar[,i]=0.5*
    # (dc_tplus_du_t * du_t_dmu + 
    # dc_tplus_du_tplus * du_tplus_dmu ) /
    # c_tplus +
    # (dc_tminus_du_tminus * du_tminus_dmu +
    # dc_tminus_du_t * du_t_dmu) / 
    # c_tminus

    ########OLD
    
    #x=margin_copula_merged_2[,"results$response"]
    #f_x=margin_copula_merged_2[,"results$margin_d"]
    #F_x=margin_copula_merged_2[,"results$margin_p"]
    #dldm_temp=margin_copula_merged_2[,paste(margin_deriv_names[i],"",sep="")]
    #d_margin_cdf1_approx=f_x+x*f_x#*dldm_temp#*dldm_temp
    #l=margin_copula_merged_2[,"mu"]#[order(dataset$time,dataset$subject),]
    #d_margin_cdf1=d_margin_cdf2=-x*exp(-x*(1/l)) * (l^-2)
    
      
    #dcdu1_temp=margin_copula_merged_2[,"dcdu1.x"] #These are correct 
    #dcdu2_temp=margin_copula_merged_2[,"dcdu2.y"] #These are correct 
    
    #dcdu1_f1=dcdu1_temp*d_margin_cdf1/margin_copula_merged_2[,"results$copula_d.x"]
    #dcdu2_f2=-dcdu2_temp*d_margin_cdf1/margin_copula_merged_2[,"results$copula_d.y"]
    
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
  
  #plot.new()
  #hist(d_margin_cdf1)
  #hist(d_margin_cdf1_approx)
  #plot(d_margin_cdf1,d_margin_cdf1_approx)
  #plot(d_margin_cdf1,d_margin_cdf1_approx)
  
  ####OK Fantastic so dcdu1 and dcdu2 are correct.
  #plot(results$num_dcdu1,dcdu1)
  #plot(results$num_dcdu2,dcdu2)
  
  #plot(dlcopdpar,as.vector(margin_copula_merged_2[,"results$num_dlcopdpar"])) #OK so this is exact

  #plot.new(); par(mfrow=c(1,3))
  #plot(dcdu1_f1,as.vector(margin_copula_merged_2[,"results$num_dlcopdpar.x"])) #OK so this is exact
  #plot(dcdu2_f2,as.vector(margin_copula_merged_2[,"results$num_dlcopdpar.y"])) #OK so this is exact
  #plot(dlcopdpar,as.vector(margin_copula_merged_2[,"results$num_dlcopdpar.x"]+margin_copula_merged_2[,"results$num_dlcopdpar.y"]))
  
  #hist(dlcopdpar)
  #hist(as.vector(margin_copula_merged_2[,"results$num_dlcopdpar.x"]))
  #hist(as.vector(margin_copula_merged_2[,"results$num_dlcopdpar.y"]))
  
  #plot(dlcopdpar,as.vector(margin_copula_merged_2[,"results$num_dlcopdpar.x"]+margin_copula_merged_2[,"results$num_dlcopdpar.y"])) #Temporary
  
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

loadDataset <- function(simOption,plot_dist=FALSE,n=100,d=3,copula.family=1,copula.link,qFUN,par.copula,par.margin) {

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
    par <- c(logit_inv(.8), logit_inv(.8), logit_inv(.8))
    par2 <- c(log_2plus_inv(2.1),log_2plus_inv(2.1),log_2plus_inv(2.1))
    
    # transform to R-vine matrix notation
    RVM <- D2RVine(order, family, par, par2)
    contour(RVM)
    
    t=d
    copsim=RVineSim(n*t,RVM)
    
    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender
    
    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {
      margin[covariates[[3]]==0,i]=qZISICHEL(copsim[,i],mu=exp(0.3+0.2*i),sigma=exp(0.3+0.2*i),nu=-0.8,tau=0.05)[covariates[[3]]==0]#Update to i*mu/sigma as needed
      margin[covariates[[3]]==1,i]=qZISICHEL(copsim[,i],mu=exp(0.3+0.2*i+0.1),sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.05)[covariates[[3]]==1]#Update to i*mu/sigma as needed
    } 
    
    response = as.data.frame(margin)
    
    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
    
  } else if (simOption==3) {
    
    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))
    
    
    if(length(par.copula)==d-1){
       par=c(copula.link$theta.linkinv(par.copula),rep(0,dd-(d-1)))
       par2=par*0
    } else {
      par <- c(copula.link$theta.linkinv(par.copula[1:(length(par.copula)/2)]), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(copula.link$zeta.linkinv(par.copula[(length(par.copula)/2+1):(length(par.copula))]),rep(0,dd-(d-1))) #+0.5*1:(d-1)  
    }

    # transform to R-vine matrix notation
    
    RVM <- D2RVine(order, family, par, par2)
    #contour(RVM)
    
    t=d
    copsim=RVineSim(n,RVM)
    
    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender
    
    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {
      
      input_list=list(p=copsim[,i],mu=exp(par.margin[1]+par.margin[2]*i),sigma=exp(0.3+0.2*i),nu=-0.8,tau=0.1)
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      input_list=list(p=copsim[,i],mu=exp(par.margin[1]+par.margin[2]*i+par.margin[3]),sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_2=do.call(qFUN,args=(input_list[args]))
      
      margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed
    } 
    
    response = as.data.frame(margin)
    
    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
    
  }
  
  
  
  if(plot_dist==TRUE) {plotDist(dataset,"ZISICHEL")}
  return(dataset)
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
  #Note doesn't print first set of starting parameters in plot - to fix
  
  stopifnegative=stopifnegative
  while ((abs(change) > stopval) & run_counter <= max_runs) { #
    if(!first_run) {start_log_lik=results$log_lik_results} else {input_par=first_input_par; phi_inner=phi}
    
    results= calc_joint_likelihood(input_par = input_par  #optim_results$par
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
    
    end_log_lik=results$log_lik_results
    
    if(first_run) {change=1} else {
      change=end_log_lik["Total"]-start_log_lik["Total"]
      log_lik_output=c(start_log_lik["Total"],end_log_lik["Total"],change,phi_inner)
      names(log_lik_output) = c("Start LogLik","End LogLik","Change","phi")
      print(log_lik_output)
    }
    
    iteration_out=newton_raphson_iteration(results
                                           ,input_par
                                           ,phi=phi_inner
                                           ,step_size=step_size
                                           ,verbose=c(verbose,verbose)
                                           ,crit_wk=crit_wk
                                           ,calc_d2=calc_d2)
    
    cbind(iteration_out[[1]],iteration_out[[2]])
    end_par=iteration_out[[2]]
    input_par=end_par
    first_run=FALSE
    end_par_matrix=rbind(end_par_matrix,end_par)
    end_loglik_matrix=rbind(end_loglik_matrix,end_log_lik)
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
    
    
    if(stopifnegative==TRUE) {
      if(nrow(end_loglik_matrix)>=3) {
        if(all(end_loglik_matrix[(length(end_loglik_matrix)-3):(length(end_loglik_matrix)-1)]-end_loglik_matrix[(length(end_loglik_matrix)-2):length(end_loglik_matrix)]>0))
          {converged_status="Negative";break}
      }
      if((change/end_log_lik["Total"])< -0.1) {converged_status="Negative";break}
    }
    
    if(abs(change) <= stopval) {converged_status="Converged"} else {converged_status="Not Converged"}
    
    phi_inner=phi*(step_adjustment^min(step_reduc_count,run_counter-1))
  }
  
  if(plot_end_only==TRUE & plot_optim==TRUE) {
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
  
  #Extract final parameters
  end_par=optim$end_par_matrix[nrow(optim$end_par_matrix),]
  
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
  names(return_list)=c("end_par","new_SE","final_results")
  return(return_list)
}
  






























