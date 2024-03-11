source("common_functions.R")

data<-generateMvtDist("NO",c(1,2,3),c(1,2,3),matrix(c(0,.5,.1,.5,0,.9,.1,.9,0),nrow=3))
plotDist(data,"NO")

data<-generateMvtDist("NO",c(1,2,3),c(1,2,3),matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),nrow=3))
plotDist(data,"NO")


data<-generateMvtDist("NO",c(1,2,3),c(1,1,1),matrix(c(0,.5,.1,.5,0,.9,.1,.9,0),nrow=3))
plotDist(data,"NO")

data<-generateMvtDist("NO",c(1,2,3),c(1,1,1),matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),nrow=3))
plotDist(data,"NO")


help(gjrm)

require(GJRM)

data_input<-data.frame(cbind(data[data[,"time"]==0,"random_variable"], data[data[,"time"]==1,"random_variable"], data[data[,"time"]==2,"random_variable"]))
colnames(data_input)<-c("V1","V2","V3")

eq.mu.1 <- formula(V1~1)
eq.mu.2 <- formula(V2~1)
eq.mu.3 <- formula(~1)
eq.theta.1 <- formula(~1)
eq.theta.2 <- formula(~1)
eq.theta.2 <- formula(~1)
eq.theta.2 <- formula(~1)
fl <- list(eq.mu.1, eq.mu.2, eq.mu.3,eq.theta.1,eq.theta.2)

t_model<-gjrm(fl,margins=c('N','N','N'),copula="N",copula2="N",model="T",data=data_input)

summary(t_model)

