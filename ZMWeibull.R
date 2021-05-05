#Reading data
library(xlsx)
AllDataN <-read.xlsx("/Users/marzieh/Desktop/Data_new.xlsx",1)   
head(AllDataN)


install.packages("DWreg")
library("DWreg")

ddwf<-function(q,beta,y) (ddw(y,q,beta))
Sddwf<-function(q,beta,y) (ddwf(q,beta,y))

#Mass probability function of zero modified Weibull
zddwf<-function(p1,q,beta,y) ifelse(y==0,p1+(1-p1)*Sddwf(q,beta,0),(1-p1)*Sddwf(q,beta,y))   


# log likelihood of the mentioned mass function
fun = function(p1,q,beta,y) zddwf(p1,q,beta,y)
fun2 = function(x, y) sum(log(fun(x[1],x[2], x[3], y)))  



d2<-dim(AllDataN)[2]

LLZMWEIBULL<-PW<-Q<-BET<-SEPW <-SEQ <-SEBET <-LLZIWEIBULL<- CIQL<- CIQU <- CIBETL <-CIBETU <-CIPWL <-CIPWU <-pvalue_Q <-pvalue_BET<- pvalue_PW <-rep(999,d2)
for(j in 1:d2){
  
  xxx<-AllDataN[,j]
  xxx<-xxx[!is.na(xxx)]
  xx<-xxx[xxx>-1] # (citation data)
  
  #optimize the log likelihood w.r.t specific initial values of parameters
  estwz = optim(c(0,.9,.8),fun2,control=list(fnscale=-1),y=xx,hessian=TRUE)   

 
  # standard error of parameters
  SEPW[j]<-round(sqrt(diag(solve(-estwz$hessian)))[1],4)
  
  SEQ[j]<-round(sqrt(diag(solve(-estwz$hessian)))[2],4)
  
  SEBET[j]<-round(sqrt(diag(solve(-estwz$hessian)))[3],4)
  
  
  # estimation of parameters
  PW[j]<-round(estwz$par[1],3)
  Q[j]<-round(estwz$par[2],2)
  BET[j]<-round(estwz$par[3],2)
  LLZMWEIBULL[j]<-round(estwz$value,2)
  
  #Confidence Interval (95%)
  CIQL[j]<-round(Q[j]-(1.96* SEQ[j]),3)
  CIQU[j]<-round(Q[j]+(1.96* SEQ[j]),3)
  
  CIBETL[j]<-round(BET[j]-(1.96*SEBET[j]),3)
  CIBETU[j]<-round(BET[j]+(1.96*SEBET[j]),3)
  
  CIPWL[j]<-round(PW[j]-(1.96*SEPW[j]),3)
  CIPWU[j]<-round(PW[j]+(1.96*SEPW[j]),3)
  
  
  #Wald Test (two-tailed P_value)
  
  pvalue_Q[j]<-round(2*pnorm(-abs((Q[j]-.9)/SEQ[j])),5)         #for Q (H0: Q=.9)    
  pvalue_BET[j]<-round(2*pnorm(-abs((BET[j]-.8)/SEBET[j])),5)   #for BET (H0: BET=.8)
  pvalue_PW[j]<-round(2*pnorm(-abs((PW[j]-0)/SEPW[j])),5)       #for PW (H0: PW=0)
  
  

} 

#The name of the subjects
Subject<-c("FoodScience","CancerResearch","Marketing","Filtration","PTChem","CSA","MSOR","GeoChem","Economics","Energy","
           CompMech","
           GlobalPC","
           Virology","
           MetalsAlloys","
           Control","
           CCICM","
           DevNeuro","
           PharmSci","
           NHEP","
           NNPP","
           HSS","
           CS","
           HealthIM")


# Displaying the outputs
DD<-data.frame(Subject,PW,SEPW,CIPWL,CIPWU,pvalue_PW,Q,SEQ, CIQL ,CIQU, pvalue_Q, BET,SEBET,CIBETL,CIBETU,pvalue_BET ,LLZMWEIBULL)
names(DD)<-c("Subject","PW","SEPW","CIPWL","CIPWU","pvalue_PW","Q","SEQ", "CIQL ","CIQU", "pvalue_Q", "BET","SEBET","CIBETL","CIBETU","pvalue_BET" ,"log ZIWEIBULL")   

print(DD)
write.csv(DD,file="Result_ZMWeibull.csv")