#Reading data
library(xlsx)
AllDataN <-read.xlsx("/Users/marzieh/Desktop/Data_new.xlsx",1)       
head(AllDataN)


#Mass probability function of zero modified hooked power law
dhpl<-function(B,alpha,y) ((B+y)^(-alpha))/(sum((B+0:1000)^(-alpha)))
zdhpl<-function(p1,B,alpha,y) ifelse(y==0,p1+(1-p1)*dhpl(B,alpha,0),(1-p1)*dhpl(B,alpha,y))     



fun = function(p1,B, alpha, y) zdhpl(p1,B, alpha,y)

# log likelihood of the mentioned mass function
fun2 = function(x, y) sum(log(fun(x[1], x[2],x[3], y)))      



d2<-dim(AllDataN)[2]

BETA<-ALPHA<- SEBETA<-SEALPHA<-SEP1<-SEB<-SEalpha<- P1<-B<-alpha<-LLZMHPL<-CIP1L<-CIP1U<-CIBL<-CIBU<-CIalphaL<-CIalphaU<-pvalue_B<-pvalue_alpha<-pvalue_p1<-rep(999,d2)
for(j in 1:d2)
  {
  xxx<-AllDataN[,j]
  xxx<-xxx[!is.na(xxx)]
  xx<-xxx[xxx>-1] # (citation data)
  
  #optimize the log likelihood w.r.t specific initial values of parameters
  esthz = optim(c(0,80,5),fun2,control=list(fnscale=-1),y=xx,hessian=TRUE)  
  
  
  # standard error of parameters
  SEP1[j]<-round(sqrt(diag(solve(-esthz$hessian)))[1],4)
  
  SEB[j]<-round(sqrt(diag(solve(-esthz$hessian)))[2],4)
  
  SEalpha[j]<-round(sqrt(diag(solve(-esthz$hessian)))[3],4)
  
  
  # estimation of parameters
  P1[j]<-round(esthz$par[1],3)
  B[j]<-round(esthz$par[2],2)
  alpha[j]<-round(esthz$par[3],2)
  LLZMHPL[j]<-round(esthz$value,2)
  
  #Confidence Interval (95%)
  CIBL[j]<-round(B[j]-(1.96* SEB[j]),3)
  CIBU[j]<-round(B[j]+(1.96* SEB[j]),3)
  
  CIalphaL[j]<-round(alpha[j]-(1.96*SEalpha[j]),3)
  CIalphaU[j]<-round(alpha[j]+(1.96*SEalpha[j]),3)
  
  CIP1L[j]<-round(P1[j]-(1.96*SEP1[j]),3)
  CIP1U[j]<-round(P1[j]+(1.96*SEP1[j]),3)
  
  
  #Wald Test (two-tailed P_value)
  pvalue_B[j]<-round(2*pnorm(-abs((B[j]-80)/SEB[j])),5)           #for B (H0: B=80)
  pvalue_alpha[j]<-round(2*pnorm(-abs((alpha[j]-5)/SEalpha[j])),5)  #for alpha (H0: alpha=5)
  pvalue_p1[j]<-round(2*pnorm(-abs((P1[j]-0)/SEP1[j])),5)            #for P1 (H0: P1=0)
  
  
}

#The names of subjects
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

# displaying the outputs
DD<-data.frame(Subject,P1,SEP1,CIP1L,CIP1U,pvalue_p1,B,SEB,CIBL,CIBU,pvalue_B,alpha,SEalpha, CIalphaL,CIalphaU,pvalue_alpha,LLZMHPL)
names(DD)<-c("Subject","P","SEP1","CI_Pl","CI_PU","p_value_P","B","SEB","CI_Bl","CI_BU","p_value_B","alpha","SEalpha","CI_alphal","CI_alphaU","p_value_alpha","loglik")   


print(DD)
write.csv(DD,file="Result_ZMHPL.csv")