library(xlsx)
AllDataN <-read.xlsx("/Users/marzieh/Desktop/Data_new.xlsx",1)    #Reading data
head(AllDataN)


check.integer <- function(t) {
  t == round(t)
}
ddlnorm=function(mu,sigma,y) 
  ifelse( (check.integer(y)&y>=1)==TRUE,(plnorm(y+1/2,mu,sigma)-plnorm(y-1/2,mu,sigma))/(1-plnorm(1/2,mu,sigma)),0)

zdlnorm<-function(p,mu,sigma,y) ifelse(y==1, p+(1-p)*ddlnorm(mu,sigma,y), (1-p)*ddlnorm(mu,sigma,y)) #density function of zero modified dis_log_normal




fun = function(p,mu, sigma, y) (zdlnorm(p,mu, sigma,y))
fun2 = function(x, y) sum(log(fun(x[1], x[2],x[3], y)))  # log likelihood of the mentioned density


d2<-dim(AllDataN)[2]

SEP<-SEMU<-SESIGMA<- P<-MU<-SIGMA<-LLZMDLN_new_density<-CIPL<-CIPU<-CIMUL<-CIMUU<-CISIGMAL<-CISIGMAU<-pvalue_mu<-pvalue_sigma<-pvalue_p<-rep(999,d2)
for(j in 1:d2){
  
  xxx<-AllDataN[,j]
  xxx<-xxx[!is.na(xxx)]
  xx<-xxx[xxx>-1]+1 # (citation data +1)

  
  estlz = optim(c(0.03,mean(log(xx)),sd(log(xx))),fun2,control=list(fnscale=-1),y=xx,hessian=TRUE)  #optimize the log likelihood w.r.t specific initial values of parameters
  
  
  # standard error of parameters
  SEP[j]<-round(sqrt(diag(solve(-estlz$hessian)))[1],4) 
  
  SEMU[j]<-round(sqrt(diag(solve(-estlz$hessian)))[2],4)
  
  SESIGMA[j]<-round(sqrt(diag(solve(-estlz$hessian)))[3],4)
  
  
  # estimation of parameters
  P[j]<-round(estlz$par[1],3)   
  MU[j]<-round(estlz$par[2],2)
  SIGMA[j]<-round(estlz$par[3],2)
  LLZMDLN_new_density[j]<-round(estlz$value,2)
  
  
  
  #Confidence Interval (95%)
  CIMUL[j]<-round(MU[j]-(1.96* SEMU[j]),3)
  CIMUU[j]<-round(MU[j]+(1.96* SEMU[j]),3)
  
  CISIGMAL[j]<-round(SIGMA[j]-(1.96*SESIGMA[j]),3)
  CISIGMAU[j]<-round(SIGMA[j]+(1.96*SESIGMA[j]),3)
  
  CIPL[j]<-round(P[j]-(1.96*SEP[j]),3)
  CIPU[j]<-round(P[j]+(1.96*SEP[j]),3)
  
  
  #Wald Test (two-tailed P_value)
  pvalue_mu[j]<-round(2*pnorm(-abs((MU[j]-2)/SEMU[j])),5)           #for MU (H0: MU=2)
  pvalue_sigma[j]<-round(2*pnorm(-abs((SIGMA[j]-1)/SESIGMA[j])),5)  #for SIGMA (H0: SIGMA=1)
  pvalue_p[j]<-round(2*pnorm(-abs((P[j]-0)/SEP[j])),5)            #for P (H0: P=0)
 
}


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
DD<-data.frame(Subject,P,SEP,CIPL,CIPU,pvalue_p,MU,SEMU,CIMUL,CIMUU,pvalue_mu,SIGMA,SESIGMA,CISIGMAL,CISIGMAU,pvalue_sigma,LLZMDLN_new_density)
names(DD)<-c("Subject","P","SEP","CI_Pl","CI_Pu","p_value_P","mu","SEmu","CI_mul","CI_muu","p_value_mu","sigma","SEsigma", "CI_sigmal","CI_sigmau","p_value_sigma","loglik")


print(DD)
write.csv(DD,file="Result_ZILN.csv")
