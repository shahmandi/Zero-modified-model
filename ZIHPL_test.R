
N<-10000
aa<-runif(۳,.1,.4)
NC<-1000

V<-rep(888,length(aa))
ddhpl_re<-function(B1,alpha,y){
  
  for(f in 1:length(aa)){
    V[f]<-sum((rep(exp(B1*aa[f]),NC)+(0:(NC-1)))^(-exp(alpha)))
  }
  
  dhpl1a<-function(B1,alpha,y) ((exp(B1*aa)+y)^(-exp(alpha)))/V
  
  
  return(dhpl1a(B1,alpha,y))
}


dzdhpl1<-function(p,B1,alpha,y) 
  ifelse(y==0,plink(p)+(plink(1-p))*ddhpl_re(B1,alpha,y),(plink(1-p))*ddhpl_re(B1,alpha,y))



B1<-2
alpha1<-3
p1<-.1
vhpl<-dzdhpl1(p1,B1*aa,alpha1,0:nn2)



nn2<-1750000

nn<-5000000
qvhpl1<-round(vhpl*nn,0)

xx<-rep(0:(nn2),qvhpl1)
x<-sample(xx,N,replace=TRUE)


alogit<-function(t) exp(t)/(1+exp((t)))
plink1<-function(t) cloglog(t)
plink22<-function(t) alogit(t)
plink3<-function(t) pnorm(t)  ## this is the probit link
plink4<-function(t) identity(t)
plink5<-function(t) cauchit(t)
plink6<-function(t) negative_loglog(t)
plink7<-function(t) loglog(t)
plink2<-function(t) (t)
plink=plink2
mulink<-function(t) exp(t)
slink<-function(t) exp(t)


#table(x)


funs= function(B,alpha,y) (ddhpl(B,alpha,y))
funsb = function(x, y) sum(log(funs(x[1], x[2],y))) 

rm(hpl) 
hpl<- optim(c(2,2),funsb,control=list(fnscale=-1),y=x,hessian=TRUE)  
hpl


##########
funs= function(B,alpha,y) (ddhpl(exp(B),exp(alpha),y))
funsb = function(x, y) sum(log(funs(x[1], x[2],y))) 

rm(hpl) 
hpl<- optim(c(1,1),funsb,control=list(fnscale=-1),y=x,hessian=TRUE)  
hpl

exp(0.6965436 )

exp(1.5243753)
exp(0.6564715)
#Fitting the zero-modified version of models


#####################################################
dzdhpl0<-function(p,B,alpha,y)
 ifelse(y==0, p+(1-p)*ddhpl(B,alpha,y),(1-p)*ddhpl(B,alpha,y)) 

fun0 = function(p,B,alpha,y) dzdhpl0(p,B,alpha,y)
fun0b = function(x, y) sum(log(fun0(x[1], x[2],x[3],y))) 
rm(zmhpl0) 
zmhpl0 <- optim(c(0,hpl$par[1],hpl$par[2]),fun0b,control=list(fnscale=-1),y=x,hessian=TRUE) 
zmhpl0


#############################NO COVARIATE#########################

dzdhpl0<-function(p,B,alpha,y)
  ifelse(y==0, plink(p)+(plink(1-p))*ddhpl(exp(B*aa),exp(alpha),y),(plink(1-p))*ddhpl(exp(B),exp(alpha),y)) 

fun0 = function(p,B,alpha,y) dzdhpl0(p,B,alpha,y)
fun0b = function(x, y) sum(log(fun0(x[1], x[2],x[3],y))) 
rm(zmhpl0) 
zmhpl0 <- optim(c(0,hpl$par[1],hpl$par[2]),fun0b,control=list(fnscale=-1),y=x,hessian=TRUE) 

zmhpl0
exp(0.7789525)
1-exp(-exp( 0.08616427))

##########################One parameter modelled###########################
NC<-1000
V<-rep(888,length(aa))
ddhpl<-function(B0,B1,alpha,y){
  
  for(f in 1:length(aa)){
    V[f]<-sum((rep(exp(B0+B1*aa[f]),NC)+(0:(NC-1)))^(-exp(alpha)))
  }
  
  dhpl1a<-function(B0,B1,alpha,y) ((exp(B0+B1*aa)+y)^(-exp(alpha)))/V
  
  
  return(dhpl1a(B0,B1,alpha,y))
}


dzdhpl1<-function(p,B0,B1,alpha,y) 
  ifelse(y==0,plink(p)+(plink(1-p))*ddhpl(B0,B1,alpha,y),(plink(1-p))*ddhpl(B0,B1,alpha,y))

fun1 = function(p,B0,B1,alpha,y) dzdhpl1(p,B0,B1,alpha,y)

fun1b = function(x, y) sum(log(fun1(x[1], x[2],x[3],x[4],y)))

rm(zmhpl1) 
zmhpl1<- optim(c(zmhpl0$par[1],zmhpl0$par[2],0,zmhpl0$par[3]),fun1b,control=list(fnscale=-1,maxit=100000),y=x,hessian=TRUE) 
zmhpl1

exp(  1.65980500 +0.01627187 *aa)
1-exp(-exp( 0.068828020))
exp(0.784853162)
#####################Two parameters modelled###################

NC<-1000
V<-rep(888,length(aa))
ddhpl_2<-function(B0,B1,alpha,y){
  
  for(f in 1:length(aa)){
    V[f]<-sum((rep(exp(B0+B1*aa[f]),NC)+(0:(NC-1)))^(-exp(alpha)))
  }
  
  dhpl1a<-function(B0,B1,alpha,y) ((exp(B0+B1*aa)+y)^(-exp(alpha)))/V
  
  
  return(dhpl1a(B0,B1,alpha,y))
}


dzdhpl2<-function(p0,p1,B0,B1,alpha,y) 
  ifelse(y==0, plink(p0+aa*p1)+plink(1-(p0+aa*p1))*ddhpl_2(B0,B1,alpha,y),plink(1-(p0+aa*p1))*ddhpl_2(B0,B1,alpha,y))


fun2 = function(p0,p1,B0,B1,alpha,y) dzdhpl2(p0,p1,B0,B1,alpha,y)
fun2b = function(x, y) sum(log(fun2(x[1], x[2],x[3],x[4],x[5],y)))


rm(zmhpl2) 
zmhpl2<-optim(c(zmhpl1$par[1],0,zmhpl1$par[2],zmhpl1$par[3],zmhpl1$par[4]),fun2b,control=list(fnscale=-1,maxit=100000),y=x,hessian=TRUE) 
zmhpl2

#rm(zmhpl2) 
#zmhpl2<-optim(c(zmhpl1$par[1],0,.6,0,zmhpl1$par[4]),fun2b,control=list(fnscale=-1,maxit=100000),y=x,hessian=TRUE) 
#zmhpl2

exp( 0.9755627587 -0.0083421404*aa)[1]
exp(0.77945957)
log(.1/.9)


#----------------------------------------------------------------------------------------
#####################Three parameters modelled###################

NC<-1000
V<-rep(888,length(aa))
ddhpl_3<-function(B0,B1,alpha0,alpha1,y){
  
  for(f in 1:length(aa)){
    V[f]<-sum((rep(exp(B0+B1*aa[f]),NC)+(0:(NC-1)))^(-exp(alpha0+alpha1*aa[f])))
  }
  
  dhpl1a<-function(B0,B1,alpha0,alpha1,y) ((exp(B0+B1*aa)+y)^(-exp(alpha0+alpha1*aa)))/V
  
  
  return(dhpl1a(B0,B1,alpha0,alpha1,y))
}

dzdhpl3<-function(p0,p1,B0,B1,alpha0,alpha1,y) 
  ifelse(y==0, plink(p0+aa*p1)+plink(1-(p0+aa*p1))*ddhpl_3(B0,B1,alpha0,alpha1,y),plink(1-(p0+aa*p1))*ddhpl_3(B0,B1,alpha0,alpha1,y))


fun3 = function(p0,p1,B0,B1,alpha0,alpha1,y) dzdhpl3(p0,p1,B0,B1,alpha0,alpha1,y)
fun3b = function(x, y) sum(log(fun3(x[1], x[2],x[3],x[4],x[5],x[6],y)))

rm(zmhpl3) 
zmhpl3= optim(c(zmhpl2$par[1],zmhpl2$par[2],zmhpl2$par[3],zmhpl2$par[4],zmhpl2$par[5],0),fun3b,control=list(fnscale=-1,maxit=100000),y=x,hessian=TRUE) 
zmhpl3

exp(0.791322506+ 0.068153304 *aa)[1]
exp(0.595464603+ 0.019654731*aa)[1]
#0.13213821 -0.01373895  *aa
