library(MASS)
library(Matrix)
library(LearnBayes)
library(MASS)
library(mvtnorm)
library(tmvtnorm)
t0=Sys.time()
library(nlme)
df<-Orthodont
dis=df$distance
y=dis
n0=4
N=27
n=N*n0

###############先验参数值
gamma=c(1,2,4,3);tau=10;lambda1=2.0001;k1=0.09;lambda2=2.003;k2=0.4
 
X0=array(0,c(4,4))
X0[,c(1,2)]=rep(1,4)
X0[,c(3,4)]=c(8,10,12,14)

X1=array(0,c(4,4))
X1[,1]=rep(1,4)
X1[,c(2,4)]=rep(0,4)
X1[,3]=c(8,10,12,14)
X00=kronecker(rep(1,16),X0)
X11=kronecker(rep(1,11),X1)
X=rbind(X00,X11)
 
Z=kronecker(diag(N),rep(1,n0)) 
l=rep(1,N*n0)
l1=rep(1,n0)
v1=kronecker(diag(N),l1%*%t(l1))
 
 
########后验分布


logpost = function(q)
{#q=propsample[1:2]=c(0,0)
mu=q[1:4]
sigma2_e=q[5]
sigma2_u=q[6]

#Sigma0=kronecker(diag(N),(sigma2_e*diag(n0)+sigma2_u*l1%*%t(l1)))
value=-(N/2*log(n0*sigma2_u+sigma2_e)+N*(n0-1)/2*log(sigma2_e)+1/2*(t(y-X%*%mu)%*%(y-X%*%mu)/sigma2_e)-
1/2*(t(y-X%*%mu)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-X%*%mu)/sigma2_e/(n0+sigma2_e/sigma2_u)))-
t(mu-gamma)%*%solve(diag(c(tau^2,tau^2,tau^2,tau^2)))%*%(mu-gamma)/2-
(lambda1+1)*log(sigma2_e)-(k1)/sigma2_e-(lambda2+1)*log(sigma2_u)-(k2)/sigma2_u
 
return(value)
} 
fit0=laplace(logpost,c(16,1.8,0.65,0.02,2.07,3.62))
mode=fit0$mode
var=fit0$var

####截断正态分布采样
rtnorm=function(n=1,mu,lo=-Inf,up=Inf,sigma ){
qnorm(runif(n,min=pnorm(lo,mean=mu,sd=sigma),
max=pnorm(up,mean=mu,sd=sigma)),
mean=mu,sd=sigma)}

##############随机游走
rwmfun<-function(logpost,start,mm)#mm迭代的链的长度
{
 pb=length(start)
 Mpar=array(0,c(mm,pb))
 b=t(start)

 accept=0
 for(i2 in 1:mm)
 { 
lb=logpost(b)
 ps1=rmnorm(1,mode[1:4],diag(c(0.004,0.004,0.004,0.004)))
ps2=rtnorm(1,2.5,0,Inf,0.2)
ps3=rtnorm(1,2.5,0,Inf,0.5)
  propsample=c(ps1,ps2,ps3)
 lbc=logpost(t(propsample))
 prob=exp(lbc-lb) 
 if(runif(1)<prob)
 {
 b=propsample#接受候选值
 accept=accept+1
 }
 Mpar[i2,]=b#也等于propsample
 }
 accept=accept
 stuff=list(par=Mpar,accept=accept)
 return(stuff)
 }
start=c(10,5,0.65,0.02,2.6,2.7)
mm=30000
rf=rwmfun(logpost,start,mm)
be0=rf$par
sigma_BE=colMeans(be0[20000:30000,5:6])
 t1=Sys.time()
  print(t1-t0)   

  


 
#####################计算LBE
t0_lbe=Sys.time()
library(nlme) 
df<-Orthodont
dis=df$distance
y=dis
n0=4
N=27
n=N*n0

###############先验参数值
gamma=c(1,2,4,3);tau=10;lambda1=2.0001;k1=0.09;lambda2=2.003;k2=0.4
#gamma=c(4,0.003);tau=3;lambda1=2.002;k1=0.09;lambda2=2.05;k2=0.1

X0=array(0,c(4,4))
X0[,c(1,2)]=rep(1,4)
X0[,c(3,4)]=c(8,10,12,14)

X1=array(0,c(4,4))
X1[,1]=rep(1,4)
X1[,c(2,4)]=rep(0,4)
X1[,3]=c(8,10,12,14)
X00=kronecker(rep(1,11),X0)
X11=kronecker(rep(1,16),X1)
X=rbind(X00,X11)
 
Z=kronecker(diag(N),rep(1,n0)) 
l=rep(1,N*n0)
l1=rep(1,n0)
v1=kronecker(diag(N),l1%*%t(l1))


#####################计算LBE
XZ=cbind(X,Z)
PXZ=XZ%*%ginv(t(XZ)%*%(XZ))%*%t(XZ)
a0=(qr(XZ)$rank-qr(X)$rank)/(n-qr(XZ)$rank)
PX=X%*%solve(t(X)%*%(X))%*%t(X)

sigma_e_an=t(y)%*%(diag(n)-PXZ)%*%y/(n-qr(XZ)$rank)
sigma_u_an=(t(y)%*%(PXZ-PX)%*%y-a0*t(y)%*%(diag(n)-PXZ)%*%y)/sum(diag((diag(n)-PX)%*%Z%*%t(Z)))
mean_mu_hat=solve(t(X) %*%(X))%*%t(X) %*%y
theta_hat=c(mean_mu_hat,sigma_e_an,sigma_u_an)
A=cbind(diag(c(0,0)),diag(2))

E_e=k1/(lambda1-1)
Cov_e=k1^2/(lambda1-1)^2/(lambda1-2)
E_u=k2/(lambda2-1)
Cov_u=k2^2/(lambda2-1)^2/(lambda2-2)
Cov_mu=diag(c(tau^2,tau^2,tau^2,tau^2))
E_mu=gamma
 
Cov1=array(c(Cov_e,0,0,Cov_u),c(2,2))
E=c(gamma,E_e,E_u)
b0=n-qr(XZ)$rank
a=a0
h=sum(diag((diag(n)-PX)%*%Z%*%t(Z)))/a
a22=2*(E_e^2+Cov_e)/b0
a23=a32=-2*(E_e^2+Cov_e)/b0/a
a33=2*(1+h/b0)*(E_e^2+Cov_e)/a^2/h+4*E_e*E_u/a/h+2*(E_u^2+Cov_u)/h
#Cov_hat=Cov+array(c(solve(t(X)%*%(X))%*%t(X)%*%Sigma%*%X%*%solve(t(X)%*%(X)),
#0,0,0,a22,a23,0,a32,a33),c(3,3))
#sigma_LBE1=A%*%(Cov%*%solve(Cov_hat)%*%theta_hat+(diag(3)-Cov%*%solve(Cov_hat))%*%E)
 
sigma_LBE2=Cov1%*%solve(Cov1+array(c(a22,a23,a32,a33),c(2,2)))%*%theta_hat[5:6]+(diag(2)-Cov1%*%solve(Cov1+array(c(a22,a23,a32,a33),c(2,2))))%*%E[5:6]
t1_lbe=Sys.time()
print(t1_lbe-t0_lbe)

dis1=sqrt(sum((sigma_LBE2-sigma_BE)^2))
dis2=sqrt(sum((theta_hat[5:6]-sigma_BE)^2))
anova=theta_hat[5:6]
lbe=sigma_LBE2
be=sigma_BE
anova
lbe
be
dis1
dis2

