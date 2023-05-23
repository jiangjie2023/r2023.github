

library(MASS)
library(Matrix)
library(LearnBayes)
library(MASS)
library(mvtnorm)
library(tmvtnorm)

##############################################beta和sigma^2相关时
n0=5
NN=seq(10,100,by=10)
m=1
sigma_TK0=sigma_LBE10=sigma_LBE20=sigma_L0=sigma_BE0<-array(0,c(m,2,length(NN)))
dis1=dis2=dis3=dis4=dis5=dis6=dis7=dis8=dis9<-array(0,c(m,length(NN)))
for(r in 1:length(NN))
{N=NN[r]
n=N*n0
 for (j in 1:m)
{
###############先验参数值
#c=1;d=5;lambda1=11;k1=80;lambda2=73;k2=600
c=-1;d=5;lambda1=10;k1=82;lambda2=70;k2=610
X=rep(1,n)
Z=kronecker(diag(N),rep(1,n0))

l=rep(1,N*n0)
l1=rep(1,n0)
v1=kronecker(diag(N),l1%*%t(l1))

set.seed(1)
mean_mu=runif(1,c,d)
sig2_e=1/rgamma(1,lambda1,k1+mean_mu)
sig2_u=1/rgamma(1,lambda2,k2+mean_mu)

meanmu=0.6
sigma_2=c(7,9)
Sigma=kronecker(diag(N),(sigma_2[1]*diag(n0)+sigma_2[2]*l1%*%t(l1)))

y=rmnorm(1,meanmu*l,Sigma)
y=t(y)
####################################lbe


XZ=cbind(X,Z)
PXZ=XZ%*%ginv(t(XZ)%*%(XZ))%*%t(XZ)
a0=(qr(XZ)$rank-qr(X)$rank)/(n-qr(XZ)$rank)
PX=X%*%solve(t(X)%*%(X))%*%t(X)

sigma_e_an=t(y)%*%(diag(n)-PXZ)%*%y/(n-qr(XZ)$rank)
sigma_u_an=(t(y)%*%(PXZ-PX)%*%y-a0*t(y)%*%(diag(n)-PXZ)%*%y)/sum(diag((diag(n)-PX)%*%Z%*%t(Z)))
mean_mu_hat=solve(t(X)%*%solve(Sigma)%*%(X))%*%t(X)%*%solve(Sigma)%*%y
theta_hat=c(mean_mu_hat,sigma_e_an,sigma_u_an)
A=cbind(rep(0,2),diag(2))
E_e=k1/(lambda1-1)
Cov_e=k1^2/(lambda1-1)^2/(lambda1-2)
E_u=k2/(lambda2-1)
Cov_u=k2^2/(lambda2-1)^2/(lambda2-2)
Cov_mu=(c-d)^2/12
E_mu=(c+d)/2
Cov12=k1/(lambda1-1)*E_mu+(E_mu^2+Cov_mu)/(lambda1-1)-E_mu*E_e
Cov13=k2/(lambda2-1)*E_mu+(E_mu^2+Cov_mu)/(lambda2-1)-E_mu*E_u
Cov21=Cov12
Cov31=Cov13
Cov=array(c(Cov_mu,Cov12,Cov13,Cov21,Cov_e,0,Cov31,0,Cov_u),c(3,3))
Cov1=array(c(Cov_e,0,0,Cov_u),c(2,2))
E=c((E_u-E_e)/2,E_e,E_u)
b0=n-qr(XZ)$rank
a=a0
h=sum(diag((diag(n)-PX)%*%Z%*%t(Z)))/a
a22=2*(E_e^2+Cov_e)/b0
a23=a32=-2*(E_e^2+Cov_e)/b0/a
a33=2*(1+h/b0)*(E_e^2+Cov_e)/a^2/h+4*E_e*E_u/a/h+2*(E_u^2+Cov_u)/h
Cov_hat=Cov+array(c(as.numeric(solve(t(X)%*%(X))%*%t(X)%*%Sigma%*%X%*%solve(t(X)%*%(X))),
0,0,0,a22,a23,0,a32,a33),c(3,3))
sigma_LBE1=A%*%(Cov%*%solve(Cov_hat)%*%theta_hat+(diag(3)-Cov%*%solve(Cov_hat))%*%E)
sigma_LBE10[j,,r]=t(sigma_LBE1)
sigma_LBE2=Cov1%*%solve(Cov1+array(c(a22,a23,a32,a33),c(2,2)))%*%theta_hat[2:3]+(diag(2)-Cov1%*%solve(Cov1+array(c(a22,a23,a32,a33),c(2,2))))%*%E[2:3]
sigma_LBE20[j,,r]=t(sigma_LBE2)
#############lindley近似


L111= -2*t(l)%*%v1%*%(l*mean_mu-y)/(n0*sig2_u+sig2_e)^3
L102=-2*n0*t(l)%*%v1%*%(l*mean_mu-y)/(n0*sig2_u+sig2_e)^3
L021=-n/(n0*sig2_u+sig2_e)^3+3*t(y-l*mean_mu)%*%v1%*%(y-l*mean_mu)/(n0*sig2_u+sig2_e)^4
L012=n0*L021
L201=t(l)%*%v1%*%l/(n0*sig2_u+sig2_e)^2
L300=0
L030=-N/(n0*sig2_u+sig2_e)^3-N*(n0-1)/sig2_e^3+3*t(y-l*mean_mu)%*%(y-l*mean_mu)/sig2_e^4-
3*t(y-l*mean_mu)%*%v1%*%(y-l*mean_mu)*(sig2_e^-4*sig2_u*(n0*sig2_u+sig2_e)^-1+sig2_e^-3*sig2_u*(n0*sig2_u+sig2_e)^-2+
sig2_e^-2*sig2_u*(n0*sig2_u+sig2_e)^-3+sig2_e^-1*sig2_u*(n0*sig2_u+sig2_e)^-4)
L003=n0*L012
L120=-2*(n*mean_mu-t(l)%*%y)*sig2_e^-3+2*t(l)%*%v1%*%(l*mean_mu-y)*(sig2_e^-3*sig2_u*(n0*sig2_u+sig2_e)^-1+
sig2_e^-2*sig2_u*(n0*sig2_u+sig2_e)^-2+sig2_e^-1*sig2_u*(n0*sig2_u+sig2_e)^-3)
L210=n/sig2_e^2-t(l)%*%v1%*%l*(sig2_e^-2*sig2_u*(n0*sig2_u+sig2_e)^-1+sig2_e^-1*sig2_u*(n0*sig2_u+sig2_e)^-2)
##信息阵
I11=n/sig2_e+t(l)%*%v1%*%l*sig2_u/sig2_e/(n0*sig2_u+sig2_e)
I22=-N/2/(n0*sig2_u+sig2_e)^2-N*(n0-1)/2/sig2_e^2+N*n0/sig2_e^2-n/sig2_e^2*sig2_u*(n0*sig2_u+sig2_e)^-1-
n/sig2_e*sig2_u*(n0*sig2_u+sig2_e)^-2
I33=N*n0^2/2/n/sig2_e^2*sig2_u*(n0*sig2_u+sig2_e)^2
I12=0
I13=0
I23=I33/n0
I21=I12
I31=I13
I32=I23

I=array(c(I11,I12,I13,I21,I22,I23,I31,I32,I33),c(3,3))
sI=solve(I)

sig2_e_ML=t(y)%*%(diag(n)-PXZ)%*%y/(n-qr(XZ)$rank)

sig2_u_ML=(t(y)%*%(PXZ-PX)%*%y-a0*t(y)%*%(diag(n)-PXZ)%*%y)/sum(diag((diag(n)-PX)%*%Z%*%t(Z)))

rho1=-1/sig2_e-1/sig2_u
rho2=-(lambda1+1)/sig2_e+(k1+mean_mu)/sig2_e^2
rho3=(lambda2+1)/sig2_u+(k2+mean_mu)/sig2_u^2


#muL=beta1mle+1/2*(l111*(2*sI[1,1]*sI[2,3]+4*sI[1,2]*sI[1,3])+l102*(sI[1,1]*sI[3,3]+2*sI[1,3]^2)+
#l021*(sI[1,3]*sI[2,2]+2*sI[1,2]*sI[2,3])+l012*(sI[1,2]*sI[3,3]+2*sI[1,3]*sI[2,3])+
#3*L201*sI[1,1]*sI[1,3])+rho1*sI[1,1]+rho2*sI[2,1]+rho3*sI[3,1]

sig2_e_L=sig2_e_ML+1/2*(L111*(2*sI[1,3]*sI[2,2]+4*sI[1,2]*sI[2,3])+L012*(sI[2,2]*sI[3,3]+2*sI[2,3]^2)+
L102*(sI[1,2]*sI[3,3]+2*sI[1,3]*sI[2,3])+L201*(sI[1,1]*sI[2,3]+2*sI[1,2]*sI[1,3])+
3*L021*sI[2,2]*sI[3,2])+rho1*sI[1,2]+rho2*sI[2,2]+rho3*sI[3,2]

sig2_u_L=sig2_u_ML+1/2*(L111*(2*sI[1,2]*sI[3,3]+4*sI[1,3]*sI[2,3])+L201*(sI[1,1]*sI[3,3]+2*sI[1,3]^2)+
L021*(sI[2,2]*sI[3,3]+2*sI[2,3]^2)+3*L012*sI[2,3]*sI[3,3]+
3*L102*sI[1,3]*sI[3,3])+rho1*sI[1,3]+rho2*sI[2,3]+rho3*sI[3,3]

sigma_L=c(sig2_e_L,sig2_u_L)
sigma_L0[j,,r]=t(sigma_L)

#############################BE


logpost = function(theta)
{
mu=theta[1]
sigma2_e=theta[2]
sigma2_u=theta[3]

value=-(N/2*log(n0*sigma2_u+sigma2_e)+N*(n0-1)/2*log(sigma2_e)+1/2*(t(y-l*mu)%*%(y-l*mu)/sigma2_e)-
1/2*(t(y-l*mu)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu)/sigma2_e/(n0+sigma2_e/sigma2_u))+
(lambda1+1)*log(sigma2_e)+(k1+mu)/sigma2_e+(lambda2+1)*log(sigma2_u)+(k2+mu)/sigma2_u)
  return(value)
}
fit0=laplace(logpost,c(-1.5,6,10))

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
 lb=logpost(b)
 
 accept=0
 for(i2 in 1:mm)
 {
 ps1=rnorm(1,mode[1],sqrt(var[1,1]))
ps2=rtnorm(1,mode[2],0,Inf,sqrt(var[2,2]))
ps3=rtnorm(1,mode[3],0,Inf,sqrt(var[3,3]))
  propsample=c(ps1,ps2,ps3)
 lbc=logpost(t(propsample))
 prob=exp(lbc-lb) 
 if(runif(1)<prob)
 {lb=lbc#这一步什么意思?,因为此时接受propsample所以b更新为propsample，则有lb=lbc.
 b=propsample#接受候选值
 accept=accept+1
 }
 Mpar[i2,]=b#也等于propsample
 }
 accept=accept/mm
 stuff=list(par=Mpar,accept=accept)
 return(stuff)
 }
start=c(-6,6,11)
mm=10000
rf=rwmfun(logpost,start,mm)
be0=rf$par

sigma_BE=colMeans(be0[5000:10000,])[2:3]
sigma_BE

sigma_BE0[j,,r]=t(sigma_BE)

######################################TK



f1 = function(theta)
{
mu=theta[1]
sigma2_e=theta[2]
sigma2_u=theta[3]

value=-(N/2*log(n0*sigma2_u+sigma2_e)+N*(n0-1)/2*log(sigma2_e)+1/2*(t(y-l*mu)%*%(y-l*mu)/sigma2_e)-
1/2*(t(y-l*mu)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu)/sigma2_e/(n0+sigma2_e/sigma2_u))+
(lambda1+1)*log(sigma2_e)+(k1+mu)/sigma2_e+(lambda2+1)*log(sigma2_u)+(k2+mu)/sigma2_u)/n
  return(value)
}
fit1=laplace(f1,c(-1.5,6,10))



f3 = function(theta)
{
mu=theta[1]
sigma2_e=theta[2]
sigma2_u=theta[3]

value=-(N/2*log(n0*sigma2_u+sigma2_e)+N*(n0-1)/2*log(sigma2_e)+1/2*(t(y-l*mu)%*%(y-l*mu)/sigma2_e)-
1/2*(t(y-l*mu)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu)/sigma2_e/(n0+sigma2_e/sigma2_u))+
(lambda1+1)*log(sigma2_e)+(k1+mu)/sigma2_e+(lambda2+1)*log(sigma2_u)+(k2+mu)/sigma2_u)/n+log(sigma2_e)/n
  return(value)
}
fit3=laplace(f3,c(-1.5,6,10))





f4 = function(theta)
{
mu=theta[1]
sigma2_e=theta[2]
sigma2_u=theta[3]

value=-(N/2*log(n0*sigma2_u+sigma2_e)+N*(n0-1)/2*log(sigma2_e)+1/2*(t(y-l*mu)%*%(y-l*mu)/sigma2_e)-
1/2*(t(y-l*mu)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu)/sigma2_e/(n0+sigma2_e/sigma2_u))+
(lambda1+1)*log(sigma2_e)+(k1+mu)/sigma2_e+(lambda2+1)*log(sigma2_u)+(k2+mu)/sigma2_u)/n+log(sigma2_u)/n
  return(value)
}
fit4=laplace(f4,c(-1.5,6,10))


e_TK=(det(fit3$var)/det(fit1$var))^(1/2)*exp(n*(f3(theta=fit3$mode )-f1(theta=fit1$mode)))

u_TK=(det(fit4$var)/det(fit1$var))^(1/2)*exp(n*(f4(theta=fit4$mode )-f1(theta=fit1$mode )))


sigma_TK=c(e_TK,u_TK)
sigma_TK0[j,,r]=t(sigma_TK)


dis1[j,r]=sqrt(sum((sigma_LBE1-sigma_BE)^2))
dis2[j,r]=sqrt(sum((sigma_L-sigma_BE)^2))
dis3[j,r]=sqrt(sum((sigma_LBE2-sigma_BE)^2))
dis4[j,r]=sqrt(sum((sigma_TK-sigma_BE)^2))

dis5[j,r]=sqrt(sum((sigma_LBE1-sigma_2)^2))
dis6[j,r]=sqrt(sum((sigma_L-sigma_2)^2))
dis7[j,r]=sqrt(sum((sigma_BE-sigma_2)^2))
dis8[j,r]=sqrt(sum((sigma_TK-sigma_2)^2))
dis9[j,r]=sqrt(sum((sigma_LBE2-sigma_2)^2))

}
}
d11=colMeans(dis1)
d22=colMeans(dis2)
d33=colMeans(dis3)
d44=colMeans(dis4)
d11
d22
d33
d44



plot(NN,d11,xlab=expression(N),ylab="Distances",ylim=c(0,3),col="red",type="o",lty=1)
lines(NN,d22,col="black",type="o",lty=2)
lines(NN,d33,col="blue",type="o",lty=3)
lines(NN,d44,col="green",type="o",lty=4)
legend(60,3,c(expression(paste("||",tilde(sigma)^"2LBE","-", hat(sigma)^"2UBE","||")),
expression(paste("||",hat(sigma)^"2Lindley","-", hat(sigma)^"2UBE","||")),
expression(paste("||",tilde(sigma)^"2LBE"[2],"-", hat(sigma)^"2UBE","||")),
expression(paste("||",hat(sigma)^"2TK","-", hat(sigma)^"2UBE","||"))),
col=c("red","black","blue","green"),
lty=c(1,2,3,4))



d55=colMeans(dis5)
d66=colMeans(dis6)
d77=colMeans(dis7)
d88=colMeans(dis8)
d99=colMeans(dis9)
d55
d66
d77
d88
d99

plot(NN,d55,xlab=expression(N),ylab="Distances",ylim=c(0,3.5),col="red",type="o",lty=1)
lines(NN,d66,col="black",type="o",lty=2)
lines(NN,d77,col="green",type="o",lty=3)
lines(NN,d88,col="blue",type="o",lty=4)
lines(NN,d99,col="orange",type="o",lty=5)
legend(70,3.5,c(expression(paste("||",tilde(sigma)^"2LBE","-", sigma^2,"||")),
expression(paste("||",hat(sigma)^"2Lindley","-", sigma^2,"||")),
expression(paste("||",hat(sigma)^"2UBE","-", sigma^2,"||")),
expression(paste("||",hat(sigma)^"2TK","-", sigma^2,"||")),
expression(paste("||",tilde(sigma)^"2LBE"[2],"-", sigma^2,"||"))),
col=c("red","black","green","blue","orange"),
lty=c(1,2,3,4,5))



 


