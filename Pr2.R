library(MASS)
library(Matrix)
library(LearnBayes)
library(MASS)
library(mvtnorm)
library(tmvtnorm)

##############################################计算的是与真实值距离的版本
t0=proc.time()###记录程序运行时间
n0=5
NN=seq(10,100,by=10)
m=1
sigma_TK0=sigma_LBE10=sigma_LBE20=sigma_L0=sigma_BE0<-array(0,c(m,2,length(NN)))
dis1=dis2=dis3=dis4=dis5=dis6=dis7=dis8<-array(0,c(m,length(NN)))
for(r in 1:length(NN))
{N=NN[r]
n=N*n0
 for (j in 1:m)
{
###############先验参数值
#gamma=-1;tau=3;lambda1=11;k1=80;lambda2=73;k2=600

gamma=-1;tau=4;lambda1=10;k1=82;lambda2=70;k2=610 

X=rep(1,n)
Z=kronecker(diag(N),rep(1,n0))

l=rep(1,N*n0)
l1=rep(1,n0)
v1=kronecker(diag(N),l1%*%t(l1))

set.seed(1)
mean_mu=rnorm(1,gamma,tau^2)
sig2_e=1/rgamma(1,lambda1,k1)
sig2_u=1/rgamma(1,lambda2,k2)
meanmu=0.6
sigma_2=c(7,9)
Sigma=kronecker(diag(N),(sigma_2[1]*diag(n0)+sigma_2[2]*l1%*%t(l1)))
y=rmnorm(1,meanmu*l,Sigma)
y=t(y)
t0_0=proc.time()
####################################lbe

t1=proc.time()
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
Cov=array(c(tau^2,0,0,0,Cov_e,0,0,0,Cov_u),c(3,3))
Cov1=array(c(Cov_e,0,0,Cov_u),c(2,2))
E=c(gamma,E_e,E_u)
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
t1_1=proc.time()

HMC = function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = matrix(rnorm(length(q),0,1), ncol = 1)  
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  for (i in 1:L)
  {
    # Make a full step for the position
    q = q + epsilon * p
    # Make a full step for the momentum, except at end of trajectory
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  # Make a half step for momentum at the end.
  p = p - epsilon * grad_U(q) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (runif(1) < exp(current_U - proposed_U + current_K - proposed_K))
  {
    return (q)  # accept
  }
  else
  {
    return (current_q)  # reject
  }
}

###############################################################
#simulate gaussian distribution using HMC
###############################################################
#potential energy
U_P = function(q)
{
mu=q[1]
sigma2_e=q[2]
sigma2_u=q[3]

value=N/2*log(n0*sigma2_u+sigma2_e)+N*(n0-1)/2*log(sigma2_e)+1/2*(t(y-l*mu)%*%(y-l*mu)/sigma2_e)-
1/2*(t(y-l*mu)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu)/sigma2_e/(n0+sigma2_e/sigma2_u))+
(mu-gamma)^2/2/tau^2+(lambda1+1)*log(sigma2_e)+k1/lambda1+(lambda2+1)*log(sigma2_u)+k2/sigma2_u
  return(value)
}

#gradient
dU = function(q)
{
mu=q[1]
sigma2_e=q[2]
sigma2_u=q[3]
K1=(as.numeric(t(l)%*%l)*mu-t(y)%*%l)/sigma2_e-(mu*t(l)-t(y))%*%(kronecker(diag(N),l1%*%t(l1)))%*%l/sigma2_e/(n0+sigma2_e/sigma2_u)-
(mu-gamma)/tau^2
K2=N/2/(n0*sigma2_u+sigma2_e)+N*(n0-1)/2/sigma2_e-t(y-l*mu)%*%(y-l*mu)/2/sigma2_e^2+
t(y-l*mu)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu)/2*(sigma2_e^-2*sigma2_u*(n0*sigma2_u+sigma2_e)^-1+sigma2_e^-1*sigma2_u*(n0*sigma2_u+sigma2_e)^-2)+
(lambda1+1)/sigma2_e-k1/sigma2_e^2
K3=N*n0/2/(n0*sigma2_u+sigma2_e)-t(y-l*mu)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu)*(n0*sigma2_u+sigma2_e)^-2/2+
(lambda2+1)/sigma2_u-k2/sigma2_u^2

K=array(c(K1,K2,K3),c(3,1))
  return(K)
}
M= 2000
q_HMC = matrix(NA, nrow = 3, ncol = M)
q_init = matrix(c(-1.5, 6,10), ncol = 1)

for (i in 1:M) 
{
  q_HMC[,i] = HMC(U = U_P, grad_U = dU, epsilon = 0.025, 
                     L = 20, current_q = q_init)
  q_init = q_HMC[,i]
}
#q_HMC[,1]
sigma_BE=rowMeans(q_HMC[2:3,])
sigma_BE0[j,,r]=t(sigma_BE)

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
I12=-(n*mean_mu-t(l)%*%y)/sig2_e^2
I13=0
I23=I33/n0
I21=I12
I31=I13
I32=I23

I=array(c(I11,I12,I13,I21,I22,I23,I31,I32,I33),c(3,3))
sI=solve(I)

sig2_e_ML=t(y)%*%(diag(n)-PXZ)%*%y/(n-qr(XZ)$rank)

sig2_u_ML=(t(y)%*%(PXZ-PX)%*%y-a0*t(y)%*%(diag(n)-PXZ)%*%y)/sum(diag((diag(n)-PX)%*%Z%*%t(Z)))

rho1=-(mean_mu-gamma)/tau^2
rho2=-(lambda1+1)/sig2_e+k1/sig2_e^2
rho3=-(lambda2+1)/sig2_u+k2/sig2_u^2


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
 

###############################################TK
 

f1 = function(theta)
{
mu_TK=theta[1]
sigma2_e_TK=theta[2]
sigma2_u_TK=theta[3]

value=-(N/2*log(n0*sigma2_u_TK+sigma2_e_TK)+N*(n0-1)/2*log(sigma2_e_TK)+1/2*(t(y-l*mu_TK)%*%(y-l*mu_TK)/sigma2_e_TK)-
1/2*(t(y-l*mu_TK)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu_TK)/sigma2_e_TK/(n0+sigma2_e_TK/sigma2_u_TK))+
(mu_TK-gamma)^2/2/tau^2+(lambda1+1)*log(sigma2_e_TK)+k1/sigma2_e_TK+(lambda2+1)*log(sigma2_u_TK)+k2/sigma2_u_TK)/n
  return(value)
}
fit1=laplace(f1,c(8,16,33))


f3 = function(theta)
{
mu_TK=theta[1]
sigma2_e_TK=theta[2]
sigma2_u_TK=theta[3]

value=-(N/2*log(n0*sigma2_u_TK+sigma2_e_TK)+N*(n0-1)/2*log(sigma2_e_TK)+1/2*(t(y-l*mu_TK)%*%(y-l*mu_TK)/sigma2_e_TK)-
1/2*(t(y-l*mu_TK)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu_TK)/sigma2_e_TK/(n0+sigma2_e_TK/sigma2_u_TK))+
(mu_TK-gamma)^2/2/tau^2+(lambda1+1)*log(sigma2_e_TK)+k1/sigma2_e_TK+(lambda2+1)*log(sigma2_u_TK)+k2/sigma2_u_TK)/n+log(sigma2_e_TK)/n
  return(value)

}
fit3=laplace(f3,c(-1.5,6,10))


f4 = function(theta)
{
mu_TK=theta[1]
sigma2_e_TK=theta[2]
sigma2_u_TK=theta[3]

value=-(N/2*log(n0*sigma2_u_TK+sigma2_e_TK)+N*(n0-1)/2*log(sigma2_e_TK)+1/2*(t(y-l*mu_TK)%*%(y-l*mu_TK)/sigma2_e_TK)-
1/2*(t(y-l*mu_TK)%*%(kronecker(diag(N),l1%*%t(l1)))%*%(y-l*mu_TK)/sigma2_e_TK/(n0+sigma2_e_TK/sigma2_u_TK))+
(mu_TK-gamma)^2/2/tau^2+(lambda1+1)*log(sigma2_e_TK)+k1/sigma2_e_TK+(lambda2+1)*log(sigma2_u_TK)+k2/sigma2_u_TK+log(sigma2_e_TK))/n+
log(sigma2_u_TK)/n
  return(value)
}
fit4=laplace(f4,c(-1.5,6,10))

 
e_TK=(det(fit3$var)/det(fit1$var))^(1/2)*exp(n*(f3(theta=fit3$mode )-f1(theta=fit1$mode )))

u_TK=(det(fit4$var)/det(fit1$var))^(1/2)*exp(n*(f4(theta=fit4$mode )-f1(theta=fit1$mode )))


sigma_TK=c(e_TK,u_TK)
sigma_TK0[j,,r]=t(sigma_TK)
 

dis1[j,r]=sqrt(sum((sigma_LBE1-sigma_BE)^2))
dis2[j,r]=sqrt(sum((sigma_L-sigma_BE)^2))
dis3[j,r]=sqrt(sum((sigma_TK-sigma_BE)^2))

dis5[j,r]=sqrt(sum((sigma_LBE1-sigma_2)^2))
dis6[j,r]=sqrt(sum((sigma_L-sigma_2)^2))
dis7[j,r]=sqrt(sum((sigma_BE-sigma_2)^2))
dis8[j,r]=sqrt(sum((sigma_TK-sigma_2)^2))

}
}
d11=colMeans(dis1)
d22=colMeans(dis2)
d33=colMeans(dis3)

d11
d22
d33


plot(NN,d11,xlab=expression(N),ylab="Distances",ylim=c(0,8),col="red",type="o",lty=1)
lines(NN,d22,col="black",type="o",lty=2)
lines(NN,d33,col="blue",type="o",lty=3)
legend(60,5.5,c(expression(paste("||",tilde(sigma)^"2LBE","-", hat(sigma)^"2UBE","||")),
expression(paste("||",hat(sigma)^"2Lindley","-", hat(sigma)^"2UBE","||")),
expression(paste("||",hat(sigma)^"2TK","-", hat(sigma)^"2UBE","||"))),
col=c("red","black","blue"),
lty=c(1,2,3))



d55=colMeans(dis5)
d66=colMeans(dis6)
d77=colMeans(dis7)
d88=colMeans(dis8)
d55
d66
d77
d88


plot(NN,d55,xlab=expression(N),ylab="Distances",ylim=c(0,8),col="red",type="o",lty=1)
lines(NN,d66,col="black",type="o",lty=2)
lines(NN,d77,col="green",type="o",lty=3)
lines(NN,d88,col="blue",type="o",lty=4)
legend(50,6.5,c(expression(paste("||",tilde(sigma)^"2LBE","-", sigma^2,"||")),
expression(paste("||",hat(sigma)^"2Lindley","-", sigma^2,"||")),
expression(paste("||",hat(sigma)^"2UBE","-", sigma^2,"||")),
expression(paste("||",hat(sigma)^"2TK","-", sigma^2,"||"))),
col=c("red","black","green","blue"),
lty=c(1,2,3,4))
 













 

