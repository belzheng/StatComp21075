## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
corr=function(n)
{
u=runif(n)
results1=(mean(u*sqrt(1-u^2))-mean(u)*mean(sqrt(1-u^2)))/(sd(u)*sd(sqrt(1-u^2)))
results1
}
corr(1000000)

## -----------------------------------------------------------------------------
rubber<-data.frame(x1=c(65,70,70,69,66,67,68,72,66,68),x2=c(45,45,48,46,50,46,47,43,47,48),x3=c(27.6,30.7,31.8,32.6,31.0,31.3,37.0,33.6,33.1,34.2))
plot(rubber)

## -----------------------------------------------------------------------------
##建立数据框
df<-data.frame(
   Name=c("张三","李四","王五","赵六","丁一"),
   Sex=c("女","男","女","男","女"),
   Age=c(14,15,16,14,15),
   Height=c(156,165,157,162,159),
   Weight=c(42.0,49.0,41.5,52.0,45.5)
);df

write.csv(df,file="foo.csv")
##将数据框转为表格

## -----------------------------------------------------------------------------
rayleigh<-function(s)
{n <- 1000
u <- runif(n)
x <- sqrt(-2*s^2*log(1-u))
hist(x, prob = TRUE, main = expression(f(x)==(x/s^2)*exp(-(x^2)/(2*s^2))))
y <- seq(0, 100, 1)
lines(y, y/s^2*exp(-(y^2)/(2*s^2)))
}
####模拟σ=1,2,3的三种情况
rayleigh(1)
rayleigh(2)
rayleigh(3)

## -----------------------------------------------------------------------------
mix<-function(n,p)
{
   u=rep(0, n)
   x=rep(0, n)
   for(i in 1:n)
   {
      u[i]=runif(1)
      if(u[i]<p){
         x[i]=rnorm(1)
         }
      else{
         x[i]=rnorm(1,3,1)
         }
   }
   hist(x,prob = TRUE,breaks=100,main = p)
   y=seq(min(x),max(x),0.01)
   lines(y,p/sqrt(2*pi)*exp(-y^2/2)+(1-p)/sqrt(2*pi)*exp(-(y-3)^2/2))
}
##模拟p1=0.75，0.4，0.5，,0.6的四种情况
mix(n = 1000,p=0.75)
mix(n = 1000,p=0.4)
mix(n = 1000,p=0.5)
mix(n = 1000,p=0.6)

## -----------------------------------------------------------------------------
##生成泊松过程和伽马分布的复合随机变量，其中泊松过程参数为lambda，伽马分布的参数为r和beta
crp=function(n,lambda,r,beta)
{  x=0
   for(i in 1:n)
      {
      n=rpois(1,10*lambda)
      y=rgamma(n,r,beta)
      x[i]=sum(y)
   }
   estimate=c(mean(x),var(x))
   theore=c(10*lambda*r/beta,10*lambda*(r/(beta^2)+(r^2)/(beta^2)))
   list(estimate=estimate,theore=theore)
}
##选择不同的参数进行模拟
crp(10000,1,1,1)
crp(10000,3,2,4)
crp(10000,5,3,5)

## -----------------------------------------------------------------------------
beta=function(q)
{
   n=length(q)
   theta.hat=c(0,n)
   for (i in 1:n) {
      g=function(u)q[i]*30*u^2*(1-u)^2
      y=runif(1e6,0,q[i])
      theta.hat[i]=mean(g(y))
  }
   theta.hat
}
q=seq(0.1,0.9,0.1)
list(MCe=round(beta(q),5),ture=pbeta(q,3,3))

## -----------------------------------------------------------------------------
raylaigh=function(x,s=sigma)
{
m=50000
u=runif(m,0,x)
T1=(x*u/s^2*exp(-u^2/(2*s^2))+x*(x-u)/s^2*exp(-(x-u)^2/(2*s^2)))/2
u1=runif(m,0,x)
T2=(x*u/s^2*exp(-u^2/(2*s^2))+x*u1/s^2*exp(-u1^2/(2*s^2)))/2
list(meanT1=mean(T1),meanT2=mean(T2),varT1=var(T1)/m,varT2=var(T2)/m,percent_reduction=(var(T2)-var(T1))/var(T2))
}
##取x=1，sigma=2
raylaigh(1,2)
##取x=2，sigma=2
raylaigh(2,2)
##取x=3，sigma=2
raylaigh(3,2)

## -----------------------------------------------------------------------------
##取f为标准正态分布的密度函数
m=1e6
x=rnorm(m)
g=x^2*(x>=1)
list(theta.hat=mean(g),vartheta.hat=var(g)/m)
##取f为指数分布的密度函数
y=rexp(m)
g1=y^2/sqrt(2*pi)*exp(y-y^2/2)*(y>=1)
list(theta1.hat=mean(g1),vartheta1.hat=var(g1)/m)
c(var(g)/var(g1))

## -----------------------------------------------------------------------------
CI=function(m,n,alpha){
   ucl=lcl=cp=numeric(m)
   for(i in 1:m){
      x=rchisq(n,2)
      ucl[i]=mean(x)+sd(x)/sqrt(n)*qt(1-alpha/2,df=n-1)
      lcl[i]=mean(x)-sd(x)/sqrt(n)*qt(1-alpha/2,df=n-1)
      cp[i]=as.numeric(ucl[i]>=2&&lcl[i]<=2)
   }
   return(mean(cp))
}
CI(m=1000,n=20,alpha=0.05)


## -----------------------------------------------------------------------------

test=function(m,n,alpha){
   p.val1=p.val2=p.val3=numeric(m)
for(i in 1:m){
   x1=rchisq(n,df=1)##x~chisq(1)
   p.val1[i]=2*(1-pt(abs(sqrt(n)*(mean(x1)-1)/sd(x1)),n-1))
   x2=runif(n,0,2)##x~u(0,2)
   p.val2[i]=2*(1-pt(abs(sqrt(n)*(mean(x2)-1)/sd(x2)),n-1))
   x3=rexp(20,1)#x~exp(1)
   p.val3[i]=2*(1-pt(abs(sqrt(n)*(mean(x3)-1)/sd(x3)),n-1))
}
list(chisq=mean(p.val1<=alpha),uniform=mean(p.val2<=alpha),exponential=mean(p.val3<=alpha))
}
test(m=10000,alpha=0.01,n=20)
test(m=10000,alpha=0.05,n=20)
test(m=10000,alpha=0.1,n=20)



## -----------------------------------------------------------------------------
##Examples 6.8
set.seed(123)
library(MASS)
alpha <- 0.05
d <- 2 #The dimension of multivariate normal variables is two-dimension.
n <- c(10,20,30,40,50,100,500) #sample size
cv <- qchisq(1-alpha,d*(d + 1)*(d + 2)/6) #crit.values for each n
msk <- function(x) {
  xbar <- mean(x)
  n <- nrow(x)
  msk <- 1/n^2*sum(diag(((x-xbar)%*%solve(cov(x))%*%t(x-xbar))^3))
  return(msk*n/6)
}

p.reject <- numeric(length(n))
m <- 10000
s <- matrix(c(1,0,0,1),2,2)
for (i in 1:length(n)) {
  sktests <- numeric(m)  #test decisions
  for (j in 1:m) {
    x <- mvrnorm(n[i],c(0,0),s)
    sktests[j] <- as.integer(msk(x) >= cv)
  }
  p.reject[i] <- mean(sktests)  #proportion rejected
}
p.reject


## -----------------------------------------------------------------------------
##Examples 6.10
library(MASS)
set.seed(123)
alpha<-0.1
d<-2
n <-30 #sample sizes 
cv <- qchisq(1-alpha,d*(d+1)*(d+2)/6) #crit. values for each n
sk <-function(x) { 
  n<-nrow(x)
  for (j in 1:d) {
    x[,j]<-x[,j]-mean(x[,j])
  }
  s<-solve(cov(x))
  b<-1/n^2*sum(diag((x%*%s%*%t(x))^3))
  return(b*n/6)
}

m <- 2500 #replicates
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05)) 
N <- length(epsilon) 
pwr <- numeric(N) #critical value for the skewness test 

for (j in 1:N) {
  e <- epsilon[j] 
  sktests <- numeric(m) 
  for (i in 1:m) {
    sig <- sample(c(1,10), replace = TRUE, size = n, prob = c(1-e, e))
    x <- matrix(nr=n,nc=d)
    for (k in 1:n) {
      sigma<-diag(rep(sig[k],d))
       x[k,]<- mvrnorm(1,rep(0,d),sigma)
    }
    sktests[i] <- as.integer(sk(x) >= cv) 
  }
  pwr[j] <- mean(sktests) 
}
#plot power vs epsilon 
plot(epsilon, pwr, type = "b", xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3) 
se <- sqrt(pwr * (1-pwr) / m) #add standard errors 
lines(epsilon, pwr+se, lty = 3) 
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
set.seed(123)

lambda_hat=eigen(cov(scor))$values
theta_hat=lambda_hat[1]/sum(lambda_hat)
B=2000 ##number of replictes

theta_star=function(data,index){
  x=data[index,]
  lambda=eigen(cov(x))$values
  theta=lambda[1]/sum(lambda)
  return(theta)
}
bootstrap_result=boot(data=cbind(scor$mec, scor$vec, scor$alg, scor$ana, scor$sta),statistic = theta_star, R = B)
theta_b <- bootstrap_result$t
bias_boot <- mean(theta_b) - theta_hat# the estimated bias of theta_hat, using bootstrap
se_boot <- sqrt(var(theta_b))# the estimated standard error (se) of theta_hat, using bootstrap

# print the answers
list(bias_boot=bias_boot,se_boot=se_boot)

## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(123)


lambda_hat=eigen(cov(scor))$values
theta_hat=lambda_hat[1]/sum(lambda_hat)
n=nrow(scor)

theta_j=numeric(n)
for(i in 1:n){
  x=scor[-i,]
  lambda=eigen(cov(x))$values
  theta_j[i]=lambda[1]/sum(lambda)
}
bias_jack=(n-1)*(mean(theta_j)-theta_hat)
se_jack=sqrt((n-1)*mean(((theta_j-mean(theta_j))^2)))

#print the answers
list(bias_jack=bias_jack,se_jack=se_jack)


## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
set.seed(123)

lambda_hat=eigen(cov(scor))$values
theta_hat=lambda_hat[1]/sum(lambda_hat)
B=2000 ##number of replictes

theta_star=function(data,index){
  x=data[index,]
  lambda=eigen(cov(x))$values
  theta=lambda[1]/sum(lambda)
  return(theta)
}
bootstrap_result=boot(data=cbind(scor$mec, scor$vec, scor$alg, scor$ana, scor$sta),statistic = theta_star, R = B)
boot.ci(bootstrap_result,conf = 0.95,type = c("perc","bca"))

## -----------------------------------------------------------------------------
library(boot)
m=1e3#number of replicates

# normal populations
ci_norm.norm=ci_norm.basic=ci_norm.perc=matrix(NA,m,2)
boot_skew=function(x,i){
  mean((x-mean(x[i]))^3)/(mean((x-mean(x[i]))^2))^(1.5)
}
for(i in 1:m){
  x=rnorm(n)
  de=boot(data=x,statistic=boot_skew, R = 1000)
  ci_norm=boot.ci(de,conf=0.95,type=c("norm","basic","perc"))
  ci_norm.norm[i,]=ci_norm$norm[2:3];ci_norm.basic[i,]=ci_norm$basic[4:5];ci_norm.perc[i,]=ci_norm$percent[4:5]
}

# $chisq$(5) distribution
ci_chisq.norm=ci_chisq.basic=ci_chisq.perc=matrix(NA,m,2)
for(i in 1:m){
  y=rchisq(n,df=5)
  de=boot(data=y,statistic=boot_skew, R = 1000)
  ci_chisq=boot.ci(de,conf=0.95,type=c("norm","basic","perc"))
  ci_chisq.norm[i,]=ci_chisq$norm[2:3];ci_chisq.basic[i,]=ci_chisq$basic[4:5];ci_chisq.perc[i,]=ci_chisq$percent[4:5]
}
chisq.skew=sqrt(8/5)#the skewness of chisquare distribution

coverage_rate=c(mean(ci_norm.norm[,1]<=0 & ci_norm.norm[,2]>=0),
                mean(ci_norm.basic[,1]<=0 & ci_norm.basic[,2]>=0),
                mean(ci_norm.perc[,1]<=0 & ci_norm.perc[,2]>=0),
                mean(ci_chisq.norm[,1]<=chisq.skew & ci_chisq.norm[,2]>=chisq.skew),
                mean(ci_chisq.basic[,1]<=chisq.skew & ci_chisq.basic[,2]>=chisq.skew),
                mean(ci_chisq.perc[,1]<=chisq.skew & ci_chisq.perc[,2]>=chisq.skew))

proportipn_miss_on_the_left=c(mean(ci_norm.norm[,1]>0),
                              mean(ci_norm.basic[,1]>0),
                              mean(ci_norm.perc[,1]>0),
                              mean(ci_chisq.norm[,1]>chisq.skew),
                              mean(ci_chisq.basic[,1]>chisq.skew),
                              mean(ci_chisq.perc[,1]>chisq.skew))
                              
  
proportipn_miss_on_the_right=c(mean(ci_norm.norm[,2]<0),
                               mean(ci_norm.basic[,2]<0),
                               mean(ci_norm.perc[,2]<0),
                               mean(ci_chisq.norm[,2]<chisq.skew),
                               mean(ci_chisq.basic[,2]<chisq.skew),
                               mean(ci_chisq.perc[,2]<chisq.skew))
cinames=c("ci_norm.norm","ci_norm.basic","ci_norm.perc","ci_chisq.norm","ci_chisq.basic","ci_chisq.perc")

df=data.frame(cinames=cinames,coverage_rate=coverage_rate,proportipn_miss_on_the_left=proportipn_miss_on_the_left,proportipn_miss_on_the_right=proportipn_miss_on_the_right)
library(knitr)

kable(df) 

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "sunflower"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)
cat(x,'\n')
cat(y,'\n')
z <- c(x, y) # pooled sample
R <- 1500;K <- 1:24;n<-length(x)
set.seed(12345)
reps <- numeric(R);t0 <- cor(x, y,method ="spearman")
for (i in 1:R) {
k <- sample(K, size = n, replace = FALSE)
x1 <- z[k]; y1 <- z[-k] #complement of x1
reps[i] <- cor(x1, y1,method ="spearman")
}
p <- mean(abs(c(t0, reps)) >= abs(t0))
round(c(p,cor.test(x,y)$p.value),3)

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(energy)
library(Ball)

Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) 
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 <= n1); i2 <- sum(block2 > n1)
(i1 + i2) / (k * n)
}

#Power comparison 
m <- 1e3; k<-1; set.seed(12345)
n1 <- n2 <- 50; R<-999; n <- n1+n2; N <- c(n1,n2);
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
p1.values <- p2.values<-p3.values<-p4.values<-matrix(NA,m,3)

#Unequal variances and equal expectations
for(i in 1:m){
x1 <- rnorm(n1,0,1.5)
y1 <- rnorm(n2,0,1)
z1 <- c(x1,y1)
p1.values[i,1]<-eqdist.nn(z1,N,k)$p.value
p1.values[i,2]<-eqdist.etest(z1,sizes=N,R=R)$p.value
p1.values[i,3]<-bd.test(x=x1,y=y1,num.permutations=R,
seed=i*12345)$p.value
}
#Unequal variances and unequal expectations
for(i in 1:m){
x2 <- rnorm(n1,0,1.5)
y2 <- rnorm(n2,1,1)
z2 <- c(x2,y2)
p2.values[i,1]<-eqdist.nn(z2,N,k)$p.value
p2.values[i,2]<-eqdist.etest(z2,sizes=N,R=R)$p.value
p2.values[i,3]<-bd.test(x=x2,y=y2,num.permutations=R,
seed=i*12345)$p.value
}
#Non-normal distributions: t distribution with 1 df (heavy-tailed distribution),bimodel distribution (mixture of two normal distributions)
rMG <- function(n){
# randomly generate n points from the Mixed Gaussian distribution
r <- runif(n, 0, 1)
x <- r
ind <- which(r < 0.3) #index for those generated from N(0,1)
x[ind] <- rnorm(length(ind), 0, 1)
x[-ind] <- rnorm(n-length(ind), 1, 0.3)
return(x)
}
for(i in 1:m){
x3 <- rt(n1,df=1)
y3 <- rMG(n2)
z3 <- c(x3,y3)
p3.values[i,1]<-eqdist.nn(z3,N,k)$p.value
p3.values[i,2]<-eqdist.etest(z3,sizes=N,R=R)$p.value
p3.values[i,3]<-bd.test(x=x3,y=y3,num.permutations=R,
seed=i*12345)$p.value
}
#Unbalanced samples
for(i in 1:m){
x4 <- rnorm(5)
y4 <- rnorm(n2)
z4 <- c(x4,y4)
p4.values[i,1]<-eqdist.nn(z4,c(5,n2),k)$p.value
p4.values[i,2]<-eqdist.etest(z4,sizes=c(5,n2),R=R)$p.value
p4.values[i,3]<-bd.test(x=x4,y=y4,num.permutations=R,
seed=i*12345)$p.value
}

alpha <- 0.05;
pow1 <-colMeans(p1.values<alpha)
pow2 <-colMeans(p2.values<alpha)
pow3 <-colMeans(p3.values<alpha)
pow4 <-colMeans(p4.values<alpha)
test_names<-c("NN","energy"," ball")
df<-data.frame(test_names=test_names,uneqva_eqme=pow1,uneqva_uneqme=pow2,Non_normal_dist=pow3,Unbalanced_samples=pow4)
library(knitr)
kable(df)

## -----------------------------------------------------------------------------
m <- 10000
cauchy.chain <- function(sigma, N, X1) {
#generates a Metropolis chain for standard Cauchy distribution
#with Normal(X[t], sigma) proposal distribution
#and starting value X1
x <- rep(0, N)
x[1] <- X1
u <- runif(N)
for (i in 2:N) {
xt <- x[i-1]
y <- rnorm(1, xt, sigma) #candidate point
r1 <- dcauchy(y) * dnorm(xt, y, sigma)
r2 <- dcauchy(xt) * dnorm(y, xt, sigma)
r <- r1 / r2
if (u[i] <= r) x[i] <- y else
x[i] <- xt
}
return(x)
}
x1=cauchy.chain(0.05,m,-10)
x2=cauchy.chain(0.5,m,-10)
x3=cauchy.chain(2,m,-10)
x4=cauchy.chain(4,m,-10)

b <- 1001 #discard the burn-in sample
z1 <- x1[b:m]
z2 <- x2[b:m]
z3 <- x3[b:m]
z4 <- x4[b:m]
deciles <- seq(0,1,0.1)
deciles_cauchy<-qcauchy(deciles)
deciles_mcmc1<-quantile(z1, deciles)
deciles_mcmc2<-quantile(z2, deciles)
deciles_mcmc3<-quantile(z3, deciles)
deciles_mcmc4<-quantile(z4, deciles)
df<-data.frame(deciles_cauchy=deciles_cauchy,sigma_0.05=deciles_mcmc1,sigma_0.5=deciles_mcmc2,sigma_2=deciles_mcmc3,sigma_4=deciles_mcmc4)
library(knitr)
kable(df)


## -----------------------------------------------------------------------------
#Gelman-Rubin method of monitoring convergence
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

cauchy.chain <- function(sigma, N, X1) {
#generates a Metropolis chain for standard Cauchy distribution
#with Normal(X[t], sigma) proposal distribution
#and starting value X1
x <- rep(0, N)
x[1] <- X1
u <- runif(N)
for (i in 2:N) {
xt <- x[i-1]
y <- rnorm(1, xt, sigma) #candidate point
r1 <- dcauchy(y) * dnorm(xt, y, sigma)
r2 <- dcauchy(xt) * dnorm(y, xt, sigma)
r <- r1 / r2
if (u[i] <= r) x[i] <- y else
x[i] <- xt
}
return(x)
}


sigma0 <- 0.5 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 15000 #length of chains
b <- 1000 #burn-in length

#choose overdispersed initial values
x0 <- c(-1, 0.5, 0.5, 1)

#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- cauchy.chain(sigma0, n, x0[i])

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n],ylim = c(0,2.5), type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)



## -----------------------------------------------------------------------------
N<-5000
burn<-1000
X <- matrix(0, N, 2)

a<-1
b<-1
n<-10

X[1, ] <- c(5, 0.5) #initialize

for (i in 2:N) {
x2 <- X[i-1, 2]
X[i, 1] <- rbinom(1, n, x2)
x1 <- X[i, 1]
X[i, 2] <- rbeta(1,x1+a, n-x1+b)
}

b <- burn + 1
x <- X[b:N, ]

plot(x, main="", cex=.5, xlab=bquote(X[1]),
ylab=bquote(X[2]), ylim=range(x[,2]))





## -----------------------------------------------------------------------------
#Gelman-Rubin method of monitoring convergence
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

cauchy.chain <- function(sigma, N, X1) {
#generates a Metropolis chain for standard Cauchy distribution
#with Normal(X[t], sigma) proposal distribution
#and starting value X1
x <- rep(0, N)
x[1] <- X1
u <- runif(N)
for (i in 2:N) {
xt <- x[i-1]
y <- rnorm(1, xt, sigma) #candidate point
r1 <- dcauchy(y) * dnorm(xt, y, sigma)
r2 <- dcauchy(xt) * dnorm(y, xt, sigma)
r <- r1 / r2
if (u[i] <= r) x[i] <- y else
x[i] <- xt
}
return(x)
}


sigma0 <- 0.5 #parameter of proposal distribution
k <- 4 #number of chains to generate
n <- 15000 #length of chains
b <- 1000 #burn-in length

#choose overdispersed initial values
x0 <- c(-1, 0.5, 0.5, 1)

#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
X[i, ] <- cauchy.chain(sigma0, n, x0[i])

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n],ylim = c(0,2.5), type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
kth_term<-function(a,k){
  a0<-as.matrix(a)
  a_norm<-norm(a0,"2")
  d<-length(a)
  tmp<-(-1)^k*a_norm^(2*k+2)/(exp(lgamma(k+1))*2^k*(2*k+1)*(2*k+2))
  kthterm<-tmp*exp(lgamma((d+1)/2)+lgamma(k+1.5)-lgamma(k+d/2+1))
  return(kthterm)
}
kth_term(c(1,2),100)
kth_term(c(1,2),400)

## -----------------------------------------------------------------------------
kth_term<-function(a,k){
  a0<-as.matrix(a)
  a_norm<-norm(a0,"2")
  d<-length(a)
  tmp<-(-1)^k*a_norm^(2*k+2)/(exp(lgamma(k+1))*2^k*(2*k+1)*(2*k+2))
  kthterm<-tmp*exp(lgamma((d+1)/2)+lgamma(k+1.5)-lgamma(k+d/2+1))
  return(kthterm)
}
mysum<-function(a,k){
  s<-0;s0<-0
  for(i in 0:k){
  s0<-kth_term(a,i)
  s<-s+s0
  }
  return(s)
}
  
mysum(c(1,2),100)
mysum(c(1,2),400)

## -----------------------------------------------------------------------------
#calculate c_k

f <- function(k, u){
  tmp <- 2*exp(lgamma((k+1)/2)-0.5*log(pi*k)-lgamma(k/2))
  inte <- tmp*(1+u^2/k)^(-(k+1)/2)
  return(inte)
}
  
equation<-function(k, a){
  up.bound.1 <-  sqrt(a^2*(k-1)/(k-a^2))
  up.bound.2 <- sqrt(a^2*k/(k+1-a^2))

  integrate.1 <- integrate(f,lower = 0, upper = up.bound.1, k = k-1)
  integrate.2 <- integrate(f,lower = 0, upper = up.bound.2 , k = k)
  
  return(integrate.1$value - integrate.2$value)
}

root <- function(k){
  a <- uniroot(equation, interval = c(0.00001, sqrt(k)-0.001), k = k)$root
  return(a)
  }



k <- 2:15
for (i in k) {
 print(root(i))

}


## -----------------------------------------------------------------------------
y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
logL <- function(lambda) {
# log-likelihood
return(-7*log(lambda)-sum(y)/lambda)
}
res <- optimize(logL,lower=0,upper=10,maximum=TRUE)
res$maximum

#the iterative formula of EM algorithm

lambda.0 <- 0.1
lambda.1 <- 0.5
while(abs(lambda.0 - lambda.1)>1e-5){
  lambda.0<- lambda.1
  lambda.1 <- lambda.0*0.3 + 0.675
}
list(MLE=res$maximum,EM_estimate=lambda.1)


## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

#####exercise 3
mpg<-mtcars$mpg
disp<-mtcars$disp
wt<-mtcars$wt
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
#using for loop
out <- vector("list", length(formulas))
for(i in seq_along(formulas)){
  out[[i]] <- rsq(lm(formulas[[i]]))
}
out#using for loop
#using lapply()
fit<-lapply(formulas, lm)
lapply(fit,rsq)#using lapply()

#####exercise 4
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
#using for loop
out1 <- vector("list", length(bootstraps))
for(i in seq_along(bootstraps)){
  out1[[i]] <- rsq(lm(mpg ~ disp,data = bootstraps[[i]]))
}
out1#using for loop
#using lapply()
fit1<-lapply(bootstraps, lm, formula=mpg ~ disp)
lapply(fit1, rsq)#using lapply()


## -----------------------------------------------------------------------------
#Compute the standard deviation of every column in a numeric data frame.
df<-data.frame(replicate(6,sample(c(1:10),10,rep=T)))
round(vapply(df,sd,FUN.VALUE=c(sd=0)),3)

#Compute the standard deviation of every numeric column in a mixed data frame.
mixdf<-data.frame(replicate(6,sample(c(1:10),10,rep=T)),x7=c("a","b","c","d","e","f","g","h","i","j"))
tmp<-vapply(mixdf,is.numeric,FUN.VALUE=c(c=0))
vapply(mixdf[,tmp],sd,FUN.VALUE=c(c=0))



## -----------------------------------------------------------------------------
library(parallel)
library(energy)
mcsapply <- function(x,f){
  num.clust <- makeCluster(5)
  f <- parSapply(num.clust,x,f)
  stopCluster(num.clust)
  return(f)
}

tl <- replicate(100,edist(matrix(rnorm(100),ncol=10),sizes=c(5,5)))
system.time(mcsapply(tl,function(x){return(x)}))
system.time(sapply(tl,function(x){return(x)}))

## -----------------------------------------------------------------------------
#R-code
bivariate.density.R <- function(a,b,n){
  N<-5000
  burn<-1000
  X <- matrix(0, N, 2)
  X[1, ] <- c(n/2, 0.5) #initialize
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, x2)
    x1 <- X[i, 1]
    X[i, 2] <- rbeta(1,x1+a, n-x1+b)
  }
  b <- burn + 1
  x <- X[b:N, ]
  return(x)
}

#C++code
library(Rcpp)
cppFunction('NumericMatrix bivariate_density_c(double a,double b,double n){
  int N=5000;
  int burn=1000;
  NumericMatrix X(N,2);
  NumericVector v={n/2,0.5};
  X(0,_)=v;
  for(int i=1;i<N;i++){
    double x2=X(i-1,1);
    X(i,0)=as<int>(rbinom(1,n,x2));
    int x1=X(i,0);
    X(i,1)=as<double>(rbeta(1,x1+a,n-x1+b));
  }
  NumericMatrix x=X(Range(burn-1,N-1), Range(0,1));
  return(x);
}')

#initialize
a<-1
b<-1
n<-10

#qqplot
X.R<-bivariate.density.R(a,b,n)
X.C<-bivariate_density_c(a,b,n)
x.r.1<-X.R[,1]
x.r.2<-X.R[,2]
x.c.1<-X.C[,1]
x.c.2<-X.C[,2]
qqplot(x.r.1,x.c.1,xlab = "Discrete variable from R",ylab = "Discrete variable from Rcpp")
qqplot(x.r.2,x.c.2,xlab = "Continuous variable from R",ylab = "Continuous variable from Rcpp")

#Campare the computation time
library(microbenchmark)
ts <- microbenchmark(bivariate.density.R=bivariate.density.R(a,b,n), bivariate_density_c=bivariate_density_c(a,b,n))
summary(ts)[,c(1,3,5,6)]


