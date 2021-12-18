#' @title PACS procedure for variable selection
#' @description This code estimates coefficient for the PACS procedure.
#' @param X scaled design matrix
#' @param y centered vector of response
#' @param lambda non‐negative tuning parameter
#' @param betawt adaptive weights, usually OLS/ridge coefficient estimates
#' @param type 1 for Adaptive PACS 2 for Adaptive Correlated PACS 3 for Threshold with Adaptive PACS 4 for Threshold with Adaptive Correlated PACS
#' @param rr correlation for Threshold PACS approaches, a value between 0 and 1(needed for type=3 and type=4)
#' @param eps criteria for convergence
#' @return estimated PACS coefficients
#' @examples 
#' \dontrun{
#' library(MASS)
#' n<-50
#' c1<-c(1,0.7,0.7,rep(0,5))
#' c2<-c(0.7,1,0.7,rep(0,5))
#' c3<-c(0.7,0.7,1,rep(0,5))
#' c4<-c(rep(0,3),1,rep(0,4))
#' c5<-c(rep(0,4),1,rep(0,3))
#' c6<-c(rep(0,5),1,rep(0,2))
#' c7<-c(rep(0,6),1,0)
#' c8<-c(rep(0,7),1)
#' s<-rbind(c1,c2,c3,c4,c5,c6,c7,c8)
#' x<-mvrnorm(n,rep(1,8),s)
#' beta<-c(2,2,2,rep(0,5))
#' eps<-rnorm(n)
#' y<-x%*%beta+eps
#' betawt<-summary(lm(y~x))$coefficients[2:9]
#' PACS(y,x,lambda=1,betawt=betawt,type=1,rr=0,eps=10^-5)
#' }
#' @import MASS
#' @import mvtnorm
#' @import Matrix
#' @importFrom stats cor cov
#' @useDynLib StatComp21075
#' @export
PACS=function(y,X,lambda,betawt,type=1,rr=0,eps=10^-5)
{
  requireNamespace(mvtnorm)
  requireNamespace(MASS)
  requireNamespace(Matrix)
  if(lambda <=0) {return(cat("ERROR: Lambda must be > 0. \n"))}
  if(rr <0) {return(cat("ERROR: RR must be >=0. \n"))}
  if(rr >1) {return(cat("ERROR: RR must be <=1. \n"))}
  if(eps <=0) {return(cat("ERROR: Eps must be > 0. \n"))}
  if(!(type %in% 1:4)){return(cat("ERROR: Type must be 1, 2, 3 or 4. \n"))}
  x=scale(x)
  y=y-mean(y)
  n=length(y)
  p=dim(x)[2]
  littleeps=10^-7
  # require p>1
  qp=p*(p-1)/2
  vc=0
  row.vec1=0
  col.vec1=0
  for(w in 1:(p-1))
  {
    row.vec1=c(row.vec1,c((vc+1):(vc+(p-w))))
    col.vec1=c(col.vec1,rep(w,length(c((vc+1):(vc+(p-w))))))
    vc=vc+(p-w)
  }
  c1=1
  c2=p-1
  row.vec2=0
  col.vec2=0
  for(w in 1:(p-1))
  {
    row.vec2=c(row.vec2,c(c1:c2))
    col.vec2=c(col.vec2,c((w+1):p))
    c1=c2+1
    c2=c2+p-w-1
  }
  dm=sparseMatrix(i=c(row.vec1[-1],row.vec2[-1]),j=c(col.vec1[-1],col.vec2[-1]),x=c(rep(1,qp),rep(-1,qp)))
  dp=sparseMatrix(i=c(row.vec1[-1],row.vec2[-1]),j=c(col.vec1[-1],col.vec2[-1]),x=c(rep(1,qp),rep(1,qp)))
  rm(c1,c2,w,vc,row.vec1,col.vec1,row.vec2,col.vec2)
  if(type==1)
  {
    ascvec=c(1/abs(betawt),1/as.vector(abs(dm%*%betawt)),1/as.vector(abs(dp%*%betawt)))
  }
  if(type==2)
  {
    cor.mat=cor(x)
    crm=1/(1-cor.mat[1,2:p])
    for(cc in 2:(p-1)){crm=c(crm,1/(1-cor.mat[cc,(cc+1):p]))}
    crp=1/(1+cor.mat[1,2:p])
    for(cc in 2:(p-1)){crp=c(crp,1/(1+cor.mat[cc,(cc+1):p]))}
    rm(cor.mat)
    ascvec=c(1/abs(betawt),crm/as.vector(abs(dm%*%betawt)),crp/as.vector(abs(dp%*%betawt)))
  }
  if(type==3)
  {
    corp=ifelse(cor(x)>rr,1,0)
    cor.mat=cor(x)
    crm=corp[1,2:p]
    for(cc in 2:(p-1)){crm=c(crm,corp[cc,(cc+1):p])}
    rr=-rr
    corm=ifelse(cor(x)<rr,1,0)
    crp=corm[1,2:p]
    for(cc in 2:(p-1)){crp=c(crp,corm[cc,(cc+1):p])}
    ascvec=c(1/abs(betawt),crm/as.vector(abs(dm%*%betawt)),crp/as.vector(abs(dp%*%betawt)))
  }
  if(type==4)
  {
    corp=ifelse(cor(x)>rr,1,0)
    cor.mat=cor(x)
    crm=corp[1,2:p]/(1-cor.mat[1,2:p])
    for(cc in 2:(p-1)){crm=c(crm,corp[cc,(cc+1):p]/(1-cor.mat[cc,(cc+1):p]))}
    rr=-rr
    corm=ifelse(cor(x)<rr,1,0)
    crp=corm[1,2:p]/(1+cor.mat[1,2:p])
    for(cc in 2:(p-1)){crp=c(crp,corm[cc,(cc+1):p]/(1+cor.mat[cc,(cc+1):p]))}
    ascvec=c(1/abs(betawt),crm/as.vector(abs(dm%*%betawt)),crp/as.vector(abs(dp%*%betawt)))
  }
  betal=betawt
  betaeps=1
  while(betaeps>eps)
  {
    betatilde=betal
    mm1=diag(ascvec[1:p]/(abs(betatilde)+littleeps))
    mm2=diag(ascvec[(p+1):(p+qp)]/(as.vector(abs(dm%*%betatilde))+littleeps))
    mm3=diag(ascvec[(p+qp+1):(p+2*qp)]/(as.vector(abs(dp%*%betatilde))+littleeps))
    betal=as.vector(solve(t(x)%*%x+lambda*(mm1+t(dm)%*%mm2%*%dm+t(dp)%*%mm3%*%dp))%*%t(x)%*%y)
    rm(mm1,mm2,mm3)
    betaeps=max(abs(betatilde-betal)/(abs(betatilde)+littleeps))
  }
  fit=NULL
  fit$coefficients=betatilde
  fit$lambda=lambda
  fit$cor=abs(rr)
  fit$weights=betawt
  fit$type=type
  fit$eps=eps
  fit$littleeps=littleeps  
  return(fit)
}

#' @title LASSO estimates 
#' @description This code estimates coefficient for the LASSO procedure.
#' @param y vector of response
#' @param x design matrix
#' @return the lasso estimator of coefficients
#' @import lars
#' @export
coef_lasso<-function(y,x){
  requireNamespace(lars)
  object1<-lars(x,y,type="lasso")
  cv_sol<-cv.lars(x,y,type="lasso",mode="step")
  fra=cv_sol$index[which.min(cv_sol$cv)]
  betahat.lasso<-object1$beta[fra,]
  return(betahat.lasso)
}


#' @title model error (ME) 
#' @description model error (ME) for prediction accuracy
#' @param beta true coefficient
#' @param betahat the estimated coeffcient from a specific method
#' @param x the design matrix
#' @return the model error
#' @importFrom stats cor cov
#' @export
ME<-function(beta,betahat,x){
  v<-cov(x)
  me<-t((betahat-beta))%*%v%*%(betahat-beta)
  return(me[1,1])
}