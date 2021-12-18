#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp to generate a chain with target joint density
//' @param a fixed number for beta distribution
//' @param b fixed number for beta distribution
//' @param n fixed number for binomial distribution
//' @return a binary Monte Carlo Markov Chain
//' @export
// [[Rcpp::export]]
NumericMatrix bivariate_density_c(double a,double b,double n){
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
}
