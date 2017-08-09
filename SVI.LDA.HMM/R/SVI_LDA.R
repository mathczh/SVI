library(Rcpp)
library(RcppArmadillo)
sourceCpp(file='src/LDA.cpp')
LDA.SVI<-function(X, # the Input data
                  K, # the number of the topics
                  n, # the number of the passes
                  alpha, # the parameter of theta
                  eta, # the parameter of beta
                  pre, # the stop criteria
                  k# the number of top words and topics

)
{
     SVI_LAD(X,K,n,alpha,eta,pre,k)
}
X=list(c('home','mother','mother','house','mother','mother','mother','mother'),c('father','dad','cat','dad','cat','cat','cat','cat','cat','cat','cat','cat','cat','cat'))
K=3;
alpha=1/K;
eta=0.01;
n=2;
pre=0.1;
t=LDA.SVI(X,K,n,alpha=1/K,eta=0.01,pre=0.1,2)
