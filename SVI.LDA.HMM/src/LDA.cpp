// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>
#include <string>
#include <map>
#include<Rinternals.h>
#include<math.h>
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]
SEXP SVI_LDA(Rcpp:: List &X,
             SEXP  KK,
             SEXP  nn,
             SEXP alphaa ,
             SEXP  etaa,
             SEXP  pree,
             SEXP  topic_lengthh
               )
{
  int K = Rcpp::as<int>(KK);
  int topic_length = Rcpp::as<int>(topic_lengthh);
  int n = Rcpp::as<int>(nn);
  double alpha = Rcpp::as<double>(alphaa);
  double eta = Rcpp::as<double>(etaa);
  double pre = Rcpp::as<double>(pree);
  std::map<std::string, int> mymap;
  std::map<std::string, int>::iterator iter;
  for(int i= 0;i< X.size();i++) //Produce map to stroe the whole vocabulary
  {
    Rcpp::StringVector document = Rcpp::as<Rcpp::StringVector>(X[i]) ;
    for(int j=0;j< document.size();j++)
    {
      mymap[""+document[j]+""]=0;
    }
  }
  int D= X.size();
  int V= mymap.size();
  Rcpp::StringVector Dictionary (V) ;
  int value =0;
  for(iter=mymap.begin();iter !=mymap.end();iter++)
  {
    mymap[""+iter->first+""]=value;
    Dictionary[iter->second]=iter->first;
    value++;
  }
  Rcpp::NumericMatrix Lambda(V,K);
  for( int i=0;i<V;i++)//initiate Lambda with rchisq
      {
         for(int j=0;j<K;j++)
           {
               Lambda(i,j)=R::rchisq(1);
           }
       }
  Rcpp::NumericMatrix gamma(K,D) ; //random initial with all 0 entry
  Rcpp::List phi(D);
  for(int i=0;i<n;i++)
  {
    for(int d=0;d<D;d++)
    {
      for(int v=0;v<K;v++)
      {
        gamma(v,d)=1;
      }
      Rcpp::NumericVector gamma1(K);
      Rcpp::StringVector Y = Rcpp::as<Rcpp::StringVector>(X[d]) ;
      Rcpp::NumericMatrix w(V,Y.size()) ;
      Rcpp::NumericMatrix mid(Y.size(),K);
      Rcpp::NumericMatrix mid_norm(Y.size(),K);
      double norm = sum(abs(gamma(_,d)-gamma1));//calculate norm as stop criteria
      while(norm>=pre)
      {
        gamma1=gamma(_,d);
        for(int j=0;j<Y.size();j++)
        {
          Rcpp::NumericVector w_temp(V);
          iter =mymap.find(""+Y[j]+"");
          int index_inmap = iter->second;
          w_temp[index_inmap]=1;
          w(_,j)= w_temp;
          mid(j,_)=exp(digamma(gamma(_,d))+digamma(Lambda(index_inmap,_))-digamma(colSums(Lambda)));
          mid_norm(j,_)=mid(j,_)/sum(mid(j,_));//normalize parameter
        }
        gamma(_,d)= alpha+colSums(mid_norm);
        norm = sum(abs(gamma(_,d)-gamma1));
      }
      phi[d]= mid_norm;
      arma::mat midnorm_temp = as<arma::mat>(mid_norm);//convert NumericMtrix into arma::mat;
      arma::mat w_temp = as<arma::mat>(w);
      arma::mat result_temp=eta+D*w_temp*midnorm_temp;
      arma::mat Lambda_temp = as<arma::mat>(Lambda);
      Lambda_temp = (1-1.0/(i*D+d+1))*Lambda_temp+1.0/(i*D+d+1)*result_temp;
      Lambda=Rcpp::wrap(Lambda_temp);
    }
  }
  arma::mat Lambda_temp = as<arma::mat>(Lambda);
  Lambda_temp.print();
  Rcpp::StringMatrix Topic(topic_length,K) ;
 for(int i=0;i<K;i++)
 {
   arma::vec topic_i = Lambda_temp.col(i);
   arma::uvec topic_index = arma::sort_index(topic_i,1);

   std::cout<<std::endl<<i<<"  "<<std::endl;
   for(int j=0;j<topic_length;j++)
   {
       std::cout<<" " <<topic_index[j]<<" "<<Dictionary[topic_index[j]];
       Topic(j,i)=Dictionary[topic_index[j]];
   }
 }

  return Rcpp::List::create(Lambda,phi,gamma,Topic,Dictionary);
}
/*** R
*/
