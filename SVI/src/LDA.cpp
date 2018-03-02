#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <string>
#include <map>
#include<Rinternals.h>
#include<math.h>
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]
SEXP SVI_LDA(Rcpp:: List &X,
             SEXP  KK,
             SEXP  max_iter_locall,
             SEXP  max_iter_globall,
             SEXP  alphaa,
             SEXP  etaa,
             SEXP  tol_locall,
             SEXP  tol_globall,
             SEXP  batch_sizee,
             SEXP  word_num
)
{
  int K = Rcpp::as<int>(KK);
  int word_length = Rcpp::as<int>(word_num);
  int batch_size = Rcpp::as<int>(batch_sizee);
  int max_iter_global = Rcpp::as<int>(max_iter_globall);
  int max_iter_local = Rcpp::as<int>(max_iter_locall);
  double alpha = Rcpp::as<double>(alphaa);
  double eta = Rcpp::as<double>(etaa);
  double tol_local = Rcpp::as<double>(tol_locall);
  double tol_global = Rcpp::as<double>(tol_globall);
  double log_likelihood = 0;
  std::map<std::string, int> mymap;
  std::map<std::string, int>::iterator iter;
  Rcpp::StringVector document ;
  for(int i= 0;i< X.size();i++) //Produce map to stroe the whole vocabulary
  {
    document = Rcpp::as<Rcpp::StringVector>(X[i]);
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
  arma::mat Lambda(V,K);
  Rcpp::NumericMatrix Lambda_temp;
  arma::mat Lambda_temp_batch(V,K);
  Lambda_temp_batch.fill(0);
  for( int i=0;i<K;i++)//initiate Lambda with rchisq
  {
      Lambda.col(i)= as<arma::vec>(rchisq(V,100));
  }
  arma::mat gamma(K,D) ; //random initial with all 0 entry
  gamma.fill(1);
  List phi(D);
  arma::mat Lambda_old(V,K);
  Lambda_old.fill(0);
  List index(D);
  Rcpp::StringVector Y;
  Rcpp::IntegerVector index_local(Y.size());
  for(int i_=0;i_<D;i_++){
    Y = Rcpp::as<Rcpp::StringVector>(X[i_]) ;
    index_local.fill(0);
    for(int j=0;j<Y.size();j++)
      {
        iter =mymap.find(""+Y[j]+"");
        index_local[j]  = iter->second;
      }
      index[i_] = index_local;
  }
  int global_iteration = 0;
  int iteration_local =0;
  int flag = 0;
  int batch_count = 0;
  while( global_iteration < max_iter_global && flag ==0){
  for(int d=0;d<D;d++){
	   if (((Rcpp::wrap( sum(abs(Lambda-Lambda_old))/((sum(abs(Lambda_old)))+0.000001) <= tol_global))|| global_iteration > max_iter_global ) &&batch_count==0){
		  flag = 1;
		  break;
	   }
	    global_iteration++;
	    Lambda_old = Lambda;
      gamma.col(d) = 1;
      arma::vec gamma1(K);
      Y = Rcpp::as<Rcpp::StringVector>(X[d]) ;
      arma::mat w(V,Y.size()) ;
      Rcpp::NumericMatrix mid(Y.size(),K);
      arma::mat mid_norm(Y.size(),K);
      mid_norm.fill(1);
      arma::mat mid_norm_old(Y.size(),K);
      index_local = index[d];
      iteration_local =0;
      while((norm(gamma.col(d)-gamma1,2)/(norm(gamma1,2)+0.00000001)>=tol_local || norm(mid_norm-mid_norm_old,2)/(norm(mid_norm_old,2)+0.00000001)>=tol_local)&&iteration_local<max_iter_local)
      {
        iteration_local++;
        mid_norm_old = mid_norm;
        gamma1=gamma.col(d);
        arma::vec w_temp(V);
        w_temp.fill(0);
        Lambda_temp = Rcpp::wrap(Lambda);
        Rcpp::NumericVector gamma_temp1 = Rcpp::wrap(gamma1);
        for(int j=0;j<Y.size();j++)
        {
          int index_inmap = index_local[j];
          w_temp[index_inmap]=1;
          w.col(j)= w_temp;
          mid(j,_)=exp(digamma(gamma_temp1)+digamma(Lambda_temp(index_inmap,_))-digamma(colSums(Lambda_temp)));
          mid(j,_)= mid(j,_)/sum(mid(j,_));
        }
        mid_norm = as<arma::mat>(mid); // mid_norm is actually the arma data type of normalized mid
        Rcpp::NumericVector gamma_d = alpha+colSums(mid);
        gamma.col(d)= as<arma::vec>(gamma_d);
        std::std::cout << "Gamma's value" << '\n';
        gamma.col(d).print();
      }
      phi[d]= mid_norm;
      batch_count++;
      Lambda_temp_batch = Lambda_temp_batch + eta+D*w*mid_norm;
      if(batch_count==batch_size){
        double lr = pow(10+global_iteration,-0.7);
        Lambda = (1-lr)*Lambda+lr/batch_size*Lambda_temp_batch;
        Lambda_temp_batch.fill(0);
        batch_count=0;
      }
    }
  }
  std::cout<<"Print without wrong 1 "<<std::endl;
  arma::mat E_beta (V,K);
  arma::mat E_theta (K,D);
  for(int i=0; i < K; i++){
    E_beta.col(i)= Lambda.col(i)/sum(Lambda.col(i));
  }
  E_beta.print();
  for(int i=0;i<D;i++){
    E_theta.col(i) = gamma.col(i)/sum(gamma.col(i));
  }
  E_theta.print();
  std::cout<<"Print without wrong 2"<<std::endl;
  Rcpp::NumericVector Each_doc_likelihood (D);
  for(int i_=0;i_< index.size();i_++){
    index_local = index[i_];
    double mysum=0;
    for(int j=0;j<index_local.size();j++)
      {
       mysum=mysum+ log(sum(E_beta.row(index_local[j])*E_theta.col(i_)));
      }
      Each_doc_likelihood[i_] = mysum;
      log_likelihood=log_likelihood+mysum;
  }
  std::cout<<"Print without wrong 3 "<<std::endl;
  Rcpp::StringMatrix Topic(word_length,K);
  arma::uvec topic_index;
  for(int i=0;i<K;i++)
  {
    topic_index = arma::sort_index(E_beta.col(i),"descend");
    for(int j=0;j<word_length;j++)
    {
      Topic(j,i)=Dictionary[topic_index[j]];
    }
  }
  std::cout<<"Print without wrong 4 "<<std::endl;
  std::cout << log_likelihood <<std::endl;
  return Rcpp::List::create(Lambda,phi,gamma,Topic,Dictionary,global_iteration,log_likelihood,Each_doc_likelihood);
}
/*** R
*/
