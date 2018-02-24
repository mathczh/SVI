#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <string>
#include <map>
#include<Rinternals.h>
//#include<lapack.h>
//#include<blas.h>
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
  for(int i= 0;i< X.size();i++) //Produce map to stroe the whole vocabulary
  {
    Rcpp::StringVector document = Rcpp::as<Rcpp::StringVector>(X[i]);
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
  arma::mat Lambda_temp_batch(V,K);
  Lambda_temp_batch.fill(0);
  for( int i=0;i<V;i++)//initiate Lambda with rchisq
  {
    for(int j=0;j<K;j++)
    {
      Lambda(i,j)=R::rchisq(100);
    }
  }
  Rcpp::NumericMatrix gamma(K,D) ; //random initial with all 0 entry
  gamma.fill(1);
  Rcpp::List phi(D);
  Rcpp::NumericMatrix Lambda_old(V,K);
  Rcpp::List index(D);
  for(int i_=0;i_<D;i_++){
    Rcpp::StringVector Y = Rcpp::as<Rcpp::StringVector>(X[i_]) ;
    Rcpp::IntegerVector index_local(Y.size());
    for(int j=0;j<Y.size();j++)
      {
        iter =mymap.find(""+Y[j]+"");
        index_local[j]  = iter->second;
      }
      index[i_] = index_local;
  }
  int global_iteration = 0;
  int flag = 0;
//  int flag_print = 0;
  int batch_count = 0;
  while( global_iteration < max_iter_global && flag ==0){
  for(int d=0;d<D;d++){
	   if ((sum(abs(Lambda-Lambda_old))/((sum(abs(Lambda_old)))+0.000001) <= tol_global || global_iteration > max_iter_global)&&batch_count==0){
		  flag = 1;
		  break;
	   }
	    global_iteration++;
	    Lambda_old = Lambda;
      for(int v=0;v<K;v++)
      {
        gamma(v,d)=1;
      }
      Rcpp::NumericVector gamma1(K);
      Rcpp::StringVector Y = Rcpp::as<Rcpp::StringVector>(X[d]) ;
      Rcpp::NumericMatrix w(V,Y.size()) ;
      Rcpp::NumericMatrix mid(Y.size(),K);
      arma::mat mid_test(Y.size(),K);
      Rcpp::NumericMatrix mid_norm(Y.size(),K);
      mid_norm.fill(1);
      Rcpp::NumericMatrix mid_norm_old(Y.size(),K);
      Rcpp::IntegerVector index_local = index[d];
      int iteration_local =0;
      while((sum(abs(gamma(_,d)-gamma1))/((sum(abs(gamma1)))+0.000001)>=tol_local
               ||sum(abs(mid_norm-mid_norm_old))/(sum(abs(mid_norm_old))+0.000001)>=tol_local)&&iteration_local<max_iter_local)
      {
        iteration_local++;
        mid_norm_old = mid_norm;
        gamma1=gamma(_,d);
        for(int j=0;j<Y.size();j++)
        {
          Rcpp::NumericVector w_temp(V);
          int index_inmap = index_local[j];
          w_temp[index_inmap]=1;
          w(_,j)= w_temp;
          mid(j,_)=exp(digamma(gamma(_,d))+digamma(Lambda(index_inmap,_))-digamma(colSums(Lambda)));
          mid_test = as<arma::mat>(mid);
          arma::mat Lambda_temp = as<arma::mat>(Lambda);
          arma::vec colsums_lambda = as<arma::vec>(colSums(Lambda));
          arma::mat gamma_d = as<arma::mat>(gamma);
          // if(mid_test.has_nan()){
          //   cout<<"mid_test.row(j).has_nan()"<<endl;
          // }
          mid_norm(j,_)=mid(j,_)/sum(mid(j,_));//normalize parameter
          arma::mat midnorm_temp = as<arma::mat>(mid_norm);//convert NumericMtrix into arma::mat;
          // if(midnorm_temp.has_nan()&&flag_print==0){
          //   cout<<"mid(j,_) "<<endl;
          //   mid_test.row(j).print();
          //   cout<<" Lambda_temp.row(index_inmap) "<<endl;
          //   Lambda_temp.row(index_inmap).print();
          //   cout<< "colSums(Lambda): " <<std::endl;
          //   colsums_lambda.print();
          //   cout<< "gamma(_,d): " <<std::endl;
          //   gamma_d.col(d).print();
          //   flag_print = 1;
          // }
        }
        gamma(_,d)= alpha+colSums(mid_norm);
      }
      phi[d]= mid_norm;
      arma::mat midnorm_temp = as<arma::mat>(mid_norm);//convert NumericMtrix into arma::mat;
      arma::mat w_temp = as<arma::mat>(w);
      // if(w_temp.has_nan()){
      //   cout<<"w_temp.has_nan"<<endl;
      // }
      // if(midnorm_temp.has_nan()){
      //   cout<<"midnorm_temp.has_nan() 2"<<endl;
      // }
      arma::mat result_temp=eta+D*w_temp*midnorm_temp;
      batch_count++;
      Lambda_temp_batch = Lambda_temp_batch + result_temp;
      // if(result_temp.has_nan()){
      //   cout<<"result_temp.has_nan"<<endl;
      // }
      // if(Lambda_temp_batch.has_nan()){
      //   cout<<"Lambda_temp_batch.has_nan"<<endl;
      // }
      if(batch_count==batch_size){
        arma::mat Lambda_temp = as<arma::mat>(Lambda);
        // if(Lambda_temp.has_nan()){
        //   cout<<"Lambda.has_nan() "<<endl;
        //   flag = 1;
        //   break;
        // }
        Lambda_temp = (1-1.0/(global_iteration*D+d+1))*Lambda_temp+1.0/(batch_size*(global_iteration*D+d+1))*Lambda_temp_batch;
        // if(Lambda_temp.has_nan()){
        //   cout<<"Lambda_temp.has_nan() "<<endl;
        //   flag =1;
        //   break;
        // }
        Lambda=Rcpp::wrap(Lambda_temp);
        // for(int i=0; i < K; i++){
        //   Lambda(_,i) = Lambda (_,i)/sum(Lambda(_,i));
        // }
        Lambda_temp_batch.fill(0);
        batch_count=0;
      }
    }
  }
  Rcpp::NumericMatrix E_beta (V,K);
  Rcpp::NumericMatrix E_theta (K,D);
  for(int i=0; i < K; i++){
    E_beta(_,i) = Lambda (_,i)/sum(Lambda(_,i));
  }
  for(int i=0;i<D;i++){
    E_theta(_,i) = gamma(_,i)/sum(gamma(_,i));
  }
  arma::mat E_beta_temp = as<arma::mat>(E_beta);
  arma::mat E_theta_temp = as<arma::mat>(E_theta);
  Rcpp::NumericVector Each_doc_likelihood (D);
  for(int i_=0;i_< index.size();i_++){
    Rcpp::IntegerVector index_local = index[i_];
    double mysum=0;
    for(int j=0;j<index_local.size();j++)
      {
       mysum=mysum+ log(sum(E_beta(index_local[j],_)*E_theta(_,i_)));
      }
      Each_doc_likelihood[i_] = mysum;
      log_likelihood=log_likelihood+mysum;
  }
  Rcpp::StringMatrix Topic(word_length,K);
  for(int i=0;i<K;i++)
  {
   // E_beta_temp.print();
    arma::vec topic_i = E_beta_temp.col(i);
    arma::uvec topic_index = arma::sort_index(topic_i,"descend");
    for(int j=0;j<word_length;j++)
    {
      Topic(j,i)=Dictionary[topic_index[j]];
    }
  }
  return Rcpp::List::create(Lambda,phi,gamma,Topic,Dictionary,global_iteration,E_beta,E_theta,log_likelihood,Each_doc_likelihood);
}
/*** R
*/
