// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>
#include <string>
#include <map>
#include<Rinternals.h>
#include<math.h>
using namespace std;
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
SEXP SVI_HMM(  Rcpp:: List &Data_sequence,
               SEXP sub_lengthh,
               SEXP  KK,
               SEXP W_AA,
               SEXP uu,
               SEXP kk,
               Rcpp:: List &sigma,
               SEXP vv,
               SEXP U_AA,
               SEXP U_phi11,
               SEXP U_phi22,
               Rcpp:: List &U_phi3,
               SEXP U_phi44,
               SEXP Dimm,
               SEXP Passs,
               SEXP pree

)
{

     int sub_length=Rcpp::as<int>(sub_lengthh);
     int Dim=Rcpp::as<int>(Dimm);
     int K=Rcpp::as<int>(KK);
     Rcpp:: NumericMatrix W_A =  Rcpp::as<Rcpp::NumericMatrix>(W_AA) ;
     Rcpp:: NumericMatrix  u  =  Rcpp::as<Rcpp::NumericMatrix>(uu) ;
     Rcpp:: NumericVector  k  =  Rcpp::as<Rcpp::NumericVector>(kk) ;
     Rcpp:: NumericVector  v  =  Rcpp::as<Rcpp::NumericVector>(vv) ;
     Rcpp:: NumericVector  U_A =  Rcpp::as<Rcpp::NumericVector>(U_AA) ;
     Rcpp:: NumericMatrix  U_phi1 =  Rcpp::as<Rcpp::NumericMatrix>(U_phi11) ;
     Rcpp:: NumericVector  U_phi2 =  Rcpp::as<Rcpp::NumericVector>(U_phi22) ;
     Rcpp:: NumericVector  U_phi4 =  Rcpp::as<Rcpp::NumericVector>(U_phi44) ;
     int Pass=Rcpp::as<int>(Passs);
     double pre =Rcpp::as<double>(pree);
     int T = Data_sequence.size();
     Rcpp::NumericMatrix P_(T,K);
     Rcpp::NumericMatrix q_x(T,K);
     Rcpp::NumericMatrix alpha_forward(T,K);
     Rcpp::NumericMatrix beta_backward(T,K);
     Rcpp::List q_transition(T);
     Rcpp::NumericMatrix A_(K,K);
     int i =0;
     arma::mat W_A_new = as<arma::mat>(W_A);
     arma::mat W_A_old(K,K);
     while(norm(W_A_new-W_A_old,"fro")/norm(W_A_new,"fro")>=pre&&i<Pass)
     {
       W_A_old = W_A_new;
       for(int j=0;j<(T-sub_length+1);j++)
       {
         for( int i_local=0;i_local<K;i_local++)
         {
           A_(_,i_local) = exp(digamma(W_A(_,i_local))-digamma(rowSums(W_A)));
         }
         for(int m =0; m<K;m++)
         {
           for(int t=j;t<=(j+sub_length-1);t++)
           {
             double pi = 3.141593;
             Rcpp:: NumericVector  Dim_vec(Dim);
             for(int i_local =0;i_local<Dim;i_local++)
             {
               Dim_vec[i_local]=0.5*(v[m]+1-(i_local+1));
             }
             arma::mat sigma_m = as<arma::mat>(sigma[m]);
             Rcpp:: NumericVector Data_t = Rcpp::as<Rcpp::NumericVector>(Data_sequence[t]) ;
              Rcpp:: NumericVector Data_sequence_t_u_m_temp = Data_t - u(_,m);
              arma::vec Data_sequence_t_u_m = as<arma::vec>(Data_sequence_t_u_m_temp);
              P_(t,m)= exp(-0.5*log(2*pi)-0.5*(sum(digamma(Dim_vec)))
                             +Dim*log(2)-log(det(sigma_m))-0.5*Dim/k[m]
                             -0.5*v[m]*sum(trans(Data_sequence_t_u_m)*inv(sigma_m)*(Data_sequence_t_u_m)));
           }
         }
           arma::mat Eq_A(K,K);
           for(int i_local=0;i_local<K;i_local++)
           {
              double my_sum = sum(W_A(i_local,_));
             for(int j_local=0;j_local<K;j_local++)
             {
               Eq_A(i_local,j_local)=W_A(i_local,j_local)/my_sum;
             }
           }
           arma::cx_vec eigval;
           arma::cx_mat eigvec;
           eig_gen(eigval,eigvec,Eq_A);
           Rcpp:: NumericVector pi_hat= Rcpp::wrap(real(eigvec.col(0))) ;
           for(int tag = j ;tag<=(j+sub_length-1);tag++)
           {
             if(tag==j)
             {
               for(int state=0; state<K;state++)
               {
                 alpha_forward(tag,state)=sum(pi_hat*A_(_,state))*P_(tag,state);
               }
               alpha_forward(tag,_)=alpha_forward(tag,_)/sum(alpha_forward(tag,_));
             }
             else
             {
               for(int state=0;state<K;state++)
               {
                 alpha_forward(tag,state)=sum(alpha_forward(tag-1,_)*A_(_,state))*P_(tag,state);
               }
               alpha_forward(tag,_)=alpha_forward(tag,_)/sum(alpha_forward(tag,_));
             }
           }
           for(int tag=(j+sub_length-1);tag>=j;tag--)
           {
             if(tag==(j+sub_length-1))
             {
               for(int state=0;state<K;state++)
               {
                 beta_backward(tag,state)=1;
               }
             }
             else
             {
               for(int state=0;state<K;state++)
               {
                 beta_backward(tag,state)= sum(A_(state,_)*beta_backward(tag+1,_)*P_(tag+1,_));
               }
                 beta_backward(tag,_)= beta_backward(tag,_)/sum(beta_backward(tag,_));
             }
           }
           for(int tag=j;tag<=(j+sub_length-1);tag++)
           {
             q_x(tag,_) = alpha_forward(tag,_)*beta_backward(tag,_);
             q_x(tag,_) = q_x(tag,_)/sum(q_x(tag,_));
           }
           for(int tag=j;tag<=(j+sub_length-1);tag++)
           {
             if(tag==j)
             {
               Rcpp::NumericMatrix temp_matrix(K,K);
               for(int j_local=0;j_local<K;j_local++)
               {
                 for(int k_local=0;k_local<K;k_local++)
                 {
                   temp_matrix(j_local,k_local)=pi_hat[j_local]*A_(j_local,k_local)*P_(tag,k_local)*beta_backward(tag,k_local);
                 }
               }
                temp_matrix=temp_matrix/sum(temp_matrix);
                q_transition(tag)=temp_matrix;
             }
             else
             {
                Rcpp::NumericMatrix temp_matrix(K,K);
               for(int j_local=0;j_local<K;j_local++)
               {
                 for(int k_local=0;k_local<K;k_local++)
                 {
                   temp_matrix(j_local,k_local)=alpha_forward(tag-1,j_local)*A_(j_local,k_local)*P_(tag,k_local)*beta_backward(tag,k_local);
                 }
               }
                temp_matrix=temp_matrix/sum(temp_matrix);
                q_transition(tag)=temp_matrix;
             }
           }
             arma::mat q_tran_sum(K,K);
             q_tran_sum.fill(0);
            for(int tag=(j+1);tag<=(j+sub_length-1);tag++)
           {
              arma::mat q_transition_tag = as<arma::mat>(q_transition[tag]);
              q_tran_sum=q_tran_sum+q_transition_tag;
           }
            Rcpp:: NumericMatrix q_x_sum1(Dim,K);
            q_x_sum1.fill(0);
           for(int k_local=0;k_local<K;k_local++)
           {
             for(int tag=j;tag<=(j+sub_length-1);tag++)
             {
               Rcpp:: NumericVector Data_t = Rcpp::as<Rcpp::NumericVector>(Data_sequence[tag]) ;
               q_x_sum1(_,k_local)=q_x_sum1(_,k_local)+q_x(tag,k_local)*Data_t;
             }
           }
            Rcpp:: NumericVector q_x_sum2(K) ;
           for(int tag=j;tag<=(j+sub_length-1);tag++)
           {
             q_x_sum2 = q_x_sum2+q_x(tag,_);
           }
           Rcpp:: NumericVector q_x_sum4 = q_x_sum2;
           Rcpp::List q_x_sum3(K);
           for(int k_local=0;k_local<K;k_local++)
           {
             arma::mat temp_matrix(Dim,Dim);
             temp_matrix.fill(0);
             for(int tag=j;tag<=(j+sub_length-1);tag++)
             {
               Rcpp:: NumericVector Data_t = Rcpp::as<Rcpp::NumericVector>(Data_sequence[tag]) ;
               arma::vec Data_sequence_t =  as<arma::vec>(Data_t);
               temp_matrix = temp_matrix+q_x(tag,k_local)*Data_sequence_t*trans(Data_sequence_t);
             }
             q_x_sum3(k_local) = Rcpp::wrap(temp_matrix) ;
           }
          // return Rcpp::List::create(q_tran_sum,q_x_sum1,q_x_sum2,q_x_sum3,q_x_sum4);
         // return Rcpp::List::create(W_A,U_A,u,U_phi1,k,U_phi2,sigma,U_phi3,v,U_phi4);
           arma::mat W_A_temp = as<arma::mat>(W_A);
           arma::mat U_A_temp = as<arma::mat>(U_A);
           W_A_temp = (1-(1.0/(i*(T-sub_length+1)+j+1)))*W_A_temp+(1.0/(i*(T-sub_length+1)+j+1))*(((T-sub_length+1.0)/(sub_length-1))*q_tran_sum+U_A_temp);
           W_A= Rcpp::wrap(W_A_temp);
           arma::mat u_temp = as<arma::mat>(u);
           arma::mat q_x_sum_temp1 = as<arma::mat>(q_x_sum1);
           arma::mat u_phi_te = as<arma::mat>(U_phi1);
           u_temp = (1-(1.0/(i*(T-sub_length+1)+j+1)))*u_temp+(1.0/(i*(T-sub_length+1)+j+1))*((T-sub_length+1.0)/sub_length*q_x_sum_temp1+u_phi_te);
           u = Rcpp::wrap(u_temp);
           k = (1-(1.0/(i*(T-sub_length+1)+j+1)))*k+(1.0/(i*(T-sub_length+1)+j+1))*((T-sub_length+1.0)/sub_length*q_x_sum2+U_phi2);
              for(int i_local=0;i_local<K;i_local++)
                  {
                   Rcpp:: NumericMatrix sigma_i = Rcpp::as<Rcpp::NumericMatrix>(sigma[i_local]) ;
                   arma::mat sigma_i_temp =  as<arma::mat>(sigma_i);
                   Rcpp:: NumericMatrix q_x_sum_temp3 = Rcpp::as<Rcpp::NumericMatrix>(q_x_sum3[i_local]) ;
                   arma::mat q_x_sum_i =  as<arma::mat>(q_x_sum_temp3);
                   Rcpp:: NumericMatrix U_phi_temp3 = Rcpp::as<Rcpp::NumericMatrix>(U_phi3[i_local]) ;
                   arma::mat U_phi_i =  as<arma::mat>(U_phi_temp3);
                   sigma_i_temp=(1-(1.0/(i*(T-sub_length+1)+j+1)))*sigma_i_temp+(1.0/(i*(T-sub_length+1)+j+1))*((T-sub_length+1.0)/sub_length*q_x_sum_i+U_phi_i);
                   sigma[i_local] = Rcpp::wrap(sigma_i_temp);
                 }
           v = (1-(1.0/(i*(T-sub_length+1)+j+1)))*v+(1.0/(i*(T-sub_length+1)+j+1))*((T-sub_length+1.0)/sub_length*q_x_sum4+U_phi4);
           }

       W_A_new = as<arma::mat>(W_A);
          i++;
     }
     return Rcpp::List::create(W_A,u,k,sigma,v);
}
/*** R
*/
