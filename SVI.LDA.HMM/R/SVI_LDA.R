# LDA.SVI<-function(X, # the Input data
#                   K, # the number of the topics
#                   n, # the number of the passes
#                   alpha, # the parameter of theta
#                   eta, # the parameter of beta
#                   pre, # the stop criteria
#                   topic_length# the number of top words and topics
#
# )
# {
#        library(Rcpp)
#        library(RcppArmadillo)
#        sourceCpp(file='src/LDA.cpp')
#        t=  SVI_LDA(X,K,n,alpha,eta,pre,topic_length)
#        names(t)=c("Lan","LSS","LL","SSS","SSSSSS")
#        return(t)
# }
#
# Tokenize<-function(Data_set, # the Input original data, which is a list data type
#                    Language = "en"# the language type you want to produce, defaut language is en
#                    )
# {
#      library(tokenizers);
#      tokenized_data <-tokenize_words(Data_set, stopwords = stopwords(Language));
#      return(tokenized_data);
# }
# Stemming <- function(Tokenized_data)
# {
#   library(tokenizers);
#    Document_num = length(Tokenized_data);
#    Stem_data <- list(Document_num);
#    for( i in 1:length(Tokenized_data))
#    {
#       Stem_data[[i]]=wordStem(unlist(Tokenized_data[i]));
#    }
#    return (Stem_data);
# }
#
#  X=list(c('home','mother','mother','house','mother',rep("cat",3)),c('father','dad','cat','dad','cat','cat','cat','cat','cat','cat'))
# K=3;
# alpha=1/K;
# eta=0.01;
# n=100;
# pre=0.1;
# t=LDA.SVI(X,K,n,alpha=1/K,eta=0.01,pre=0.01,topic_length= 3)
# doc_a = "Brocolli is good to eat. My brother likes to eat good brocolli, but not my mother."
# doc_b = "My mother spends a lot of time driving my brother around to baseball practice."
# doc_c = "Some health experts suggest that driving may cause increased tension and blood pressure."
# doc_d = "I often feel pressure to perform well at school, but my mother never seems to drive my brother to do better."
# doc_e = "Health professionals say that brocolli is good for your health."
# data <- list(doc_a,doc_b,doc_c,doc_d,doc_e)

Forward <- function( observed_data,#Observed data
                     pi_hat,      # initial distribution
                     A_, # transition matrix
                     P_,
                     s_begin,
                     sublength

)
{


}


SVI.HMM <- function( Data_sequence, # Observed data
                     sub_length =2, # subchain length
                     K,#Hidden_state_number
                     alpha, # initial w_A parameter
                     u0, # initial W_phi
                     k0,
                     sigma0,
                     v0,
                     U_A, # Hyperparameter
                     U_phi1, # Hyperparameter A
                     U_phi2,
                     U_phi3,
                     U_phi4 ,
                     Pre, # stop cre
                     Pass
)
{
   Dim=length(Data_sequence[[1]])# data_dimention
   u = matrix(rep(0,Dim*K),Dim) ;
   k = rep(k0,K);
   v = rep(v0,K);
   for( i in 1:K)
   {
     u[,i]=u0;
   }
   sigma = list(K)
   for(i in 1:K)
   {
     sigma[[i]]= sigma0;
   }
   T = length(Data_sequence)
   W_A = matrix(rep(alpha,K*K),K);
   P_ = matrix(rep(1/K,T*K),T);
   q_x = matrix(rep(1/K,T*K),T);
   alpha_forward=matrix(rep(1/K,T*K),T);
   beta_backward= matrix(rep(1/K,T*K),T);
   q_transition = array(rep(0,T*K*K),c(T,K,K)) # index transition matrix according to beta_index
   for(i in 1:Pass )
   {
     for(j in 1:(T-sub_length+1))
     {
       A_ = exp(digamma(W_A)-digamma(rowSums(W_A)));#compute A_ to use forward_backward algorithm
       for(m in 1:K)
       {
         for(t in j:(j+sub_length-1))
         {
             P_[t,m]= exp(-0.5*logb(2*pi)-0.5*(sum(digamma(0.5*(v[m]+1-1:Dim))))
                          +Dim*logb(2)-logb(det(sigma[[m]]))-0.5*Dim/k[m]
                          -0.5*v[m]*t(Data_sequence[[t]]-u[,m])%*%solve(sigma[[m]])%*%(Data_sequence[[t]]-u[,m]));
         }
       }
       Eq_A = W_A/rowSums(W_A)#Get the expectation of A for q(A)
       eigen_value= eigen(Eq_A)
       pi_hat = eigen_value$vectors[,1]
       for(tag in j:(j+sub_length-1))
       {
           if(tag==j)
           {
             for(state in 1:K)
             {
               alpha_forward[tag,state]=crossprod(pi_hat,A_[,state])*P_[tag,state]
             }
           }
           else
           {
             for(state in 1:K)
             {
               alpha_forward[tag,state]=sum(alpha_forward[tag-1,]*A_[,state])*P_[tag,state]
             }

           }
       }
       for(tag in (j+sub_length-1):j)
       {
            if(tag==j+sub_length-1)
            {
              beta_backward [tag,]=1
            }
            else
            {
              for(state in 1:K)
              {
                beta_backward[tag,state]=sum(A_[state,]*beta_backward[tag+1,]*P_[tag+1,])
              }
            }
       }
       beta_original = rep(0,K)
       for(state in 1:K)
       {
         beta_original[state]=sum(A_[state,]*beta_backward[j,]*P_[j,])
       }
       for(tag in j:(j+sub_length-1))
       {
         q_x[tag,]=alpha_forward[tag,]*beta_backward[tag,]
         q_x[tag,]= q_x[tag,]/sum(q_x[tag,])
       }
       for(tag in j:(j+sub_length-1))
       {
         if(tag==j)
         {
           for(j_local in 1:K)
           {
             for(k_local in 1:K)
             {
               q_transition[tag,j_local,k_local]=pi_hat[j_local]*A_
             }
           }

         }
       }


       print(A_)
       print(P_)
       print(pi_hat)
       print(alpha_forward)
       print(beta_backward)
       print(beta_original)
     }

   }

}
data <-list(c(1:3),c(3:1),c(1:3),c(3:1),c(3:1),c(3:1))
K=2
Dim=3
alpha=2
u0 = rep(2,3)
k0=1
v0=6
sigma0 = diag(rep(1,3))
U_A = 2
U_phi1 = k0*u0
U_phi2 = k0
U_phi3 = sigma0+k0*u0%*%t(u0)
U_phi4 = v0+2+3
result= SVI.HMM(data,sub_length = 6,K=2,alpha, u0, k0, sigma0, v0, U_A, U_phi1,U_phi2, U_phi3, U_phi4, Pre=0.2,Pass=1)





