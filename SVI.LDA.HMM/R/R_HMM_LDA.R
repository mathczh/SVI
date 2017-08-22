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
SVI.HMM <- function( Data_sequence, # Observed data
                     sub_length , # subchain length
                     K,             #Hidden_state_number
                     W_A,         # initial parameter,which controls transition matrix
                     u,            # initial mean value
                     k,            # initial parameter, which controls Gassuian Distribution
                     sigma,        # initial covariance matrix
                     v,             # initial parameter, which controls Gassuian Distribution
                     U_A,    # Hyperparameter
                     U_phi1, # Hyperparameter
                     U_phi2, # Hyperperameter the same notation with the paper
                     U_phi3, # Hyperperameter the same notation with the paper
                     U_phi4, # Hyperperameter the same notation with the paper
                     Pass,    # toal pass
                     pre
)
{
  Dim=length(Data_sequence[[1]])  # data_dimention
  T = length(Data_sequence)    # data squence length T
  P_ = matrix(rep(1/K,T*K),T);   #  Have the same meaning with paper in S10
  q_x = matrix(rep(1/K,T*K),T);  # transition matrix, In S14 of the paper
  alpha_forward=matrix(rep(1/K,T*K),T); # the same meaning with that in forward and backward algorithm
  beta_backward= matrix(rep(1/K,T*K),T);
  q_transition = array(rep(0,T*K*K),c(T,K,K)) # Equal to S15 in the paper
  i = 0;
  W_A_old=0
  while(norm(W_A-W_A_old,"F")/norm(W_A)>=pre&&i<Pass)
  {
    W_A_old=W_A
    i=i+1
    for(j in 1:(T-sub_length+1))  # sample each subchain
    {
      A_ = exp(digamma(W_A)-digamma(rowSums(W_A)));# the same with S11
      trick_sum = rep(0,K)
        for(t in j:(j+sub_length-1))
        {

          for(m in 1:K)
          {
              trick_sum[m]=(-0.5*Dim*logb(pi)+0.5*(sum(digamma(0.5*(v[m]+1-1:Dim))))-0.5*logb(det(sigma[[m]]))-0.5*Dim/k[m]-0.5*v[m]*t(Data_sequence[[t]]-u[,m])%*%solve(sigma[[m]])%*%(Data_sequence[[t]]-u[,m]));
          }
          shift_center = max(trick_sum)
          for(m in 1:K)
          {
            P_[t,m]= exp(trick_sum[m]-shift_center-log(sum(exp(trick_sum-shift_center))));
          }
         }

       Eq_A = W_A/rowSums(W_A)#Get the expectation of A with q(A) density function
      eigen_value= eigen(Eq_A)
      pi_hat = eigen_value$vectors[,1] # corresponding to biggest eigen_value, initiate original distribution in forward algorithm
       for(tag in j:(j+sub_length-1))  # tag means time
      {
        if(tag==j)
        {
          for(state in 1:K)
          {
            alpha_forward[tag,state]=crossprod(pi_hat,A_[,state])*P_[tag,state] # S12 in the paper
          }
          alpha_forward[tag,]= alpha_forward[tag,]/sum(alpha_forward[tag,])
        }
        else
        {
          for(state in 1:K)
          {
            alpha_forward[tag,state]=sum(alpha_forward[tag-1,]*A_[,state])*P_[tag,state]
            # recursion
          }
          alpha_forward[tag,]= alpha_forward[tag,]/sum(alpha_forward[tag,])
        }
       }
     # write.table(alpha_forward,"alpha_forward.txt",append = T)
      for(tag in (j+sub_length-1):j)
      {
        if(tag==j+sub_length-1)
        {
          beta_backward [tag,]=1 # initiate
        }
        else
        {
          for(state in 1:K)
          {
            beta_backward[tag,state]=sum(A_[state,]*beta_backward[tag+1,]*P_[tag+1,]) #S13 ,but I think paper is wrong,and please check again
          }
          beta_backward[tag,]=beta_backward[tag,]/sum(beta_backward[tag,])
        }
      }
    #  write.table(beta_backward,"beta.txt",append = T)
      beta_original = rep(0,K) # this corresponds to beta_0, even it is no use for following caculation
      for(state in 1:K)
      {
        beta_original[state]=sum(A_[state,]*beta_backward[j,]*P_[j,])
      }
      for(tag in j:(j+sub_length-1))  # S14 in this paper, we only need time from 1, that is from j in subchain situation.
      {                              #Please check my idea accuracy, I think we don't need calculate t=0
        q_x[tag,]=alpha_forward[tag,]*beta_backward[tag,]
        q_x[tag,]= q_x[tag,]/sum(q_x[tag,]) # normalize
      }
    #  write.table(q_x,"q_x.txt",append = T)
    #  write.table(A_,"A_.txt",append = T)
    #  write.table(P_,"P_.txt",append = T)
      for(tag in j:(j+sub_length-1))
      {
        if(tag==j) # we need pi_hat to calculate S15 in the first iteration
        {
          for(j_local in 1:K)
          {
            for(k_local in 1:K)
            { #q_transition according to index of Xt. And S15 is wrong in my opinion. I substitute P,A with P_,A_. Check please
              q_transition[tag,j_local,k_local]=pi_hat[j_local]*A_[j_local,k_local]*P_[tag,k_local]*beta_backward[tag,k_local];
            }
          }
          q_transition[tag,,]=q_transition[tag,,]/sum(q_transition[tag,,])

        }
        else
        {
          for(j_local in 1:K)
          {
            for(k_local in 1:K)
            {
              q_transition[tag,j_local,k_local]=alpha_forward[tag-1,j_local]*A_[j_local,k_local]*P_[tag,k_local]*beta_backward[tag,k_local];
            }
          }
          q_transition[tag,,]=q_transition[tag,,]/sum(q_transition[tag,,]) # normalize
        }
      }
      q_tran_sum =  matrix(rep(0,K*K),K); #To calculate the first sum in S8 to update global phi_A
      for(tag in (j+1):(j+sub_length-1))
      {
        q_tran_sum=q_tran_sum+q_transition[tag,,];
      }
      q_x_sum1 = matrix(rep(0,K*Dim),Dim)#To calculate the second sum in S8, and paper is wrong in notation
      for(k_local in 1:K)
      {
        for(tag in j:(j+sub_length-1))
        {
          q_x_sum1[,k_local]=q_x_sum1[,k_local]+q_x[tag,k_local]*Data_sequence[[tag]]
        }
      }
      q_x_sum2 = rep(0,K);
      for(tag in j:(j+sub_length-1)) #To calculate the third sum in S8
      {
        q_x_sum2=q_x_sum2+q_x[tag,]
      }
      q_x_sum4 = q_x_sum2
      q_x_sum3 = array(rep(0,K*Dim*Dim),c(K,Dim,Dim))
      for(k_local in 1:K) # To calculate the fourth sum in S8
      {
        for(tag in j:(j+sub_length-1))
        {
          q_x_sum3[k_local,,]=q_x_sum3[k_local,,]+q_x[tag,k_local]*Data_sequence[[tag]]%*%t(Data_sequence[[tag]])
        }
      }
      # 调和数列
      W_A = ((1-((i-1)*(T-sub_length+1)+j)^(-1))*W_A+((i-1)*(T-sub_length+1)+j)^(-1)*(((T-sub_length+1)/(sub_length-1))*q_tran_sum+U_A))
      u = (1-((i-1)*(T-sub_length+1)+j)^(-1))*u+((i-1)*(T-sub_length+1)+j)^(-1)*((T-sub_length+1)/sub_length*q_x_sum1+U_phi1)
      k = (1-((i-1)*(T-sub_length+1)+j)^(-1))*k+((i-1)*(T-sub_length+1)+j)^(-1)*((T-sub_length+1)/sub_length*q_x_sum2+U_phi2)
      for(i_local in 1:K)
      {
        sigma[[i_local]]=(1-((i-1)*(T-sub_length+1)+j)^(-1))*sigma[[i_local]]+((i-1)*(T-sub_length+1)+j)^(-1)*((T-sub_length+1)/sub_length*q_x_sum3[i_local,,]+U_phi3[[i_local]])
      }
      v = (1-((i-1)*(T-sub_length+1)+j)^(-1))*v+((i-1)*(T-sub_length+1)+j)^(-1)*((T-sub_length+1)/sub_length*q_x_sum4+U_phi4)
    }
  }
    result <- list(W_A)
   result[[1]]=  W_A
   result[[2]]=  u
   result[[3]] = k
   result[[4]] =sigma
   result[[5]] = v
  return(result)
}
data <-list(c(1.8,2.0))
sigma0 = diag(rep(100,2))
for(i in 1:100)
{
  if(i%%6==0||i%%6==1)
  {
    data[[i]]=rnorm(2);
  }
  else
  {
    data[[i]]= rnorm(2,c(200,200),10);
  }
}
K=2
Dim=length(data[[1]])# Gassiuan Component Number
u0 = 1:Dim # initial mean value
k0=1             # initial parameter
v0=6             # initial parameter
 # initial covariance with identity matrix
U_A = matrix(rep(1/K,K*K),K)        # Hyperperameter the same notation with the paper
U_phi1 = matrix(rep(0,Dim*K),Dim) # Hyperperameter the same notation with the paper
for(i in 1:K)
{
  U_phi1[,i]=2*k0*u0;
}
U_phi2 = rep(2*k0,K)    # Hyperperameter the same notation with the paper
U_phi3 = list(K)
for(i in 1:K)
{
  U_phi3[[i]]=  sigma0+2*k0*u0%*%t(u0);           # initial every component with the same covariance
}
U_phi4 = rep(v0+2+2,K)  # Hyperperameter the same notation with the paper
sub_length = 5   # sample chain length  # data_dimention
u = matrix(rep(0,Dim*K),Dim) ;  # claim mean value of K components.
k = 1:K;                  # claim k, which controls Gussuian Distribution
v = (3+Dim):(2+Dim+K);                  # claim v, which controls Gussuian Distribution
for( i in 1:K)
{
  u[,i]=u0;                    # initial every component with the same mean value
}
sigma = list(K)
for(i in 1:K)
{
  sigma[[i]]= sigma0;         # initial every component with the same covariance
}
W_A = matrix(1:(K*K),K);
result = SVI.HMM(data,sub_length =100,K=2,W_A, u, k, sigma, v, U_A, U_phi1,U_phi2, U_phi3, U_phi4, Pass=1000, pre=0.01)
