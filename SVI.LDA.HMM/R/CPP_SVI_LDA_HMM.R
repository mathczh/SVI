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
      library(Rcpp)
      library(RcppArmadillo)
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
  sourceCpp(file='src/HMM.cpp')
  result= SVI_HMM(Data_sequence,sub_length,K,W_A,u,k,sigma,v,U_A,U_phi1,U_phi2,U_phi3,U_phi4,Dim,Pass,pre)
  return(result)
}
data <-list(c(1.8,2.0),c(1,2),c(1.1,1.6),c(1.7,1.7),c(1,2),c(1.1,1.6),c(1,2),c(1.1,1.6),c(1.7,1.7),c(1,2),c(1.1,1.6),c(2,4),c(3,3))
K=2
Dim=length(data[[1]])# Gassiuan Component Number
alpha=1/K        #
u0 = 1:Dim # initial mean value
k0=1             # initial parameter
v0=6             # initial parameter
sigma0 = diag(rep(1,2)) # initial covariance with identity matrix
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
  u[,i]=i*u0;                    # initial every component with the same mean value
}
sigma = list(K)
for(i in 1:K)
{
  sigma[[i]]= i*sigma0;         # initial every component with the same covariance
}
W_A = matrix(1:(K*K),K);
result = SVI.HMM(data,sub_length =13,K=2,W_A, u, k, sigma, v, U_A, U_phi1,U_phi2, U_phi3, U_phi4, Pass=100, pre=0.01)
