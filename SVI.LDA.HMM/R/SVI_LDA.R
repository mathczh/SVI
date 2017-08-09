LDA.SVI<-function(X, # the Input data
                  K, # the number of the topics
                  n, # the number of the passes
                  alpha, # the parameter of theta
                  eta, # the parameter of beta
                  pre, # the stop criteria
                  topic_length# the number of top words and topics

)
{
       library(Rcpp)
       library(RcppArmadillo)
       library(tokenizers)
       sourceCpp(file='src/LDA.cpp')
       return(SVI_LDA(X,K,n,alpha,eta,pre,topic_length))
}

Tokenize<-function(Data_set, # the Input original data, which is a list data type
                   Language = "en"# the language type you want to produce, defaut language is en
                   )
{
     library(tokenizers);
     tokenized_data <-tokenize_words(Data_set, stopwords = stopwords(Language));
     return(tokenized_data);
}

X=list(c('home','mother','mother','house','mother','mother','mother','mother'),c('father','dad','cat','dad','cat','cat','cat','cat','cat','cat','cat','cat','cat','cat'))
K=3;
alpha=1/K;
eta=0.01;
n=2;
pre=0.1;
t=LDA.SVI(X,K,n,alpha=1/K,eta=0.01,pre=0.1,2)
doc_a = "Brocolli is good to eat. My brother likes to eat good brocolli, but not my mother."
doc_b = "My mother spends a lot of time driving my brother around to baseball practice."
doc_c = "Some health experts suggest that driving may cause increased tension and blood pressure."
doc_d = "I often feel pressure to perform well at school, but my mother never seems to drive my brother to do better."
doc_e = "Health professionals say that brocolli is good for your health."
data <- list(doc_a,doc_b,doc_c,doc_d,doc_e)
