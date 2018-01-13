#' Solve the hyperparameters in LDA
#'
#' @param X The Input data, should be a list
#' @param K The number of the topics
#' @param n The number of passes we run, also the maximum number of iteration
#' @param alpha The prior of theta, distribution of vocabulary over different topics
#' @param eta The prior of beta, distribution of topics in each position over different document.
#' @param pre The stop criteria
#' @param word_length The number of top words
#'
#' @author Zhehui Chen, Shiyang Li, Xingguo Li, Tuo Zhao
#'
#' @references Matthew D. Hoffman, David M. Blei, Chong Wang, John Paisley. Stochastic Variational Inference. \emph{Journal of Machine Learning Research}, 2013\cr
#' @description In this function,we require users to input their data by a list, that is, each document is one vector as a component in the list. First, we remove the stop words in each document to make those top words more meaningful. Then we build a dictionary of these document. Finally, we implement the LDA with SVI. Since we initialize the latent parameters randomly, you may not get the same result even you have the same input.
#' @return All the hyperparameters in the model. \item{Lambda}{The hyperparameter for each words of the dictionary in each topic.}\item{phi}{The hyperparameter for different topics of each position in each document.}\item{gamma}{The hyperparameter for topics proportion in each document.}\item{Topic}{The top words in each topic.}\item{Dictionary}{The dictionary for all documents except the stop words.}
#' @examples #####################################################################################
#' ## test data
#' doc_a = "Brocolli is good to eat. My brother likes to eat good brocolli, but not my mother."
#' doc_b = "My mother spends a lot of time driving my brother around to baseball practice."
#' doc_c = "Some health experts suggest that driving may cause increased tension and blood pressure."
#' doc_d = "I often feel pressure to perform well at school, but my mother never seems to drive my brother to do better."
#' doc_e = "Health professionals say that brocolli is good for your health."
#' data <- list(doc_a,doc_b,doc_c,doc_d,doc_e)
#' ## Initial prior parameters
#' K=3;
#' alpha=1/K;
#' eta=0.01;
#' n=1;
#' pre=0.01;
#' word_length=3;
#' ## run the function
#' t=LDA.SVI(data,K,n,alpha,eta,pre,word_length)
LDA.SVI<-function(X, # the Input data
                  K, # the number of the topics
                  n, # the number of the passes
                  alpha, # the parameter of theta
                  eta, # the parameter of beta
                  pre, # the stop criteria
                  word_length # the number of top words
)
{
  library(Rcpp)
  library(RcppArmadillo)
  sourceCpp(file='src/LDA.cpp')
  t=  SVI_LDA(X,K,n,alpha,eta,pre,word_length);
  names(t) <- c("Lambda","phi","gamma","Topic","Dictionary","Epoch")
  K_topic_name <- c("Topic 1")
  for(i in  2:K)
  {
    K_topic_name <-c(K_topic_name,paste0("Topic ",toString(i)))
  }
  Top_word_name <- c("Top word 1")
  for(i in 2:word_length)
  {
    Top_word_name <-c(Top_word_name,paste0("Top word ",toString(i)))
  }
  rownames(t$Topic) <- Top_word_name
  colnames(t$Topic) <- K_topic_name
  return(t)
}
Tokenize<-function(Data_set, # the Input original data, which is a list data type
                   Language = "en"# the language type you want to produce, defaut language is en
){
  
  tokenized_data <-tokenize_words(Data_set, stopwords = stopwords(Language));
  return(tokenized_data);
}
remove_HTML_markup <- function(s) tryCatch({
  doc <- htmlTreeParse(paste("<!DOCTYPE html>", s),
                       asText = TRUE, trim = FALSE) 
  xmlValue(xmlRoot(doc)) 
}, error = function(s) s)
library(corpus.JSS.papers)
library(tokenizers)
library(SnowballC)
library("tm")
library("XML")
data("JSS_papers")
LDA_document <- JSS_papers[,"description"]
LDA_document <- Corpus(VectorSource(sapply(LDA_document,
                                     remove_HTML_markup)))
Sys.setlocale("LC_COLLATE", "C")
LDA_document <- tm_map(LDA_document, content_transformer(tolower))
LDA_document <- tm_map(LDA_document,stripWhitespace)
LDA_document <- tm_map(LDA_document,  removePunctuation)
LDA_document <- tm_map(LDA_document,  removeWords, c("the","and",stopwords("english")))
LDA_document <- tm_map(LDA_document,  stemDocument)
LDA_document <- tm_map(LDA_document,  removeNumbers)
LDA_document <- as.list(LDA_document)
# you may use content_transformer(FUN) function or whatever else to custom your content transforms. More details are in "tm" package.
LDA_document<- Tokenize(LDA_document)
print("Finish Tokenize")
#print(LDA_document)
print("Run Start: ")
time_start <- Sys.time()
K=30
alpha=50/K;
eta=0.1;
n=100;
pre=0.0001;
top_word = 5
result=LDA.SVI(LDA_document,K,n,alpha,eta,pre, top_word)
print("Run End: ")
time_end <- Sys.time()
print("Run total time: ")
print(time_end-time_start)
file_out_name <- "file/log2.txt"
write.table(result$Topic,file_out_name)
cat("Epoch", file = file_out_name,append = T, sep = "\n")
cat(result$Epoch,file = file_out_name,append = T, sep = "\n")
cat("Run time",file = file_out_name,append = T, sep = "\n")
cat("Run Start time",file = file_out_name,append = T, sep = "\n")
cat(toString(time_start),file = file_out_name,append = T, sep = "\n")
cat("Run End time",file = file_out_name,append = T, sep = "\n")
cat(toString(time_end),file = file_out_name,append = T, sep = "\n")
cat("Run Total time",file = file_out_name,append = T, sep = "\n")
cat(time_end-time_start,file = file_out_name,append = T, sep = "\n")
