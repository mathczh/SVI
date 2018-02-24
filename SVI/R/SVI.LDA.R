LDA.SVI<-function(Documents , # the Input data
                  K,# the number of the topics
                  max_iter_local, 
                  max_iter_global, # the number of the passes
                  alpha, # the parameter of theta
                  eta, # the parameter of beta
                  tol_local, # the stop criteria
                  tol_global,
                  batch_size, # batch size
                  word_num  # the number of top words
)
{
  t=  SVI_LDA(Documents,K,max_iter_local,max_iter_global,alpha,eta,tol_local,tol_global,batch_size, word_num);
  names(t) <- c("Lambda","phi","gamma","Topic","Dictionary","global_iteration","E_beta","E_theta","log_likelihood","Each_doc_likelihood")
  K_topic_name <- c("Topic 1")
  for(i in  2:K)
  {
    K_topic_name <-c(K_topic_name,paste0("Topic ",toString(i)))
  }
  Top_word_name <- c("Top word 1")
  for(i in 2:word_num)
  {
    Top_word_name <-c(Top_word_name,paste0("Top word ",toString(i)))
  }
  rownames(t$Topic) <- Top_word_name
  colnames(t$Topic) <- K_topic_name
  return(t)
}
dtm2svi.ldaformat <- function(dtm){
  LDA_document <- dtm2ldaformat(dtm)
  Documents_matrix <- LDA_document$documents
  Dictionary <- LDA_document$vocab
  Documents <- list()
  for(i in 1:length(Documents_matrix)){
    document= Documents_matrix[[i]];
    document_word <- NULL;
    for(j in 1:ncol(document)){
      index <- document[1,j]
      rep_times <-document[2,j]
      word <- Dictionary[index]
      document_word <- c(document_word,rep(word,rep_times))
    }
    Documents[[i]] <- document_word
  }
  return (Documents)
}
lda2svi.ldaformat <- function(documents,
                        vocab){
  Documents_matrix <- documents
  Dictionary <- vocab
  Documents <- list()
  for(i in 1:length(Documents_matrix)){
    document= Documents_matrix[[i]];
    document_word <- NULL;
    for(j in 1:ncol(document)){
      index <- document[1,j]
      rep_times <-document[2,j]
      word <- Dictionary[index]
      document_word <- c(document_word,rep(word,rep_times))
    }
    Documents[[i]] <- document_word
  }
  return (Documents)
}
library(corpus.JSS.papers)
library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
sourceCpp(file='src/LDA.cpp')
data("JSS_papers")
#JSS_papers <- JSS_papers[JSS_papers[,"date"] < "2010-08-05",]
JSS_papers <- JSS_papers[sapply(JSS_papers[, "description"],Encoding) == "unknown",]
library("tm")
library("XML")
library("topicmodels")
remove_HTML_markup <- function(s) tryCatch({
  doc <- htmlTreeParse(paste("<!DOCTYPE html>", s),
                       asText = TRUE, trim = FALSE) 
  xmlValue(xmlRoot(doc)) 
}, error = function(s) s)
sequ_documents <- 1:length(JSS_papers[, "description"])
corpus <- Corpus(VectorSource(sapply(JSS_papers[, "description"],
                                     remove_HTML_markup)))
Sys.setlocale("LC_COLLATE", "C")
JSS_dtm <- DocumentTermMatrix(corpus,
                              control = list(stemming = TRUE, stopwords = TRUE, minWordLength = 3,
                                             removeNumbers = TRUE, removePunctuation = TRUE))
dim(JSS_dtm)
library("slam")
summary(col_sums(JSS_dtm))
term_tfidf <- tapply(JSS_dtm$v/row_sums(JSS_dtm)[JSS_dtm$i], JSS_dtm$j, mean) *log2(nDocs(JSS_dtm)/col_sums(JSS_dtm > 0))
summary(term_tfidf)
JSS_dtm <- JSS_dtm[,term_tfidf >= 0.1]
JSS_dtm <- JSS_dtm[row_sums(JSS_dtm) > 0,]
Documents <- dtm2svi.ldaformat(JSS_dtm)
topic_num = seq(2,100,5)
repete_times <-1
log_likelihood <- matrix(rep(0,repete_times*length(topic_num)), ncol = repete_times)
time_cost <- matrix(rep(0,repete_times*length(topic_num)), ncol = repete_times)
topicmodel_loglik <- matrix(rep(0,repete_times*length(topic_num)), ncol = repete_times)
topicmodel_time_cost <- matrix(rep(0,repete_times*length(topic_num)), ncol = repete_times)
for( rep_i in 1 : repete_times)
{
  rep_j = 0
  for( K in topic_num){
    rep_j = rep_j+1
    print(K)
    print("Run Start: ")
    time_start <- unclass(as.POSIXct(Sys.time())) 
    alpha=50/K;
    eta=0.1;
    max_iter_global=1000;
    max_iter_local=1000;
    tol_local=0.00001;
    tol_global= 0.00001;
    top_word_num = 5
    batch_size = 1
    result=LDA.SVI(Documents,K,max_iter_local,max_iter_global,alpha,eta,tol_local,tol_global,batch_size,top_word_num)
    print("Run End: ")
    time_end <- unclass(as.POSIXct(Sys.time())) 
    log_likelihood[rep_j,rep_i] <- result$log_likelihood
    print("Epoch Num: ")
    print(result$global_iteration)
    print("log_likelihood")
    print(result$log_likelihood)
    print("Run total time: ")
    print(time_end-time_start)
    time_cost[rep_j,rep_i] <- (time_end-time_start)
  }
  rep_j = 0
  SEED <- 2010
  for (k in topic_num)
  {
    rep_j = rep_j+1
    print("Topicmodel Run Start: ")
    time_start <- unclass(as.POSIXct(Sys.time()))
    VEM = LDA(JSS_dtm, k = k, control = list(seed = SEED))
    time_end <- unclass(as.POSIXct(Sys.time()))
    print("Topicmodel Run End: ")
    print(k)
    print(logLik(VEM))
    topicmodel_loglik[rep_j,rep_i] <- logLik(VEM)
    print("Run total time: ")
    print(time_end-time_start)
    topicmodel_time_cost[rep_j,rep_i] <- (time_end-time_start)
  }
}
log_likelihood <- rowSums(log_likelihood)/repete_times
time_cost <- rowSums(time_cost)/repete_times
topicmodel_loglik <- rowSums(topicmodel_loglik)/repete_times
topicmodel_time_cost <- rowSums(topicmodel_time_cost)/repete_times

model_data_loglik <- data.frame(topic_k = topic_num, topicmodel_loglik = topicmodel_loglik, SVI_LDA_loglik =log_likelihood)
model_data_time <- data.frame(topic_k = topic_num, topicmodel_time = topicmodel_time_cost, svi_lda_time =time_cost)
model_data_loglik= melt(model_data_loglik,id = "topic_k")
model_data_time= melt(model_data_time,id = "topic_k")
names(model_data_loglik)[names(model_data_loglik) =="value"] = "log_likehood"
names(model_data_time)[names(model_data_time)=="value"] = "runtime"
filename <- paste0("Doc_1to",as.character(length(sequ_documents)),"_tol_global=",as.character(tol_global),"_repete_times=",as.character(repete_times) ,"_batch_size=",as.character(batch_size),"tol_local=",as.character(tol_local),"max_global_iteration=",as.character(max_iter_global),"max_local_iteration=",as.character(max_iter_local),".pdf")
p <- ggplot(data = model_data_loglik ,aes(x = topic_k, y = log_likehood, colour=variable)) + geom_point(size = 2)+geom_line()
ggsave(p,file=filename, width = 12, height = 4)
filename <- paste0("time_cost_",filename)
plot_time <- ggplot(data= model_data_time, aes(x = topic_k, y = runtime)) + geom_point(size = 2)+geom_line(aes(color=variable))
ggsave(plot_time, file=filename, width=12, height=4)
