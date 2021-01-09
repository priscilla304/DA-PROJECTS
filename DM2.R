library(rtweet)
library(ggplot2)
library(dplyr)
library(tidytext)
library(igraph)
library(ggraph)
library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")
library(wordcloud)
library(plotrix)
library(dendextend)
library(ggplot2)
library(ggthemes)
library(RWeka)
library(qdapRegex)
library("twitteR")
h1b_tweets <- search_tweets(q = "#H1B", n = 500,
                            lang = "en",
                            include_rts = FALSE)


#tweets <- searchTwitter(searchString = "data science", n=100, lang="en", since="2017-03-01", until="2017-04-07" , resultType = "popular")
lck_tweets <- search_tweets(q = "#pandemic#usa", n = 300,
                            lang = "en",
                            include_rts = FALSE)


#df_corpus <- VCorpus(h1b_tweets$text)


#h1b <- VectorSource(h1b_tweets$text)
#h1b_corpus <- VCorpus(h1b)
#h1b_corpus

h1b_tweets$text <- qdapRegex::rm_twitter_url(
  h1b_tweets$text,
  replacement = " ",
  clean = TRUE
)


h1b_tweets$text <- qdapRegex::rm_url(
  h1b_tweets$text,
  replacement = " ",
  clean = TRUE
)

h1b_tweets$text <- qdapRegex::rm_hash(
  h1b_tweets$text,
  replacement = " ",
  clean = TRUE
)

h1b_tweets$text <- qdapRegex::rm_tag(
  h1b_tweets$text,
  replacement = " ",
  clean = TRUE
)

h1b_tweets$text <- qdapRegex::rm_emoticon(
  h1b_tweets$text,
  replacement = " ",
  clean = TRUE
)

h1b_tweets$text <- qdapRegex::rm_email(
  h1b_tweets$text,
  replacement = " ",
  clean = TRUE
)

h1b_tweets$text <- qdapRegex::rm_between(
  h1b_tweets$text,
  left = "<",
  right = ">",
  replacement = " ",
  clean = TRUE
)


h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "â???o",
  replacement = " "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "â???T",
  replacement = "'"
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "ã",
  replacement = "a"
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "s",
  replacement = "s"
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "â",
  replacement = "a"
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "¿",
  replacement = "?"
)


#h1b_tweets$text


h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\\+",
  replacement = " plus "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "=",
  replacement = " equals "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "±",
  replacement = " plus minus "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\u009d",
  replacement = " "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "¼",
  replacement = " quarter "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "½",
  replacement = " half "
)

h1b_tweets$text


h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\U0001f615",
  replacement = " "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern =  "³",
  replacement = " 3 "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "º",
  replacement = " degree "
)


h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\U0001f1e8",
  replacement = " "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\U0001f1e6",
  replacement = " "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\U0001f1fa",
  replacement = " "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\U0001f1f8",
  replacement = " "
)

#h1b_tweets$text


h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = " i'm ",
  replacement = " i am "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "'re ",
  replacement = " are "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "'t ",
  replacement = " not "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "'ve ",
  replacement = " have "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "'ll ",
  replacement = " will "
)


h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = " doesn't ",
  replacement = " does not "
)

h1b_tweets$text <- qdapRegex::rm_phone(
  h1b_tweets$text,
  clean = TRUE
)

h1b_tweets$text <- qdapRegex::rm_zip(
  h1b_tweets$text,
  clean = TRUE
)

h1b_tweets$text <- qdapRegex::rm_time(
  h1b_tweets$text,
  replacement = " ",
  clean = TRUE
)


h1b_tweets$text <- qdapRegex::rm_date(
  h1b_tweets$text,
  replacement = " ",
  clean = TRUE
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "â???T",
  replacement = "'"
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "â",
  replacement = " "
)
h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "???",
  replacement = " "
)
h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "T",
  replacement = " "
)

h1b_tweets$text <- stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\\W",
  replacement = " "
)


#stop words

h1b_tweets$text <- tm::removeWords(
  x = h1b_tweets$text,
  words = tm::stopwords(kind = "SMART")
)

h1b_tweets$text <- tm::removeWords(
  x = h1b_tweets$text,
  words = tm::stopwords(kind = "english")
)

#install.packages("qdapDictionaries")
library(qdapDictionaries)
h1b_tweets$text <- tm::removeWords(
  x = h1b_tweets$text,
  words = qdapDictionaries::Top200Words
)
#h1b_tweets$text

h1b_tweets$text <- trimws(stringr::str_replace_all(
  string = h1b_tweets$text,
  pattern = "\\s+",
  replacement = " "
))

h1b_tweets$text <- gsub(
  x = h1b_tweets$text,
  pattern = "https|http|ita|amp",
  replacement = " "
)
h1b_tweets$text <- gsub(
  x = h1b_tweets$text,
  pattern = "american",
  replacement = "america"
)
h1b_tweets$text <- gsub(
  x = h1b_tweets$text,
  pattern = "thata",
  replacement = "that"
)
h1b_tweets$text <- gsub(x = h1b_tweets$text,pattern = "realdonaldtrump",replacement = "donald trump")
h1b_tweets$text <- gsub(
  x = h1b_tweets$text,
  pattern = "trumpa",
  replacement = "trump"
)
h1b_tweets$text <- gsub(
  x = h1b_tweets$text,
  pattern = "presidenti",
  replacement = "presid"
)

h1b_tweets$text <- gsub(
  x = h1b_tweets$text,
  pattern = "kill",
  replacement = ""
)


strsplit_tweets <- strsplit(h1b_tweets$text," ")
dictionary_tweets <- sort(unique(unlist(strsplit_tweets)))
strsplit_tweets <- lapply(
  X = strsplit_tweets,
  FUN = tm::stemDocument
)

strsplit_tweets <- lapply(
  X = strsplit_tweets,
  FUN = paste,
  collapse = " "
)


h1b_tweets$text <- unlist(strsplit_tweets)





docs <- Corpus(VectorSource(h1b_tweets$text))
#docs <- tm_map(docs, removeWords, c("kill", "american")) 
doc_mat <- TermDocumentMatrix(docs)

m <- as.matrix(doc_mat)

v <- sort(rowSums(m), decreasing = TRUE)

d_data <- data.frame(word = names(v), freq = v)

head(d_data, 5)

wordcloud(words = d_data$word, 
          freq = d_data$freq, 
          min.freq = 1,
          max.words = 100,
          random.order = FALSE,
          rot.per = 0.0, 
          colors = brewer.pal(4, "Set1"))


require(devtools)
library("wordcloud2")
wordcloud2(data = d_data)

barplot(d_data[1:25,]$freq, 
        las = 2, 
        names.arg = d_data[1:25,]$word,
        col ="orange", 
        main ="Most Frequent Words from H1B keyword",
        ylab = "Word Count")

h1b <- VectorSource(h1b_tweets$text)
h1b_corpus <- VCorpus(h1b)
hco <- DocumentTermMatrix(h1b_corpus)

m  <- as.matrix(hco)

distMatrix <- dist(m, method="euclidean")


groups <- hclust(distMatrix,method="ward.D")
plot(clustering.hierarchical, cex=0.9, hang=-1)
rect.hclust(groups, k=5)


#TOPIC AND SENTIMENT FOLLOWS

#TOPIC MODELLING
Corpus_tweets <- tm::VCorpus(tm::VectorSource(h1b_tweets$text))
dtm_tweets <- tm::DocumentTermMatrix(Corpus_tweets)
dtm_tweets <- tm::removeSparseTerms(
  x = dtm_tweets,
  sparse = 0.95
)
tm::inspect(dtm_tweets)


rowTotals <- apply(dtm_tweets , 1, sum)
tweets <- h1b_tweets[rowTotals> 0, ]
dtm_tweets <- dtm_tweets[rowTotals> 0, ]
library(MASS)
#install.packages("topicmodels")
library(topicmodels)
lda_5 <- topicmodels::LDA(dtm_tweets, k = 5)
top_10terms <- topicmodels::terms(lda_5,10)
top_10terms


cat<-modeltools::posterior(lda_5)$terms
cat
ggplot(cat, aes(x = american, y = worker)) + geom_point()
library(ggplot2)
#install.packages("GmAMisc")
library(GmAMisc)
help(ppdPlot)
n<-data.frame(cat)
scatterplotMatrix(n[1:10])
help(ppdPlot)


M_topics <- modeltools::posterior(lda_5)$topics
head(M_topics)
tail(M_topics)
cor(M_topics)
topicmodels::terms(lda_5, threshold=0.0075)
topicmodels::logLik(lda_5)
topicmodels::perplexity(lda_5)
list_LDA <- lapply(
  X = 2:20,
  FUN = function(x) topicmodels::LDA(dtm_tweets, k = x)
)
v_logLik <- sapply(
  X = list_LDA,
  FUN = topicmodels::logLik
)
plot(2:20,v_logLik,type = "l",main = "Log likelihood of LDA models")


v_perplexity <- sapply(
  X = list_LDA,
  FUN = topicmodels::perplexity
)
plot(2:20,v_perplexity,type = "l",main = "Log likelihood of LDA models")


#sentiment
#install.packages("syuzhet")
library("syuzhet")
sentiment <- syuzhet::get_nrc_sentiment(h1b_tweets$text)
head(sentiment)

#EMOT

# Emotions for each tweet using NRC dictionary
emotions <- get_nrc_sentiment(h1b_tweets$text)
emo_bar = colSums(emotions)
emo_sum = data.frame(count=emo_bar, emotion=names(emo_bar))
emo_sum$emotion = factor(emo_sum$emotion, levels=emo_sum$emotion[order(emo_sum$count, decreasing = TRUE)])



# Visualize the emotions from NRC sentiments
library(plotly)
p <- plot_ly(emo_sum, x=~emotion, y=~count, type="bar", color=~emotion) %>%
  layout(xaxis=list(title=""), showlegend=FALSE,
         title="Emotion Type for hashtag: #H1B")
p
#api_create(p,filename="Sentimentanalysis")

wordcloud_tweet = c(
  paste(h1b_tweets$text[emotions$anger > 0], collapse=" "),
  paste(h1b_tweets$text[emotions$anticipation > 0], collapse=" "),
  paste(h1b_tweets$text[emotions$disgust > 0], collapse=" "),
  paste(h1b_tweets$text[emotions$fear > 0], collapse=" "),
  paste(h1b_tweets$text[emotions$joy > 0], collapse=" "),
  paste(h1b_tweets$text[emotions$sadness > 0], collapse=" "),
  paste(h1b_tweets$text[emotions$surprise > 0], collapse=" "),
  paste(h1b_tweets$text[emotions$trust > 0], collapse=" ")
)

# create corpus
corpus = Corpus(VectorSource(wordcloud_tweet))

# remove punctuation, convert every word in lower case and remove stop words

corpus = tm_map(corpus, tolower)
corpus = tm_map(corpus, removePunctuation)
corpus = tm_map(corpus, removeWords, c(stopwords("english")))
corpus = tm_map(corpus, stemDocument)

# create document term matrix

tdm = TermDocumentMatrix(corpus)

# convert as matrix
tdm = as.matrix(tdm)
tdmnew <- tdm[nchar(rownames(tdm)) < 11,]

# column name binding
colnames(tdm) = c('anger', 'anticipation', 'disgust', 'fear', 'joy', 'sadness', 'surprise', 'trust')
colnames(tdmnew) <- colnames(tdm)
par(mfrow=c(1,1))
comparison.cloud(tdmnew, random.order=FALSE,
                 colors = c("#00B2FF", "red", "#FF0099", "#6600CC", "green", "orange", "blue", "brown"),
                 title.size=1, max.words=100, scale=c(2.5, 0.4),rot.per=0.4)


#transpose
td<-data.frame(t(emotions))
#The function rowSums computes column sums across rows for each level of a grouping variable.
td_new <- data.frame(rowSums(td[2:253]))
#Transformation and cleaning
names(td_new)[1] <- "count"
td_new <- cbind("sentiment" = rownames(td_new), td_new)
rownames(td_new) <- NULL
td_new2<-td_new[1:8,]
#Plot One - count of words associated with each sentiment
quickplot(sentiment, data=td_new2, weight=count, geom="bar", fill=sentiment, ylab="count")+ggtitle("Survey sentiments")

tdm <- tm::DocumentTermMatrix(h1b_corpus) 
tdm.tfidf <- tm::weightTfIdf(tdm)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra)
tdm.tfidf <- tm::removeSparseTerms(tdm.tfidf, 0.999) 
tfidf.matrix <- as.matrix(tdm.tfidf) 
# Cosine distance matrix (useful for specific clustering algorithms) 
dist.matrix = proxy::dist(tfidf.matrix, method = "cosine")
truth.K = 5
clustering.kmeans <- kmeans(tfidf.matrix, truth.K) 
clustering.hierarchical <- hclust(dist.matrix, method = "ward.D2") 
clustering.dbscan <- dbscan::hdbscan(dist.matrix, minPts = 10)

master.cluster <- clustering.kmeans$cluster 
slave.hierarchical <- cutree(clustering.hierarchical, k = truth.K) 
slave.dbscan <- clustering.dbscan$cluster 
stacked.clustering <- rep(NA, length(master.cluster))  
names(stacked.clustering) <- 1:length(master.cluster)
for (cluster in unique(master.cluster)) { 
  indexes = which(master.cluster == cluster, arr.ind = TRUE) 
  slave1.votes <- table(slave.hierarchical[indexes]) 
  slave1.maxcount <- names(slave1.votes)[which.max(slave1.votes)]   
  slave1.indexes = which(slave.hierarchical == slave1.maxcount, arr.ind = TRUE) 
  slave2.votes <- table(slave.dbscan[indexes]) 
  slave2.maxcount <- names(slave2.votes)[which.max(slave2.votes)]   
  stacked.clustering[indexes] <- slave2.maxcount 
}

points <- cmdscale(dist.matrix, k = 2) 
palette <- colorspace::diverge_hcl(truth.K) # Creating a color palette 
previous.par <- par(mfrow=c(2,2), mar = rep(1.5, 4)) 

plot(points, main = 'K-Means clustering', col = as.factor(master.cluster), 
     mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') 
plot(points, main = 'Hierarchical clustering', col = as.factor(slave.hierarchical), 
     mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),  
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') 
plot(points, main = 'Density-based clustering', col = as.factor(slave.dbscan), 
     mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') 
plot(points, main = 'Stacked clustering', col = as.factor(stacked.clustering), 
     mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') 
par(previous.par) # recovering the original plot space paramete


######################

#Taking tweets lockdown specific

# Combine both corpora: all_tweets
all_h1b <- paste(h1b_tweets$text, collapse = "")
all_lck <- paste(lck_tweets$text, collapse = "")
all_tweets <- c(all_h1b, all_lck)

# clean all_tweets
all_tweets <- VectorSource(all_tweets)
all_corpus <- VCorpus(all_tweets)

all_tweets <- qdapRegex::rm_twitter_url(
  all_tweets,
  replacement = " ",
  clean = TRUE
)

all_tweets <- qdapRegex::rm_url(
  all_tweets,
  replacement = " ",
  clean = TRUE
)

all_tweets <- qdapRegex::rm_hash(
  all_tweets,
  replacement = " ",
  clean = TRUE
)

all_tweets <- qdapRegex::rm_tag(
  all_tweets,
  replacement = " ",
  clean = TRUE
)

all_tweets <- qdapRegex::rm_emoticon(
  all_tweets,
  replacement = " ",
  clean = TRUE
)

all_tweets <- qdapRegex::rm_email(
  all_tweets,
  replacement = " ",
  clean = TRUE
)

all_tweets <- qdapRegex::rm_between(
  all_tweets,
  left = "<",
  right = ">",
  replacement = " ",
  clean = TRUE
)


all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "â???o",
  replacement = " "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "â???T",
  replacement = "'"
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "ã",
  replacement = "a"
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "s",
  replacement = "s"
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "â",
  replacement = "a"
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "¿",
  replacement = "?"
)


#all_tweets


all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "\\+",
  replacement = " plus "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "=",
  replacement = " equals "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "±",
  replacement = " plus minus "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "\u009d",
  replacement = " "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "¼",
  replacement = " quarter "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "½",
  replacement = " half "
)

#all_tweets


all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "\U0001f615",
  replacement = " "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern =  "³",
  replacement = " 3 "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "º",
  replacement = " degree "
)



all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "\U0001f1e8",
  replacement = " "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "\U0001f1e6",
  replacement = " "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "\U0001f1fa",
  replacement = " "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "\U0001f1f8",
  replacement = " "
)

#all_tweets


all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = " i'm ",
  replacement = " i am "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "'re ",
  replacement = " are "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "'t ",
  replacement = " not "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "'ve ",
  replacement = " have "
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "'ll ",
  replacement = " will "
)


all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = " doesn't ",
  replacement = " does not "
)

all_tweets <- qdapRegex::rm_phone(
  all_tweets,
  clean = TRUE
)

all_tweets <- qdapRegex::rm_zip(
  all_tweets,
  clean = TRUE
)

all_tweets <- qdapRegex::rm_time(
  all_tweets,
  replacement = " ",
  clean = TRUE
)


all_tweets <- qdapRegex::rm_date(
  all_tweets,
  replacement = " ",
  clean = TRUE
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "â???T",
  replacement = "'"
)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "â",
  replacement = " "
)
all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "???",
  replacement = " "
)
all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "T",
  replacement = " "
)



all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "[:punct:]",
  replacement = " "
)

#all_tweets <- stringr::str_replace_all(
  #string = all_tweets,
  #pattern = "[:digit:]",
  #replacement = " "
#)

all_tweets <- stringr::str_replace_all(
  string = all_tweets,
  pattern = "\\W",
  replacement = " "
)

#all_tweets

#stop words

all_tweets <- tm::removeWords(
  x = all_tweets,
  words = tm::stopwords(kind = "SMART")
)

all_tweets <- tm::removeWords(
  x = all_tweets,
  words = tm::stopwords(kind = "english")
)

#install.packages("qdapDictionaries")
library(qdapDictionaries)
all_tweets <- tm::removeWords(
  x = all_tweets,
  words = qdapDictionaries::Top200Words
)
#all_tweets

all_tweets <- trimws(stringr::str_replace_all(
  string = all_tweets,
  pattern = "\\s+",
  replacement = " "
))


strsplit_tweets2 <- strsplit(all_tweets," ")
dictionary_tweets2 <- sort(unique(unlist(strsplit_tweets)))
strsplit_tweets2 <- lapply(
  X = strsplit_tweets,
  FUN = tm::stemDocument
)

strsplit_tweets2 <- lapply(
  X = strsplit_tweets,
  FUN = paste,
  collapse = " "
)

all_tweets <- unlist(strsplit_tweets)
#all_tweets

#h1b_tweets$text<-all_tweets
#h1b_tweets$text
write.csv(
  x = all_tweets$text,
  file = "cleaned_tweets.csv",
  row.names = FALSE
)

file_path_tweets2 <- c(
  "cleaned_tweets.csv"
)

list_tweets2 <- lapply(
  X = file_path_tweets,
  FUN = read.csv,
  colClasses = "character"
)

names(list_tweets) <- file_path_tweets


#####################################

# Add new stop words to clean_corpus()
clean_corpus2 <- function(corpus){
  corpus <- tm_map(corpus, removePunctuation)
  corpus <- tm_map(corpus, stripWhitespace)
  corpus <- tm_map(corpus, removeNumbers)
  corpus <- tm_map(corpus, content_transformer(tolower))
  corpus <- tm_map(corpus, removeWords, 
                   c(stopwords("en"), "amp", "foreign", "train", "freez", "don"))
  return(corpus)
}


all_clean2 <- clean_corpus2(all_corpus)
all_tdm2 <- TermDocumentMatrix(all_clean2)
all_m <- as.matrix(all_tdm2)
commonality.cloud(all_m, 
                  colors = "steelblue1",
                  max.words = 100)


###



allc <- VectorSource(all_tweets)
allc_corpus <- VCorpus(allc)
#h1b_corpus <- tm_map(h1b_corpus, removeWords, c(stopwords("en"), "amp", "web", "false","true","app","web","twitter","photo"))
allc_corpus

alhc <- DocumentTermMatrix(allc_corpus)
#htm<-tm::DocumentTermMatrix(hco)

# Convert coffee_dtm to a matrix:
all_mat <- as.matrix(alhc)

# Print the dimensions
dim(all_mat)

# Create a TDM from clean_corp: tdm
tdm2 <- TermDocumentMatrix(allc_corpus)

# Print tdm data
#tdm

tdm_mat2 <- as.matrix(tdm)
dim(tdm_mat2)


# Create a matrix: 
cm2 <- as.matrix(tdm2)

# Calculate the rowSums: term_frequency
term_frequency2 <- rowSums(cm2)


#Dissimilar

# Clean the corpus
all_clean2 <- clean_corpus(all_corpus)

# Create all_tdm
all_tdm2 <- TermDocumentMatrix(all_clean)

# Give the columns distinct names
colnames(all_tdm2) <- c("h1b", "lockdown")

# Create all_m
all_m <- as.matrix(all_tdm2)

# Create comparison cloud
comparison.cloud(all_m,
                 colors = c("orange", "blue"),
                 scale=c(3, .2),
                 max.words = 50)


# Identify terms shared by both documents
common_words <- subset(
  all_m,
  all_m[, 1] > 0 & all_m[, 2] > 0
)

tail(common_words)


# calc common words and difference
difference <- abs(common_words[, 1] - common_words[, 2])
common_words <- cbind(common_words, difference)
common_words <- common_words[order(common_words[, 3],
                                   decreasing = T), ]
head(common_words)




# Create a TDM from clean_corp: tdm
tdm2 <- TermDocumentMatrix(all_corpus)

# Print tdm data
#tdm

tdm_mat2 <- as.matrix(tdm2)
dim(tdm_mat2)


# Create a matrix: 
#cm2 <- as.matrix(tdm2)

# Calculate the rowSums: term_frequency
term_frequency <- rowSums(tdm_mat2)

# Sort term_frequency in descending order
term_frequency <- sort(term_frequency, decreasing = T)

# View the top 10 most common words
term_frequency[1:10]

#common words bar plot

#barplot(term_frequency[1:10], col = "orange", las = 2)


DocumentTermMatrix_tweets <- tm::removeSparseTerms(
  tdm2,
  0.995
)

term_frequency2 <- data.frame(
  Term = colnames(tdm_mat2),
  Frequency = colSums(tdm_mat2),
  stringsAsFactors = FALSE
)

#Frequent terms
term_frequency2 <- data.frame(
  Term = colnames(tdm_mat2),
  Frequency = colSums(tdm_mat2),
  stringsAsFactors = FALSE
)
term_frequency2 <- term_frequency[order(term_frequency$Frequency),]
tail(term_frequency)


colnames(tdm_mat2) <- c("h1b", "lockdown")
###Frequency Word cloud
comparison.cloud(tdm_mat2,
                 colors = c("orange", "blue"),
                 max.words = 50)




wordcloud_tweet = c(
  paste(all_tweets[emotions$anger > 0], collapse=" "),
  paste(all_tweets[emotions$anticipation > 0], collapse=" "),
  paste(all_tweets[emotions$disgust > 0], collapse=" "),
  paste(all_tweets[emotions$fear > 0], collapse=" "),
  paste(all_tweets[emotions$joy > 0], collapse=" "),
  paste(all_tweets[emotions$sadness > 0], collapse=" "),
  paste(all_tweets[emotions$surprise > 0], collapse=" "),
  paste(all_tweets[emotions$trust > 0], collapse=" ")
)

# create corpus
corpus = Corpus(VectorSource(wordcloud_tweet))

# remove punctuation, convert every word in lower case and remove stop words

corpus = tm_map(corpus, tolower)
corpus = tm_map(corpus, removePunctuation)
corpus = tm_map(corpus, removeWords, c(stopwords("english")))
corpus = tm_map(corpus, stemDocument)

# create document term matrix

tdm2 = TermDocumentMatrix(corpus)

# convert as matrix
tdm2 = as.matrix(tdm2)
tdmnew2 <- tdm2[nchar(rownames(tdm2)) < 11,]

# column name binding
colnames(tdm2) = c('anger', 'anticipation', 'disgust', 'fear', 'joy', 'sadness', 'surprise', 'trust')
colnames(tdmnew2) <- colnames(tdm2)
comparison.cloud(tdmnew2, random.order=FALSE,
                 colors = c("#00B2FF", "red", "#FF0099", "#6600CC", "green", "orange", "blue", "brown"),
                 title.size=1, max.words=250, scale=c(2.5, 0.4),rot.per=0.4)



#######CLUSTERING#############
tdm <- tm::DocumentTermMatrix(alhc) 
tdm.tfidf <- tm::weightTfIdf(alhc)

library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra)
tdm.tfidf <- tm::removeSparseTerms(tdm.tfidf, 0.999) 
tfidf.matrix <- as.matrix(tdm.tfidf) 
# Cosine distance matrix (useful for specific clustering algorithms) 
dist.matrix = proxy::dist(tfidf.matrix, method = "cosine")
truth.K = 5
clustering.kmeans <- kmeans(tfidf.matrix, truth.K) 
clustering.hierarchical <- hclust(dist.matrix, method = "ward.D2") 
clustering.dbscan <- dbscan::hdbscan(dist.matrix, minPts = 10)

master.cluster <- clustering.kmeans$cluster 
slave.hierarchical <- cutree(clustering.hierarchical, k = truth.K) 
slave.dbscan <- clustering.dbscan$cluster 
stacked.clustering <- rep(NA, length(master.cluster))  
names(stacked.clustering) <- 1:length(master.cluster)
for (cluster in unique(master.cluster)) { 
  indexes = which(master.cluster == cluster, arr.ind = TRUE) 
  slave1.votes <- table(slave.hierarchical[indexes]) 
  slave1.maxcount <- names(slave1.votes)[which.max(slave1.votes)]   
  slave1.indexes = which(slave.hierarchical == slave1.maxcount, arr.ind = TRUE) 
  slave2.votes <- table(slave.dbscan[indexes]) 
  slave2.maxcount <- names(slave2.votes)[which.max(slave2.votes)]   
  stacked.clustering[indexes] <- slave2.maxcount 
}

points <- cmdscale(dist.matrix, k = 5) 
palette <- colorspace::diverge_hcl(truth.K) # Creating a color palette 
previous.par <- par(mfrow=c(2,2), mar = rep(1.5, 4)) 

plot(points, main = 'K-Means clustering', col = as.factor(master.cluster), 
     mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') 
plot(points, main = 'Hierarchical clustering', col = as.factor(slave.hierarchical), 
     mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),  
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') 
plot(points, main = 'Density-based clustering', col = as.factor(slave.dbscan), 
     mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') 
plot(points, main = 'Stacked clustering', col = as.factor(stacked.clustering), 
     mai = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), 
     xaxt = 'n', yaxt = 'n', xlab = '', ylab = '') 
par(previous.par) # recovering the original plot space parameters

#########CLUSTERING ENDS##########