library(tseries)
library(factoextra)
library(devtools)
library(umap)
library(dplyr)
library(stringr)
library(randomForest)
library(caret)
##################ML and visualization script for the practical course
#Development of a risk score to analyze SARS-CoV-2 genomes
#winter term 24/25
#Friederike Wohlfarth and Ibraim Ibraimi
setwd("C:/Users/Nutzer/Documents/covid2")
mat <- read.matrix('outtest.txt',header=FALSE)
mat <- matrix(mat, ncol = 443, nrow=46,byrow = TRUE)
col <- read.delim('lineages.txt', header=FALSE)

#data.pca <- princomp(t(mat))
res.pca <- prcomp((mat))
summary(res.pca)


fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
fviz_pca_ind(res.pca, habillage=col$V2, 
             palette = "Dark2",
             label="none",
          #   addEllipses=TRUE, ellipse.level=0.95
             )
fviz_pca_var(res.pca, col.var = "black")
pc <- prcomp(mat)
summary(res.pca)


#umap

iris.umap <- umap::umap(mat)
plot.iris(iris.umap, as.factor(col$V2))
#now kmeans
#find no of clusters
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(mat, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot


k2 <- kmeans(mat, centers = 4, nstart = 25)
fviz_cluster(k2, data = mat)
k2$cluster<-col$V2

fviz_cluster(k2, data = mat)

#now add risc scores
#risc scores
risc<- read.delim('escape.csv',sep=",")


risc2<-risc %>%
  group_by(site) %>%
  summarise(score = sum(escape))

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
risc2$score_scaled<-scale_values(risc2$score)
risc2$obergrenze<-risc2$site*3+21563
risc2$untergrenze<-(risc2$site*3+21563)-2
mut<-read.delim('list_of_unique_mutations.txt',sep=",", header=FALSE)
#process new data matrix

for (i in 1:46){
  for (j in 1:443){
    if (mat[i,j]==1){
      for (k in 1:187){
        if(risc2$untergrenze[k] <= as.numeric(unlist(str_extract_all(mut[j,], "\\d+"))[1]) && as.numeric(unlist(str_extract_all(mut[j,], "\\d+"))[1])<=risc2$obergrenze[k]){
          mat[i,j]=as.numeric(risc2$score[k])
        }
        
      }
      
    }
  } 
}

res.pca <- prcomp(mat)
summary(res.pca)


fviz_pca_biplot(res.pca, repel = FALSE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
fviz_pca_ind(res.pca, habillage=col$V2, 
             palette = "Dark2",
             label="none",
             #   addEllipses=TRUE, ellipse.level=0.95
)
fviz_pca_var(res.pca, col.var = "black")

pc <- prcomp(mat)
summary(res.pca)


#umap

iris.umap <- umap::umap(mat)
plot.iris(iris.umap, as.factor(col$V2))
#now kmeans
#find no of clusters
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(mat, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot


k2 <- kmeans(mat, centers = 4, nstart = 25)
fviz_cluster(k2, data = mat)
k2$cluster<-col$V2

fviz_cluster(k2, data = mat)


#now make row sums
target<-rowSums(mat)
quantile(target)
qnt <- quantile(target,seq(0,1,.25))

quant<-cut(target,unique(qnt),include.lowest=TRUE)
k2$cluster<-quant
fviz_cluster(k2, data = mat)



#random forest
set.seed(222)
colnames(mat)<-as.character(mut$V1)
colnames(mat)<-as.character(1:443)
colnames(mat)<-NULL

ind <- sample(2, nrow(mat), replace = TRUE, prob = c(0.7, 0.3))
train <- mat[ind==1,]
test <- mat[ind==2,]


rf <- randomForest(quant[ind==1]~., data=train, proximity=TRUE) 
print(rf)

p1 <- predict(rf, train)
confusionMatrix(p1, quant[ind==1])

#on test data
p2 <- predict(rf, test)
confusionMatrix(p2, quant[ind==2
                          ])

plot(rf)
t <- tuneRF(train, quant[ind==1],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")
#Variable Importance

varImpPlot(rf,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf)

#######now the same thing for spike matrix

mat <- read.matrix('outtest2.txt',header=FALSE)
mat <- matrix(mat, ncol = 443, nrow=46,byrow = TRUE)
col <- read.delim('lineages.txt', header=FALSE)

#data.pca <- princomp(t(mat))
res.pca <- prcomp((mat))
summary(res.pca)


fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
fviz_pca_ind(res.pca, habillage=col$V2, 
             palette = "Dark2",
             label="none",
             #   addEllipses=TRUE, ellipse.level=0.95
)
fviz_pca_var(res.pca, col.var = "black")
pc <- prcomp(mat)
summary(res.pca)


#umap

iris.umap <- umap::umap(mat)
plot.iris(iris.umap, as.factor(col$V2))
#now kmeans
#find no of clusters
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(mat, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot


k2 <- kmeans(mat[,-(which(colSums(mat)==0))], centers = 5, nstart = 25)
fviz_cluster(k2, data = mat[,-(which(colSums(mat)==0))])
k2$cluster<-col$V2

fviz_cluster(k2, data = mat[,-(which(colSums(mat)==0))])

#now add risc scores
#risc scores
risc<- read.delim('escape.csv',sep=",")


risc2<-risc %>%
  group_by(site) %>%
  summarise(score = sum(escape))

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
risc2$score_scaled<-scale_values(risc2$score)
risc2$obergrenze<-risc2$site*3+21563
risc2$untergrenze<-(risc2$site*3+21563)-2
write.csv(risc2,"score_new.csv", row.names = FALSE)

mut<-read.delim('list_of_unique_mutations.txt',sep=",", header=FALSE)
#process new data matrix

for (i in 1:46){
  for (j in 1:443){
    if (mat[i,j]==1){
      for (k in 1:187){
        if(risc2$untergrenze[k] <= as.numeric(unlist(str_extract_all(mut[j,], "\\d+"))[1]) && as.numeric(unlist(str_extract_all(mut[j,], "\\d+"))[1])<=risc2$obergrenze[k]){
          mat[i,j]=as.numeric(risc2$score[k])
        }
        
      }
      
    }
  } 
}

res.pca <- prcomp(mat)
summary(res.pca)


fviz_pca_biplot(res.pca, repel = FALSE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)
fviz_pca_ind(res.pca, habillage=col$V2, 
             palette = "Dark2",
             label="none",
             #   addEllipses=TRUE, ellipse.level=0.95
)
fviz_pca_var(res.pca, col.var = "black")

pc <- prcomp(mat)
summary(res.pca)


#umap

iris.umap <- umap::umap(mat)
plot.iris(iris.umap, as.factor(col$V2))
#now kmeans
#find no of clusters
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(mat, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot


k2 <- kmeans(mat[,-(which(colSums(mat)==0))], centers = 3, nstart = 25)
fviz_cluster(k2, data = mat[,-(which(colSums(mat)==0))])
k2$cluster<-col$V2

fviz_cluster(k2, data = mat[,-(which(colSums(mat)==0))])


#now make row sums
target<-rowSums(mat)
quantile(target)
qnt <- quantile(target,seq(0,1,.25))

quant<-cut(target,unique(qnt),include.lowest=TRUE)
k2$cluster<-quant
fviz_cluster(k2, data = mat[,-(which(colSums(mat)==0))])



#random forest
set.seed(222)
colnames(mat)<-as.character(mut$V1)
colnames(mat)<-as.character(1:443)
colnames(mat)<-NULL

ind <- sample(2, nrow(mat), replace = TRUE, prob = c(0.7, 0.3))
train <- mat[ind==1,]
test <- mat[ind==2,]


rf <- randomForest(quant[ind==1]~., data=train, proximity=TRUE,mtry=84) 
print(rf)

p1 <- predict(rf, train)
confusionMatrix(p1, quant[ind==1])

#on test data
p2 <- predict(rf, test)
confusionMatrix(p2, quant[ind==2
])

plot(rf)
t <- tuneRF(train, quant[ind==1],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")
#Variable Importance

varImpPlot(rf,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
importance(rf)


#statistics
library(ggplot2)
dat <- data.frame(read.delim('stats_seq.txt', sep=" ",header=FALSE))
dat<-dat[order(-dat$V2),]
ggplot(dat, aes(x = reorder(V1,-V2), y = V2)) + 
  geom_col(fill = "#0099f9",width=0.5)+
  labs(x="mutations",y="frequency of genomes")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
        axis.text.x=element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank())#+theme_bw() #+
  
dat <- data.frame(read.delim('stats_seq_spike.txt', sep=" ",header=FALSE))
dat<-dat[order(-dat$V2),]
ggplot(dat, aes(x = reorder(V1,-V2), y = V2)) + 
  geom_col(fill = "#0099f9",width=0.5)+
  labs(x="mutations",y="frequency of genomes")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    axis.text.x=element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+



dat <- data.frame(read.delim('stats_mut.txt', sep=" ",header=FALSE))
dat<-dat[order(-dat$V2),]
ggplot(dat, aes(x = reorder(V1,-V2), y = V2)) + 
  geom_col(fill = "#0099f9",width=0.5)+
  labs(x="genomes",y="frequency of mutations")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    axis.text.x=element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+

dat <- data.frame(read.delim('stats_mut_spike.txt', sep=" ",header=FALSE))
dat<-dat[order(-dat$V2),]
ggplot(dat, aes(x = reorder(V1,-V2), y = V2)) + 
  geom_col(fill = "#0099f9",width=0.5)+
  labs(x="genomes",y="frequency of mutations")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    axis.text.x=element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+


#ns
dat <- data.frame(read.delim('numberofns.txt', sep="\t",header=FALSE))
dat<-dat[order(-dat$V2),]
ggplot(dat, aes(x = reorder(V1,-V2), y = V2)) + 
  geom_col(fill = "#0099f9",width=0.5)+
  labs(x="genomes",y="frequency of N's")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    axis.text.x=element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+

ggplot(dat, aes( y = V2)) + 
  geom_boxplot(fill = "#0099f9",width=0.5)+
  labs(x="",y="number of Ns")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    axis.text.x=element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+


ggplot(dat, aes( y = V2)) + 
  geom_boxplot(fill = "#0099f9",width=0.5)+
  labs(x="",y="number of Ns")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    axis.text.x=element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+



dat <- data.frame(read.delim('stats_del.txt', sep=" ",header=FALSE))
dat<-dat[order(-dat$V2),]
ggplot(dat, aes(x = reorder(V1,-V2), y = V2)) + 
  geom_col(fill = "#0099f9",width=0.5)+
  labs(x="genomes",y="frequency of deletions")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    axis.text.x=element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+

dat <- data.frame(read.delim('stats_seq_del.txt', sep=" ",header=FALSE))
dat<-dat[order(-dat$V2),]
ggplot(dat, aes(x = reorder(V1,-V2), y = V2)) + 
  geom_col(fill = "#0099f9",width=0.5)+
  labs(x="deletion",y="frequency of genomes")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    axis.text.x=element_blank(),panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+

ggplot(data.frame(col$V2), aes(x=col$V2)) +
  geom_bar(fill = "#0099f9",width=0.5)+  labs(x="lineage",y="frequency of genomes")+
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+


risc2<-data.frame(group=quant, score=target,lineage=col$V2)
ggplot(risc2, aes( y = score,x=group)) + 
  geom_boxplot(fill = "#0099f9",width=0.5,outlier.shape=NA)+
  labs(x="",y="risk score")+
 # stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(width = 0.2)) +
  theme(#axis.title.x=paste("mutations"),axis.title.y=paste("frequency"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x=element_blank())#+theme_bw() #+
