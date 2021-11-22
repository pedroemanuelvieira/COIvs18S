library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(magrittr)
library(dplyr)
library(purrr)

df<-read.csv(file = "Matriz_Metabarcoding_new_COIvs18S_edited.csv", sep=";", dec=".")


rownames(df)<-df$Species

df <- df[,-1]
df2<-df

#Remove non-variant cases
df2<-df2[ , which(apply(df2, 2, var) != 0)]

#To check the distance between samples
distance <- get_dist(df2)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))




#check the different clustering patterns
k2 <- kmeans(df2, centers = 2, nstart = 25)
k3 <- kmeans(df2, centers = 3, nstart = 25)
k4 <- kmeans(df2, centers = 4, nstart = 25)
k5 <- kmeans(df2, centers = 5, nstart = 25)
k6 <- kmeans(df2, centers = 6, nstart = 25)
k7 <- kmeans(df2, centers = 7, nstart = 25)


# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df2) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = df2) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df2) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 5")
p5 <- fviz_cluster(k6, geom = "point",  data = df2) + ggtitle("k = 6")
p6 <- fviz_cluster(k7, geom = "point",  data = df2) + ggtitle("k = 6")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)


#check the optimal number of clusters - elbow

set.seed(123)

fviz_nbclust(df2, kmeans, method = "wss")




# function to compute total within-cluster sum of square - longer version
set.seed(123)

wss <- function(k) {
  kmeans(df2, k, nstart = 10 )$tot.withinss
}

k.values <- 1:10

wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(df2, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(df2))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:10

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")
fviz_nbclust(df2, kmeans, method = "silhouette")



# Compute k-means clustering with best fitting k (in this case k = 2)
#k2 <- kmeans(df2, centers = 2, nstart = 25)
#better save as landscape
fviz_cluster(k2, data = df2)+theme_bw()


# Hierarchical clustering
# ++++++++++++++++++++++++
# Use hcut() which compute hclust and cut the tree
hc.cut <- hcut(df2, k = 2, hc_method = "complete")
# Visualize dendrogram
fviz_dend(hc.cut, show_labels = TRUE, rect = TRUE)
# Visualize cluster (just another way to plot the k-means)
fviz_cluster(hc.cut, ellipse.type = "convex")
