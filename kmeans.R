# basic kmeans
custom.kmeans <- function(dataFile, minK, maxK, method="euclidean", bootstrap=100, seed=1234){
    library(fpc)
    library(cluster)
    library(amap)

    # loading data
    data <- read.table(dataFile, row.names=1, header=T)
    zero <- data[rowSums(data) == 0,]
    data <- data[rowSums(data) > 0,]

    # compute distance
    data.dist <- Dist(data, method=method)

    all.sil <- numeric(0)
    for(k in minK:maxK){
        if(method == "euclidean"){
            kmean.result <- clusterboot(data, B=bootstrap, bootmethod="boot",
                                        clustermethod=kmeansCBI, count=FALSE, krange=k, seed=seed)
        } else {
            kmean.result <- clusterboot(data.dist, B=bootstrap, bootmethod="boot", distances=TRUE,
                                        clustermethod=kmeansCBI, count=FALSE, krange=k, seed=seed)
        }

        # useful kmean result statistics
        # bootmean    = values close to 1 indicate stable clutsers
        # bootrecover = values of recovered members in a cluster
        # bootbrd     = values of dissolved members in a cluster
        k.stats <- cbind(1:k, kmean.result$bootmean, kmean.result$result$result$size,
                         round(kmean.result$bootrecover/bootstrap, 4),
                         round(kmean.result$bootbrd/bootstrap, 4))

        # compute silhouette coefficient
        sil <- silhouette(kmean.result$partition, dist(data))
        all.sil <- cbind(all.sil, sil[,3])

        # kmean data result statistics
        # (1) feature name (2) kmean group (3) neighbor group (4) silhouette coefficient
        data.stats <- cbind(rownames(data), kmean.result$partition, sil[,2], sil[,3])

        write.table(data.stats, paste("K",k,"_cluster.tab",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
        write.table(k.stats, paste("K",k,"_stats.tab",sep=""), quote=F, sep="\t", row.names=F, col.names=F)
    }
    # draw box plot of silhouette coefficient
    postscript("silhouette.eps", onefile=F, horizontal=F, paper="special", width=10.5, height=4.5)
    boxplot(all.sil, col="skyblue", names=minK:maxK, ylab="Silhouette width", xlab="K", boxwex=0.5)
    dev.off()
}