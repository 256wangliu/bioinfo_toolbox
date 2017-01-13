# The script searches the optimal k(s) for k-means and
# outputs the clustering results with the optimal k(s).
# It also draws the diagnostic plot of R squared, AIC,
# and BIC values.


# load libraries in outside of the function
library(flexclust)
library(parallel)

# make it reproducible 
set.seed(144)

opt_kmeans <- function(data_file, kRange=2:30, distFun=distCor, centFun=centMean,
                       name_list=NA, thread_num=16){

    data <- read.table(data_file, row.names=1, as.is=TRUE)

    # if want to filter with names
    if(!is.na(name_list)){
        target_id <- as.character(read.table(name_list, as.is=TRUE)[,1])
        data <- data[target_id,]
    }

    # function for different ks
    kmFam <- kccaFamily(dist=distFun, cent=centFun)
    FUN <- function(k_i){
        res_i <- kcca(data, k_i, family=kmFam)
        list(k=res_i@k, centers=res_i@centers, cluster=res_i@cluster,
             totaldist=res_i@totaldist, distsum=info(res_i, "distsum"),
             cldist=res_i@cldist)
    }

    # parallelization and run
    cl <- makeCluster(thread_num, type="PSOCK")
    clusterCall(cl, function() require("flexclust"))
    res <- parLapply(cl, kRange, FUN)

    # get total distance
    unsorted_k <- unlist(lapply(res, function(x) x$k))
    unsorted_distsum <- unlist(lapply(res, function(x) x$distsum))
    totaldist <- mean(unique(unlist(lapply(res, function(x) x$totaldist))))

    index <- sort(unsorted_k, index.return=TRUE)
    k <- c(1, index$x)
    dist_sum <-c(totaldist, unsorted_distsum[index$ix])

    # compute R squared, AIC, BIC
    rsq <- 1 - (dist_sum*(nrow(data)-1)) / (dist_sum[1]*(nrow(data)-k))    
    aic <- dist_sum + 2 * ncol(data) * k
    bic <- dist_sum + log(nrow(data)) * ncol(data) * k

    # R squared values: estimate optimal k
    vrsq <- (rsq-min(rsq)) / (max(rsq)-min(rsq))
    vk <- (k - min(k)) / (max(k)-min(k))
    dsq <- (vk)^2 + (vrsq - 1)^2
    rsq_opt_k <- k[which.min(dsq)]

    # AIC and BIC: estimate optimal k
    aic_opt_k <- k[which.min(aic)]
    bic_opt_k <- k[which.min(bic)]

    # plot R squared, AIC, and BIC
    prefix <- unlist(strsplit(basename(data_file), "\\."))[1]
    pdf(paste(prefix,"_qc.pdf", sep=""), onefile=FALSE, paper="special", height=3.8, width=10)
    par(mfrow=c(1,3))
    plot(k, rsq, main=paste("R squared (k = ",rsq_opt_k,")", sep=""),
          ylab="R squared", pch=20, cex=1.5, cex.lab=1.2, cex.axis=1.2)
    points(rsq_opt_k, rsq[which.min(dsq)], col="red", pch=20, cex=2)

    plot(k, aic, main=paste("AIC (k = ",aic_opt_k,")", sep=""),
          ylab="AIC", pch=20, cex=1.5, cex.lab=1.2, cex.axis=1.2)
    points(aic_opt_k, aic[which.min(aic)], col="red", pch=20, cex=2)

    plot(k, bic, main=paste("BIC (k = ",bic_opt_k,")", sep=""),
          ylab="BIC", pch=20, cex=1.5, cex.lab=1.2, cex.axis=1.2)
    points(bic_opt_k, bic[which.min(bic)], col="red", pch=20, cex=2)
    dev.off()

    # output clustering results: labels and centers
    opt_k <- unique(c(rsq_opt_k, aic_opt_k, bic_opt_k))
    for(opt_k_i in opt_k){
        i <- which(unsorted_k == opt_k_i)

        label_file <- paste(prefix,"_k",opt_k_i,"_label.txt", sep="")
        out_data <- cbind(res[[i]]$cluster, res[[i]]$cldist[,1])
        write.table(out_data, label_file, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)

        center_file <- paste(prefix,"_k",opt_k_i,"_center.txt", sep="")
        out_data <- res[[i]]$centers
        rownames(out_data) <- 1:opt_k_i
        write.table(out_data, center_file, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
    }
}
