fold_split <- function(K=3, # default for now
                       index=NULL){
  data_indices <- data.frame(cbind(id=1:length(index),
                        zip_id=index))
  ziplist <- unique(index)
  sapply(ziplist, function(x){
    temp_data <- data_indices[!is.na(match(data_indices[,2],x)),1]
    n <- length(temp_data)
    kn <- floor(n/K)
    setn <- temp_data
    id <- list()
    i <- 1
    while(i<K){
      idx <- sample(setn,kn, replace=F)
      id[[i]] <- idx
      new.setn = setn[! setn %in% idx]
      setn = new.setn
      i <- i+1
    }
    id[[K]] <- setn
    id
  },simplify = F)
}

# Not Run
# split <- test_train_ind(K=5,index=ind)
# fold_1 <- unlist(lapply(split, function(x) x[[1]]))
# fold_2 <- unlist(lapply(split, function(x) x[[2]]))
# fold_3 <- unlist(lapply(split, function(x) x[[3]]))
# fold_4 <- unlist(lapply(split, function(x) x[[4]]))
# fold_5 <- unlist(lapply(split, function(x) x[[5]]))


