dist <- function(x,y){
  #return(max(abs(x-y)))
  return(sqrt(mean((x-y)^2)))
  #return(sum(abs(x-y))/length(x))
}
#-------------------------------------------------------------------------------------------------------
db_star <- function(part,C,traj){
  centroids <- t(sapply(1:C,function(k){apply(traj[part==k,],2,mean)}))
  delta <- sapply(1:C , 
                  function(k){mean(sapply(which(part==k) , 
                                          function(i){dist(traj[i,],centroids[k,])}))})
  sum <- rep(0,C)
  for(c in 1:C){
    sum[c] <- delta[c] + max(delta[-c])
    temp <- sapply(1:C,function(l){dist(centroids[c,],centroids[l,])})
    sum[c] <- sum[c]/min(temp[-c])
  }
  return(mean(sum))
}
#-------------------------------------------------------------------------------------------------------
COP <- function(part,C,traj){
  centroids <- t(sapply(1:C,function(k){apply(traj[part==k,],2,mean)}))
  delta <- sapply(1:C , 
                  function(k){mean(sapply(which(part==k) , 
                                          function(i){dist(traj[i,],centroids[k,])}))})
  summ<-rep(0,C)
  for(c in 1:C){
    summ[c] <- delta[c]*sum(part==c)
    summ[c] <- summ[c]/min(sapply(which(part!=c) , 
                                  function(i){max(sapply(which(part==c) ,                                                                     function(j){dist(traj[i,],traj[j,])}))}))
  }
  return(sum(summ)/nrow(traj))
}
#-------------------------------------------------------------------------------------------------------
sym33 <- function(part,C,traj){
  centroids <- t(sapply(1:C,function(k){apply(traj[part==k,],2,mean)}))
  nom <- min(sapply(1:C,function(k){min(sapply((1:C)[-k],function(l){
    sum(sapply(which(part==l),function(i){sum(sapply(which(part==k),function(j){
      dist(traj[j,],traj[i,])
    }))}))/sum(part==l)/sum(part==k)}))}))
  denom <- max(sapply(1:C , 
                      function(k){sum(sapply(which(part==k) , 
                                             function(i){sort(sapply(which(part==k) , 
                                                                     function(j){dist(2*centroids[k,]-traj[i,],traj[j,])}))[1:2]}))/sum(part==k)}))
  return(nom/denom)
}
#-------------------------------------------------------------------------------------------------------
gdi33 <- function(part,C,traj){
  centroids <- t(sapply(1:C,function(k){apply(traj[part==k,],2,mean)}))
  nom <- min(sapply(1:C,function(k){min(sapply((1:C)[-k],function(l){
    sum(sapply(which(part==l),function(i){sum(sapply(which(part==k),function(j){
      dist(traj[j,],traj[i,])
    }))}))/sum(part==l)/sum(part==k)}))}))
  denom <- max(sapply(1:C , 
                      function(k){2*sum(sapply(which(part==k) , 
                                               function(i){dist(centroids[k,],traj[i,])}))/sum(part==k)}))
  return(nom/denom)
}
#-------------------------------------------------------------------------------------------------------
SDbw <- function(part, C, traj) {
  centroids <-
    t(sapply(1:C, function(k) {
      apply(traj[part == k, ], 2, mean)
    }))
  total_sd <- sqrt(sum(apply(traj, 2, function(v){var(v)*(length(v)-1)/length(v)}) ^ 2))
  cluster_sd <-
    sapply(1:C, function(k) {
      sqrt(sum(apply(traj[part == k, ], 2, function(v){var(v)*(length(v)-1)/length(v)}) ^ 2))
    })
  sd_c <- sqrt(sum(cluster_sd)) / C
  
  part_1 <- mean(cluster_sd) / total_sd
  
  part_2 <- sum(sapply(1:C , function(k) {
    sum(sapply((1:C)[-k] ,
               function(l) {
                 sum(sapply(which(part %in% c(l, k)) , function(i) {
                   ifelse(dist(traj[i, ], (centroids[k, ] + centroids[l, ]) / 2) <= sd_c, 1, 0)
                 })) / max(c(sum(sapply(which(part == k), function(i) {
                   ifelse(dist(traj[i, ], centroids[k, ]) <= sd_c, 1, 0)
                 })) , sum(sapply(which(part == l), function(j) {
                   ifelse(dist(traj[j, ], centroids[l, ]) <= sd_c, 1, 0)
                 }))))
                 
               }))
  }))
  
  return(part_1 + part_2 / (C * (C - 1)))
}