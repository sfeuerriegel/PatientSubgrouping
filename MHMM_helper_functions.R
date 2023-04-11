#Function to compute both information criteria 
#and lppd for MHMM models to choose the appropriate number of Clusters

model_kpi <- function(ll_train = NULL ,ll_val=NULL){
  # ll_train is a NxS matrix which stores the log-likelihood for each observation n in 1,...,N and each samples s in 1,...,S from the training set
  # ll_test is a NxS matrix which stores the log-likelihood for each observation n in 1,...,N and each samples s in 1,...,S from the test set
  waic_res <- waic(ll_train)
  loo_res<-loo(ll_train)
  lpd_train <- sum(log(apply(exp(ll_train),2,mean)))
  
  
  if(is.null(ll_val)){
    lpd_val <- NA
  } else{
    lpd_val <- sum(log(apply(exp(ll_val),2,mean)))
  }
  results <- data.frame( WAIC = waic_res$waic , 
                         LOOIC = loo_res$looic , 
                         LL_Train = lpd_train , 
                         LL_Val = lpd_val)
  return(results)
}
#-----------------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------
viterbi_uni_basic <- function(Y_pain , l_pi , l_tpm , l_theta){
  back_track <- best_logp <- matrix(rep(0,T*S),nrow=T)
  hidden_seq <- rep(0,T)
  for(s in 1:S){
    best_logp[1,s] <- l_pi[s] + ifelse(Y_pain[1] == 0 , 0 , l_theta[s,Y_pain[1]])
  }
  for(t in 2:T){
    for(s in 1:S){
      best_logp[t,s] <- -Inf;
      for(s_from in 1:S){
        logp <- best_logp[t-1,s_from] + l_tpm[s_from,s] + ifelse(Y_pain[t] == 0 , 0 , 
                                                                 l_theta[s,Y_pain[t]])
        if(logp > best_logp[t,s]){
          back_track[t,s] <- s_from
          best_logp[t,s] <- logp
        }
      }
    }
  }
  hidden_seq[T] <- which.max(best_logp[T,])
  for(t in 1:(T-1)){
    hidden_seq[T-t] <- back_track[T-t+1,hidden_seq[T-t+1]]
  }
  return(hidden_seq)
}
#------------------------------------------------------------------------------------------------------
spike_count <- function(y,lag,threshold,influence) {
  y<-y[!is.na(y)]
  if(length(y)<lag+1){
    return(0)
  }
  signals <- rep(0,length(y))
  filteredY <- y[1:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[1:lag])
  stdFilter[lag] <- sd(y[1:lag])
  for (i in (lag+1):length(y)){
    if (y[i]-avgFilter[i-1] > threshold*stdFilter[i-1]) {
      signals[i] <- 1;
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i])
    stdFilter[i] <- sd(filteredY[(i-lag):i])
  }
  return(sum(sapply(2:length(y),function(k){ifelse(signals[k]==1 & signals[k-1]==0,1,0)})))
}
