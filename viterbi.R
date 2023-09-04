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
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
argmax_uni_basic <- function(Y_pain , l_pi , l_tpm , l_theta){
  alpha_vec <- beta_vec <- matrix(rep(0,T*S),nrow=T)
  hidden_seq <- rep(0,T)
  for(s in 1:S){
    alpha_vec[1,s] <- l_pi[s] + ifelse(Y_pain[1] == 0 , 0 , l_theta[s,Y_pain[1]])
    beta_vec[T,s] <- 0
  }
  for(t in 2:T){
    for(s in 1:S){
      acc <- beta_temp <- rep(0,S)
      for(s_from in 1:S){
        acc[s_from] <- alpha_vec[t-1,s_from] + l_tpm[s_from,s]
      }
      for(s_to in 1:S){
        beta_temp[s_to] <- l_tpm[s,s_to]+ ifelse(Y_pain[T-t+2] == 0 , 0 , l_theta[s_to,Y_pain[T-t+2]]) + beta_vec[T-t+2,s_to]
      }
      alpha_vec[t,s] = rethinking::log_sum_exp(acc) + ifelse(Y_pain[t] == 0 , 0 , l_theta[s,Y_pain[t]])
      beta_vec[T-t+1,s] = rethinking::log_sum_exp(beta_temp)
    }
  }
  for(t in 1:T){
    hidden_seq[t] <- which.max(alpha_vec[t,] + beta_vec[t,])
  }
  return(hidden_seq)
}
#-----------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------
viterbi_uni_cov_hier_trend <- function(Y_pain , cov , l_pi , l_theta , beta , alpha){
  back_track <- best_logp <- matrix(rep(0,T*S),nrow=T)
  tpm <- l_tpm <- matrix(rep(0,S*S),nrow=S)
  hidden_seq <- rep(0,T)
  k <- 1
  for(s in 1:S){
    best_logp[1,s] <- l_pi[s] + ifelse(Y_pain[1] == 0 , 0 , l_theta[s,Y_pain[1]])
    for(s_from in 1:S){
      if(s_from != s){
        tpm[s_from,s] <- beta[k,]%*%cov
        k<-k+1
      }
    }
  }
  for(t in 2:T){
    k <- 1
    for(s_from in 1:S){
      l_tpm[s_from,] <- tpm[s_from,]-rethinking::log_sum_exp(tpm[s_from,])
    }
    for(s in 1:S){
      best_logp[t,s] <- -Inf;
      logp<-rep(0,S)
      for(s_from in 1:S){
        logp[s_from] <- best_logp[t-1,s_from] + l_tpm[s_from,s]
      }
      best_logp[t,s] = max(logp) + ifelse(Y_pain[t] == 0 , 0 , l_theta[s,Y_pain[t]])
      back_track[t,s] = which.max(logp)
      for(s_from in 1:S){
        if(s_from != s){
          tpm[s_from,s] <-tpm[s_from,s] + (log(t)-log(t-1))*alpha[k]  
          k<-k+1
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


hidden_seq_alt <- function(Y_pain , cov , l_pi , l_theta , beta , alpha){
  alpha_vec <- beta_vec <- matrix(rep(0,T*S),nrow=T)
  tpm <- l_tpm <- matrix(rep(0,S*S),nrow=S)
  hidden_seq <- rep(0,T)
  k <- 1
  for(s in 1:S){
    alpha_vec[1,s] <- l_pi[s] + ifelse(Y_pain[1] == 0 , 0 , l_theta[s,Y_pain[1]])
    beta_vec[T,s] <- 0
    for(s_from in 1:S){
      if(s_from != s){
        tpm[s_from,s] <- beta[k,]%*%cov
        k<-k+1
      }
    }
  }
  for(t in 2:T){
    k <- 1
    for(s_from in 1:S){
      l_tpm[s_from,] <- tpm[s_from,]-rethinking::log_sum_exp(tpm[s_from,])
    }
    for(s in 1:S){
      acc <- beta_temp <- rep(0,S)
      for(s_from in 1:S){
        acc[s_from] <- alpha_vec[t-1,s_from] + l_tpm[s_from,s]
      }
      for(s_to in 1:S){
        beta_temp[s_to] <- l_tpm[s,s_to]+ ifelse(Y_pain[T-t+2] == 0 , 0 , l_theta[s_to,Y_pain[T-t+2]]) + beta_vec[T-t+2,s_to]
      }
      alpha_vec[t,s] = rethinking::log_sum_exp(acc) + ifelse(Y_pain[t] == 0 , 0 , l_theta[s,Y_pain[t]])
      beta_vec[T-t+1,s] = rethinking::log_sum_exp(beta_temp)
      for(s_from in 1:S){
        if(s_from != s){
          tpm[s_from,s] <-tpm[s_from,s] + (log(t)-log(t-1))*alpha[k]  
          k<-k+1
        }
      }
    }
  }
  for(t in 1:T){
    hidden_seq[t] <- which.max(alpha_vec[t,] + beta_vec[t,])
  }
  return(hidden_seq)
}

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

vitervi_mvt_basic <- function(Y_pain , Y_limit , l_pi , l_tpm , l_limit_prob , l_pain_prob , rho){
  back_track <- best_logp <- matrix(rep(0,T*S),nrow=T)
  hidden_seq <- rep(0,T)
  for(s in 1:S){
    bounds <- findBounds(Y_limit[1], Y_pain[1] , l_limit_prob[s,] , l_pain_prob[s,])
    best_logp[1,s] <-  l_pi[s] + cop_ll(bounds[[1]] , bounds[[2]] , rho[1])
  }
  for(t in 2:T){
    for(s in 1:S){
      best_logp[t,s] <- -Inf;
      logp<-rep(0,S)
      for(s_from in 1:S){
        logp[s_from] <- best_logp[t-1,s_from] + l_tpm[s_from,s]
      }
      back_track[t,s] <- which.max(logp)
      bounds <- findBounds(Y_limit[t], Y_pain[t] , l_limit_prob[s,] , l_pain_prob[s,])
      best_logp[t,s] <- max(logp) + cop_ll(bounds[[1]] , bounds[[2]] , rho[t])
      }
  }

  hidden_seq[T] <- which.max(best_logp[T,])
  for(t in 1:(T-1)){
    hidden_seq[T-t] <- back_track[T-t+1,hidden_seq[T-t+1]]
  }
  return(hidden_seq)  
}


max_state_mvt_basic <- function(Y_pain , Y_limit , l_pi , l_tpm , l_limit_prob , l_pain_prob , rho){
  alpha_vec <- beta_vec <- matrix(rep(0,T*S),nrow=T)
  hidden_seq <- rep(0,T)
  for(s in 1:S){
    bounds <- findBounds(Y_limit[1], Y_pain[1] , l_limit_prob[s,] , l_pain_prob[s,])
    alpha_vec[1,s] <- l_pi[s] + cop_ll(bounds[[1]] , bounds[[2]] , rho[1])
    beta_vec[T,s] <- 0
  }
  for(t in 2:T){
    for(s in 1:S){
      acc <- beta_temp <- rep(0,S)
      for(s_from in 1:S){
        acc[s_from] <- alpha_vec[t-1,s_from] + l_tpm[s_from,s]
      }
      for(s_to in 1:S){
        bounds <- findBounds(Y_limit[T-t+2], Y_pain[T-t+2] , l_limit_prob[s_to,] , l_pain_prob[s_to,])
        beta_temp[s_to] <- l_tpm[s,s_to]+ cop_ll(bounds[[1]] , bounds[[2]] , rho[T-t+2]) + beta_vec[T-t+2,s_to]
      }
      bounds <- findBounds(Y_limit[t], Y_pain[t] , l_limit_prob[s,] , l_pain_prob[s,])
      alpha_vec[t,s] <- rethinking::log_sum_exp(acc) + cop_ll(bounds[[1]] , bounds[[2]] , rho[t])
      beta_vec[T-t+1,s] <- rethinking::log_sum_exp(beta_temp)
    }
  }
  for(t in 1:T){
    hidden_seq[t] <- which.max(alpha_vec[t,] + beta_vec[t,])
  }
  return(hidden_seq)
}
