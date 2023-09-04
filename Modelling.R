library(rstan)
library(tidyverse)
rstan_options(auto_write = TRUE)
options(mc.cores = 8)
#Preparing stan data

load(file="Preprocessed.MHMM.Data.Rdata")


QR <- qr(risk_factors[train_idx,])

data_stan <- list(S = NA,
                  K = NA,
                  T = T,
                  P = 11,
                  L = 8,
                  N_train = N,
                  N_test = nrow(risk_factors) - N,
                  nCovs = ncol(risk_factors),
                  Y_pain_train = Y_pain[train_idx,],
                  Y_pain_test = Y_pain[-train_idx,],
                  Y_dis_train = Y_limit[train_idx,],
                  Y_dis_test = Y_limit[-train_idx,],
                  Q = qr.Q(QR)*sqrt(N-1),
                  R_mat = qr.R(QR)/sqrt(N-1),
                  X_test = risk_factors[-train_idx,])

extract_int_para <-function(chain_id,map,S,K,nCovs){
  pi<-lambda_pain<-lambda_dis<-inflation_raw<-array(0,dim=c(K,S))
  tpm<-array(0,dim=c(K,S,S))
  alpha<-array(0,dim=c(K-1))
  theta_raw<-rep(0,K)
  beta_tilde <- array(0,dim=c(nCovs,K-1))
  for(k in 1:K){
    theta_raw[k] <- map$par[paste0("theta_raw[",k,"]")]
    for(s in 1:S){
      pi[k,s] <- map$par[paste0("pi[",k,",",s,"]")]
      lambda_dis[k,s] <- map$par[paste0("lambda_dis[",k,",",s,"]")]
      lambda_pain[k,s] <- map$par[paste0("lambda_pain[",k,",",s,"]")]
      inflation_raw[k,s] <- map$par[paste0("inflation_raw[",k,",",s,"]")]
      for(s_from in 1:S){
        tpm[k,s_from,s] <- map$par[paste0("tpm[",k,",",s_from,",",s,"]")]
      }
    }
  }
  for(k in 1:(K-1)){
    alpha[k] <- map$par[paste0("alpha[",k,"]")]
    for(m in 1:nCovs){
      beta_tilde[m,k] <- map$par[paste0("beta_tilde[",m,",",k,"]")]
    }
  }
  return(list(pi = pi,
              lambda_pain = lambda_pain,
              lambda_dis = lambda_dis,
              inflation_raw = inflation_raw,
              tpm = tpm,
              alpha = alpha,
              beta_tilde = beta_tilde))
}

data_stan$S<-3
model <- rstan::stan_model(file="StanCode/mhmm_copula.stan")
for(cluster in 7:9){
  print("---------------------------------")
  print(cluster)
  print("---------------------------------")
  data_stan$K<-cluster
  loop=TRUE
  while(loop){
    mhmm_map <- optimizing(model, data = data_stan,verbose =TRUE)
    if(mhmm_map$return_code == 0){
      loop = FALSE
    }
  }
  init_ll <- lapply(1:1, function(id) extract_int_para(chain_id = id , map = mhmm_map , data_stan$S,data_stan$K,data_stan$nCovs))
  mhmm_model <- stan(file="StanCode/mhmm_copula.stan", 
                  chains = 1, 
                  init = init_ll,
                  iter = 1000, 
                  data = data_stan,
                  refresh = 10,
                  control = list(max_treedepth = 8))
  save(mhmm_model, file=paste0("ModelFits/mhmm_3_",cluster,".R"))
}
rm(mhmm_model,model,loop,cluster)
#################################################################################################################
#State Check
#################################################################################################################
for(state in 1:4){
  print(state)
  data_stan$S<-state
  mhmm_model <- stan(file="StanCode/hmm_state.stan", 
                     chains = 1,
                     iter = 1000, 
                     data = data_stan,
                     refresh = 10,
                     control = list(max_treedepth = 8))
  save(mhmm_model, file = paste0("ModelFits/hmmm_",state,".R"))
}
rm(state,mhmm_model)
#################################################################################################################
#Model with C=8 S=1 als robustness check
#################################################################################################################

data_stan$S<-1
data_stan$K<-8
model <- rstan::stan_model(file="StanCode/mhmm_copula.stan")
loop=TRUE
  while(loop){
    mhmm_map <- optimizing(model,init="0", data = data_stan,verbose =TRUE, algorithm = "BFGS", as_vector =FALSE)
    if(mhmm_map$return_code == 0){
      loop = FALSE
    }
  }
  init_ll <- lapply(1:1, function(id) extract_int_para(chain_id = id , map = mhmm_map , data_stan$S,data_stan$K,data_stan$nCovs))
  mhmm_model <- stan(file="StanCode/mhmm_copula.stan", 
                     chains = 1, 
                     init = "0",
                     iter = 1000, 
                     data = data_stan,
                     refresh = 10,
                     control = list(max_treedepth = 8))
  save(mhmm_model, file="ModelFits/mhmm_1_8.R")
 
