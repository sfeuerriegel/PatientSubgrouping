//------------------------------------------------------------------------------------------------------------------
//Bivariate model
//Pain margin as truncated poisson
//Disability margin as Inflated truncated poisson model
//Different copula models
//------------------------------------------------------------------------------------------------------------------
functions{
//Copula function where different copula approaches are implemented
//-----------------------------------------------------------------------------
//Independence
//return log_u + log_v;
//-----------------------------------------------------------------------------
//Gumbel
//real R;
//  if(log_u == 0.0 || log_v == 0.0){
//    return log_u +log_v;
//  } 
//  else{
//    R = log_sum_exp(theta*log(-log_u),theta*log(-log_v));
//    return -exp(R/theta);
//  }
//-----------------------------------------------------------------------------
//Galambos
//real R;
//  if(log_u == 0.0 || log_v == 0.0 || theta == 0.0){
//    return log_u +log_v;
//  } 
//  else{
//    R = log_sum_exp(-theta*log(-log_u),-theta*log(-log_v));  
//   return log_u + log_v + exp(-R/theta);
//  }
//-----------------------------------------------------------------------------
//BB7
//real log_nu_inv_u;
//  real log_nu_inv_v;
//  real inv_rho = inv(rho);
//  real inv_theta = inv(theta);
//  real R;
//  if(log_u == 0.0 || log_v == 0.0){
//    return log_u + log_v;
//  } else{
//    log_nu_inv_u = log_diff_exp(-rho*log1m_exp(theta*log1m_exp(log_u)),0.0);
//    log_nu_inv_v = log_diff_exp(-rho*log1m_exp(theta*log1m_exp(log_v)),0.0);
//    R = log_sum_exp(log_nu_inv_v,log_nu_inv_u);
//    return log1m_exp(inv_theta*log1m_exp(-inv_rho*log1p_exp(R)));
//  }
//}
//-----------------------------------------------------------------------------
//Tawn type 2 (Joe et al. p.)
//  real log_copula(real log_u , real log_v , real theta , real a_1 , real a_2){
//    real A;
//    real t;
//    if(log_u == 0.0 && log_v ==0.0){
//      return 0.0;
//    } else{
//      t= log_v /(log_v+ log_u);
//      if(t == 1){
//        A = 1 - a_1 + a_2;
//      } else if(t== 0){
//        A = 1- a_2 + a_1;
//      } else{
//        A = (1-a_1)*t + (1-a_2)*(1-t) + pow(pow(a_1*(1-t),theta) + pow(a_2*t,theta),1.0/theta);
//      }
//      return A*(log_u+log_v);
//    }
//  }
//-----------------------------------------------------------------------------
//BB1
//  real log_nu_inv_u;
//  real log_nu_inv_v;
//  real inv_rho = inv(rho);
//  real inv_theta = inv(theta);
//  real R;
//  if(log_u == 0.0 || log_v == 0.0){
//    return log_u + log_v;
//  } else{
//    log_nu_inv_u = rho*log_diff_exp(-theta*log_u,0.0);
//    log_nu_inv_v = rho*log_diff_exp(-theta*log_v,0.0);
//    R = log_sum_exp(log_nu_inv_v,log_nu_inv_u);
//    return -inv_theta*log1p_exp(R*inv_rho);
//  }
  
real log_copula(real log_u , real log_v, real theta ){
real R;
if(log_u == 0.0 || log_v == 0.0){
  return log_u + log_v;
} else{
  R = log_sum_exp(theta*log(-log1m_exp(log_u)), theta*log(-log1m_exp(log_v)));
  return log_diff_exp(log_sum_exp(log_sum_exp(log_u,log_v),-exp(R/theta)),0.0);
}
}
real cop_ll(int Y_1, int Y_2 , vector l_p_1, vector l_p_2 , real theta){
  real LL;
  vector[2] temp;
  if(l_p_1[Y_1+1] == negative_infinity() || l_p_2[Y_2+1] == negative_infinity()){
    return negative_infinity();
  }else{
    if(Y_1 == 0 || Y_2 == 0){
      if(Y_1 == 0 && Y_2 == 0){
        LL = log_copula(l_p_1[1] , l_p_2[1] , theta );  
      } else if(Y_1 == 0 && Y_2 != 0){
        temp[1] = log_copula(l_p_1[1] , l_p_2[Y_2+1] , theta );
        temp[2] = log_copula(l_p_1[1] , l_p_2[Y_2] , theta );
        LL = (temp[1] <= temp[2] ? negative_infinity() : log_diff_exp(temp[1],temp[2]));
      } else{
        temp[1] = log_copula(l_p_1[Y_1+1] , l_p_2[1] , theta );
        temp[2] = log_copula(l_p_1[Y_1] , l_p_2[1] , theta);
        LL = (temp[1] <= temp[2] ? negative_infinity() : log_diff_exp(temp[1],temp[2]));
      }
    } else{
      temp[1] = log_sum_exp(log_copula(l_p_1[Y_1+1] , l_p_2[Y_2+1], theta ) , 
                            log_copula(l_p_1[Y_1] , l_p_2[Y_2], theta ));
      temp[2] = log_sum_exp(log_copula(l_p_1[Y_1+1] , l_p_2[Y_2], theta) , 
                            log_copula(l_p_1[Y_1],l_p_2[Y_2+1], theta));
      LL = (temp[1] <= temp[2] ? negative_infinity() : log_diff_exp(temp[1],temp[2]));
    }
    return LL;
  }
}
real compute_log_lik(int[] Y_pain , int[] Y_dis , int P , int L ,  int T , int S , vector l_pi , vector[] l_tpm , vector[] l_prob_1 , vector[] l_prob_2 , real theta){
  real forw[T,S];
  for(s in 1:S){
    if(Y_pain[1] == P || Y_dis[1] == L){
      forw[1,s] =  l_pi[s]; 
    }else{
      forw[1,s] =  l_pi[s] + cop_ll(Y_pain[1],Y_dis[1],l_prob_1[s],l_prob_2[s] ,theta);
    }
  }
  for(t in 2:T){
    for(s in 1:S){
      vector[S] acc;
      for(s_from in 1:S){
        acc[s_from] = forw[t-1,s_from] + l_tpm[s_from,s];
      }
      if(Y_pain[t] == P || Y_dis[t] == L){
        forw[t,s] = log_sum_exp(acc);  
      }else{
        forw[t,s] = log_sum_exp(acc) + cop_ll(Y_pain[t],Y_dis[t],l_prob_1[s],l_prob_2[s] ,theta);
      }
    }
  }
  return log_sum_exp(forw[T]);
}
int[] viterbi(int[] Y_pain , int[] Y_dis , int P , int L ,  
              int T , int S , vector l_pi , vector[] l_tpm , 
              vector[] l_prob_1 , vector[] l_prob_2 , real theta){
    int hidden_seq[T];
    {
      int back_track[T,S];
      real best_logp[T,S];
      real log_p_max;
      for(s in 1:S){
        if(Y_pain[1] == P || Y_dis[1] == L){
          best_logp[1,s] = l_pi[s];  
        }else{
          best_logp[1,s] = l_pi[s] + cop_ll(Y_pain[1],Y_dis[1],l_prob_1[s],l_prob_2[s] ,theta);
        }
      }
      for(t in 2:T){
        for(s in 1:S){
          best_logp[t,s] = negative_infinity();
          for(s_from in 1:S){
            real logp;
            if(Y_pain[t] == P || Y_dis[t] == L){
              logp = best_logp[t-1,s_from] + l_tpm[s_from,s];  
            }else{
              logp = best_logp[t-1,s_from] + l_tpm[s_from,s] + cop_ll(Y_pain[t],Y_dis[t],l_prob_1[s],l_prob_2[s] ,theta);
            }
            if(logp > best_logp[t,s]){
              back_track[t,s] = s_from;
              best_logp[t,s] = logp;
            }
          }
        }
      }
      log_p_max = max(best_logp[T]);
      for(s in 1:S){
        if(best_logp[T,s] == log_p_max){
          hidden_seq[T] = s;
        } 
      }
      for(t in 1:(T-1)){
        hidden_seq[T-t] = back_track[T-t+1,hidden_seq[T-t+1]];
      }
    }
    return hidden_seq;
  }

}

data {
  int<lower=1> S; // number of states
  int<lower=1> K; // number of clusters
  int<lower=1> T; // length of timeseries
  int<lower=1> P; // number of emission levels - pain
  int<lower=1> L; // number of emission levels - disability
  int<lower=1> N_train; // number of patients in training set
  int<lower=1> N_test; // number of patients in test set
  int<lower=1> nCovs;
  int<lower=0,upper=P> Y_pain_train[N_train,T]; // pain level of patient n at time t level;equals 0 if original data was missing
  int<lower=0,upper=P> Y_pain_test[N_test,T]; // pain level of patient n at time t level;equals 0 if original data was missing
  int<lower=0,upper=L> Y_dis_train[N_train,T]; // disability level of patient n at time t level;equals 0 if original data was missing
  int<lower=0,upper=L> Y_dis_test[N_test,T]; // disability level of patient n at time t level;equals 0 if original data was missing
  matrix[N_train,nCovs] Q; // thinned Q matrix of Covariates for training
  matrix[nCovs,nCovs] R_mat; // thinned R matrix corresponding to Q matrix
  matrix[N_test,nCovs] X_test;
}
transformed data{
  vector[S] dis_mean;
  vector[S] pain_mean;
  if(S>1){
    for(s in 1:S){
      dis_mean[s] = 0 + (5.0/(S-1))*(s-1);
      pain_mean[s] = 0 + (5.0/(S-1))*(s-1);
    }
  } else{
    dis_mean[1] = 3.0;
    pain_mean[1] = 3.0;
  }
}
parameters {
  simplex[S] pi[K] ; //initial state distribution for each patient
  simplex[S] tpm[K,S] ;  //Transition probability matrix for each patient, time independent
  vector<lower=0>[S] lambda_pain[K]; // state specific coefficients for pain lvl distribution
  vector<lower=0>[S] lambda_dis[K];
  vector<lower=0>[K] theta_raw;
  vector[S] inflation_raw[K];
  matrix[nCovs,K-1] beta_tilde;//covariate slopes
  row_vector[K-1] alpha;
}
transformed parameters{
  vector[N_train] log_lik;
  matrix[N_train,K] prior_cluster_train;
  matrix[N_train,K] posterior_cluster_train;
  {
    // log transforms of parameters to avoid duplicative evaluations for each patient
    vector[S] l_pi[K]=log(pi);
    vector[S] l_tpm[K,S]=log(tpm);
    vector[P] l_cum_pain[K,S];
    vector[L] l_cum_dis[K,S];
    vector[K] theta = 1.0 +theta_raw;
    for(k in 1:K){
      for(s in 1:S){
          for(p in 0:(P-1)){
        //First emission vector stores log cummulative probability for a truncated poisson
            l_cum_pain[k,s,p+1] = (p == P-1 ? 0.0 : (poisson_lcdf(p|lambda_pain[k,s]) - 
                                              poisson_lcdf(P-1|lambda_pain[k,s])));
          }                                              
            for(d in 0:(L-1)){
            //Second emission vector stores log cummulative probability for a truncated poisson inflated at L-1
              l_cum_dis[k,s,d+1] = (d == L-1 ? 0.0 : log1m_inv_logit(inflation_raw[k,s]) + poisson_lcdf(d|lambda_dis[k,s]) - poisson_lcdf(L-2|lambda_dis[k,s]));
            }
      }
    }
    prior_cluster_train[,1] = rep_vector(0.0,N_train);
    for(n in 1:N_train){
      prior_cluster_train[n,2:K] = Q[n] * beta_tilde + alpha;
      prior_cluster_train[n] = to_row_vector(log_softmax(to_vector(prior_cluster_train[n])));
      for(k in 1:K){
        posterior_cluster_train[n,k] =  compute_log_lik(Y_pain_train[n], Y_dis_train[n] , P , L , T , S , 
                                                             l_pi[k] , l_tpm[k]  , l_cum_pain[k] , l_cum_dis[k] , theta[k]);
    }
    posterior_cluster_train[n] = posterior_cluster_train[n] + prior_cluster_train[n];
    
    log_lik[n] = log_sum_exp(posterior_cluster_train[n]);
  }
}
}
model {
  for(k in 1:K){
    pi[k] ~ dirichlet(rep_vector(1,S));
    inflation_raw[k] ~ normal(0,2);
    theta_raw[k] ~ normal(0,5);
    for(s in 1:S){
      tpm[k,s] ~ dirichlet(rep_vector(1,S));
      lambda_pain[k,s] ~ normal(pain_mean[s],1);
      lambda_dis[k,s] ~ normal(dis_mean[s],1);
    }
  }
  alpha ~ normal(0,5);
  to_vector(beta_tilde) ~ normal(0,1);
  target += sum(log_lik);
}
generated quantities{
    matrix[nCovs,K-1] beta = (R_mat\beta_tilde);
    matrix[N_test,K] prior_cluster_test;
    matrix[N_test,K] posterior_cluster_test;
    vector[N_test] log_lik_test;
    {
    vector[S] l_pi[K]=log(pi);
    vector[S] l_tpm[K,S]=log(tpm);
    vector[P] l_cum_pain[K,S];
    vector[L] l_cum_dis[K,S];
    vector[K] theta = 1.0 +theta_raw;
    for(k in 1:K){
      for(s in 1:S){
          for(p in 0:(P-1)){
        //First emission vector stores log cummulative probability for a truncated poisson
            l_cum_pain[k,s,p+1] = (p == P-1 ? 0.0 : (poisson_lcdf(p|lambda_pain[k,s]) - 
                                              poisson_lcdf(P-1|lambda_pain[k,s])));
          }                                              
            for(d in 0:(L-1)){
            //Second emission vector stores log cummulative probability for a truncated poisson inflated at L-1
              l_cum_dis[k,s,d+1] = (d == L-1 ? 0.0 : log1m_inv_logit(inflation_raw[k,s]) + poisson_lcdf(d|lambda_dis[k,s]) - poisson_lcdf(L-2|lambda_dis[k,s]));
            }
      }
    }
    prior_cluster_test[,1] = rep_vector(0.0,N_test);
    for(n in 1:N_test){
      prior_cluster_test[n,2:K] = X_test[n] * beta + alpha;
      prior_cluster_test[n] = to_row_vector(log_softmax(to_vector(prior_cluster_test[n])));
      for(k in 1:K){
        posterior_cluster_test[n,k] =  compute_log_lik(Y_pain_test[n], Y_dis_test[n] , P , L , T , S , 
                                                             l_pi[k] , l_tpm[k]  , l_cum_pain[k] , l_cum_dis[k] , theta[k]);
    }
    posterior_cluster_test[n] = posterior_cluster_test[n] + prior_cluster_test[n];
    
    log_lik_test[n] = log_sum_exp(posterior_cluster_test[n]);
  }  
  }
}
   
   

