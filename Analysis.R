library(ggplot2)
library(rstan)
library(xtable)
library(rethinking)
library(dplyr)
library(tidyr)
library(loo)
library(clusterCrit)
library(flexclust)
###################################################################################################################################
#Some descriptive statistics
###################################################################################################################################
#Correlation between pain/disability at 3m
load(file="Preprocessed.MHMM.Data.Rdata")
disability <- lbp.data_filtered[,c("rmprop2u","rmprop3","rmprop12m")]
pain <- lbp.data_filtered[,c("vasl2u","vasl3m","vasl12m")]
#Correlation  2 weeks  0.5470066***
cor.test(lbp.data_filtered$vasl2u,lbp.data_filtered$nrs2,method = "pearson",use="pairwise")
#Correlation  3 months 0.664574***
cor.test(lbp.data_filtered$vasl3m,lbp.data_filtered$nrs12,method = "pearson",use="pairwise")
#Correlation 12 months 0.5972775***
cor.test(lbp.data_filtered$vasl12m,lbp.data_filtered$nrs52,method = "pearson",use="pairwise")
#Duration of LBP
# 0-2w   | 2-4w   | 1-3m   | >3m
# 65.25% | 13.69% | 10.68% | 13.09% 
prop.table(table(lbp.data_filtered$dlva0))
#Number of non-missing responses per patient 47.82
(sum(Y_pain>=0)-sum(Y_pain==11))/nrow(Y_pain)
#Number of patients with no missing replies 378
sum(apply(Y_pain,1,function(v){sum(v==11)})==0)

rm(list=ls())

###################################################################################################################################
#Determination of correct number of latent states
###################################################################################################################################
STATES <- 1:4
results_df <- data.frame(LOOIC = rep(0,2*length(STATES)), 
                         WAIC = rep(0,2*length(STATES)), 
                         LPD_TRAIN = rep(0,2*length(STATES)), 
                         LPD_VAL = rep(0,2*length(STATES)))
for(S in STATES){
  load(file=paste0("ModelFits/hmmm_",S,".R"))
  samples<-rstan::extract(mhmm_model,c("log_lik","log_lik_test"))
  t_res <-loo(samples$log_lik)
  results_df$LOOIC[2*(S-1)+1]<-t_res$looic
  results_df$LOOIC[2*(S-1)+2]<-t_res$se_looic
  t_res <-waic(samples$log_lik)
  results_df$WAIC[2*(S-1)+1]<-t_res$waic
  results_df$WAIC[2*(S-1)+2]<-t_res$se_waic
  results_df$LPD_TRAIN[2*(S-1)+1]<--2*sum(log(apply(exp(samples$log_lik),2,mean)))
  results_df$LPD_TRAIN[2*(S-1)+2]<- 2*sqrt(425*var(log(apply(exp(samples$log_lik),2,mean))))
  results_df$LPD_VAL[2*(S-1)+1]<--2*sum(log(apply(exp(samples$log_lik_test),2,mean)))
  results_df$LPD_VAL[2*(S-1)+2]<- 2*sqrt(422*var(log(apply(exp(samples$log_lik_test),2,mean))))
  rm(mhmm_model,samples,t_res)
  gc()
}
results_df<-matrix(formatC(unlist(results_df), 
                           digits = 2, 
                           format = "f"), 
                   nrow = 2*length(STATES))
results_df[seq(from=2,to=8,by=2),] <- paste0("(",results_df[seq(from=2,to=8,by=2),],")")
first_col <- rep("",2*length(STATES))
first_col[seq(1,7,by=2)] <- as.character(STATES)
print.xtable(xtable(cbind(first_col,results_df)) , 
             sanitize.text.function = identity , 
             include.rownames = FALSE , 
             include.colnames = FALSE ,
             type = "latex" , 
             hline.after = NULL,
             booktabs = TRUE , 
             floating = FALSE ,
             only.contents = TRUE,
             file="number_of_states.tex")
rm(list=ls())
###################################################################################################################################
#Determination of correct numbers of Clusters
###################################################################################################################################


S<- 3
CLUSTERS <- 4:9
results_df <- data.frame(LOOIC = rep(0,2*length(CLUSTERS)), 
                         WAIC = rep(0,2*length(CLUSTERS)), 
                         LPD_TRAIN = rep(0,2*length(CLUSTERS)), 
                         LPD_VAL = rep(0,2*length(CLUSTERS)))


for(C in CLUSTERS){
  load(file=paste0("ModelFits/mhmm_",S,"_",C,".R"))
  samples<-rstan::extract(mhmm_model,c("log_lik","log_lik_test"))
  t_res <-loo(samples$log_lik)
  results_df$LOOIC[2*(C-4)+1]<-t_res$looic
  results_df$LOOIC[2*(C-4)+2]<-t_res$se_looic
  t_res <-waic(samples$log_lik)
  results_df$WAIC[2*(C-4)+1]<-t_res$waic
  results_df$WAIC[2*(C-4)+2]<-t_res$se_waic
  results_df$LPD_TRAIN[2*(C-4)+1]<--2*sum(log(apply(exp(samples$log_lik),2,mean)))
  results_df$LPD_TRAIN[2*(C-4)+2]<- 2*sqrt(425*var(log(apply(exp(samples$log_lik),2,mean))))
  results_df$LPD_VAL[2*(C-4)+1]<--2*sum(log(apply(exp(samples$log_lik_test),2,mean)))
  results_df$LPD_VAL[2*(C-4)+2]<- 2*sqrt(422*var(log(apply(exp(samples$log_lik_test),2,mean))))
  rm(mhmm_model,samples,t_res)
  gc()
}
results_df<-matrix(formatC(unlist(results_df), 
                           digits = 2, 
                           format = "f"), 
                   nrow = 2*length(CLUSTERS))
results_df[seq(from=2,to=12,by=2),] <- paste0("(",results_df[seq(from=2,to=12,by=2),],")")
first_col <- rep("",2*length(CLUSTERS))
first_col[seq(1,11,by=2)] <- as.character(CLUSTERS)
print.xtable(xtable(cbind(first_col,results_df)) , 
             sanitize.text.function = identity , 
             include.rownames = FALSE , 
             include.colnames = FALSE ,
             type = "latex" , 
             hline.after = NULL,
             booktabs = TRUE , 
             floating = FALSE ,
             only.contents = TRUE,
             file="number_of_clusters.tex")
rm(list=ls())

###################################################################################################################################
#Comparison of posterior MHMM clustering to LCA and k-means clustering
###################################################################################################################################
load(file="Preprocessed.MHMM.Data.Rdata")
source("MHMM_helper_functions.R")

C <- 8
S <- 3
load(paste0("ModelFits/mhmm_",S,"_",C,".R"))
samples <- rstan::extract(mhmm_model,c("posterior_cluster_test","log_lik_test"))
for(c in 1:C){
  samples$posterior_cluster_test[,,c] <- samples$posterior_cluster_test[,,c] - samples$log_lik_test
}
map_idx_test <- which.max(apply(samples$log_lik_test,1,sum))
mhmm_cluster_test<-apply(samples$posterior_cluster_test[map_idx_test,,],1,which.max)
rm(samples,mhmm_model,c,map_idx_test,PreProcess)





##Comparison to other clustering techniques


part_SS <- as.integer(lbp.data_filtered$pp_single_modal[-train_idx])
part_TS <- as.integer(lbp.data_filtered$pp_2stage_modal[-train_idx])
max_criteria <- c("Silhouette", "Calinski_Harabasz" , "GDI33" , "GDI43" , "GDI53")
min_criteria <- c("Davies_Bouldin")
set.seed(2208)
kmeans_cl <- cclust(risk_factors[train_idx,],8)
kmeans_clus <- predict(kmeans_cl,risk_factors[-train_idx,])
crit<-list()
for(k in 1:2){
obs <- switch(k, 
              matrix(as.numeric(Y_pain[-train_idx,]),nrow=422),
              matrix(as.numeric(Y_limit[-train_idx,]),nrow=422))
non_NA<-ifelse(apply(obs,1,function(v){sum(v==11)})==0,TRUE,FALSE)
obs <- obs[non_NA,]
cluster<-list(mhmm_cluster_test[non_NA] , part_TS[non_NA] , part_SS[non_NA] , kmeans_clus[non_NA])
crit[[k]]<-as.data.frame(sapply(1:length(cluster) , 
                           function(x){as.numeric(intCriteria(obs , 
                                                              cluster[[x]] , 
                                                              max_criteria))}))
crit[[k]] <- rbind(crit[[k]] , 
              sapply(1:length(cluster) , 
                     function(x){as.numeric(sym33(cluster[[x]] , 
                                                  length(table(cluster[[x]])) , 
                                                  obs))}))
crit[[k]]<-rbind(crit[[k]] , 
            sapply(1:length(cluster) , 
                   function(x){as.numeric(intCriteria(obs , 
                                                      cluster[[x]] , 
                                                      min_criteria))}))              
crit[[k]] <- rbind(crit[[k]] , 
              sapply(1:length(cluster) , 
                     function(x){as.numeric(db_star(cluster[[x]] , 
                                                    length(table(cluster[[x]])) , 
                                                    obs))}))
crit[[k]] <- rbind(crit[[k]] , 
              sapply(1:length(cluster) , 
                     function(x){as.numeric(COP(cluster[[x]] , 
                                                length(table(cluster[[x]])) , 
                                                obs))}))
crit[[k]] <- rbind(crit[[k]] , 
              sapply(1:length(cluster) , 
                     function(x){as.numeric(SDbw(cluster[[x]] , 
                                                 length(table(cluster[[x]])) , 
                                                 obs))}))

rownames(crit[[k]]) <-c(max_criteria,"Sym33",min_criteria,"DB*","COP","SD_bw")
colnames(crit[[k]]) <- c("MHMM with posterior assignment","LCA (two-stage)","LCA (single-stage)","$K$-means clustering")
}
results_df <- rbind(t(crit[[1]]),t(crit[[2]]))
results_df <- formatC(results_df,digits = 2,format = "f")
print.xtable(xtable(results_df) , 
             sanitize.text.function = identity , 
             include.rownames = TRUE , 
             include.colnames = FALSE ,
             type = "latex" , 
             hline.after = c(2) , 
             booktabs = TRUE , 
             floating = FALSE ,
             only.contents = TRUE,
             file="CVI_comparison.tex")

rm(list=ls())
###################################################################################################################################
#KDE for prior and posterior cluster probabilities 10%,50% and 90% quantiles added
###################################################################################################################################
C<- 8
S<-3
load(paste0("ModelFits/mhmm_",S,"_",C,".R"))
samples <- rstan::extract(mhmm_model,c("posterior_cluster_test","prior_cluster_test","log_lik_test"))
for(c in 1:C){
  samples$posterior_cluster_test[,,c] <- samples$posterior_cluster_test[,,c] - samples$log_lik_test
}
map_idx_test <- which.max(apply(samples$log_lik_test,1,sum))

df <- data.frame(Prob = 100*c(apply(exp(samples$prior_cluster_test[map_idx_test,,]),1,max),apply(exp(samples$posterior_cluster_test[map_idx_test,,]),1,max)),
                          Type = rep(c("Offline","Online: 52 weeks"),each=422))

df2 <- df %>% group_by(Type) %>% summarize(Q10 = quantile(Prob,0.1),Q50 = quantile(Prob,0.5),Q90 = quantile(Prob,0.9))

  ggplot() + 
  geom_density(data=df, mapping = aes(x=Prob),size=1.25) + 
  geom_hline(yintercept=0,color="white",size=1.25)+  
  geom_vline(data = df2, mapping = aes(xintercept =Q10),size=1.25,linetype="dotted") +
  geom_vline(data = df2, mapping = aes(xintercept =Q50),size=1.25,linetype="dotted") +  
  geom_vline(data = df2, mapping = aes(xintercept =Q90),size=1.25,linetype="dotted") +
  facet_wrap(~Type, scales = "free_y") + 
  theme_bw() +
  xlab("Maximum cluster probability") + 
  ylab("Density") +
  theme(text = element_text(size=20)) + 
  scale_x_continuous(expand = c(0, 10)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 20),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))


rm(list=ls())

###################################################################################################################################
#Cost Analysis
###################################################################################################################################
C<- 8
S<- 3
load(paste0("ModelFits/mhmm_",S,"_",C,".R"))
samples <- rstan::extract(mhmm_model,c("posterior_cluster_test","prior_cluster_test","log_lik_test"))
for(c in 1:C){
  samples$posterior_cluster_test[,,c] <- samples$posterior_cluster_test[,,c] - samples$log_lik_test
}
rm(c)
#Posterior cluster after observing the complete trajectory == TRUE Clustering
mhmm_cluster_posterior<-apply(apply(exp(samples$posterior_cluster_test),c(2,3),median),1,which.max)

cluster_probs<-list()
#Prior clustering without any observationy Y
cluster_probs[[1]]<-samples$prior_cluster_test
rm(samples,mhmm_model)
#Online clustering after 4 and 26 weeks
t_vec<-c(5,10,15)
load(file="ModelFits/time_cluster_prob.R")
for(k in t_vec){
  cluster_probs[[length(cluster_probs)+1]] <- time_dependent_cluster[,k,,]
}
rm(time_dependent_cluster)
#Calculating different decision boundaries at 50,55,60,65,70,75,80,85 % cluster probability
cut_offs<-seq(from = 0.5, to = 0.85, by = 0.05)
###############################################################################################
#Currently
###############################################################################################
df <- data.frame(x = numeric(), y= numeric(), Probs = character(), Type=character())
alpha <- seq(from=0,to=1,by=0.01)
arrow_start <- rep(0,length(cluster_probs))
beta_at_85 <- rep(0,length(cluster_probs))
for(k in 1:length(cluster_probs)){
  clustering <-apply(apply(exp(cluster_probs[[k]]),c(2,3),median),1,which.max)
  probabilities <- apply(apply(exp(cluster_probs[[k]]),c(2,3),median),1,max)
  boundaries<-array(0,dim=c(length(cut_offs),length(alpha)))
  for(i in 1:length(cut_offs)){
    patient_idx <- which(probabilities>=cut_offs[i])
    x<-sum(diag(table(factor(clustering[patient_idx],levels = 1:8), 
                      factor(mhmm_cluster_posterior[patient_idx],levels=1:8))))/length(patient_idx)
    boundaries[i,]<-(1-alpha*x)/(1-x)-1
  }
  arrow_start[k] <- (1-x)/x
  beta_at_85[k] <- (1-0.85*x)/(1-x)-1
  df_poly <- data.frame(x=rep(alpha,length(cut_offs)),
                        y=as.vector(t(boundaries)), 
                        Probs=factor(as.character(100*rep(cut_offs,each=length(alpha)))),
                        Type = rep(switch(k,"Offline","Online: 5 weeks", "Online: 10 weeks","Online: 15 weeks"),length(alpha)*length(cut_offs)))
  df <-rbind(df,df_poly)
  rm(df_poly)
}
dat_arrow <- data.frame(
  Type  = c("Offline","Online: 5 weeks", "Online: 10 weeks","Online: 15 weeks") ,
  x1    = rep(0.5,length(cluster_probs)),
  x2    = rep(0.05,length(cluster_probs)),
  y1    = 0.5*arrow_start^-1,
  y2    = (0.5*arrow_start^(-1))-0.45*arrow_start)

df %>% ggplot(aes(x=x,y=y,colour=Probs)) + 
    geom_line(size=1.25) + 
    scale_color_grey(breaks=c("50","55","60","65","70","75","80","85"), 
                   labels=c("50","55","60","65","70","75","80","85"), 
                   name=TeX("Threshold $\\tilde{\\omega}$ (in %)"),
                   start = 0.1, end = 0.9) +
  
  xlab(TeX("$\\gamma_T$")) + 
  ylab(TeX("$\\gamma_F$")) +
  # geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
  #              data = dat_arrow,
  #              color = "black",
  #              size = 1.1,
  #              arrow = arrow(length = unit(0.025, "npc"), type = "closed")) +
    facet_rep_wrap(~Type,nrow = 2,repeat.tick.labels = FALSE) + 
    scale_x_continuous(expand = c(0, 0),labels = function(x) sprintf("%.1f", x)) + 
    scale_y_continuous(expand = c(0, 0)) + 
    theme_bw() + 
    theme(legend.position="top",
          text = element_text(size=20),
          axis.ticks = element_line(size=1.25),
          axis.title = element_text(size=25),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.spacing = unit(1.8, "lines"),
          panel.border = element_rect(colour = "black"),
          plot.margin=unit(c(0,0.4,0,0.2),"inches"),
          legend.direction = "horizontal") + 
    guides(color = guide_legend(nrow = 1)) #+
  #annotate("text",x=0.15,y=0.6,label="Increase in\ncost effectiveness",size=4.5)


rm(list=ls())


###############################################################################################
#Old
###############################################################################################
df <- data.frame(x = numeric(), y= numeric(), Probs = character(), Type=character())
beta_opt <- data.frame(0,dim=c(length(cluster_probs),length(cut_offs)))
for(k in 1:length(cluster_probs)){
  clustering <-apply(apply(exp(cluster_probs[[k]]),c(2,3),median),1,which.max)
  probabilities <- apply(apply(exp(cluster_probs[[k]]),c(2,3),median),1,max)
  decision_boundary<-rep(0,length(cut_offs))
  for(i in 1:length(cut_offs)){
    patient_idx <- which(probabilities>=cut_offs[i])
    x<-sum(diag(table(factor(clustering[patient_idx],levels = 1:8), 
                      factor(mhmm_cluster_posterior[patient_idx],levels=1:8))))/length(patient_idx)
    decision_boundary[i]<-1/(1-x)-1
    beta_opt[k,i] <- (1-0.85*x)/(1-x)-1
  }
  df_poly <- data.frame(x=rep(c(0,0,1),length(cut_offs)),
                        y=rep(c(1,0,0),length(cut_offs)), 
                        Probs=factor(as.character(100*rep((cut_offs),each=3))),
                        Type = rep(switch(k,"Offline","Online: 5 weeks", "Online: 10 weeks","Online: 15 weeks"),3*length(cut_offs)))
  df_poly$x[seq(from=3*length(cut_offs),to=3,by=-3)] <- decision_boundary
  df <-rbind(df,df_poly)
  rm(df_poly)
}

df %>% ggplot(aes(x=x,y=y,fill=Probs)) + 
  geom_polygon() + 
  facet_wrap(~Type,nrow = 1) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme_bw() + 
  scale_fill_grey(breaks=c("85","80","75","70","65","60","55","50"), 
                  labels=c("50","55","60","65","70","75","80","85"), 
                  name=TeX("Threshold $\\tilde{\\omega}$ (in %)"),
                  start = 0.9, end = 0.1)+
  xlab(TeX("$\\gamma_F$")) + 
  ylab(TeX("$\\gamma_T$")) + 
  theme(legend.position="top",
        text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.direction = "horizontal") + 
  guides(fill = guide_legend(nrow = 1))
rm(list=ls())
########################################################################################################################
#Convergence of cluster membership for different probabilities over time
########################################################################################################################
C<- 8
S<- 3
load(paste0("ModelFits/mhmm_",S,"_",C,".R"))
samples <- rstan::extract(mhmm_model,c("posterior_cluster_test","prior_cluster_test","log_lik_test"))
for(c in 1:C){
  samples$posterior_cluster_test[,,c] <- samples$posterior_cluster_test[,,c] - samples$log_lik_test
}
rm(c)
#Posterior cluster after observing the complete trajectory == TRUE Clustering
mhmm_cluster_posterior<-apply(apply(exp(samples$posterior_cluster_test),c(2,3),median),1,which.max)

cluster_probs<-list()
#Prior clustering without any observationy Y
rm(samples,mhmm_model)
#Online clustering after1 through 52 weeks
load(file="ModelFits/time_cluster_prob.R")
for(k in 1:51){
  cluster_probs[[length(cluster_probs)+1]] <- time_dependent_cluster[,k,,]
}
rm(time_dependent_cluster)
cut_offs<-seq(from = 0.5, to = 0.85, by = 0.05)
accuracy<-array(0,dim=c(length(cut_offs),51))
overall_accuracy <- rep(0,51)
for(k in 1:length(cluster_probs)){
  clustering <-apply(apply(exp(cluster_probs[[k]]),c(2,3),median),1,which.max)
  overall_accuracy[k] <- sum(diag(table(factor(clustering,levels = 1:8), 
                                        factor(mhmm_cluster_posterior,levels = 1:8))))/422
  probabilities <- apply(apply(exp(cluster_probs[[k]]),c(2,3),median),1,max)
  for(i in 1:length(cut_offs)){
    patient_idx <- which(probabilities>=cut_offs[i])
    accuracy[i,k]<-sum(diag(table(factor(clustering[patient_idx],levels = 1:8), 
                      factor(mhmm_cluster_posterior[patient_idx],levels = 1:8))))/length(patient_idx)
  }
}
df <- data.frame(Probs = factor(as.character(100*rep((cut_offs),each=51))),
                 Time  = rep(1:51,8),
                 Accuracy = 100*as.vector(t(accuracy)))
df %>% ggplot(aes(x=Time,y=Accuracy,color=Probs)) + 
  geom_line(size=1.5) +
  scale_color_grey(name=TeX("Threshold $\\tilde{\\omega}$ (in %)")) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  ylab("Accuracy of subgroup assignment (in %)") + 
  xlab("Time (in Weeks)")+
  theme_bw() + 
  theme(legend.position="top",
        text = element_text(size=20),
        legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 1))
rm(list=ls())
########################################################################################################################
#Computation of clustering for different t (via t_vec)
#DO NOT RUN unless new run needed (takes very long)
########################################################################################################################
load("Preprocessed.MHMM.Data.Rdata")
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
rm(lbp.data_filtered,PreProcess,QR,risk_factors,Y_dur,Y_limit,Y_pain,N,T,train_idx)
N_samples <- 500
t_vec <- 1:51
time_dependent_cluster <- array(0,dim=c(N_samples,length(t_vec),422,8))
C<- 8
S<- 3
load(paste0("ModelFits/mhmm_",S,"_",C,".R"))
samples <- rstan::extract(mhmm_model,c("prior_cluster_test"))
samples_parameters <- rstan::extract(mhmm_model,c("pi","tpm","lambda_dis","inflation_raw","lambda_pain","theta_raw"))
expose_stan_functions(mhmm_model)
l_cum_pain <- array(0,dim=c(C,S,11))
l_cum_dis <- array(0,dim=c(C,S,8))
for(i in 1:N_samples){
  print(i)
  for(k in 1:C){
    for(s in 1:S){
      for(p in 0:10){
        l_cum_pain[k,s,p+1] <- ifelse(p==10 , 
                                  0.0 , 
                                  ppois(p, 
                                        lambda = samples_parameters$lambda_pain[i,k,s] , 
                                        log.p = TRUE) - ppois(10, 
                                                            lambda = samples_parameters$lambda_pain[i,k,s], 
                                                            log.p = TRUE))  
      }
      for(d in 0:7){
        l_cum_dis[k,s,d+1] <- ifelse(d == 7, 
                                     0.0, log1p(-inv_logit(samples_parameters$inflation_raw[i,k,s])) + 
                                     ppois(d, 
                                           lambda = samples_parameters$lambda_dis[i,k,s] , 
                                           log.p = TRUE) - ppois(6, 
                                                               lambda = samples_parameters$lambda_dis[i,k,s], 
                                                               log.p = TRUE))
      }
    }
  }
  for(t in 1:length(t_vec)){
    for(n in 1:422){
      for(k in 1:C){
        time_dependent_cluster[i,t,n,k] <- samples$prior_cluster_test[i,n,k] +  compute_log_lik(Y_pain = data_stan$Y_pain_test[n,1:t_vec[t]],
                                                                                                Y_dis = data_stan$Y_dis_test[n,1:t_vec[t]],
                                                                                                P = 11,
                                                                                                L = 8,
                                                                                                S = 3,
                                                                                                T = t_vec[t],
                                                                                                l_pi = log(samples_parameters$pi[i,k,]),
                                                                                                l_tpm = list(log(samples_parameters$tpm[i,k,1,]),
                                                                                                             log(samples_parameters$tpm[i,k,2,]),
                                                                                                             log(samples_parameters$tpm[i,k,3,])),
                                                                                                theta = 1+samples_parameters$theta_raw[i,k],
                                                                                                l_prob_1 = list(l_cum_pain[k,1,],
                                                                                                                l_cum_pain[k,2,],
                                                                                                                l_cum_pain[k,3,]),
                                                                                                l_prob_2 = list(l_cum_dis[k,1,],
                                                                                                                l_cum_dis[k,2,],
                                                                                                                l_cum_dis[k,3,]))
      }
      time_dependent_cluster[i,t,n,] <- time_dependent_cluster[i,t,n,]-log_sum_exp(time_dependent_cluster[i,t,n,])
    }
  }
}

#################################################################################################################################
#Cluster interpretation - Comparison of dynamics
#################################################################################################################################
##Closer inspection of different dynamics in each cluster
load("Preprocessed.MHMM.Data.Rdata")
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
rm(lbp.data_filtered,PreProcess,QR,risk_factors,Y_dur,Y_limit,Y_pain,N,T,train_idx)
C<- 8
S<- 3
load(paste0("ModelFits/mhmm_",S,"_",C,".R"))
samples <- rstan::extract(mhmm_model,c("posterior_cluster_test","log_lik_test"))
for(c in 1:C){
  samples$posterior_cluster_test[,,c] <- samples$posterior_cluster_test[,,c] - samples$log_lik_test
}
map_idx_test <- which.max(apply(samples$log_lik_test,1,sum))
mhmm_cluster_test<-apply(samples$posterior_cluster_test[map_idx_test,,],1,which.max)
rm(samples,c)
samples_parameters <- rstan::extract(mhmm_model,c("pi","tpm","lambda_dis","inflation_raw","lambda_pain","theta_raw"))
expose_stan_functions(mhmm_model)
rm(mhmm_model)
l_cum_pain <- array(0,dim=c(S,11))
l_cum_dis <- array(0,dim=c(S,8))


hidden_seq <- array(0,dim=c(422,52))
for(n in 1:422){
  c <- mhmm_cluster_test[n]
  for(s in 1:S){
    for(p in 0:10){
      l_cum_pain[s,p+1] <- ifelse(p==10 , 
                                    0.0 , 
                                    ppois(p, 
                                          lambda = samples_parameters$lambda_pain[map_idx_test,c,s] , 
                                          log.p = TRUE) - ppois(10, 
                                                                lambda = samples_parameters$lambda_pain[map_idx_test,c,s], 
                                                                log.p = TRUE))  
    }
    for(d in 0:7){
      l_cum_dis[s,d+1] <- ifelse(d == 7, 
                                   0.0, log1p(-inv_logit(samples_parameters$inflation_raw[map_idx_test,c,s])) + 
                                     ppois(d, 
                                           lambda = samples_parameters$lambda_dis[map_idx_test,c,s] , 
                                           log.p = TRUE) - ppois(6, 
                                                                 lambda = samples_parameters$lambda_dis[map_idx_test,c,s], 
                                                                 log.p = TRUE))
    }
  }
  hidden_seq[n,] <- viterbi(Y_pain = data_stan$Y_pain_test[n,],
                            Y_dis = data_stan$Y_dis_test[n,],
                            S = 3,
                            T = 52,
                            P = 11,
                            L = 8,
                            l_pi = log(samples_parameters$pi[map_idx_test,c,]),
                            l_tpm = list(log(samples_parameters$tpm[map_idx_test,c,1,]),
                                         log(samples_parameters$tpm[map_idx_test,c,2,]),
                                         log(samples_parameters$tpm[map_idx_test,c,3,])),
                            theta = 1 + samples_parameters$theta_raw[map_idx_test,c],
                            l_prob_1 = list(l_cum_pain[1,],
                                            l_cum_pain[2,],
                                            l_cum_pain[3,]),
                            l_prob_2 = list(l_cum_dis[1,],
                                            l_cum_dis[2,],
                                            l_cum_dis[3,]))
}

rm(s,c,map_idx_test,samples_parameters,l_cum_pain,l_cum_dis,data_stan,d,p,n)

#Recoding of the hidden states: Inspection of posterior mean of emission parameters reveals
# That all cluster but 5,2 and 1 follow the normal ordering, i.e., 1=Stable, 2=Unstable, 3=Acute
# Cluster 5,2, and 1 are recoded accordingly

hidden_seq[mhmm_cluster_test==5,] <- ifelse(hidden_seq[mhmm_cluster_test==5,] == 1 , 
                                         2 , 
                                         ifelse(hidden_seq[mhmm_cluster_test==5,]==2, 
                                                3, 
                                                1))
hidden_seq[mhmm_cluster_test==2,] <- ifelse(hidden_seq[mhmm_cluster_test==2,] == 1 , 
                                         2 , 
                                         ifelse(hidden_seq[mhmm_cluster_test==2,]==2, 
                                                1, 
                                                3))
hidden_seq[mhmm_cluster_test==1,] <- ifelse(hidden_seq[mhmm_cluster_test==1,] == 1 , 
                                         3 , 
                                         ifelse(hidden_seq[mhmm_cluster_test==1,]==2, 
                                                2, 
                                                1))



shares<-array(0,dim=c(52,3))
df <- data.frame(Cluster = character(),
                 Time = numeric(),
                 State = character(),
                 Share = numeric())
for(c in 1:C){
  for(t in 1:52){
    shares[t,] <- prop.table(table(factor(hidden_seq[mhmm_cluster_test==c,t],levels = 1:3)))
  }
  df <- rbind(df,data.frame(Cluster = rep(paste("Subgroup",c),156),
                            Time = rep(1:52,3),
                            State = rep(paste("State",3:1),each=52),
                            Share = as.vector(shares)))
}
rm(shares)

df  %>% ggplot(aes(x = Time, y = 100*Share , fill=factor(State))) + 
  geom_area(alpha = 0.6) + 
  facet_wrap(~Cluster,nrow=2) + 
  xlab("Time") + 
  ylab("Proportion of patients") +
  scale_x_continuous(expand = c(0, 0),limits = c(1,52)) + 
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + 
  scale_fill_grey(breaks=c("State 1","State 2","State 3"), 
                  labels=c("Acute","Unstable","Stable"), 
                  name="Patient health state: ") +
  theme(legend.position="top",
        text = element_text(size=20), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
  


rm(list=ls())

######################################################################################################################
#Cluster interpretation - Cluster descriptives
######################################################################################################################
load(file="Preprocessed.MHMM.Data.Rdata")

C <- 8
S <- 3
load(paste0("ModelFits/mhmm_",S,"_",C,".R"))
samples <- rstan::extract(mhmm_model,c("posterior_cluster_test","log_lik_test"))
for(c in 1:C){
  samples$posterior_cluster_test[,,c] <- samples$posterior_cluster_test[,,c] - samples$log_lik_test
}
map_idx_test <- which.max(apply(samples$log_lik_test,1,sum))
mhmm_cluster_test<-apply(samples$posterior_cluster_test[map_idx_test,,],1,which.max)
rm(samples,mhmm_model,c,map_idx_test,PreProcess,Y_dur)
raw_data_test <- lbp.data_filtered[-train_idx,]
rm(lbp.data_filtered,risk_factors)


test_pain <- Y_pain[-train_idx,]
test_limit <- Y_limit[-train_idx,]
test_pain[test_pain==11]<-NA
test_limit[test_limit==8]<-NA
Activity <- ifelse(as.numeric(raw_data_test$bfbe0) %in%c(3,4),1,ifelse(as.numeric(raw_data_test$bfbe0) %in%c(1,2),0,NA))
description_df <- data.frame(Cluster = paste("Cluster",1:8) ,  
                             #Description = c("Mild fluctuating back pain (short-term but fluctuating recovery, small improvements in pain, yet no clear tendency regarding patient outcomes)",
                             #                "Severe back pain (highly unstable phases involving recurrent and extent pain, as well as disability, with mostly females)", 
                             #                "Recovery (with mostly stable phases and hardly any activity limitations)", "Relapsing back pain (largely for females with marginal improvements, though low pain levels)",
                             #                "Early improving back pain (initial recovery but fluctuating condition with continuous unstable phases for young patients with low BMI and rare physical work)",
                             #                "Relapsing back pain (low pain and disability levels, yet with frequent unstable phases for largely females)",
                             #                "Persistent severe back pain (affecting mostly males with high BMI)",
                             #                "Mild episodic (prevalent recovery from a largely work-induced condition for males, resulting in very stable phases, low pain and modest disability level)") ,
                             Patients = sapply(1:C,function(k){sum(mhmm_cluster_test==k)}),
                             Male = 100*sapply(1:C,function(k){mean(as.numeric(raw_data_test$bsex0[mhmm_cluster_test==k])-1)}) , 
                             Age = sapply(1:C,function(k){mean(raw_data_test$age[mhmm_cluster_test==k])}) , 
                             BMI = sapply(1:C,function(k){mean(raw_data_test$bmi[mhmm_cluster_test==k],na.rm = TRUE)}) , 
                             Physical_Work = 100*sapply(1:C,function(k){mean(Activity[mhmm_cluster_test==k],na.rm = TRUE)}), 
                             Average_pain = sapply(1:C , 
                                                   function(k){mean(test_pain[mhmm_cluster_test==k],na.rm = TRUE)}) , 
                             SD_pain = sapply(1:C , 
                                              function(k){sd(test_pain[mhmm_cluster_test==k],na.rm = TRUE)}) , 
                             Average_disability = sapply(1:C , 
                                                         function(k){mean(test_limit[mhmm_cluster_test==k , ],na.rm = TRUE)}) , 
                             SD_disability = sapply(1:C , 
                                                    function(k){sd(test_limit[mhmm_cluster_test==k , ],na.rm = TRUE)}) ,
                             Pain_12m = sapply(1:C , 
                                               function(k){mean(raw_data_test[mhmm_cluster_test==k,"vasl12m"] , 
                                                                na.rm = TRUE)}) , 
                             Disability_12m = sapply(1:C , 
                                                     function(k){mean(raw_data_test[mhmm_cluster_test==k,"rmprop12m"] , 
                                                                      na.rm = TRUE)}))




print.xtable(xtable(description_df,digits = 2) , 
             sanitize.text.function = identity , 
             include.rownames = FALSE , 
             include.colnames = FALSE ,
             type = "latex" , 
             hline.after = NULL ,
             booktabs = FALSE , 
             floating = FALSE ,
             only.contents = TRUE,
             file="cluster_interpretation.tex")

rm(list=ls())



#OLD

#As an initial analysis, we compare the different mean transition matrices from the clusters, namely cluster 2 and cluster 7



tpm_samples <- rstan::extract(basic_model,"tpm")$tpm
tpm_5 <- tpm_3 <- array(0,dim=c(6,3))
for(i in 1:3){
  for(j in 1:3){
    tpm_5[(2*i-1),j] <- round(100*mean(tpm_samples[,5,i,j]),2)
    tpm_5[2*i,j]<-paste0("(" , 
                         paste(round(100*quantile(tpm_samples[,5,i,j],c(0.025,0.975)),2) , 
                               collapse = ","),")")
    
    tpm_3[(2*i-1),j] <- round(100*mean(tpm_samples[,3,i,j]),2)
    tpm_3[2*i,j]<-paste0("(" , 
                         paste(round(100*quantile(tpm_samples[,3,i,j],c(0.025,0.975)),2) , 
                               collapse = ","),")")
  }
}
rm(i,j,tpm_samples)
print.xtable(xtable(cbind(c("Stable","" , "Unstable" ,"", "Acute",""),tpm_3,tpm_5)) , 
             sanitize.text.function = identity , 
             include.rownames = FALSE , 
             include.colnames = FALSE ,
             type = "latex" , 
             hline.after = NULL ,
             booktabs = FALSE , 
             floating = FALSE ,
             only.contents = TRUE,
             file="tpm_comparison.tex")
rm(tpm_3,tpm_5)







