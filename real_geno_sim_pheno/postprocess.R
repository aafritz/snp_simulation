setwd("/home/amelie/nas/PKU_Inter99/merged_kasper/data/snp_simulation/real_geno_sim_pheno/")
rm(list=ls())

res_all <- fread("res_all.txt")


res_all2 <- res_all %>% 
  add_column(OR_freq=exp(res_all$beta_freq), .after="sd_freq") %>% 
  add_column(OR_bayes=exp(res_all$beta_int_mb), .after="sd_int_mb") %>% 
  add_column(prob=pnorm(0.1, abs(res_all$beta_int_mb), res_all$sd_int_mb, lower.tail=FALSE), .before="b_int") %>% 
  add_column(OR_true=exp(res_all$b_int), .after="b_int")

res_all3 <- res_all2 %>% 
  add_column(int_freq=ifelse(res_all2$OR_freq>1.1 | res_all2$OR_freq<0.9, 1, 0), .after="OR_true") %>% 
  add_column(int_bayes=ifelse(res_all2$OR_bayes>1.1 | res_all2$OR_bayes<0.9, 1, 0), .after="OR_true") %>% 
  add_column(int_true=ifelse(res_all2$OR_true>1.1 | res_all2$OR_true<0.9, 1, 0), .after="OR_true") %>% 
  order(decreasing=T)

sum(res_all3$int_true)
sum(res_all3$int_bayes)
sum(res_all3$int_freq)

res_all3 %>% 
  filter(prob>0.95)

res_all3 %>% 
  filter(int_true==1)
