setwd("/home/amelie/nas/PKU_Inter99/merged_kasper/data/snp_simulation/sim_snps_1/")
rm(list=ls())
library(data.table)
library(broom)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(parallel)

set.seed(123)
N_SNPs <- 10
N_ind <- 10000
N_int <- 5

# sample beta g, sex and interaction with no effect
b_g <- matrix(rnorm(N_SNPs, 0, 0.01), 1, N_SNPs)
b_g <- rnorm(N_SNPs, 0, 0.01)
b_sex <- rnorm(1, 0, 0.01)
b_int <- rnorm(N_SNPs, 0, 0.01)

# sample N_int positions for SNPs with effects from all SNPs and exchange their beta (effect size)
pos_int <- sample(N_SNPs, N_int)
b_int[pos_int] <- rnorm(N_int, 0.5, 0.01)

# intercept
b_0 <- 0

MAF <- runif(N_SNPs*N_ind, 0.1, 0.49)

g <- matrix(rbinom(N_SNPs*N_ind, 2, MAF), N_ind, N_SNPs)

sex <- rbinom(N_ind, 1, 0.5)

model <- stan_model(file = "model_anders_nocovar_student_t.stan")  
res_all <- NULL
i <- 1
r <- mclapply(1:N_SNPs, function(i) {

  int <- sex*g[,i]
  logit <- b_0 + b_sex*sex + b_g[i]*g[,i] + b_int[i]*int
  prob <- exp(logit)/(exp(logit)+1)
  pheno <- rbinom(N_ind, 1, prob)
  
  X <- cbind(sex, g[,i], int)
  
  data <- list(N = N_ind, # number of samples (phenotypes)
               pheno = pheno, # phenotype
               M = 3, # SNP, sex & interaction
               X = X) # predictor matrix
  
  ##--## sampling from stan model
  #system.time(
  fit <- sampling(
    model,  # Stan program
    data = data,    # named list of data
    chains = 4,             # number of Markov chains
    warmup = 100,          # number of warmup iterations per chain
    iter = 1000,            # total number of iterations per chain
    cores = 4,              # number of cores (could use one per chain)
    refresh = 0             # no progress shown
  )
  fit
  fit_summary <- summary(fit)
  
  b_sex <- round(as.numeric(fit_summary$summary[2,"mean"]), 2)
  sd_sex <- round(as.numeric(fit_summary$summary[2,"sd"]), 2)
  n_eff_sex <- round(as.numeric(fit_summary$summary[2,"n_eff"]), 2)
  rhat_sex <- round(as.numeric(fit_summary$summary[2,"Rhat"]), 2)
  
  b_SNP <- round(as.numeric(fit_summary$summary[3,"mean"]), 2)
  sd_SNP<- round(as.numeric(fit_summary$summary[3,"sd"]), 2)
  n_eff_SNP <- round(as.numeric(fit_summary$summary[3,"n_eff"]), 2)
  rhat_SNP <- round(as.numeric(fit_summary$summary[3,"Rhat"]), 2)
  
  b_int <- round(as.numeric(fit_summary$summary[4,"mean"]), 2)
  sd_int <- round(as.numeric(fit_summary$summary[4,"sd"]), 2)
  n_eff_int <- round(as.numeric(fit_summary$summary[4,"n_eff"]), 2)
  rhat_int <- round(as.numeric(fit_summary$summary[4,"Rhat"]), 2)
  
  print(c(i, b_sex, sd_sex, n_eff_sex, rhat_sex, b_SNP, sd_SNP, n_eff_SNP, rhat_SNP, 
          b_int, sd_int, n_eff_int, rhat_int))
  
}, mc.cores=10)

result_df <- data.frame(matrix(unlist(r), ncol=13, byrow=TRUE))
colnames(result_df) <- c("snp", "b_sex", "sd_sex", "n_eff_sex", "rhat_sex", 
                         "b_SNP", "sd_SNP", "n_eff_SNP", "rhat_SNP", 
                         "b_int", "sd_int", "n_eff_int", "rhat_int")

data.frame(snp=result_df$snp, b_int_real = b_int, b_int_inf = result_df$b_int) %>% 
  write.table("results_int_student_t.txt", row.names = F, quote = F)

write.table(result_df, file = "results_student_t.txt", row.names = F, quote = F)


