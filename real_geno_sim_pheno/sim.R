setwd("/home/amelie/nas/PKU_Inter99/merged_kasper/data/snp_simulation/real_geno_sim_pheno/")
rm(list=ls())
library(data.table)
library(broom)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(parallel)

geno <- BEDMatrix::BEDMatrix("/home/amelie/nas/PKU_Inter99/merged_kasper/data/clumped_data/clumped_data/PKU_Inter99_clumped.bed", simple_names=TRUE) %>% 
  as.matrix()
fam <- fread("/home/amelie/nas/PKU_Inter99/merged_kasper/data/clumped_data/clumped_data/PKU_Inter99_clumped.fam")

identical(rownames(geno), fam$V2)

N_ind <- nrow(geno)
N_SNPs <- ncol(geno)
N_int <- 10
N_m_ef <- 100

# sample beta g, sex and interaction with no effect
b_g <- rnorm(N_SNPs, 0, 0.01)
b_sex <- rnorm(1, 0, 0.01)
b_int <- rnorm(N_SNPs, 0, 0.01)

# sample N_m_eff positions for SNPs with main effects from all SNPs and exchange their beta (effect size)
pos_main <- sample(N_SNPs, N_m_ef, replace = FALSE)
m_dir <- rbinom(N_m_ef,1 , 0.5) 
b_g[pos_main] <- ifelse(m_dir==0, rnorm(N_m_ef, -0.5, 0.01), rnorm(N_m_ef, 0.5, 0.01))

# sample N_int positions for SNPs with effects from all SNPs and exchange their beta (effect size)
pos_int <- sample(N_SNPs, N_int, replace = FALSE)
m_dir <- rbinom(N_int, 1, 0.5)
b_int[pos_int] <- ifelse(m_dir==0, rnorm(N_int, -0.5, 0.01), rnorm(N_int, 0.5, 0.01))

# intercept
b_0 <- 0

sex <- fam$V5

# replace NA with a sampled genotype
N_NA <- geno[is.na(geno)] %>% 
  length()
geno[is.na(geno)] <- rbinom(N_NA, 2, runif(N_NA, 0.10, 0.49)) 

# calculate interactions
int <- sex*geno
identical(int[,1], sex*geno[,1])

#Combine into one genotype matrix
X <- cbind(geno, sex, int)
betas <- c(b_g, b_sex, b_int)

logit <- b_0 + X%*%betas

prob <- exp(logit)/(exp(logit)+1)

pheno <- rbinom(N_ind, 1, prob)

d_test <- data.frame(pheno = pheno, geno, sex = sex)
write.table(d_test, file = "d_test.txt")

res <- list()
ii <- 1
#for(ii in colnames(d_test[,2:(ncol(d_test)-1)])) {
for(ii in 1:N_SNPs) {
  d_cur <- d_test %>% select(pheno, (ii+1), sex)
  snp <- colnames(d_test %>% select(ii+1))
  colnames(d_cur)[2] <- c("SNP")
  
  m <- glm(pheno ~ SNP + sex + SNP*sex, data = d_cur, family = binomial()) %>% tidy() %>% 
    dplyr::filter(term == "SNP:sex") %>% 
    mutate(index = ii)
  m$snp <- snp
  res[[ii]] <- m
}


results <- bind_rows(res) %>% 
  select(snp, beta_freq=estimate, sd_freq=std.error, p.value) 

#snp <- 3
model <- stan_model(file = "/nas/users/amelie/PKU_Inter99/merged_kasper/model/model_anders_nocovar.stan")  
system.time(
  r <- mclapply(1:N_SNPs, function(snp) {
  
    X <- d_test %>% select(pheno, "sex", (snp+1))
    
    X <- cbind(X, 0)
    snp_name <- colnames(X)[3]
    colnames(X)[c(3,4)] <- c("g_SNP", "int") 
    
    X[,"g_SNP"] <- X[,"g_SNP"] - mean(X[,"g_SNP"])
    X[,"sex"] <- X[,"sex"] - mean(X[,"sex"]) #sex
    X[,"int"] <- X[,"g_SNP"] * X[,"sex"] #interaction
    
    
    data <- list(N = nrow(X), # number of samples (phenotypes)
                 pheno = X[,"pheno"], # phenotype
                 M = 3, # sex, snp, interaction
                 #P = 10, # number of covariates
                 X = X[,2:4]) # predictor matrix
    
    fit <- sampling(
      model,  # Stan program
      data = data,    # named list of data
      chains = 4,             # number of Markov chains
      warmup = 100,          # number of warmup iterations per chain
      iter = 1000,            # total number of iterations per chain
      cores = 4,              # number of cores (could use one per chain)
      refresh = 0             # no progress shown
    )
    
    fit_summary <- summary(fit)$summary[4,]
    
    beta_int <- round(as.numeric(fit_summary["mean"]), 3)
    sd_int <- round(as.numeric(fit_summary["sd"]), 3)
    
    print(c(snp_name, beta_int, sd_int))
  }, mc.cores=25)
)

res <- as.data.frame(matrix(unlist(r), ncol=3, byrow=TRUE))
colnames(res) <- c("snp", "beta_int_mb", "sd_int_mb")
write.table(res, "res_bayes.txt", quote = FALSE, row.names = FALSE)

r_bayes <- fread("res_bayes.txt")

results2 <- results %>% 
  left_join(r_bayes, by=c("snp"="snp")) %>% 
  cbind(b_int)

write.table(results2, "res_all.txt", row.names = F, quote = F)


