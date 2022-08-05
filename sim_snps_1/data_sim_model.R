setwd("/home/amelie/nas/PKU_Inter99/merged_kasper/data/snp_simulation/sim_snps_1/")
rm(list=ls())
library(data.table)
library(broom)
library(tidyverse)

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

res_all <- NULL
for (i in 1:N_SNPs) {
  logit <- b_0 + b_sex*sex + b_g[i]*g[,i] + b_int[i]*(sex*g[,i])
  prob <- exp(logit)/(exp(logit)+1)
  pheno <- rbinom(N_ind, 1, prob)
  
  fit <- glm(pheno ~ g[,i] + sex + g[,i]*sex, family = binomial) %>% 
    tidy() 
  fit_t <- fit[,c("term","estimate")]
  res <- data.frame(true = c(b_0, b_g[i], b_sex, b_int[i]))
  #res <- data.frame(true = c(b_0, b_g[i]))
  fit_t <- cbind(fit_t, res, i)
  fit_t
  res_all <- rbind(res_all, fit_t)
}
filter(res_all, term=="g[, i]:sex")



