setwd("/home/amelie/nas/PKU_Inter99/merged_kasper/data/sim_snps_1/")
rm(list=ls())
library(data.table)
library(broom)
library(tidyverse)

set.seed(123)
N_SNPs <- 10
N_ind <- 5
N_int <- 5

b_main <- matrix(rnorm(N_SNPs*N_ind, 0, 0.01), N_ind, N_SNPs)
b_sex <- rnorm(1, 0, 0.01)
b_int <- matrix(rnorm(N_SNPs, 0, 0.01), N_ind, N_SNPs)

pos_int <- sample(N_SNPs, N_int)
b_int[pos_int] <- rnorm(N_int, 0.5, 0.01)

b_0 <- 0

MAF <- runif(N_SNPs*N_ind, 0.1, 0.49)

g <- matrix(rbinom(N_SNPs*N_ind, 2, MAF), N_ind, N_SNPs)

sex <- rbinom(N_ind, 1, 0.5)

logit <- b_0 + b_main*g + b_sex*sex + b_int*g*sex
prob <- exp(logit)/(exp(logit)+1)
pheno <- matrix(rbinom(N_ind*N_SNPs, 1, prob), N_ind, N_SNPs)

i <- 1
ii <- 1
for (i in 1:N_ind) {
  for (ii in 1:N_SNPs) {
    data_cur <- list(pheno_cur = pheno[i, ii],
                 g_cur = g[i, ii],
                 sex_cur = sex[i])
    pheno_cur <- pheno[i, ii]
    g_cur <- g[i, ii]
    sex_cur <- sex[i]
    glm(pheno_cur ~ g_cur + sex_cur + g_cur*sex_cur, family = binomial) %>% 
      tidy()
  }
}

y <- 1
x <- 1

glm(y ~ x, family=binomial) %>% 
  tidy()



set.seed(42)
df <- data.frame(target = sample(0:1, 10, TRUE),
                 F1 = rnorm(10), F2 = rnorm(10), F3 = rnorm(10))


predictors <- c("F1", "F2", "F3")    
fits <- sapply(predictors,
               function(x) {
                 tmp <- try(coef(summary(glm(as.formula(paste("target", x, sep = "~")),
                                             family=binomial, data = df)))[2, ], TRUE)
                 if (class(tmp) == "try-error") NULL else tmp})

glm(df$target[1] ~ df$F1[1] + df$F2[1] + df$F1[1]*df$F2[1]) %>% 
  tidy()

glm(df$target ~ df$F1 + df$F2 + df$F1*df$F2) %>% 
  tidy()
