rm(list=ls())
path_WD <- "/home/amelie/nas/PKU_Inter99/merged_kasper/data/snp_simulation/"

setwd(path_WD)
library(data.table)
library(broom)
library(tidyverse)

N_ind <- 60 # number of individuals
N_SNPs <- 1000 # number of SNPs to simulate
N_int <- 10 # number of SNPs that have an interaction with sex
N_m_ef <- 100 # number of SNPs that have a main effect

MAF <- runif(N_SNPs*N_ind, 0.1, 0.49)

geno <- matrix(rbinom(N_SNPs*N_ind, 2, MAF), N_ind, N_SNPs)

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

sex <- rbinom(N_ind, 1, 0.5)

# calculate interactions
int <- sex*geno
identical(int[,1], sex*geno[,1])

#Combine into one genotype matrix
X <- cbind(geno, sex, int)
betas <- c(b_g, b_sex, b_int)

logit <- b_0 + X%*%betas

prob <- exp(logit)/(exp(logit)+1)

pheno <- rbinom(N_ind, 1, prob)
betas_df <- data.frame(b_g = b_g, b_sex = b_sex, b_int = b_int, b_0 = b_0)
d_test <- data.frame(pheno = pheno, geno, sex = sex)

write.table(betas_df, "betas_df.txt", quote = F, row.names = F)
write.table(d_test, file = "d_test.txt", quote = F)