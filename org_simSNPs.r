library(data.table)
library(tidyverse)
library(broom)
rm(list=ls())

N <- 2000
SNPs <- 10
N_interactions_with_effect <- 5

# sample genotype matrix
SNPmat <- matrix(data = 0, nrow = N, ncol = SNPs)

for(ii in 1:SNPs){
  # sample MAF between 0.1 & 0.49
  MAF <- runif(1, 0.10, 0.49)
  # sample genotype for current SNP for all ind with MAF (0,1,2)
  cur_SNP <- rbinom(N, 2, MAF)
  # enter genotype value to SNPmat
  SNPmat[,ii] <- cur_SNP
}
# SNPmat is the "normal" genotype matrix

#Simulate main effects for # of SNPs, sample from normal distribution
main_effects <- rnorm(SNPs, mean = 0, sd = 0.01)

# number combinations of all possible SNP combinations
combinations <- combn(SNPs, 2) 
# sample 5 (or N) SNPs with effect
effect_idx <- sample(dim(combinations)[2], N_interactions_with_effect, replace = FALSE)

# create "empty" interaction matrix
SNPmat_int <- matrix(data = 0, nrow = N, ncol = dim(combinations)[2])

# fill interaction matrix
for(ii in 1:dim(combinations)[2]) {
  
  idx_1 <- combinations[,ii][1]
  idx_2 <- combinations[,ii][2]
  
  SNP1 <- SNPmat[,idx_1]
  SNP2 <- SNPmat[,idx_2]
  
  SNPmat_int[,ii] <- SNP1*SNP2
}

beta_int <- rep(0, dim(combinations)[2])
# sets b to 0.5 if interaction and to 0 if no interaction
for(ii in 1:dim(combinations)[2]){
  
  if(ii %in% effect_idx){
    
    beta_int[ii] <- rnorm(1, mean = 0.5, sd = 0.01)
    
  } else {
    
    beta_int[ii] <- rnorm(1, mean = 0, sd = 0.01)
    
  }
  
}


#Combine into one genotype matrix

X <- cbind(SNPmat, SNPmat_int)
betas <- c(main_effects, beta_int)

intercept <- -1

logit <- intercept + X%*%betas

prob <- exp(logit)/(exp(logit)+1)

pheno <- rep(0, N)

for(ii in 1:length(pheno)){
  
  pheno[ii] <-  rbinom(1, 1, prob = prob[ii])
  
}


d_test <- data.frame(pheno = pheno, SNPs = SNPmat)
write.csv(d_test, file = "d_test.csv")

class(d_test)

res <- list()

for(ii in 1:dim(combinations)[2]) {
  
  idx_1 <- combinations[,ii][1] + 1
  idx_2 <- combinations[,ii][2] + 1
  
  d_cur <- d_test %>% select(pheno, idx_1, idx_2)
  colnames(d_cur)[c(2,3)] <- c("SNP1", "SNP2")
  
  m <- glm(pheno ~ SNP1*SNP2, data = d_cur, family = binomial) %>% tidy() %>% 
    filter(term == "SNP1:SNP2") %>% 
    mutate(index = ii)
  
  res[[ii]] <- m
}


results <- bind_rows(res) %>% 
  mutate(real_int = beta_int)

