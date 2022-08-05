setwd("/home/amelie/nas/PKU_Inter99/merged_kasper/data/snp_simulation/sim_snps_1/")
library(data.table)
library(broom)
library(tidyverse)

rm(list=ls())
N_SNPs <- 10000
list <- NULL
for (b_0 in c(-5, 0, 5)) {
  for (g in c(0, 1, 2)) {
    print(c(b_0, g))
    
    case_vector <- NULL
    for (i in 1:500) {
      cases <- 0
      for (i in 1:N_SNPs) {
        #g <- 2 #(0/1/2) 0 effect alleles, 1 effect allele, 2 effect alleles
        sex <- 0 #(0/1) 0 male, 1 female
        
        #b_0 <- 5 # determines case/control ratio
        b_main <- rnorm(1, mean = 0, sd = 0.01) #no effect
        b_sex <- rnorm(1, mean = 0, sd = 0.01) #no effect
        b_int <- rnorm(1, mean = 0.5, sd = 0.01) # effect
        
        logit <- b_0 + b_main*g + b_sex*sex + b_int*g*sex
        prob <- exp(logit)/(exp(logit)+1)
        pheno <- rbinom(1, 1, prob)
        cases <- cases + pheno
      }
      case_vector <- append(case_vector, cases)
    }
    case_vector
    res <- (c("gender: ", sex, "genotype: ", g, "intercept: ", b_0, "average number of cases: ", mean(case_vector)))
    list <- append(list, res)
    
    case_vector_female <- NULL
    for (i in 1:500) {
      cases <- 0
      for (i in 1:N_SNPs) {
        #g <- 2 #(0/1/2) 0 effect alleles, 1 effect allele, 2 effect alleles
        sex <- 1 #(0/1) 0 male, 1 female
        
        #b_0 <- 5 # determines case/control ratio
        b_main <- rnorm(1, mean = 0, sd = 0.01) #no effect
        b_sex <- rnorm(1, mean = 0, sd = 0.01) #no effect
        b_int <- rnorm(1, mean = 0.5, sd = 0.01) # effect
        
        logit <- b_0 + b_main*g + b_sex*sex + b_int*g*sex
        prob <- exp(logit)/(exp(logit)+1)
        pheno <- rbinom(1, 1, prob)
        cases <- cases + pheno
      }
      case_vector_female <- append(case_vector_female, cases)
    }
    case_vector_female
    res <- (c("gender: ", sex, "genotype: ", g, "intercept: ", b_0, "average number of cases: ", mean(case_vector_female)))
    list <- append(list, res)
    
  }
}

res2 <- as.data.frame(matrix(unlist(list), nrow=18, byrow=TRUE))
write.table(res2, "intercept_testing_int_effect.txt", quote = F, row.names = F)
