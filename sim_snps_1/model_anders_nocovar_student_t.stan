data {
int<lower=1> N; // number of samples (phenotypes)
int<lower=1> M; // 3; sex, SNP, interaction
//int<lower=1> P; // number of covariates
matrix[N, M] X; // predictor matrix
int<lower=0, upper=1> pheno[N]; // phenotype
}

parameters {
real alpha; 
vector[M] beta;
}

model {
alpha  ~ normal(0,1);
beta[1:M] ~ student_t(1, 0, 0.1);
//double_exponential(0, 0.1); // SNP
//beta[M+1:P] ~ normal(0,1);

pheno ~ bernoulli_logit_glm(X, alpha, beta);
}
