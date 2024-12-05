
##### FUNCTIONS
sim_gene <- function(n,alpha,sigma,Li, sigma_multiplier=1){

  # Simulate a single gene with following parameters:
  # n:                number of cells
  # alpha:            alpha parameter
  # sigma:            sigma parameter
  # Li:               library size
  # sigma_multiplier: factor to increase variance

  # Returns expression values for single gene in n cells

  # initialize Y
  Y <- rep(NA,n)

  # calculate latent expression using gamma distribution
  # ensure that simulation has mean within 10% of target
  # perc_diff = 1
  # while(perc_diff > 0.05){
  lambdaj <- rgamma(n,shape = alpha/sigma_multiplier,scale = sigma*sigma_multiplier)
  # perc_diff = (alpha*sigma - mean(lambdaj))/alpha*sigma
  # }


  #lambdaj <- rgamma(n,shape = alpha/sigma_multiplier,scale = sigma*sigma_multiplier)

  # for each cell use the latent expression and library size to simulate expression
  for(i in 1:n){

    Y[i] <- rpois(1,Li*lambdaj[i])
  }

  return(Y)
}

sigma_sq_est <- function(df,n_cells,q,N_c){
  sigma_sq_term1 = (df^2 - df*(1-q))/(N_c^2*q^2 - N_c*q*(1-q))
  sigma_sq_term1 = colSums(sigma_sq_term1)/n_cells

  sigma_sq_term2 = df/(N_c*q)
  sigma_sq_term2 = (colSums(sigma_sq_term2)/n_cells)^2

  sigma_sq = sigma_sq_term1 - sigma_sq_term2


  return(sigma_sq)
}

##### PARAMETERS
avg_li = 4064.29
mu <- 1
phi <- 0.8917347
sigma <- mu*phi/avg_li
alpha <- 1/phi
n_sims = 1000
B = 1000
n_cells = 1000
n_indv = 17
q = 0.4

library_size = 4064.29
N_c = library_size/q
varying_factor <- seq(1, 2, length.out = 100)

##### RUN SIMULATIONS
n_sig <- rep(0,(length(varying_factor)-1))
for(j in 1:n_sims){
  gene_sim_matrix = c()
  num_indv_grp = n_cells*n_indv

  for(i in 1:length(varying_factor)){
    gene_sim <- sim_gene(num_indv_grp,alpha,sigma,avg_li,sigma_multiplier=varying_factor[i])
    gene_sim_matrix <- cbind(gene_sim_matrix,gene_sim)

  }
  colnames(gene_sim_matrix) <- make.names(varying_factor)

  sigma_sq_obs = sigma_sq_est(gene_sim_matrix,num_indv_grp,q,N_c)
  delta_t = sigma_sq_obs[2:length(sigma_sq_obs)] - sigma_sq_obs[1]


  asl_geq <- rep(0, (length(varying_factor)-1))
  asl_lt <- rep(0, (length(varying_factor)-1))
  for(i in 1:B){
    sample_idx <- sample(1:num_indv_grp,size = num_indv_grp,replace=TRUE)
    sample_df <- gene_sim_matrix[sample_idx,]

    sigma_sq_samp = sigma_sq_est(sample_df,num_indv_grp,q,N_c)
    delta_t_boot = sigma_sq_samp[2:length(sigma_sq_samp)] - sigma_sq_samp[1]
    delta_t_star = delta_t_boot - delta_t

    asl_geq = asl_geq + as.numeric(delta_t_star >= delta_t)
    asl_lt = asl_lt + as.numeric(delta_t_star <= delta_t)

  }


  asl <-c()
  for(i in 1:length(delta_t)){
    if(delta_t[i] >=0){
      asl <- c(asl,asl_geq[i])
    }
    else{
      asl <- c(asl,asl_lt[i])
    }
  }

  asl <- asl/B

  n_sig = n_sig + as.numeric(asl<0.05)
}


##### SAVE DATA
power_vec <- n_sig/n_sims
names(power_vec) <- varying_factor[2:length(varying_factor)]

path = "/project2/gilad/awchen55/differentialDispersion/data/simulations/power_analysis"

save(power_vec, file = file.path(path,paste("power_analysis_100ef_indv_",n_indv,".RData",sep="")))




