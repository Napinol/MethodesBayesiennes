school = function(nchain, data, init, prop.sd){
  
  ## Initialisation avec les donnees
  LRT = data$LRT
  VR_1 = data$VR.1.
  VR_2 = data$VR.2.
  Gender = data$Gender
  School_gender_1 = data$School_gender.1.
  School_gender_2 = data$School_gender.2.
  School_denom_1 = data$School_denom.1.
  School_denom_2 = data$School_denom.2.
  School_denom_3 = data$School_denom.3.
  Y = data$Y
  
  X = cbind(1, LRT, VR_1,
            LRT^2, VR_2, Gender, School_gender_1, School_gender_2,
            School_denom_1, School_denom_2, School_denom_3)
  ## Initialisation de la precision pour les lois a priori
  tau_0 = 10^-4
  
  
  ## Chaine pour les alpha, beta et gamma
  chain.norm = matrix(NA, nchain+1, 14)
  colnames(chain.norm) = c("Alpha_1","Alpha_2","Alpha_3", 
                      "Beta_1","Beta_2","Beta_3","Beta_4","Beta_5","Beta_6","Beta_7","Beta_8",
                      "Gamma_1","Gamma_2","Gamma_3")
  chain.norm[1,] = init[1:14]
  
  ## Chaine pour les autres parametres theta et phi dans tau
  chain.tau = matrix(NA, nchain+1, 2)
  colnames(chain.tau) = c("Theta", "Phi")
  chain.tau[1,] = init[15:16]
  acc.rates = rep(0, 2)
  
  ## Chaine pour la wishart (sigma-1)
  sigma_chain = array(NA, dim = c(nchain+1,3,3))
  R = matrix(c(0.1, 0.005, 0.005, 
               0.005, 0.01, 0.005, 
               0.005, 0.005, 0.01), nrow = 3)
  R_inv = solve(R)
  sigma_chain[1,,] = rWishart(1, 3, R_inv)[,,1]
  
  ## Parametres de la loi a priori de gamma
  gamma_mean = c(0,0,0)
  gamma_sd = matrix(c(1, 0, 0,
                      0, 1, 0,
                      0, 0, 1), 3, 3, byrow=TRUE)
  
  
  for (iter in (1:nchain)){
    cat(iter)
    
    ## Current
    current.norm = chain.norm[iter,]
    current.tau = chain.tau[iter,]
    sigma = sigma_chain[iter,,] # Mise a jour de sigma
    inv_sigma = solve(sigma_chain[iter,,])
    
    # MAj des tau_ij et mu_ij
    tau = exp(current.tau[1]-current.tau[2]*LRT) # Mise à jour de tau_ij
    mu = X %*% current.norm[1:11] # Mise à jour de mu_ij
    
    ## MAJ des beta
    for (i in 1:8){
      update_mean = sum(X[,i+3]*tau*Y) / (tau_0 + sum(X[,i+3]^2*tau))
      update_sd = tau_0 + sum(X[,i]^2*tau)
      beta = rnorm(1, update_mean, update_sd)
      
      chain.norm[iter+1, i+3] = beta
    }
    
    
    ## MAJ de gamma
    inv_sigma = solve(sigma_chain[iter,,])
    inv_sigma_gamma = solve(gamma_sd)
    mu_gamma = c(0, 0, 0)    #loi a priori de gamma, normale de dim 3 non informative
    alpha = chain.norm[iter, 1:3]
    
    update_moy = solve(inv_sigma_gamma + 3 * inv_sigma) %*% 
      (inv_sigma_gamma %*% mu_gamma + 3 * inv_sigma %*% alpha)
    update_var = solve(inv_sigma_gamma + 3 * inv_sigma)
    
    gamma = MASS::mvrnorm(n=1,  mu=update_moy,  Sigma=update_var)
    
    # Update des gamma dans la chaine et dans le vecteur courant pour calculer alpha
    chain.norm[iter+1,12:14] = gamma
    current.norm[12:14] = gamma
    
    
    ## Mise a jour des alpha
    gamma = chain.norm[iter,12:14]
    beta = chain.norm[iter, 4:11]
    
    
    # probleme dans mise a jour de la moyenne, update_mean a directement des valeurs tres grandes
    update_mean = inv_sigma %*% gamma + sum(t(Y - X[,4:11] %*% beta) %*% X[,1:3])
    update_sd = solve(inv_sigma + t(X[,1:3]) %*% X[,1:3])
    
    alpha = MASS::mvrnorm(n=1,  mu=update_mean,  Sigma=update_sd)
    
    chain.norm[iter+1,1:3] = alpha
    
    
    
    ## MAJ pour theta avec M-H
    theta = chain.tau[iter, 1]
    phi = chain.tau[iter, 2]
    theta.prop = rnorm(1, theta, prop.sd[1])
    
    tau_ij = exp(theta + phi*LRT)
    tau_ij.prop = exp(theta.prop + phi*LRT)
    
    top = dnorm(theta.prop, 0, prop.sd[1], log = TRUE) + sum(dnorm(Y, mu, tau_ij.prop,log = TRUE))
    bottom = dnorm(1, theta, prop.sd[1], log = TRUE) + sum(dnorm(Y, mu, tau_ij, log = TRUE))
    
    acc.prob = exp(top - bottom)
    
    if (runif(1) < acc.prob){
      chain.tau[iter+1, 1] = theta.prop
      acc.rates[1] = acc.rates[1] + 1
    }
    
    chain.tau[iter+1,1] = theta
    
    
    
    ## MAJ pour phi avec M_H
    theta = chain.tau[iter, 1]
    phi = chain.tau[iter, 2]
    phi.prop = rnorm(1, phi, prop.sd[2])

    tau_ij = exp(theta + phi*LRT)
    tau_ij.prop = exp(theta + phi.prop*LRT)
    
    top = dnorm(phi.prop, 0, prop.sd[2], log = TRUE) + sum(dnorm(Y, mu, tau_ij.prop,log = TRUE))
    bottom = dnorm(1, phi, prop.sd[2], log = TRUE) + sum(dnorm(Y, mu, tau_ij, log = TRUE))
    
    acc.prob = exp(top - bottom)
    
    if (runif(1) < acc.prob){
      chain.tau[iter+1, 2] = phi.prop
      acc.rates[2] = acc.rates[2] + 1
    }
    
    chain.tau[iter+1, 2] = phi
    
    
    
    ## MAJ de sigma-1
    df = 6
    update_matrix = solve(R_inv + sum((current.norm[1:3] - current.norm[12:14])
                                         %*% t(current.norm[1:3] - current.norm[12:14]))) 
    # probleme dans l'inversion de sigma_1, car la matrice devient singuliere
    sigma_1 = rWishart(1, df=df, Sigma=update_matrix)
    
    sigma_chain[iter+1,,] = sigma_1[,,1]
    
    
  }
  
  return(list(chain.norm = chain.norm, chain.tau = chain.tau, acc.rates = acc.rates / nchain))
  
}
# Init : "Alpha_1","Alpha_2","Alpha_3", "Beta_1","Beta_2","Beta_3","Beta_4","Beta_5","Beta_6",
# "Beta_7","Beta_8", "Gamma_1","Gamma_2","Gamma_3", "theta", "phi"
init = rep(1, 16)
prop.sd = c(1, 1)
res = school(1000, data, init, prop.sd)
res

