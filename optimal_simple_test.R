library(tidyverse)

# Calculate optimal simple PAPs
# Based on deterministic cutoff rules for subsets of 1...n

pap_binary_data_simple = function(etaJ, # prob of observing each component
                               minp = .2, size = .1, # null hypothesis and size constraint
                               alpha = 1, beta = 1) { # parameters of Beta prior for theta
    
    n = length(etaJ)

    
    eta_sorted = sort(etaJ, index.return = T, decreasing = T)
    
    optimal_power = 0
    for (i in 1:n){ # consider the i components with highest values of etaJ
        # Find the size-controlling cutoff (for test based on no of successes >= z)
        achieved_size = 0
        for (z in seq(i,0, by=-1)){
            p = choose(i,z) * minp^z   * (1 - minp)^(1-z)
            if (achieved_size + p >size) {
                kappa =  (size - achieved_size) / p # Probability of rejecting for successes = z
                break
            }
            achieved_size = achieved_size + p
        }
        
        # calculate power
        # all subsets of 1:i
            J = map(0:(2 ^ i - 1), ~ as.numeric(intToBits(.x))[1:i]) 
            # bernoulli distributions for the components of J, with different probabilities
            P_J = map_dbl(J, ~ prod(.x * eta_sorted$x[1:i] + (1 - .x) * (1 - eta_sorted$x[1:i])))
            X = map(0:(2 ^ i - 1), ~ as.numeric(intToBits(.x))[1:i])
            P_X = map_dbl(X, ~ beta(sum(.x) + alpha, i - sum(.x) + beta) / beta(alpha, beta))
            power = 0
            for (a in 1:length(J)) {
                for (b in 1:length(X)){
                    if (sum(J[[a]]*X[[b]]) > z) {power = power + P_J[a] * P_X[b]}
                    else if (sum(J[[a]]*X[[b]]) == z) {power = power + kappa * P_J[a] * P_X[b]}
                }
            }
        
        if (power > optimal_power){
            optimal_power = power
            istar = i
            zstar = z
            kappastar = kappa
        }
        # print(c(i, power, P, expected_power, kappa, z))
    }
    # browser()
    list(test = tibble(Parameter = c("Components", "Cutoff for sum", "Rejection prob at the margin"),
                       Value = c(paste0(eta_sorted$ix[1:istar],collapse = ","), zstar, round(kappastar, digits = 2))),
         expected_power = optimal_power)
}





pap_normal_data_simple = function(etaJ, 
                           mu0, Sigma0,
                           mu, Sigma,
                           size = .1) {
    n = length(etaJ)
    J = map(0:(2 ^ n - 1), ~ as.numeric(intToBits(.x))[1:n]) 
   
    optimal_power = 0
    for (j in J){
        # Probability of observing at least the components j
        p_j_plus = prod(prod(j * etaJ + (1 - j)))
        
        # Calculate the distribution of the test statistic Z (sum of components in j) under the null and the alternative
        mu0_Z = sum(mu0*j)
        Sigma0_Z = j %*% Sigma0 %*% j
        mu_Z = sum(mu*j)
        Sigma_Z = j %*% Sigma %*% j
        
        # Critical value for Z
        z = qnorm(1-size, mu0_Z, sqrt(Sigma0_Z))
        # expected power
        power = p_j_plus*(1-pnorm(z, mu_Z, sqrt(Sigma_Z)))

        # print(c(j, z, power))
        if (power > optimal_power){
            optimal_power = power
            jstar = j
            zstar = z
        }
    }
    
    list(test = tibble(Parameter = c("Components", "Cutoff for sum"),
                       Value = c(paste0((1:n)[as.logical(jstar)],collapse = ","), round(zstar, digits = 3))),
         expected_power = optimal_power)
}








