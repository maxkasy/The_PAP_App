library(tidyverse)
library(lpSolve) # For finding the optimal matching
library(mvtnorm)
library(hrbrthemes)
library(latex2exp)

setwd(here::here())

optimal_test = function(X, #list of values for x
                        P_X, #vector mapping elements of X into prior marginal probabilities
                        P_J, # vector mapping J (binary interpreted as integer) into probability of observing these components of X
                        P0_X, # vector mapping elements of X into probabilities under the null
                        size) # size constraint
{
    n = length(X[[1]]) # number of components of X
    N = length(X) # number of possible values for X
    
    dim = N + N * 2 ^ n # dimension of parameter vector (t, actreduced)
    
    objective = c(rep(0, N),
                  rep(P_X, 2 ^ n) * rep(P_J, each = N))
    
    # Next we construct the constraint matrices
    # In dense matrix representation, with columns representing rows, columns, and values of a standard matrix
    dense_size_constraint = matrix(c(rep(1,N), 1:N, P0_X),N,3)

    # [0,1] support constraint for the solution. Rows starting after the size constraint.
    dense_support_constraint = matrix(c(2:(N+1),1:N,rep(1,N)), N, 3)

    # Initializing monotonicity constraint    
    dense_monotonicity_constraint = matrix(0, 2 * (N ^ 2) * (2 ^ n), 3)
    i = 0 # row in dense constraint matrix
    j = N+1 # row in original constraint matrix. Starting at N, after the size constraints.
    for (w_int in 0:(2 ^ n - 1)) {
        w = as.logical(intToBits(w_int))[1:n]
        for (x_int in 1:N) {
            for (xprime_int in 1:N) {
                if (all(X[[x_int]][w] == X[[xprime_int]][w])) {
                    i = i + 1
                    j = j + 1
                    dense_monotonicity_constraint[i,1] = j
                    dense_monotonicity_constraint[i,2] = xprime_int
                    dense_monotonicity_constraint[i,3] = -1
                    i = i + 1
                    dense_monotonicity_constraint[i,1] = j
                    dense_monotonicity_constraint[i,2] = N + x_int + N * w_int
                    dense_monotonicity_constraint[i,3] = 1
                }
            }
        }
    }
    # Drop unused rows
    dense_monotonicity_constraint = dense_monotonicity_constraint[1:i, ]
    
    # combine constraints
    dense_constr_mat = rbind(dense_size_constraint,
                             dense_support_constraint,
                             dense_monotonicity_constraint)
    cap = c(size, rep(1, N), rep(0, j - 1 - N))
   
    solve = lp(
        direction = "max",
        objective.in = objective,
        dense.const = dense_constr_mat,
        const.dir = rep("<=", j),
        const.rhs = cap
    )
    
    # return only the solution for the full data test
    t = solve$solution[1:N]
    names(t) = map(X, ~ paste0(round(.x, digits=2), collapse = ","))
    t
}






binary_data = function(etaJ, # prob of observing each component
                       minp = .2, size = .1, # null hypothesis and size constraint
                       alpha = 1, beta = 1) { # parameters of Beta prior for theta
    n = length(etaJ)
    J = map(0:(2 ^ n - 1), ~ as.numeric(intToBits(.x))[1:n]) 
    # bernoulli distributions for the components of J, with different probabilities
    P_J = map_dbl(J, ~ prod(.x * etaJ + (1 - .x) * (1 - etaJ)))
    
    X = map(0:(2 ^ n - 1), ~ as.numeric(intToBits(.x))[1:n])
    # bernoulli distributions for X under the null
    P0_X = map_dbl(X, ~ prod(.x * minp + (1 - .x) * (1 - minp)))
    # beta-binomial distribution for sum(x), for beta prior over theta
    P_X = map_dbl(X, ~ beta(sum(.x) + alpha, n - sum(.x) + beta) / beta(alpha, beta))

    t = optimal_test(X, P_X, P_J, P0_X, size)

    # Formatting output
    X_tibble =
        do.call(rbind, X) |> 
        as_tibble(.name_repair = "unique") 
    names(X_tibble) = paste0("X", 1:n)
    X_tibble|> 
        mutate(t = t)
}

normal_data = function(etaJ, 
                       mu0, Sigma0,
                       mu, Sigma,
                       size = .1,
                       steps = 10) {
    n = length(etaJ)
    J = map(0:(2 ^ n - 1), ~ as.numeric(intToBits(.x))[1:n]) 
    # bernoulli distributions for the components of J, with different probabilities
    P_J = map_dbl(J, ~ prod(.x * etaJ + (1 - .x) * (1 - etaJ)))
  
    #Create list of all combinations of possible values for the steps, for n componens
    # QUESTION: Can we choose this more intelligently, e.g. dropping lower values where we would not reject anyway?    
    pstep_vec = seq(1/steps,1,by=1/steps)
    X_component_quantiles = map(1:n, ~ qnorm(pstep_vec, mean = mu0[.x], sd = sqrt(Sigma0[.x,.x])))

    X_steps = expand_grid(!!!rep(list(1:length(pstep_vec)), n)) |> 
        t() |> as.data.frame() |> as.list()
    # Create list of corresponding quantiles for the P0 distribution    
    X_q_fun = function(step) { # helper function to return quantile values for a vector
        xq = rep(0,n)
        for (i in 1:n) {
            if (step[i] > 0) {xq[i] = X_component_quantiles[[i]][step[i]]}
            else {xq[i] = -Inf}
        }
        xq
    }
    X = map(X_steps, X_q_fun)

    # normal distributions for X under the null
    P0_X = map_dbl(X_steps, 
                   ~ pmvnorm(lower = X_q_fun(.x - 1),
                             upper = X_q_fun(.x),
                             mean = mu0, 
                             sigma = Sigma0))
    
    # normal distributions for X, integrating over the prior
    P_X = map_dbl(X_steps, 
                  ~ pmvnorm(lower = X_q_fun(.x - 1),
                            upper = X_q_fun(.x),
                            mean = mu, 
                            sigma = Sigma))

    t = optimal_test(X, P_X, P_J, P0_X, size)

    # Formatting output
    X_upper = do.call(rbind, 
                       X) |>  
        as_tibble(.name_repair = "unique")
    names(X_upper) = paste0("X_upper", 1:n)
    X_lower = do.call(rbind, 
                      map(X_steps, ~ X_q_fun(.x-1))) |>  
        as_tibble(.name_repair = "unique")
    names(X_lower) = paste0("X_lower", 1:n)
    bind_cols(X_upper, X_lower) |> 
        mutate(t = t)
}

