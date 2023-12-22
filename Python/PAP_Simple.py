import numpy as np
import pandas as pd
import itertools as it
from scipy.special import comb, beta
from scipy.stats import norm

def pap_binary_data_simple(input_binary, size = .05):
    """From a  specification of null-hypothesis and prior,
    calculate the optimal simple (cutoff) PAP,
    for binary data X."""
    
    etaJ = input_binary["etaJ"]
    minp = input_binary["minp"]
    a = input_binary["alpha"]
    b = input_binary["beta"]    
    n = len(etaJ)

    eta_indices = np.argsort(etaJ)[::-1]
    eta_sorted = etaJ[eta_indices]    

    optimal_power = 0
    # consider a test based on the i components with the highest values of etaJ
    for i in range(1, n+1): 
        # Find the size-controlling cutoff
        achieved_size = 0
        for z in range(i, -1, -1):
            p = comb(i, z) * minp**z * (1 - minp)**(i-z)
            if achieved_size + p >= size:
                kappa = (size - achieved_size) / p
                break
            achieved_size += p
        
        # calculate expected power        
        J = [np.array(j) for j in it.product([False, True], repeat=i)]
        P_J = [np.prod(eta_sorted[0:i][j]) * np.prod(1-eta_sorted[0:i][~j]) for j in J]
        X = [np.array(j) for j in it.product([0, 1], repeat=i)]
        P_X = [beta(sum(x) + a, i - sum(x) + b) / beta(a, b) for x in X]        

        power = 0
        for a in range(len(J)):
            for b in range(len(X)):
                j = J[a]
                P_j = P_J[a]
                x = X[b]
                P_x = P_X[b]

                if sum(j*x) > z:
                    power += P_j * P_x
                elif sum(j*x) == z:
                    power += kappa * P_j * P_x        

        if power > optimal_power:
            optimal_power = power
            istar = i
            zstar = z
            kappastar = kappa

    return {'test': pd.DataFrame({
                'Parameter': ['Components', 'Cutoff for sum', 'Rejection prob at the margin'],
                'Value': [', '.join(map(str, 1+eta_indices[0:istar])), zstar, round(kappastar, 2)]
            }),
        'power': optimal_power
    }


def pap_normal_data_simple(input_normal, size=0.05):
    """From a  specification of null-hypothesis and prior,
    calculate the optimal simple (cutoff) PAP,
    for normal data X."""
    
    etaJ = input_normal["etaJ"] 
    mu0 = input_normal["mu0"]
    Sigma0 = input_normal["Sigma0"] 
    mu = input_normal["mu"]
    Sigma = input_normal["Sigma"]
    n = len(etaJ)
    
    J = [np.array(j) for j in it.product([False, True], repeat=n)]

    optimal_power = 0
    for j in J[1:]:
        # Probability of observing at least the components j
        p_j_plus = np.prod(etaJ[j])

        # Distribution of the test statistic Z (sum of components in j) under the null and the alternative
        mu0_Z = np.sum(mu0[j])
        Sigma0_Z = np.sum(Sigma0[j,j]) 
        mu_Z = np.sum(mu[j])
        Sigma_Z = np.sum(Sigma[j,j]) 

        # Critical value for Z
        z = norm.ppf(1 - size, mu0_Z, np.sqrt(Sigma0_Z))
        # expected power
        power = p_j_plus * (1 - norm.cdf(z, mu_Z, np.sqrt(Sigma_Z)))

        if power > optimal_power:
            optimal_power = power
            jstar = j
            zstar = z
            
    return {'test': pd.DataFrame({'Parameter': ["Components", "Cutoff for sum"],
                                  'Value': [', '.join(map(str, 1 + np.nonzero(jstar)[0])), round(zstar, 3)]}),
            'power': optimal_power}