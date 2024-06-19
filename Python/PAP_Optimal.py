import numpy as np
import pandas as pd
import cvxpy as cp
import itertools as it
from scipy.special import beta
from scipy.stats import mvn, norm
from scipy.sparse import csr_matrix, lil_matrix
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

def powerset(iterable):
    """Helper function to create powerset from iterable"""
    s = list(iterable)
    return [list(j) for j in it.chain.from_iterable(it.combinations(s, r) for r in range(len(s)+1))]

def optimal_test(test_args, size = .05):
    """Find the optimal test, which maximizes expected power,
    subject to the constraints of size control and implementability.
    test_args is a list which needs to include the matrix X 
    (with rows for possible realizations), and dictionaries P_J, P_X, and P0_X
    (with the probabilities of these realizations)."""

    X = test_args["X"]
    P_X = test_args["P_X"]
    P0_X = test_args["P0_X"]
    P_J = test_args["P_J"]
    
    N,n = X.shape
    K = list(range(n))
    JJ = powerset(K)

    # Calculate distribution of X_J, J from distribution of X and J
    P_pi = {(tuple(x[j]),tuple(j)):0 for x in X for j in JJ}
    for x in X:
        for j in JJ:
            P_pi[(tuple(x[j]), tuple(j))] += P_X[tuple(x)] * P_J[tuple(j)]

    dimb = len(P_pi)
    dim = N + dimb # dimension of t + dimension of b 

    # Define objective and constraints
    objective = np.concatenate((np.zeros(N), 
                            np.array(list(P_pi.values()))))

    size_constraint = np.concatenate((np.array(list(P0_X.values())), 
                                  np.zeros(dimb)))



    # helper dictionaries to pin down indices in constraint matrix
    t_index = {tuple(x): index for index, x in enumerate(X)}
    b_index = {key: N + index for index, key in enumerate(P_pi)}
    
    # Construct monotonicity constraints by iteration
    monotonicity_constraint = lil_matrix((N * 2 ** n, dim))   
    r = 0 # row in constraint matrix  
    for x in X:
        for j in JJ:
            # add one constraint to the constraint matrix
            monotonicity_constraint[r,t_index[tuple(x)]] = -1
            monotonicity_constraint[r,b_index[(tuple(x[j]),tuple(j))]] = 1
            r+=1

    # Define objective and constraint for cvxpy
    tb = cp.Variable(dim)
    cp_objective = cp.Maximize(objective @ tb)
    cp_constraints = [size_constraint @ tb <= size,
                      tb <= 1, tb >=0,
                      monotonicity_constraint @ tb <= 0]
    
    cp_problem = cp.Problem(cp_objective, cp_constraints)    
    # Solve the problem
    result = cp_problem.solve(solver = cp.CLARABEL)

    # Format output
    t = np.hstack((X,tb.value[:N].round(3).reshape((N,1))))
    t = pd.DataFrame(t, 
                  columns=[f"X{num}" for num in range(1,X.shape[1]+1)] + ["t"])
    if "X_lower" in test_args: # Add lower bounds for cells to the output, for normal data
        Xl = pd.DataFrame(test_args["X_lower"], 
                  columns=[f"Xl{num}" for num in range(1,X.shape[1]+1)])
        t = pd.concat([Xl, t], axis = 1)
    return {"tb": tb.value.round(3), 
            "power": result, 
            "t": t}
    
def pap_binary_data(input_binary):
    """From a  specification of null-hypothesis and prior,
    prepare the test_args required for optimal_test,
    for binary data X."""
    
    etaJ = input_binary["etaJ"]
    minp = input_binary["minp"]
    a = input_binary["alpha"]
    b = input_binary["beta"]    
    n = len(etaJ)

    K = list(range(n))
    JJ = powerset(K)
    # Distribution of J
    P_J = {tuple(j):np.prod(etaJ[list(j)]) * np.prod(1-etaJ[list(set(K) - set(j))]) for j in JJ}

    # Construct matrix with binary entries in rows for all possible X vectors
    X = np.array(list(it.product([0, 1], repeat=n)))

    # Distribution of X under the null and averaged over the prior
    P0_X = {tuple(x):np.prod(x * minp + (1 - x) * (1 - minp)) for x in X}
    P_X = {tuple(x):beta(sum(x) + a, n - sum(x) + b) / beta(a, b) for x in X}

    return {"X": X, "P_J": P_J, "P0_X": P0_X, "P_X": P_X}


def pap_normal_data(input_normal, steps = 40):
    """From a  specification of null-hypothesis and prior,
    prepare the test_args required for optimal_test,
    for normal data X."""
    
    etaJ = input_normal["etaJ"] 
    mu0 = input_normal["mu0"]
    Sigma0 = input_normal["Sigma0"] 
    mu = input_normal["mu"]
    Sigma = input_normal["Sigma"]
    n = len(etaJ)
    
    K = list(range(n))
    JJ = powerset(K)
    # Distribution of J
    P_J = {tuple(j):np.prod(etaJ[list(j)]) * np.prod(1-etaJ[list(set(K) - set(j))]) for j in JJ}
    
    # Vector of normal quantiles for discretization
    quantile_vec = np.linspace(0,1,steps+1)
    quantile_vec[0] = 1e-10 # to avoid infinite values
    quantile_vec[steps] = 1-1e-10 # to avoid infinite values
    norm_step_vec = norm.ppf(quantile_vec)
    
    X_component_quantiles = np.array([mu0[i] + np.sqrt(Sigma0[i,i]) * norm_step_vec for i in range(n)]).transpose()
    # all index combinations for the different components of X
    indices = np.array(list(it.product(range(1,len(norm_step_vec)), repeat = n)))
    X_upper = X_component_quantiles[indices,  1:2].squeeze()
    X_lower = X_component_quantiles[indices-1,1:2].squeeze()

    P_X  = {tuple(X_upper[r]):mvn.mvnun(X_upper[r], X_lower[r], mu, Sigma)[0]   for r in range(len(X_upper))}
    P0_X = {tuple(X_upper[r]):mvn.mvnun(X_upper[r], X_lower[r], mu0, Sigma0)[0] for r in range(len(X_upper))}
    
    return {"X": X_upper, "X_lower": X_lower, "P_J": P_J, "P0_X": P0_X, "P_X": P_X}



def plot_normal_pap(t, steps = 40, max_coord = 2.5):
    """Plotting the optimal PAP for the normal case, when n = 2"""

    fig, ax = plt.subplots()
    X1 = np.append(np.unique(t["Xl1"]), max_coord)
    X2 = np.append(np.unique(t["Xl2"]), max_coord)
    values = np.reshape(t['t'], (steps, steps))

    ax.pcolormesh(X1, X2, values, cmap='Blues')

    ax.set_xlim([-max_coord,max_coord])
    ax.set_ylim([-max_coord,max_coord])
    ax.set_aspect('equal')

    plt.close(fig)

    return fig