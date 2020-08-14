import numpy as np
import cvxpy as cp

def compute_antibodies(fitness_matrix):
    n = len(fitness_matrix[0])
    c = cp.Variable(n, boolean = True)
    incidence_matrix = np.sign(fitness_matrix)
    
    objective = cp.Minimize(-1*cp.sum(np.sign(incidence_matrix @ c)) - cp.matmul(cp.sum(fitness_matrix, axis=0), c) + cp.norm1(c))
    problem = cp.problem(objective)
    x = problem.solve()
    return x.value
