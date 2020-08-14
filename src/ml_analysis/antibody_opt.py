import numpy as np
import cvxpy as cp

def compute_antibodies(fitness_matrix, k):
    n = len(fitness_matrix[0])
    num_unmanaged = math.floor(n*k)
    c = cp.Variable(n, boolean = True)
    incidence_matrix = np.sign(fitness_matrix).clip(min=0)
    constraints = [cp.sum(c) >= 1, cp.sum_smallest(incidence_matrix@c, num_unmanaged+1) >= 1]
    objective = cp.Minimize(cp.norm1(c)-cp.matmul(cp.sum(fitness_matrix, axis=0), c))
    problem = cp.Problem(objective, constraints)
    result = problem.solve()
    return c.value
