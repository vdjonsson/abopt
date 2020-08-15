import numpy as np
import cvxpy as cp
import pandas as pd

def compute_antibodies(fitness_matrix, k, lmbda):
    m,n = fitness_matrix.shape
    num_unmanaged = int(m*k)
    c = cp.Variable(n, boolean = True)
    incidence_matrix = np.sign(fitness_matrix).clip(min=0)
    constraints = [cp.sum(c) >= 1, cp.sum_smallest((incidence_matrix@c), num_unmanaged+1) >= 1]
    objective = cp.Minimize(lmbda*cp.norm1(c)-cp.matmul(cp.sum(fitness_matrix, axis=0), c))
    problem = cp.Problem(objective, constraints)
    result = problem.solve()
    return c.value

filepath = '../data/'
filename = 'antibody_matrix'
k = 0.15
lmbda = 580

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0, index_col=0)
fitness_matrix = df.values
results = compute_antibodies(fitness_matrix, k, lmbda)
print(df.columns[results.astype(bool)])
