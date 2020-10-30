import numpy as np
import cvxpy as cp
import pandas as pd

def compute_antibodies(fitness_matrix, k, lmbda, virus_matrix):
    m,n = fitness_matrix.shape
    num_unmanaged = int(m*k)
    c = cp.Variable(n, boolean = True)
    incidence_matrix = np.sign(fitness_matrix).clip(min=0)
    virus_incidence_matrix = np.sign(virus_matrix).clip(min=0)
    constraints = [cp.sum(c) >= 1, cp.sum_smallest((incidence_matrix@c), num_unmanaged+1) >= 1]
    #objective = cp.Minimize(lmbda*cp.norm1(c)-cp.matmul(cp.sum(fitness_matrix, axis=0), c))
    objective = cp.Minimize(lmbda*cp.norm1(c)-cp.matmul(cp.sum(fitness_matrix, axis=0), c)-incidence_matrix@c@virus_incidence_matrix)
    problem = cp.Problem(objective, constraints)
    result = problem.solve()
    return c.value

filepath = '../data/'
filename = 'rbd_ab_fitness'
virus_filename = 'rbd_ace2_fitness'
k = 0.15
lmbda = 580

df = pd.read_csv(filepath+filename+'.csv', sep=',', header=0, index_col=0)
virus_df = pd.read_csv(filepath+virus_filename+'.csv', sep=',', header=0, index_col=0)
virus_df = virus_df[~df.isnull().any(axis=1).values]
df = df.dropna()
virus_matrix = virus_df.values
fitness_matrix = df.values
results = compute_antibodies(fitness_matrix, k, lmbda, virus_matrix)
choices = df.columns[results.astype(bool)]

with open(filepath+'opt_results_new.txt', 'w') as writer:
    for choice in choices:
        writer.write(choice+'\n')
