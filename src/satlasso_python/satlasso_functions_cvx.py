import numpy as np
import random
import statistics
import math
import cvxpy as cp

def partition_data(data, num_partitions):
    size = math.floor(len(data)/num_partitions)
    random.shuffle(data)
    partitions = []
    for i in range(0, num_partitions):
        if i==num_partitions-1:
            partitions.append(data[i*size:len(data)])
        else:
            partitions.append(data[i*size:(i+1)*size])
    return tuple(partitions)

def separate_data(X, y):
    # X - data
    # y - labels
    mode = statistics.mode(y)
    saturated_indices = [i for i in range(len(y)) if y[i] == mode]
    unsaturated_indices = [i for i in range(len(y)) if y[i] != mode]
    yu = [y[i] for i in unsaturated_indices]
    ys = [y[i] for i in saturated_indices]
    Xu = [X[i] for i in unsaturated_indices]
    Xs = [X[i] for i in saturated_indices]
    return (Xu, Xs, yu, ys)

def objective_function(beta, Xu, Xs, yu, ys, lmbdas):
    m = len(Xu)
    n = len(Xu)+len(Xs)
    return lmbdas[0]*(1/m)*cp.norm2(np.array(yu) - np.array(Xu) @ beta)**2+lmbdas[1]*(1/m)*cp.norm1(beta)+lmbdas[2]*cp.max(cp.hstack([np.array(ys) - np.array(Xs) @ beta, 0]))

def satlasso(Xu, Xs, yu, ys, lmbdas):
    # Xu - unsaturated data
    # Xs - saturated data
    # lmbdas - lambda values (lambda1, lambda2, lambda3)
    # yu - unsaturated labels
    # ys - saturated labels
    Xu_with_bias = np.hstack((Xu, [[1] for i in range(0, len(Xu))])).tolist()
    Xs_with_bias = np.hstack((Xs, [[1] for i in range(0, len(Xs))])).tolist()
    
    beta = cp.Variable(len(Xu_with_bias[0]))
    problem = cp.Problem(cp.Minimize(objective_function(beta, Xu_with_bias, Xs_with_bias, yu, ys, lmbdas)))
    solvers = [cp.ECOS, cp.SCS, cp.CVXOPT]
    for solver_choice in solvers:
        try:
            problem.solve(solver = solver_choice)
            print('Used solver: '+solver_choice)
            break
        except SolverError:
            continue
    return beta.value.tolist()
    
def cross_validation(Xu, Xs, yu, ys, lmbda1_vals, lmbda2_vals, lmbda3_vals, num_folds):
    # Xu - unsaturated data
    # Xs - saturated data
    # lmbda1_vals - potential values for lambda1
    # lmbda2_vals - potential values for lambda2
    # lmbda3_vals - potential values for lambda3
    # yu - unsaturated labels
    # ys - saturated labels
    X = np.vstack((Xu, Xs)).tolist()
    y = yu + ys
    data = np.hstack((X, [[y[i]] for i in range(0, len(y))])).tolist()
    partitions = partition_data(data, num_folds)
    
    lmbda_combns = np.array(np.meshgrid(lmbda1_vals, lmbda2_vals, lmbda3_vals)).T.reshape(-1,3)
    mses = []
    for i in range(0, len(lmbda_combns)):
        print('Trying lambda1 = '+str(lmbda_combns[i][0])+', lambda2 = '+str(lmbda_combns[i][1])+', lambda3 = '+str(lmbda_combns[i][2]))
        sses = []
        for j in range(0, num_folds):
            test_data = partitions[j]
            test_data_X = np.array(test_data)[:,:-1].tolist()
            test_data_y = np.array(test_data)[:,-1].tolist()
            
            train_data = np.vstack(list(partitions[n] for n in range(0,len(partitions)) if n!=j)).tolist()
            train_data_X = np.array(train_data)[:,:-1].tolist()
            train_data_y = np.array(train_data)[:,-1].tolist()
            train_data_Xu = [train_data_X[i] for i in np.array(np.all((np.array(train_data_X)[:,None,:]==np.array(Xu)[None,:,:]),axis=-1).nonzero()).tolist()[0]]
            train_data_Xs = [train_data_X[i] for i in np.array(np.all((np.array(train_data_X)[:,None,:]==np.array(Xs)[None,:,:]),axis=-1).nonzero()).tolist()[0]]
            train_data_yu = [train_data_y[i] for i in np.array(np.all((np.array(train_data_X)[:,None,:]==np.array(Xu)[None,:,:]),axis=-1).nonzero()).tolist()[0]]
            train_data_ys = [train_data_y[i] for i in np.array(np.all((np.array(train_data_X)[:,None,:]==np.array(Xs)[None,:,:]),axis=-1).nonzero()).tolist()[0]]
            beta = satlasso(train_data_Xu, train_data_Xs, train_data_yu, train_data_ys, lmbda_combns[i])
            predicted_ys = np.add(np.matmul(test_data_X,beta[0:len(beta)-1]), [beta[len(beta)-1]]*len(test_data_X))
            error = np.sum(np.square(np.add(test_data_y, [x*-1 for x in predicted_ys])))
            sses.append(error)
        mses.append(statistics.mean(sses))
        
    lmbdas = lmbda_combns[mses.index(min(mses))]
    return lmbdas

def satlasso_solve(X, y, lmbda1_vals, lmbda2_vals, lmbda3_vals, num_folds):
    # X - data
    # y - labels
    # lmbda1_vals - potential values for lambda1
    # lmbda2_vals - potential values for lambda2
    # lmbda3_vals - potential values for lambda3
    # num_folds - number of folds for k-fold cross validation
    Xu, Xs, yu, ys = separate_data(X,y)
    optimal_lmbdas = cross_validation(Xu, Xs, yu, ys, lmbda1_vals, lmbda2_vals, lmbda3_vals, num_folds)
    print(optimal_lmbdas)
    beta = satlasso(Xu, Xs, yu, ys, optimal_lmbdas)
    return beta
