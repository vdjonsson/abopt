import numpy as np
import random
import statistics
import math
from scipy.optimize import minimize

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

def satlasso_loss_function(beta, Xu, Xs, yu, ys, lmbdas):
    m = len(Xu)
    n = len(Xu)+len(Xs)
    loss = lmbdas[0]*(1/m)*np.linalg.norm(np.add(yu, -1*np.matmul(Xu,beta)), ord=2) + lmbdas[1]*(1/m)*np.linalg.norm(beta, ord=1) + lmbdas[2]*np.max(np.concatenate((np.add(ys, -1*np.matmul(Xs,beta)), [0])))
    return loss

def compute_gradient_satlasso_loss(beta, Xu, Xs, yu, ys, lmbdas):
    m = len(Xu)
    if np.all(np.add(ys, -1*np.matmul(Xs,beta)) < 0):
        sat_grad = np.zeros(len(beta))
    else:
        max_element = np.argmax(np.add(ys, -1*np.matmul(Xs,beta)))
        sat_grad = -1*np.array(Xs[max_element])
        
    l1_grad = np.sign(beta)
    l1_grad[l1_grad == 0] = np.random.uniform(low=-1,high=1, size=len(l1_grad[l1_grad == 0]))
    
    grad = lmbdas[0]*(1/m)*1/(np.linalg.norm(np.add(yu, -1*np.matmul(Xu,beta)), ord=2))*np.sum(np.add(yu, -1*np.matmul(Xu,beta)))*np.sum(-1*np.array(Xu), axis=0)+lmbdas[1]*(1/m)*l1_grad+lmbdas[2]*sat_grad
    return grad
    
def gradient_descent(initial_beta, Xu, Xs, yu, ys, alpha, lmbdas, tol):
    while(True):
        beta = np.add(initial_beta, -1*alpha*compute_gradient_satlasso_loss(initial_beta, Xu, Xs, yu, ys, lmbdas)).tolist()
        if(np.linalg.norm(np.add(beta, -1*np.array(initial_beta)), ord=2) < tol):
            break
        initial_beta = beta
    return beta

def satlasso(Xu, Xs, yu, ys, alpha, lmbdas, tol):
    # Xu - unsaturated data
    # Xs - saturated data
    # lmbdas - lambda values (lambda1, lambda2, lambda3)
    # yu - unsaturated labels
    # ys - saturated labels
    Xu_with_bias = np.hstack((Xu, [[1] for i in range(0, len(Xu))])).tolist()
    Xs_with_bias = np.hstack((Xs, [[1] for i in range(0, len(Xs))])).tolist()
    
    initial_beta = [0]*len(Xu_with_bias[0])
    beta = gradient_descent(initial_beta, Xu_with_bias, Xs_with_bias, yu, ys, alpha, lmbdas, tol)
    # min_obj = minimize(satlasso_loss_function, initial_beta, args=(Xu_with_bias, Xs_with_bias, yu, ys, lmbdas))
    # beta = min_obj.x
    return beta

def cross_validation(Xu, Xs, yu, ys, alphas, lmbda1_vals, lmbda2_vals, lmbda3_vals, num_folds, tol):
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
    
    param_combns = np.array(np.meshgrid(alphas, lmbda1_vals, lmbda2_vals, lmbda3_vals)).T.reshape(-1,4)
    mses = []
    for i in range(0, len(param_combns)):
        alpha = param_combns[i][0]
        lmbda_combn = [param_combns[i][1], param_combns[i][2], param_combns[i][3]]
        print('Trying alpha = '+str(alpha)+', lambda1 = '+str(lmbda_combn[0])+', lambda2 = '+str(lmbda_combn[1])+', lambda3 = '+str(lmbda_combn[2]))
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
            beta = satlasso(train_data_Xu, train_data_Xs, train_data_yu, train_data_ys, alpha, lmbda_combn, tol)
            predicted_ys = np.add(np.matmul(test_data_X,beta[0:len(beta)-1]), [beta[len(beta)-1]]*len(test_data_X))
            error = np.sum(np.square(np.add(test_data_y, [x*-1 for x in predicted_ys])))
            sses.append(error)
        mses.append(statistics.mean(sses))
        
    params = param_combns[mses.index(min(mses))]
    return (params[0], [params[1], params[2], params[3]])

def satlasso_solve(X, y, alphas, lmbda1_vals, lmbda2_vals, lmbda3_vals, num_folds, tol):
    # X - data
    # y - labels
    # lmbda1_vals - potential values for lambda1
    # lmbda2_vals - potential values for lambda2
    # lmbda3_vals - potential values for lambda3
    # num_folds - number of folds for k-fold cross validation
    Xu, Xs, yu, ys = separate_data(X,y)
    optimal_alpha, optimal_lmbdas = cross_validation(Xu, Xs, yu, ys, alphas, lmbda1_vals, lmbda2_vals, lmbda3_vals, num_folds, tol)
    beta = satlasso(Xu, Xs, yu, ys, optimal_alpha, optimal_lmbdas)
    return beta
