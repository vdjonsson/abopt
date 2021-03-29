import statistics
import numpy as np
import cvxpy as cp
from collections.abc import Iterable

from sklearn.linear_model._base import LinearModel, RegressorMixin, MultiOutputMixin
from sklearn.linear_model._base import _preprocess_data
from sklearn.utils.validation import check_X_y, check_is_fitted, check_array, _deprecate_positional_args
from sklearn.utils.extmath import safe_sparse_dot
from sklearn.model_selection import KFold

@_deprecate_positional_args
def sat_separate_data(X, y, saturation = 'max'):
    """Used to separate data X and labels y into saturated data and unsaturated data.
    
    Data separated based on value of parameter 'saturation'.
        
    Parameters
    ----------
    X : ndarray of (n_samples, n_features)
        Data
    y : ndarray of shape (n_samples,) or \
        (n_samples, n_targets)
        Target. Will be cast to X's dtype if necessary
    saturation : float or string, {'mode', 'max'}, default='max'
        If float, regressors X with label y equal to float used as saturated data.
        If string in {'mode', 'max'}, use regressors X with label y equal to mode or maximum of labels y, respectively, as saturated data.
        
    Returns
    ----------
    Xu : ndarray of (n_unsaturated_samples, n_features)
        Unsaturated data according to saturated parameter
    Xs : ndarray of (n_saturated_samples, n_features)
        Saturated data according to saturated parameter
    yu : ndarray of shape (n_unsaturated_samples,) or \
        (n_unsaturated_samples, n_targets)
        Target of unsaturated data.
    ys : ndarray of shape (n_saturated_samples,) or \
        (n_saturated_samples, n_targets)
        Target of saturated data.
    saturation_val : float
        value chosen to determine saturation
    """
    if isinstance(saturation, float):
        if saturation not in y:
            error_msg = 'Saturation value passed : ' + str(saturation) + ' not found in y array.'
            raise ValueError(error_msg)
        saturation_val = saturation
        
    elif isinstance(saturation, str):
        if saturation not in ['mode', 'max']:
            # Raise ValueError if saturation arg not in accepted args
            raise ValueError('Saturation value must be in {"max", "mode"} or a float.')
        if saturation == 'max':
            saturation_val = max(y)
        elif saturation == 'mode':
            saturation_val = statistics.mode(y)
    else:
        # Raise ValueError if saturation arg not in accepted args
        raise ValueError('Saturation value must be in {"max", "mode"} or a float.')
    
    # Set unsaturated and saturated data based on saturated value
    Xu = X[y != saturation_val]
    Xs = X[y == saturation_val]
    yu = y[y != saturation_val]
    ys = y[y == saturation_val]
    
    return Xu, Xs, yu, ys, saturation_val

def objective_function(coefs, Xu, Xs, yu, ys, lambda_1, lambda_2, lambda_3):
    """Objective function for SatLasso method.
    
    The function returned:
    
    lambda_1 * ||y_u - X_uw||^2_2 + lambda_2 * ||w||_1 + lambda_3 * max(max(y_s-X_sw), 0)
        
    Parameters
    ----------
    Xu : ndarray of (n_unsaturated_samples, n_features)
        Unsaturated data according to saturated parameter
    Xs : ndarray of (n_saturated_samples, n_features)
        Saturated data according to saturated parameter
    yu : ndarray of shape (n_unsaturated_samples,) or \
        (n_unsaturated_samples, n_targets)
        Target of unsaturated data.
    ys : ndarray of shape (n_saturated_samples,) or \
        (n_saturated_samples, n_targets)
        Target of saturated data.
    lambda_1 : float
        Constant that multiplies the least squares loss.
    lambda_2 : float
        Constant that multiplies the L1 term.
    lambda_3 : float
        Constant that multiplies the penalty on saturated data.
    """
    
    # Convert to numpy arrays
    Xu = np.asarray(Xu)
    Xs = np.asarray(Xs)
    yu = np.asarray(yu)
    ys = np.asarray(ys)
    
    m = len(Xu)
    # Compute and return objective function
    # Check unsaturated and saturated data for empty arrays
    if yu.size > 0 and ys.size > 0:
        # Compute and return objective function with unsaturated and saturated loss and l1-regularization
        return lambda_1*(1/m)*cp.norm2(yu - Xu @ coefs)**2+lambda_2*(1/m)*cp.norm1(coefs)+lambda_3*cp.max(cp.hstack([ys - Xs @ coefs, 0]))
    elif yu.size > 0 and ys.size == 0:
        # Compute and return objective function with unsaturated penalty only and l1-regularization
        return lambda_1*(1/m)*cp.norm2(yu - Xu @ coefs)**2+lambda_2*(1/m)*cp.norm1(coefs)
    elif yu.size == 0 and ys.size > 0:
        # Compute and return objective function with saturated penalty only and l1-regularization
        return lambda_2*(1/m)*cp.norm1(coefs)+lambda_3*cp.max(cp.hstack([ys - Xs @ coefs, 0]))
    else:
        # If unsaturated and saturated data both empty, raise error
        raise ValueError('Encountered empty y: ', yu, ys)

@_deprecate_positional_args
def satlasso_cvxopt(Xu, Xs, yu, ys, lambda_1, lambda_2, lambda_3):
    """Compute optimal coefficient vector for saturated lasso problem using convex optimization.
        
    Parameters
    ----------
    Xu : ndarray of (n_unsaturated_samples, n_features)
        Unsaturated data according to saturated parameter
    Xs : ndarray of (n_saturated_samples, n_features)
        Saturated data according to saturated parameter
    yu : ndarray of shape (n_unsaturated_samples,) or \
        (n_unsaturated_samples, n_targets)
        Target of unsaturated data.
    ys : ndarray of shape (n_saturated_samples,) or \
        (n_saturated_samples, n_targets)
        Target of saturated data.
    lambda_1 : float
        Constant that multiplies the least squares loss.
    lambda_2 : float
        Constant that multiplies the L1 term.
    lambda_3 : float
        Constant that multiplies the penalty on saturated data.
    
    Returns
    ----------
    coefs : ndarray of shape (n_features,) or (n_targets, n_features)
        parameter vector found by cvxpy
    """
    
    coefs = cp.Variable(len(Xu[0]))
    problem = cp.Problem(cp.Minimize(objective_function(coefs, Xu, Xs, yu, ys, lambda_1, lambda_2, lambda_3)))
    solvers = [cp.ECOS, cp.SCS, cp.CVXOPT]
    for solver_choice in solvers:
        try:
            problem.solve(solver = solver_choice)
            break
        except cp.error.SolverError:
            continue
    return np.asarray(coefs.value)

class SatLasso(MultiOutputMixin, RegressorMixin, LinearModel):
    """Linear Model trained with L1 prior as regularizer (aka Lasso) and penalty on underestimated saturated data
    
    The optimization objective for SatLasso is::
    
        lambda_1 * ||y_u - X_uw||^2_2 + lambda_2 * ||w||_1 + lambda_3 * max(max(y_s-X_sw), 0)
        
    Parameters
    ----------
    lambda_1 : float, default=1.0
        Constant that multiplies the least squares loss. Defaults to 1.0.
    lambda_2 : float, default=1.0
        Constant that multiplies the L1 term. Defaults to 1.0.
        ``lambda_2 = 0`` is equivalent to an ordinary least square,
        with penalty on underestimated saturated data.
    lambda_3 : float, default=1.0
        Constant that multiplies the penalty on saturated data. Default to 1.0.
    saturation : float or string, {'mode', 'max'}, default='max'
        If float, regressors X with label y equal to float used as saturated data.
        If string in {'mode', 'max'}, use regressors X with label y equal to mode or maximum of labels y, respectively, as saturated data.
    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to False, no intercept will be used in calculations
        (i.e. data is expected to be centered).
    normalize : bool, default=False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.
    copy_X : bool, default=True
        If ``True``, X will be copied; else, it may be overwritten.
        
    Attributes
    ----------
    coef_ : ndarray of shape (n_features,) or (n_targets, n_features)
        parameter vector (w in the cost function formula)
    intercept_ : float or ndarray of shape (n_targets,)
        independent term in decision function.
    saturation_val_ : float
        value used to determine saturation
    """
    def __init__(self, lambda_1 = 1.0, lambda_2 = 1.0, lambda_3 = 1.0, fit_intercept = True, saturation = 'max', normalize = False, copy_X = True): # potentially add warm_start; combine lambdas?
        self.lambda_1 = lambda_1
        self.lambda_2 = lambda_2
        self.lambda_3 = lambda_3
        self.fit_intercept = fit_intercept
        self.saturation = saturation
        self.normalize = normalize
        self.copy_X = copy_X
    
    method = staticmethod(satlasso_cvxopt)
    
    def fit(self, X, y, check_input=True):
        """Fit model with convex optimization.
        Parameters
        ----------
        X : ndarray of (n_samples, n_features)
            Data
        y : ndarray of shape (n_samples,) or \
            (n_samples, n_targets)
            Target. Will be cast to X's dtype if necessary
        check_input : bool, default=True
            Allow to bypass several input checking.
            Don't use this parameter unless you know what you do.
            
        Returns
        -------
        self : object
            Returns self.
        """
        
        # Copy X if copy X is True
        if self.copy_X:
            X = X.copy()
        
        # Convert to numpy arrays
        X = np.asarray(X)
        y = np.asarray(y)
        
        # Check for correct input shape
        if check_input:
            X, y = check_X_y(X, y, accept_sparse=False)
        
        # Normalize with mean centering and l2-norm if normalize
        #   set to true, and fit_intercept set to true
        if self.normalize and self.fit_intercept:
            # Already copied so do not need to copy again
            X, y, X_offset, y_offset, X_scale = _preprocess_data(X, y, fit_intercept=self.fit_intercept, normalize=self.normalize, copy=False)
            if isinstance(self.saturation, float):
                self.saturation = self.saturation - y_offset
                
        # Separate unsaturated data from saturated data
        Xu, Xs, yu, ys, self.saturation_val_ = sat_separate_data(X, y, saturation = self.saturation)
        
        # Use convex optimization to solve for minimized objective function
        if self.fit_intercept and not self.normalize:
            # Add a bias variable to each data point if fit intercept = True
            Xu_with_bias = np.hstack((Xu, [[1] for i in range(0, len(Xu))]))
            Xs_with_bias = np.hstack((Xs, [[1] for i in range(0, len(Xs))]))
            
            # Use convex optimization to solve for coefficients
            coefs = self.method(Xu_with_bias, Xs_with_bias, yu, ys, self.lambda_1, self.lambda_2, self.lambda_3)
            self.coef_ = coefs[:-1]
            self.intercept_ = coefs[-1]
        else:
            # Use convex optimization to solve for coefficients
            self.coef_ = self.method(Xu, Xs, yu, ys, self.lambda_1, self.lambda_2, self.lambda_3)
            self.intercept_ = 0.
            
        # Set intercept and rescale coefficient if data was normalized
        if self.normalize and self.fit_intercept:
            self._set_intercept(X_offset, y_offset, X_scale)
            self.saturation_val_ = self.saturation_val_ + y_offset
            
        self.is_fitted_ = True
        return self
    
    def _decision_function(self, X):
        """Decision function of the linear model
        Parameters
        ----------
        X : numpy array of shape (n_samples, n_features)
        
        Returns
        -------
        T : ndarray of shape (n_samples,)
            The predicted decision function
        """
        
        # Check if fit has been called
        check_is_fitted(self, 'is_fitted_')
        
        # Check input
        X = check_array(X)
        return safe_sparse_dot(X, self.coef_.T, dense_output=True) + self.intercept_

@_deprecate_positional_args
def satlassoCV_cvxopt(Xu, Xs, yu, ys, lambda_1s, lambda_2s, lambda_3s, cv=5):
    """Compute optimal coefficient vector for saturated lasso problem using convex optimization.
        
    Parameters
    ----------
    Xu : ndarray of (n_unsaturated_samples, n_features)
        Unsaturated data according to saturated parameter
    Xs : ndarray of (n_saturated_samples, n_features)
        Saturated data according to saturated parameter
    yu : ndarray of shape (n_unsaturated_samples,) or \
        (n_unsaturated_samples, n_targets)
        Target of unsaturated data.
    ys : ndarray of shape (n_saturated_samples,) or \
        (n_saturated_samples, n_targets)
        Target of saturated data.
    lambda_1s : ndarray, default=None
        List of lambda_1s where to compute the models.
        If ``None`` lambda_1s are set automatically
    lambda_2s : ndarray, default=None
        List of lambda_2s where to compute the models.
        If ``None`` lambda_2s are set automatically
    lambda_3s : ndarray, default=None
        List of lambda_3s where to compute the models.
        If ``None`` lambda_3s are set automatically
    
    Returns
    ----------
    lambda_1 : float
        The amount of penalization on unsaturated data chosen by cross validation
    lambda_2 : float
        The amount of l1-norm penalization chosen by cross validation
    lambda_3 : float
        The amount of penalization on saturated data chosen by cross validation
    mse_dict : python dict
            keys : tuple of form (lambda_1, lambda_2, lambda_3)
            values : mean square error for values of lambda_1, lambda_2, lambda_3
        mean square error for the test set on each fold, varying lambda_1, lambda_2, lambda_3
    """
    
    # TO DO: only have one problem
    
    # Concatenate X and y arrays in order to split for KFold cross validation
    X = np.vstack((Xu, Xs))
    y = np.hstack((yu, ys))
    
    # Create iterable object to split training and test indices
    if isinstance(cv, int):
        # Check that cv does not exceed size of data
        if cv > len(X):
            raise ValueError('Cannot have number of splits cv=' + str(cv) + 'greater than the number of samples: n_samples='+ str(len(X)))
        # Create KFold object for iteration if int provided
        kfold = KFold(n_splits=cv, shuffle=True, random_state=0)
        _ = kfold.get_n_splits(X)
        cv_iter = list(kfold.split(X))
    elif isinstance(cv, Iterable):
        # Use iterable if provided
        cv_iter = list(cv)
    else:
        # Raise ValueError if cv not of accepted type
        raise ValueError('Expected cv as an integer, or an iterable.')
    
    # Iterate over possible lambda values and keep track of MSEs
    lambda_combns = np.array(np.meshgrid(lambda_1s, lambda_2s, lambda_3s)).T.reshape(-1,3)
    mses = {}
    for i in range(0, len(lambda_combns)):
        sses = []
        lambda_1, lambda_2, lambda_3 = lambda_combns[i]
        for train_index, test_index in cv_iter:
            # Split data
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]
            
            # Separate saturated and unsaturated data
            Xu_train = X_train[np.argwhere(y_train != np.unique(ys).item()).flatten()]
            Xs_train = X_train[np.argwhere(y_train == np.unique(ys).item()).flatten()]
            yu_train = y_train[np.argwhere(y_train != np.unique(ys).item()).flatten()]
            ys_train = y_train[np.argwhere(y_train == np.unique(ys).item()).flatten()]
            
            # Calculate optimal coefficient
            coefs = satlasso_cvxopt(Xu_train, Xs_train, yu_train, ys_train, lambda_1, lambda_2, lambda_3)
            
            # Calculate error on test data set
            y_predicted = safe_sparse_dot(X_test, coefs.T, dense_output=True)
            error = np.sum(np.square(y_test-y_predicted))
            sses.append(error)
        mses[tuple(lambda_combns[i])] = statistics.mean(sses)
        
    # Retrieve optimal lambda values from mses dictionary
    lambda_1, lambda_2, lambda_3 = min(mses, key = mses.get)
    return lambda_1, lambda_2, lambda_3, mses

class SatLassoCV(MultiOutputMixin, RegressorMixin, LinearModel):
    """Lasso linear model with iterative fitting along a regularization path.
    
    See glossary entry for :term:`cross-validation estimator`.
    
    The best model is selected by cross-validation.
    
    The optimization objective for SatLasso is::
    
        lambda_1 * ||y_u - X_uw||^2_2 + lambda_2 * ||w||_1 + lambda_3 * max(max(y_s-X_sw), 0)
        
    Parameters
    ----------
    n_lambdas : int, default=10
        Number of lambda_1 values, lambda_2 values, lambda_3 values along the regularization path
    lambda_1s : ndarray, default=None
        List of lambda_1s where to compute the models.
        If ``None`` lambda_1s are set automatically
    lambda_2s : ndarray, default=None
        List of lambda_2s where to compute the models.
        If ``None`` lambda_2s are set automatically
    lambda_3s : ndarray, default=None
        List of lambda_3s where to compute the models.
        If ``None`` lambda_3s are set automatically
    fit_intercept : bool, default=True
        whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).
    saturation : float or string, {'mode', 'max'}, default='max'
        If float, regressors X with label y equal to float used as saturated data.
        If string in {'mode', 'max'}, use regressors X with label y equal to mode or maximum of labels y, respectively, as saturated data.
    normalize : bool, default=False
        This parameter is ignored when ``fit_intercept`` is set to False.
        If True, the regressors X will be normalized before regression by
        subtracting the mean and dividing by the l2-norm.
        If you wish to standardize, please use
        :class:`sklearn.preprocessing.StandardScaler` before calling ``fit``
        on an estimator with ``normalize=False``.
    copy_X : bool, default=True
        If ``True``, X will be copied; else, it may be overwritten.
    cv : int, or iterable, default=5
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
        - int, to specify the number of folds.
        - An iterable yielding (train, test) splits as arrays of indices.
        
    Attributes
    ----------
    lambda_1_ : float
        The amount of penalization on unsaturated data chosen by cross validation
    lambda_2_ : float
        The amount of l1-norm penalization chosen by cross validation
    lambda_3_ : float
        The amount of penalization on saturated data chosen by cross validation
    mse_dict_ : python dict
            keys : tuple of form (lambda_1, lambda_2, lambda_3)
            values : mean square error for values of lambda_1, lambda_2, lambda_3
        mean square error for the test set on each fold, varying lambda_1, lambda_2, lambda_3
    coef_ : ndarray of shape (n_features,) or (n_targets, n_features)
        parameter vector (w in the cost function formula)
    intercept_ : float or ndarray of shape (n_targets,)
        independent term in decision function.
    lambdas_1s_ : ndarray of shape (n_lambdas_1s,)
        The grid of lambdas_1s used for fitting
    lambdas_2s_ : ndarray of shape (n_lambdas_2s,)
        The grid of lambdas_2s used for fitting
    lambdas_3s_ : ndarray of shape (n_lambdas_3s,)
        The grid of lambdas_3s used for fitting
    saturation_val_ : float
        value used to determine saturation
    """
    
    def __init__(self, n_lambdas = 10, lambda_1s = None, lambda_2s = None, lambda_3s = None, fit_intercept = True, saturation = 'max', normalize = False, copy_X = True, cv = 5):
        self.n_lambdas=n_lambdas
        self.lambda_1s = lambda_1s
        self.lambda_2s = lambda_2s
        self.lambda_3s = lambda_3s
        self.fit_intercept = fit_intercept
        self.saturation = saturation
        self.normalize = normalize
        self.copy_X = copy_X
        self.cv = cv

    cvmethod = staticmethod(satlassoCV_cvxopt)
    method = staticmethod(satlasso_cvxopt)
    
    def fit(self, X, y, check_input=True):
        """Fit model with convex optimization.
        Parameters
        ----------
        X : ndarray of (n_samples, n_features)
            Data
        y : ndarray of shape (n_samples,) or \
            (n_samples, n_targets)
            Target. Will be cast to X's dtype if necessary
        check_input : bool, default=True
            Allow to bypass several input checking.
            Don't use this parameter unless you know what you do.
            
        Returns
        -------
        self : object
            Returns self.
        """
        
        # Copy X if copy X is True
        if self.copy_X:
            X = X.copy()
        
        # Convert to numpy arrays
        X = np.asarray(X)
        y = np.asarray(y)
        
        # Check for correct input shape
        if check_input:
            X, y = check_X_y(X, y, accept_sparse=False)
        
        # Check for lambda values and update if None
        self.lambda_1s_ = self.lambda_1s
        self.lambda_2s_ = self.lambda_2s
        self.lambda_3s_ = self.lambda_3s
        
        if self.lambda_1s_ is None:
            self.lambda_1s_ = np.linspace(start=1, stop=10, num=self.n_lambdas_1s) # TO DO
        
        if self.lambda_2s_ is None:
            self.lambda_2s_ = np.linspace(start=1, stop=10, num=self.n_lambdas_2s) # TO DO
            
        if self.lambda_3s_ is None:
            self.lambda_3s_ = np.linspace(start=1, stop=10, num=self.n_lambdas_3s) # TO DO
        
        # Normalize with mean centering and l2-norm if normalize
        #   set to true, and fit_intercept set to true
        if self.normalize and self.fit_intercept:
            # Already copied so do not need to copy again
            X, y, X_offset, y_offset, X_scale = _preprocess_data(X, y, fit_intercept=self.fit_intercept, normalize=self.normalize, copy=False)
            if isinstance(self.saturation, float):
                self.saturation = self.saturation - y_offset
        
        # Separate unsaturated data from saturated data
        Xu, Xs, yu, ys, self.saturation_val_ = sat_separate_data(X, y, saturation = self.saturation)
        
        # Use convex optimization to solve for minimized objective function
        if self.fit_intercept and not self.normalize:
            # Add a bias variable to each data point if fit intercept = True
            Xu_with_bias = np.hstack((Xu, [[1] for i in range(0, len(Xu))]))
            Xs_with_bias = np.hstack((Xs, [[1] for i in range(0, len(Xs))]))
            
            # Use convex optimization to solve for coefficients
            self.lambda_1_, self.lambda_2_, self.lambda_3_, self.mse_dict_ = self.cvmethod(Xu_with_bias, Xs_with_bias, yu, ys, self.lambda_1s_, self.lambda_2s_, self.lambda_3s_, cv = self.cv)
            coefs = self.method(Xu_with_bias, Xs_with_bias, yu, ys, self.lambda_1_, self.lambda_2_, self.lambda_3_)
            self.coef_ = coefs[:-1]
            self.intercept_ = coefs[-1]
        else:
            # Use convex optimization to solve for coefficients
            self.lambda_1_, self.lambda_2_, self.lambda_3_, self.mse_dict_ = self.cvmethod(Xu, Xs, yu, ys, self.lambda_1s_, self.lambda_2s_, self.lambda_3s_, cv = self.cv)
            self.coef_ = self.method(Xu, Xs, yu, ys, self.lambda_1_, self.lambda_2_, self.lambda_3_)
            self.intercept_ = 0.
            
        # Set intercept and rescale coefficient if data was normalized
        if self.normalize and self.fit_intercept:
            self._set_intercept(X_offset, y_offset, X_scale)
            self.saturation_val_ = self.saturation_val_ + y_offset
            
        self.is_fitted_ = True
        return self
    
    def _decision_function(self, X):
        """Decision function of the linear model
        Parameters
        ----------
        X : numpy array of shape (n_samples, n_features)
        
        Returns
        -------
        T : ndarray of shape (n_samples,)
            The predicted decision function
        """
        
        # Check if fit has been called
        check_is_fitted(self, 'is_fitted_')
        
        # Check input
        X = check_array(X)
        return safe_sparse_dot(X, self.coef_.T, dense_output=True) + self.intercept_
