import math
import numpy as np
import scipy.optimize

_tolerance = 0.02
_delta = 0.0001
_maxiters = 20
_maxdistance = 3


def kimura_distance(N):
    '''
    This method is not tied to a rate matrix. It is basically a classic
    model-free heuristic. It is implemented here primarily to provide a
    starting distance for ml_distance_estimate.

    Note: written for protein alignments
    '''
    n_sum = np.sum(N)
    changed_fraction = 1 -  (np.sum(np.diag(N)) / n_sum)
    adjusted = changed_fraction + 0.2 * changed_fraction**2
    if adjusted > 0.854:    # "basically infinite"
        adjusted = 0.854    # Safe value for the logarithm
    d = - math.log(1-adjusted)
    return d


def log_likelihood_at_t(model, N, t):
    '''
    Return the loglikelihood of observing a pair of sequence that yield a
    replacement count matrix N.

    Use that:
        log L(t) = log \prod_{i,j} p_{i,j}^{N_{i,j}}(t) = \sum N_{i,j} log(p_{i,j}(t))

    '''
    if t>0:
        P_t = model.get_replacement_probs(t)
        log_likelihood = np.sum(np.multiply(N, np.log(P_t)))
        return log_likelihood
    else:
        off_diagonal=np.sum(N) - np.sum(np.diag(N))
        if off_diagonal > 0:
            return 0.0
        else:
            return 1.0


def ml_distance_estimate(model, N, start_dist=None):
    '''
    Maximum likelihood estimation of evolutionary distance between the two
    sequences, and represented by the
    20x20 matrix N that contains the observed replacements (and
    non-replacements).

    You can give a starting value for the ML distance estimation with the
    start_dist parameter.
    '''
    return ml_distance_estimate_scipy(model, N, start_dist)


def ml_distance_estimate_scipy(model, N, start_dist=None):
    '''
    Same as ml_distance_estimate, but using the "bounded" method in SciPy (really Brent's,
    apparently), and Kimura's heuristic in a first step to determine suitable bounds.

    '''
    d_kimura = kimura_distance(N)
    if d_kimura < 0.1:
        bounds = (_delta, 0.2)
    elif d_kimura < 0.25:
        bounds = (_delta, 0.5)
    elif d_kimura < 0.5:
        bounds = (0.1, 0.75)
    elif d_kimura < 1.0:
        bounds = (0.25, 1.5)
    else:
        bounds = (1.0, _maxdistance)

    bracket = (0.75 * d_kimura, 1.25 * d_kimura)
    method='Bounded'
    ml_func = lambda x: - log_likelihood_at_t(model, N, x)
    res = scipy.optimize.minimize_scalar(ml_func, bounds=bounds, bracket=bracket, method=method, options={'xatol':_tolerance, 'maxiter':_maxiters})
    return res.x



def _experimental_ml_distance_estimate_scipy(model, N, start_dist=None):
    '''
    Same as ml_distance_estimate, but using the "bounded" method in SciPy (really Brent's,
    apparently), and Kimura's heuristic in a first step to determine suitable bounds.

    '''
    d_kimura = kimura_distance(N)
    if d_kimura < 0.1:
        bounds = (_delta, 0.2)
    elif d_kimura < 0.25:
        bounds = (_delta, 0.5)
    elif d_kimura < 0.5:
        bounds = (0.1, 0.75)
    elif d_kimura < 1.0:
        bounds = (0.25, 1.5)
    else:
        bounds = (1.0, _maxdistance)

    bracket = (0.75 * d_kimura, 1.25 * d_kimura)
    method='Bounded'
    ml_func = lambda x: - log_likelihood_at_t(model, N, x)
    res = scipy.optimize.minimize_scalar(ml_func, bounds=bounds, bracket=bracket, method=method, options={'xatol':_tolerance, 'maxiter':_maxiters})
    return res.x, res.nfev, res.nit # estimate, evaluations of function, and number of iterations


