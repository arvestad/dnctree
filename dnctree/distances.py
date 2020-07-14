from modelmatcher.models import RateMatrix
import math
import numpy as np
import scipy.optimize
import sys

_tolerance = 0.001
_delta = 0.0001
_maxiters = 50
_maxdistance = 10


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


def ml_distance_estimate_arve(model, N, start_dist=None):
    '''
    Maximum likelihood estimation of evolutionary distance between the two
    sequences, and represented by the
    20x20 matrix N that contains the observed replacements (and
    non-replacements).

    We apply Newton-Raphson estimation to find the evolutionary distance
    where the derivative of the likelihood of N is zero (or close enough).
    Note that Newton-Raphson works with derivatives of functions, so there are
    two derivatives mentioned here: the likelihood derivative and double derivative.
    '''
    if start_dist is None:
        d = kimura_distance(N)
    else:
        d = start_dist
    delta = _delta           # The little step we take when estimating derivative

    # These important contants should of course be stored more prominently, or
    # even be programmatically changable.

    if d < _tolerance:             # Too close to zero. Not good for computational reasons, bump it up:
        d = _tolerance

    L_derivative = likelihood_derivative_at_t(model, N, d)
    for iteration in range(_maxiters): # Limit number of iterations
        print(f'Iter {iteration}: d = {d}, L_der = {L_derivative}')
        if abs(L_derivative) < _tolerance: # Derivative is basically zero.
            print(f'Returning bc L_derivative < _tolerance: {L_derivative} < {RateMatrix._tolerance}')
            return d                   # This is the _normal_ exit point!
        new_L_derivative = likelihood_derivative_at_t(model, N, d + delta)
        deriv = (new_L_derivative - L_derivative) / delta
        d = d - L_derivative / deriv # Newton-Raphson update

        # Three possible termination conditions
        if d < _tolerance:
            print(f'Returning bc d < _tolerance: {d} < {RateMatrix._tolerance}')
            return d
        if d > _maxdistance:    # Basically infinity! Don't go further
            print(f'Returning bc d > _maxdistance: {d} > {_maxdistance}')
            return 5
        if abs(L_derivative) < abs(new_L_derivative): # Bad sign, this derivative should not increase. It might be an effect of working with small numbers, so jsut return.
            print(f'Returning bc abs(L_derivative) < abs(new_L_deriv): {abs(L_derivative)} < {abs(new_L_derivative)}')
            return d

        # Prepare for next iteration
        L_derivative = new_L_derivative

    # After max allowed iterations, we return
    return d


def ml_distance_estimate_scipy(model, N, start_dist=None):
    '''
    Same as ml_distance_estimate, but using the "bounded" method in SciPy (really Brent's, apparently).
    '''
    bounds = (_delta, _maxdistance)
    ml_func = lambda x: - log_likelihood_at_t(model, N, x)
    res = scipy.optimize.minimize_scalar(ml_func, bounds=bounds, method='Bounded', options={'xatol':_tolerance, 'maxiter':_maxiters})
    return res.x

def ml_distance_estimate_secant(model, N, start_dist=None):
    '''
    Same as ml_distance_estimate, but using the secant method in SciPy.
    '''
    if start_dist is None:
        d = kimura_distance(N)
    else:
        d = start_dist

    d = max(d, _tolerance)

    ml_func = lambda x: likelihood_derivative_at_t(model, N, x)
    res = scipy.optimize.newton(ml_func, d, tol=_tolerance, maxiter=_maxiters)
    #secant method is used because we do not provide a function for the derivative
    return res


def likelihood_derivative_at_t(model, N, t):
    '''
    Derivative of the loglikelihood. Helper for ml_distance_estimate.

    log L(t) = log \prod_{i,j} p_{i,j}^{N_{i,j}}(t) = \sum N_{i,j} log(p_{i,j}(t))
    d/dt log L(t) = \sum_{i,j} N_{i,j} \cdot \frac{1}{p_{i,j}(t)} \cdot p_{i,j}'t()
    Which means that it becomes
       log L(t) = \sum N .* P'(t) ./ P(t)
    where
       P(t) = e^{Qt}
       P'(t) = Q e^{Qt}

    Use the eigen decomposition to exponentiate Q.
    '''
    P_t = model.get_replacement_probs(t)

    # Compute  QP(t)
    direction = np.matmul(P_t, model.Q)

    likelihood_derivative = np.sum(np.divide(np.multiply(N, direction), P_t))
    return likelihood_derivative


def __poisson_distance(t1, t2, s1, s2, supress_warnings=False):
    '''
    Older and worse alternative than kimura_distance().
    Currently not up-to-date! This function will not work in the new framework
    '''
    n_chars = 0
    n_diffs = 0
    for c1, c2 in zip(s1, s2):
        if c1=='-' or c2=='-':
            continue
        if c1 != c2:
            n_diffs += 1
        n_chars += 1

    if n_chars == 0:
        if not supress_warnings:
            print(f'Warning: No shared characters between {t1} and {t2}, so cannot estimate their distance. Defaulting to distance=2.5.', file=sys.stderr)
        return 2.5
    if n_diffs == n_chars:
        if not supress_warnings:
            print(f'Warning: Degenerate sequence pair: {t1} and {t2}. All shared characters are different. Defaulting to distance=2.5.', file=sys.stderr)
        return 2.5

    distance = - math.log(1 - n_diffs/n_chars)
    return distance
