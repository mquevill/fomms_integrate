"""
This function implements 1d  Monte Carlo integration
"""

import numpy as np

def monte_1d(x, f, trials):
   """
    Compute a 1D definite integral
    Parameters
    ----------
    f : function
        User defined function.
    x : numpy array
        Integration domain.
    trials : integer
        Total number of generated random samples.
    Returns
    -------
    I : float
        Integration result.
   """
   a = x[0]
   b = x[1]
   d = (b - a) *  np.random.rand(1, trials) + a
   y = f(d)
   return (b-a) * np.sum(y) / trials

def monte_2d(f, v, domain, trials):
    """
    Compute area of circle using Monte Carlo sampling
    Parameters
    ----------
    f : function
        User defined function.
    domain : numpy array
        Integration domain.
    v : function
        Integration volume .
    Returns
    -------
    I : float
        Integration result.
    """

    area = np.prod(domain[1] - domain[0])
    matrix = np.sqrt(area) * np.random.rand(trials, 2) + np.min(domain)
    v_eval = v(matrix)
    idx = np.where(v_eval < 1.0)
    points_inside = v_eval[idx]
    k = points_inside.size
    return area * k / trials
