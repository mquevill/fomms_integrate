"""
This file contains the implementation of the Newton-Cotes rules
"""

import numpy as np

def trapz(x, f):
    """ 
    Compute a 1D definite integral using the trapezoidal rule
    Parameters
    ----------
    f : function
        User defined function.
    x : numpy array
        Integration domain.
    Returns
    -------
    I : float
        Integration result.
    """   
    a = x[0]
    b = x[1]
    ya = f(a)
    yb = f(b)
    I = (b-a) * (ya + yb) / 2
    return I

def simpson(x, f):
    """ 
    Compute a 1D definite integral using the Simpson's rule.
    Parameters
    ----------
    f : function
        User defined function.
    x : numpy array
        Integration domain.
    Returns
    -------
    I : float
        Integration result.
    """   
    a = x[0]
    b = x[1]
    ya = f(a)
    yb = f((a+b)/2)
    yc = f(b)
    I = (b-a) * (ya + 4 * yb + yc) / 6
    return I
