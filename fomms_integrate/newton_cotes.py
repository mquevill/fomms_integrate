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

def simpson3_8(x, f):
    """ 
    Compute a 1D definite integral using the 3/8 Simpson's rule.
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
    yb = f((2*a+  b)/3)
    yc = f((  a+2*b)/3)
    yd = f(b)
    I = (b-a) * (ya + 3 * (yb + yc) + yd) / 8
    return I

def boole(x, f):
    """ 
    Compute a 1D definite integral using the 3/8 Simpson's rule.
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
    yb = f((3*a+  b)/4)
    yc = f((  a+  b)/2)
    yd = f((  a+3*b)/4)
    ye = f(b)
    I = (b-a) * (7 * (ya + ye) + 32 * (yb + yd) + 12 * yc) * 2 / 45
    return I

