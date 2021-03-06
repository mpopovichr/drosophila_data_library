import numpy as np
from scipy.optimize import fmin_bfgs

Nfeval= 1

def rosen(X):
    return (1. - X[0])**2 + 100.0*(X[1]-X[0]**2)**2 + (1. - X[1])**2 + 100.0*(X[2]-X[1]**2)**2

def callbackF(Xi):
    global Nfeval
    print '{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}'.format(Nfeval, Xi[0], Xi[1], Xi[2], rosen(Xi))
    Nfeval += 1

print  '{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s}'.format('Iter', ' X1', ' X2', ' X3', 'f(X)')
x0 = np.array([1.1, 1.1, 1.1], dtype=np.double)
[xopt, fopt, gopt, Bopt, func_calls, grad_calls, warnflg] = fmin_bfgs(rosen, x0, callback=callbackF, maxiter=2000, full_output=True, retall=False)

import numpy as np
import scipy
import pymorph
import mahotas
from scipy import ndimage
import cv2
import matplotlib.pyplot as plt

dna= mahotas.imread('/Users/mpopovic/Documents/Work/Projects/drosophila_wing_analysis/drosophila_data_library/test/dna.jpeg')