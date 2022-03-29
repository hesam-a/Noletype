# program for multi-variable optimization

import subprocess
import numpy as np
from scipy.optimize import minimize
from subMM import *

def MSE(A,B):
  A = np.array(A)
  B = np.array(B)
  return np.square(A-B).mean()

def wMSE(A,B,W):
  A = np.array(A)
  B = np.array(B)
  W = np.array(W)
  return (W*np.square(A-B)).mean()

def L1_reg(A,B,lambdaa,weight):
  A = np.array(A)
  B = np.array(B)
  return lambdaa*np.sum(weight*np.abs(A-B))

def L2_reg(A,B,lambdaa,weight):
  A = np.array(A)
  B = np.array(B)
  return lambdaa*np.sum(weight*np.square(A-B))

def costFUNC(params):
  QM = np.loadtxt('QM-energy.dat',usecols=(-1,), unpack=True)
  MM = getEnergy(params, 'filelist')
  weight = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0]) 
  refprm = np.loadtxt("p0.txt", usecols=(-1,)) 
  cost = 0.5*MSE(QM,MM) +L2_reg(params,refprm, 0.8, weight) 
  print(" func %15.8f "%(cost))
  return cost

def main():
  x0 = np.loadtxt("p0.txt",usecols=(-1,))
  myBounds = [[1, 2.5], [0.25, 0.39], [2.0, 4.5], [0.6, 0.99], [2.0, 4.5], [0.6, 0.99]]
  myOptions = {'eps': 1.30234e-4, 'maxiter': 150, 'ftol': 1e-6, 'disp': True, 'iprint': 1}
  ret = minimize(costFUNC, x0, method='SLSQP', jac=None, tol=1e-4, bounds=myBounds, options = myOptions)
  np.savetxt("p1.txt", ret.x,fmt='%15.10f')
  subprocess.run("paste p0.txt p1.txt >temp && mv temp p0.txt", shell=True)

if __name__ == "__main__":
    main()
