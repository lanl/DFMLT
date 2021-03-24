"""
© (or copyright) 2021. Triad National Security, LLC. All rights reserved.  This
program was produced under U.S. Government contract 89233218CNA000001 for Los
Alamos National Laboratory (LANL), which is operated by Triad National
Security, LLC for the U.S.  Department of Energy/National Nuclear Security
Administration. All rights in the program are reserved by Triad National
Security, LLC, and the U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting on its
behalf a nonexclusive, paid-up, irrevocable worldwide license in this material
to reproduce, prepare derivative works, distribute copies to the public,
perform publicly and display publicly, and to permit others to do so.
"""

from pylab import *
import sys
import argparse
from datetime import datetime
from scipy.integrate import quad
#import numpy.random

parser = argparse.ArgumentParser(
  description='Command line inputs for 1D hohlraum problem')
parser.add_argument('x', metavar='x', type=float,
  help='Position at which to evaluate distribution function')
parser.add_argument('t', metavar='t', type=float,
  help='Time at which to evaluate distribution function')
parser.add_argument('N', metavar='N', type=int,
  help='Number of samples to generate from distribution function')

args = parser.parse_args()

x = args.x
t = args.t
if x > t:
  print('ERROR x > t means I = 0!')
  sys.exit()
N = args.N
theta_max = arccos(x/t)

theta = zeros(N)
phi = zeros(N)
I = zeros(N)

def get_I(x, t, theta):
  if theta < theta_max:
    return 1.
  else:
    return 0.

# Randomly sample positions on sphere
for n in range(N):
  theta[n] = arccos(1. - 2.*random())
  phi[n] = 2.*pi*random()
  I[n] = get_I(x, t, theta[n])

def func(th, i, j):
  return

# Evaluate eddington tensor exactly
def get_n(i, theta, phi):
  if i == 0:
    return cos(theta)
  elif i == 1:
    return sin(theta)*cos(phi)
  elif i == 2:
    return sin(theta)*sin(phi)
  else:
    print("ERROR i != 0, 1, 2!")
    sys.exit()

def get_p():
  def int_p(theta):
    if theta < theta_max:
      return sin(theta)
    else:
      return 0
  return 2.*pi*quad(int_p, 0, pi)[0]
def get_pij(i, j):
  def int_pij_th(theta, phi, i, j):
    if theta < theta_max:
      return sin(theta)*get_n(i, theta, phi)*get_n(j, theta, phi)
    else:
      return 0
  def int_pij_phi(phi, i, j):
    return quad(int_pij_th, 0, pi, args=(phi,i,j))[0]
  return quad(int_pij_phi, 0, 2.*pi, args=(i, j))[0]

pij = zeros([3,3])
for i in range(3):
  for j in range(3):
    pij[i,j] = get_pij(i,j)/get_p()
    # Clean up
    if (fabs(pij[i,j]) < 1.e-7):
      pij[i,j] = 0.;

f = open("1d_hohlraum_%.2g_" % x + "%.2g" % t + "_%i" % N + ".txt", "w")
now = datetime.now()
f.write("# 1D hohlraum problem " + now.strftime("%Y/%m/%d %H:%M:%S") + "\n")
f.write("x %e" % x + "\n")
f.write("t %e" % t + "\n")
f.write("N %i" % N + "\n")
f.write("# theta phi I\n")
for n in range(N):
  f.write("%e " % theta[n] + "%e " % phi[n] + "%e " % I[n] + "\n")
f.write("# Exact (~1e-8) pressure tensor\n")
for i in range(3):
  f.write("%e " % pij[i,0] + "%e " % pij[i,1] + "%e " % pij[i,2] + "\n")
f.close()
