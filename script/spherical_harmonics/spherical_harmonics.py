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
from scipy.integrate import quad, dblquad
from scipy.special import lpmn
from math import factorial
#import numpy.random
import numpy as np

parser = argparse.ArgumentParser(
  description='Command line inputs for 1D hohlraum problem')
parser.add_argument('order', metavar='order', type=int,
  help='Order of spherical harmonics to include')
parser.add_argument('N', metavar='Nsamples', type=int,
  help='Number of samples to generate from distribution function')
parser.add_argument('--amp', dest='amp', type=float,
  help='Amplitude of perturbative (inhomogeneous) component')
parser.add_argument('--plot', dest='plot', action='store_true')
parser.set_defaults(plot=False)

args = parser.parse_args()

plot_f = args.plot
order = args.order
if order < 0:
  print("ERROR order < 0!")
  sys.exit()
lmax = args.order
mmax = lmax
nplotpoint = 50
Nsamples = args.N
user_amp = 1.
if args.amp is not None:
  user_amp = args.amp

# Real spherical harmonics (scipy bug in sph_harm?)
def get_ylm(l, m, theta, phi):
  if l == 0:
    return 1./sqrt(2.*pi)
  elif l == 1:
    prefac = 1./2.*sqrt(3./pi)
    if m == -1:
      return prefac*sin(phi)*sin(theta)
    elif m == 0:
      return prefac*cos(theta)
    elif m == 1:
      return prefac*cos(phi)*sin(theta)
  elif l == 2:
    prefac = 1./4.*sqrt(15./pi)
    if m == -2:
      return prefac*sin(2*phi)*sin(theta)**2
    elif m == -1:
      return prefac*sin(phi)*sin(2.*theta)
    elif m == 0:
      return 1./8.*sqrt(5./pi)*(1. + 3.*cos(2.*theta))
    elif m == 1:
      return prefac*cos(phi)*sin(2.*theta)
    elif m == 2:
      return prefac*cos(2.*phi)*sin(theta)**2
  elif l == 3:
    if m == -3:
      return 1./4.*sqrt(35./(2.*pi))*sin(3*phi)*sin(theta)**3
    elif m == -2:
      return 1./4.*sqrt(105./pi)*cos(theta)*sin(2.*phi)*sin(theta)**2
    elif m == -1:
      return 1./8.*sqrt(21./(2.*pi))*(3. + 5.*cos(2.*theta))*sin(phi)*sin(theta)
    elif m == 0:
      return 1./16*sqrt(7./pi)*(3.*cos(theta) + 5.*cos(3.*theta))
    elif m == 1:
      return 1./8.*sqrt(21./(2.*pi))*(3. + 5.*cos(2.*theta))*cos(phi)*sin(theta)
    elif m == 2:
      return 1./4.*sqrt(105./pi)*cos(theta)*cos(2.*phi)*sin(theta)**2
    elif m == 3:
      return 1./4.*sqrt(35./(2.*pi))*cos(3*phi)*sin(theta)**3

  print("ERROR l = %i" % l + " m = %i" % m + " not supported!")
  sys.exit()

from scipy.special import sph_harm
def get_I(theta, phi, amp):
  ans = 0.
  for l in range(lmax+1):
    for m in range(-l, l+1):
      ans += amp[l,m+l]*get_ylm(l, m, theta, phi)
  return ans

# Realize a (probably) non-zero everywhere distribution function.
# Start with [l,m]=[0,0] max = 1. Then sample an amplitude for the inhomogeneous
# (l >= 1) part between 0 and 0.5. Then sample another amplitude factor for each
# l >= 1 harmonic between -1 and 1. This
sampled = False
while sampled == False:
  amp = zeros([lmax+1, 2*(lmax+1)+1])
  pamp = user_amp*random()/(lmax**2/2 + 1)
  amp[0,0] = 1.
  for l in range(1, lmax+1):
    for m in range(-l, l+1):
      amp[l,m+l] = pamp*2.*(random()-0.5)
  theta_p = linspace(0, pi, nplotpoint)
  phi_p = linspace(0, 2*pi, nplotpoint)
  theta_p, phi_p = meshgrid(theta_p, phi_p)

  Icolors = get_I(theta_p, phi_p, amp)
  Imax, Imin = Icolors.max(), Icolors.min()
  if Imin > 0:
    sampled = True
  else:
    print("Imin < 0! Resampling!")

if plot_f:
  Icolors = Icolors/Imax
  import matplotlib.pyplot as plt
  from matplotlib import cm, colors
  from mpl_toolkits.mplot3d import Axes3D
  x = sin(theta_p) * cos(phi_p)
  y = sin(theta_p) * sin(phi_p)
  z = cos(theta_p)
  fig = plt.figure(figsize=plt.figaspect(1.))
  ax = fig.add_subplot(1,1,1,projection='3d')
  ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.afmhot(Icolors))
  plt.show()

# Randomly sample positions on sphere
theta = zeros(Nsamples)
phi = zeros(Nsamples)
I = zeros(Nsamples)
for n in range(Nsamples):
  theta[n] = arccos(1. - 2.*random())
  phi[n] = 2.*pi*random()
  I[n] = get_I(theta[n], phi[n], amp)

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
  def int_p_th(theta, phi):
    return (sin(theta)*
           get_I(theta, phi, amp))
  def int_p_phi(phi):
    return quad(int_p_th, 0, pi, args=(phi))[0]
  return quad(int_p_phi, 0, 2.*pi)[0]
def get_pij(i, j):
  def int_pij_th(theta, phi, i, j):
    return (sin(theta)*get_n(i, theta, phi)*get_n(j, theta, phi)*
           get_I(theta, phi, amp))
  def int_pij_phi(phi, i, j):
    return quad(int_pij_th, 0, pi, args=(phi,i,j))[0]
  return quad(int_pij_phi, 0, 2.*pi, args=(i, j))[0]

pij = zeros([3,3])
p = get_p();
for i in range(3):
  for j in range(3):
    pij[i,j] = get_pij(i,j)#/p
    if fabs(pij[i,j]) < 1.e-14:
      pij[i,j] = 0.
    if p > 1.e-14:
      pij[i,j] /= p

print("Pressure tensor:")
print(pij)

print("sh_%i" % order + "_%i" % Nsamples + ".txt")
f = open("sh_%i" % order + "_%i" % Nsamples + ".txt", "w")
now = datetime.now()
f.write("# Spherical harmonics problem " + now.strftime("%Y/%m/%d %H:%M:%S") + "\n")
f.write("order %i" % order + "\n")
f.write("N %i" % Nsamples + "\n")
f.write("# theta phi I\n")
for n in range(Nsamples):
  f.write("%e " % theta[n] + "%e " % phi[n] + "%e " % I[n] + "\n")
f.write("# Exact (~1e-8) pressure tensor\n")
for i in range(3):
  f.write("%e " % pij[i,0] + "%e " % pij[i,1] + "%e " % pij[i,2] + "\n")
f.close()
