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
import os
from subprocess import call

Os = [0, 1, 2, 3]
Ns = [8, 16, 32, 64, 128, 256, 512]
Nrealizations = 2

for iO in range(len(Os)):
  for iN in range(len(Ns)):
    for n in range(Nrealizations):
      prename = "sh_%i" % Os[iO] + "_%i" % Ns[iN] + ".txt"
      name = "sh_%i" % Os[iO] + "_%i" % Ns[iN] + "_%i" % n + ".txt"
      print("Creating " + name + "... ", end="")
      call(["python", "spherical_harmonics.py", str(Os[iO]), str(Ns[iN])])
      os.rename(prename, "../../data/spherical_harmonics/" + name)
      print("done")
