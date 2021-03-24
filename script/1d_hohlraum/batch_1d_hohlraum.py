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

xs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
ts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
Ns = [8, 16, 32, 64, 128, 256, 512]

for ix in range(len(xs)):
  for it in range(len(ts)):
    for iN in range(len(Ns)):
      x = xs[ix]
      t = ts[it]
      N = Ns[iN]
      if x < t:
        print(x, t, N)
        name = "1d_hohlraum_%.2g_" % x + "%.2g" % t + "_%i" % N + ".txt"
        call(["python", "1d_hohlraum.py", str(x), str(t), str(N)])
        os.rename(name, "../../data/1d_hohlraum/" + name)
