Collection of scripts and data for training machine learning algorithms to
reconstruct pressure tensors given sampled distribution functions.

Contains scripts for generating distribution functions, randomly sampling those
distribution functions uniformly in solid angle, and accurately integrating the
associated pressure tensors. Sampled intensities and integrated pressure tensors
are output in a uniform datafile format. Distribution functions are taken from
radiation transport test problems in the literature or randomly generated by
summing spherical harmonics. Current methods are:

 - 1d_hohlraum: The 1D hohlraum test problem from Ryan, B.R. & Dolence, J. C. 2020, ApJ,
     891, 118 is an analytic transport problem in which radiation is emitted from a
     hot boundary condition at x = 0 into an initially empty region x > 0 in plane-
     parallel symmetry. This problem contains a discontinuity in moment space due to
     the lack of radiation propagating in the -x direction. At late time, although the
     radiation is highly nonthermal (a half-spherical shell in momentum space), the
     pressure tensor is exactly the thermal value, diag(1/3).
 - spherical_harmonics: This script generates a random distribution of spherical harmonics
     up to order 3. The realization is chosen such that the distribution function can be
     substantially anisotropic but is (almost) guaranteed to be nonthermal. Note that
     only l=0 and l=2 modes contribute to the radiation pressure tensor; l=1 and l=3
     modes add structure to the distribution function that should be ignored by a
     method to evaluate the radiation pressure tensor.

Questions/comments: brryan@lanl.gov
                    jdolence@lanl.gov
