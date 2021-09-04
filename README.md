# stability-receptivity-channel-flow
This repository contains some codes to generate the figures in the paper "analysis of fluid systems: stability,receptivity,sensitivity".

  Code 'resolvent_norm.m' plots Fig.12 (Couette flow) in the paper. 
  
  Code 'transient_energy.m' plots Fig.4 (Poiseuille flow) in the paper.
  
  Code 'H2_norm.m' plots contour of the H2 norm of the resolvent as functions of streamwise and spanwise wavenumbers for Poiseuille flow.
  
  Code 'eig_pseudo.m' plots Fig.6 (Couette flow) in the paper.

Four functions 'cheb.m', 'clenCurt.m', 'oss_operator.m' and 'eigenvalue_selection.m' are needed for the main functions. 

Variables (Re, kx, ky, type of flow, number of eigenvalues) can be modified at the start of each main function code. 
