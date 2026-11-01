These programs are used for analysis of histogram data from grand canonical Monte Carlo simulations.

Patch.cpp is used to determine the weights when combining histograms through the method of Ferrenberg and Swendsen.  The original code used a simple successive substitution method, which was stable, but very slow.  A later version combined successive substitution with a DIIS solver, but DIIS could fall into a limit cycle and never converge when the weights span many orders of magnitude. This version solves for ln(weight) by Newton-Raphson on the convex objective of Shirts and Chodera (J. Chem. Phys. 129, 124105 (2008)), whose minimum is the Ferrenberg-Swendsen solution, with a line search and a fixed-point fallback step far from the solution. Convergence requires both the fixed-point residual to fall below the tolerance and the Newton step in every ln(weight) to fall below 1e-6. phase.cpp uses the same solver to refine the weights read from input_fsp2.dat. 

"Histograms" are generated from GOMC in as a list of discrete state points.  This is essentially "histogram-free" reweighting.  There is no loss in fidelity since we are not binning data into a histogram.  

Once the weights have been determined, one can run phase.cpp to determine the phase behavior or PVT behavior of the fluid. 
