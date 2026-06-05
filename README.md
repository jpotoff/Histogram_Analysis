These programs are used for analysis of histogram data from grand canonical Monte Carlo simulations.

Patch.cpp is used to determine the weights when combining histograms through the method of Ferrenberg and Swendsen.  
"Histograms" are generated from GOMC in as a list of discrete state points.  

Once the weights have been determined, one can run phase.cpp to determine the phase behavior or PVT behavior of the fluid. 
