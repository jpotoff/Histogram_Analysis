These programs are used for analysis of histogram data from grand canonical Monte Carlo simulations.

Patch.cpp is used to determine the weights when combining histograms through the method of Ferrenberg and Swendsen.  The original code used a simple successive substitution method, which was stable, but very slow.  This version uses successive substition combined with a DIIS solver and a revised order of operations, which provides approximately 20 times the performance of the original code I wrote 30 years ago as a graduate student. 

"Histograms" are generated from GOMC in as a list of discrete state points.  This is essentially "histogram-free" reweighting.  There is no loss in fidelity since we are not binning data into a histogram.  

Once the weights have been determined, one can run phase.cpp to determine the phase behavior or PVT behavior of the fluid. 
