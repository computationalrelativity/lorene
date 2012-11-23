#################### PHYSICAL PARAMETERS ######################################
1		m : boson mass
0		k : coefficient in the azimuthal dependence of Phi
0.2   	rphi_c : central value of the real part of the scalar field
0.1		iphi_c : central value of the imaginary part of the scalar field
#################### COMPUTATIONAL PARAMETERS #################################
100      mer_max : maximum number of steps
1.e-7   precis : threshold on Phi relative change for ending the computation
0.3    	thres_adapt : threhold on ?? for the mapping adaptation
1800    mer_triax : step at which the 3-D perturbation is switched on
1.e-3   ampli_triax : relative amplitude of the 3-D perturbation
0.5     relax : relaxation factor in the main iteration 
4       mermax_poisson : maximum number of steps in Map_et::poisson
1.5     relax_poisson :  relaxation factor in Map_et::poisson
1.e-14  precis_adapt : precision in Map_et::adapt
1       graph : 1 = graphical outputs during the computation 
#################### MULTI-GRID PARAMETERS ###################################
3	nz : total number of domains
1	nzadapt : number of domains of where the mapping adaptation will be done.
5	nt: number of points in theta (the same in each domain)
1	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
9	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)  	
9	1.	<-   nr	  &   min(r)  in domain 1
9	2.	<-   nr   &   min(r)  in domain 2

