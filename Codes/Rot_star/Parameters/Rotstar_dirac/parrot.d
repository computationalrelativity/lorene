#################### PHYSICAL PARAMETERS ######################################
1     Relativity parameter: 1 = relativistic computation , 0 = Newtonian
0.228    ent_c : central enthalpy [c^2]
550.     freq_si : rotation frequency [Hz]
1.    fact_omega : 1.01 = search for the Keplerian frequency, 1. = otherwise.
1.60    Requested baryon mass [M_sol] (effective only if mer_mass > mer_max)
#################### COMPUTATIONAL PARAMETERS #################################
100      mer_max : maximum number of steps
1.e-9   precis : threshold on the enthalpy relative change for ending the computation
10      mer_rot : step at which the rotation is switched on
550.    freq_ini_si : initial rotation frequency [Hz] (switched on at mer = mer_rot)
10      mer_change_omega : step at which f is changed to reach freq_si
20	mer_fix_omega : step at which f must have reached 
1      delta_mer_kep : number of steps after mer_fix_omega to search for Kepler.
2000    mer_mass : step from which the baryon mass is forced to converge (if negative, variation of Omega)
0.5     aexp_mass : exponent for the increase factor of the central enthalpy
0.5     relax : relaxation factor in the main iteration 
1       graph : 1 = graphical outputs during the computation 
#################### MULTI-GRID PARAMETERS ###################################
4	nz : total number of domains
2	nzet : number of domains inside the star
2	nzadapt : number of domains of where the mapping adaptation will be done.
17	nt: number of points in theta (the same in each domain)
1	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
25	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)  	
33	0.5	<-   nr	  &   min(r)  in domain 1
25	1.	<-   nr   &   min(r)  in domain 2
17	2.	<-   nr   &   min(r)  in domain 3
0.08	enthalpy defining boundary between domains 0 and 1
2	order of spectral filtering (for rotstar_dirac). 0 means no filtering.
25	mer_hij : step at which the equation for h^{ij} is started being solved (put ~1000 for CFC)
