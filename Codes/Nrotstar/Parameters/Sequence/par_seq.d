#################### PHYSICAL PARAMETERS ######################################
1    Relativity parameter: 1 = relativistic computation , 0 = Newtonian
0.3   entc_min : initial central enthalpy [c^2]
0.5   entc_max : initial central enthalpy [c^2]
300. 	freq_min_si : initial central rotation frequency [Hz]
300. 	freq_max_si : final central rotation frequency [Hz]
10	Number of configurations in the sequence
2.3     Requested baryon mass [M_sol] (effective only if mer_mass < mer_max)
#################### COMPUTATIONAL PARAMETERS #################################
200     mer_max : maximum number of steps
1.e-7  precis : threshold on the enthalpy relative change for ending the computation
0.1     thres_adapt : threhold on (dH/dr_eq)/dH/dr_pole) for the mapping adaptation
1000       mer_mass : step from which the baryon mass is forced to converge
0.5     aexp_mass : exponent for the increase factor of the central enthalpy or frequency
0.5     relax : relaxation factor in the main iteration
4       mermax_poisson : maximum number of steps in Map_et::poisson
1.5     relax_poisson :  relaxation factor in Map_et::poisson
1.e-14  precis_adapt : precision in Map_et::adapt
0       graph : 1 = graphical outputs during the computation
#################### MULTI-GRID PARAMETERS ###################################
3	nz : total number of domains
1	nzet : number of domains inside the star
1	nzadapt : number of domains of where the mapping adaptation will be done.
17	nt: number of points in theta (the same in each domain)
1	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
33	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)  	
17	1.	<-   nr	  &   min(r)  in domain 1
17	2.	<-   nr   &   min(r)  in domain 2
0.1	enthalpy defining boundary between domains 0 and 1
