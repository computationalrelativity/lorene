#################### PHYSICAL PARAMETERS ######################################
0.       freq_si : rotation frequency [Hz]
1.       fact_omega : 1.01 = search for the Keplerian frequency, 1. = otherwise.
0.1	 gravastar core density
#################### COMPUTATIONAL PARAMETERS #################################
100      mer_max : maximum number of steps
1.e-7    precis : threshold on the enthalpy relative change for ending the computation
10       mer_rot : step at which the rotation is switched on
0.       freq_ini_si : initial rotation frequency [Hz] (switched on at mer = mer_rot)
10       mer_change_omega : step at which f is changed to reach freq_si
20       mer_fix_omega : step at which f must have reached freq_si
1        delta_mer_kep : number of steps after mer_fix_omega to search for Kepler.
0.3      thres_adapt : threhold on (dH/dr_eq)/dH/dr_pole) for the mapping adaptation
0.5      relax : relaxation factor in the main iteration 
4        mermax_poisson : maximum number of steps in Map_et::poisson
1.5      relax_poisson :  relaxation factor in Map_et::poisson
1.e-14   precis_adapt : precision in Map_et::adapt
0        graph : 1 = graphical outputs during the computation 
#################### MULTI-GRID PARAMETERS ###################################
4	 nz : total number of domains
2	 nzet : number of domains inside the star
2	 nzadapt : number of domains where the mapping adaptation will be done
5	 nt: number of points in theta (the same in each domain)
1	 np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
17	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)  	
17	1.	<-   nr	  &   min(r)  in domain 1
17	2.	<-   nr   &   min(r)  in domain 2
17	3.	<-   nr   &   min(r)  in domain 3
################### ENTHALPY AT GRID BOUNDARIES ##############################
0.1	 enthalpy defining boundary between domains 0 and 1 inside gravastar
