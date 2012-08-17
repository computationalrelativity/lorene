#################### PHYSICAL PARAMETERS ######################################
1     Relativity parameter: 1 = relativistic computation , 0 = Newtonian
0.1155520093   ent_c : central enthalpy [c^2]
0.     freq_si : rotation frequency [Hz]
1.    fact_omega : 1.01 = search for the Keplerian frequency, 1. = otherwise.
1.60    Requested baryon mass [M_sol] (effective only if mer_mass > mer_max)
#################### MAGNETIC PARAMETERS ######################################
1.0e-6  m_max deformation amplitude M_max * 6.673e-11 (see Gabler et al 2012)
1.	a0_j  deformation amplitude coefficient a_0 
0.      a1_j  deformation amplitude coefficient a_1
0.      a2_j  deformation amplitude coefficient a_2
0.      c0_j  deformation amplitude coefficient a_c
1.      rc0_j current shape coefficient c
0.      b0_j  force-free mixed toroidal-poloidal field coeff. b_0
0.      b1_j  force-free mixed toroidal-poloidal field coeff. b_1
0. 	b2_j  force-free mixed toroidal-poloidal field coeff. b_2
0.5     relax_mag
1	initial_j
#################### COMPUTATIONAL PARAMETERS #################################
2000      mer_max : maximum number of steps
1.e-6   precis : threshold on the enthalpy relative change for ending the computation
10      mer_rot : step at which the rotation is switched on
0.    freq_ini_si : initial rotation frequency [Hz] (switched on at mer = mer_rot)
10      mer_change_omega : step at which f is changed to reach freq_si
20      mer_fix_omega : step at which f must have reached freq_si
1      delta_mer_kep : number of steps after mer_fix_omega to search for Kepler.
0.2    thres_adapt : threhold on (dH/dr_eq)/dH/dr_pole) for the mapping adaptation
5000000    mer_mass : step from which the baryon mass is forced to converge (if negative, variation of Omega)
0.5     aexp_mass : exponent for the increase factor of the central enthalpy
0.5     relax : relaxation factor in the main iteration 
4       mermax_poisson : maximum number of steps in Map_et::poisson
1.5     relax_poisson :  relaxation factor in Map_et::poisson
1.e-15  precis_adapt : precision in Map_et::adapt
1       graph : 1 = graphical outputs during the computation 
#################### MULTI-GRID PARAMETERS ###################################
5	nz : total number of domains
1	nzet : number of domains inside the star
1	nzadapt : number of domains of where the mapping adaptation will be done.
17	nt: number of points in theta (the same in each domain)
1	np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
33	0.	<-   nr	  &   min(r)  in domain 0  (nucleus)  	
33	1.	<-   nr	  &   min(r)  in domain 1
33	2.	<-   nr   &   min(r)  in domain 2
33      4.
33      8.
0.1	enthalpy defining boundary between domains 0 and 1
