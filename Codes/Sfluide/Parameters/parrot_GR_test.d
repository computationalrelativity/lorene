#################### PHYSICAL PARAMETERS ######################################
1       Relativity parameter: 1 = relativistic computation , 0 = Newtonian
0.22    ent_c : central enthalpy [c^2]
0.227   ent2_c : second central enthalpy [c^2]
800.      freq_si : rotation frequency [Hz]
-700.    freq2_si : second rotation frequency [Hz]
#################### COMPUTATIONAL PARAMETERS #################################
200     mer_max : maximum number of steps
1.e-7  precis : threshold on the enthalpy relative change for ending the computation
10      mer_rot : step at which the rotation is switched on
200.      freq_ini_si : initial rotation frequency [Hz] (switched on at mer = mer_rot)
0.      freq2ini_si : initial second rotation frequency [Hz] (switched on at mer = mer_rot)
10      mer_change_omega : step at which f is changed to reach freq_si
20      mer_fix_omega : step at which f must have reached freq_si
0.5     relax : relaxation factor in the main iteration 
6       mermax_poisson : maximum number of steps in Map_et::poisson
1.5     relax_poisson :  relaxation factor in Map_et::poisson
1       graph : 1 = graphical outputs during the computation 
#################### MULTI-GRID PARAMETERS ###################################
3       nz : total number of domains
1       nzet : number of domains inside the star
33       nt: number of points in theta (the same in each domain)
1       np: number of points in phi   (the same in each domain)
# Number of points in r and (initial) inner boundary of each domain:
65      0.      <-   nr   &   min(r)  in domain 0  (nucleus)    
17      1.      <-   nr   &   min(r)  in domain 1
9       2.      <-   nr   &   min(r)  in domain 2
0.2     enthalpy defining boundary between domains 0 and 1
