# Parameters for the binary equilibrium computation by coal
###################################################################
ini.d
0.8	fact_separ : Factor by which the initial separation is multiplied
1.55	mbar_voulue[0] : Baryon mass required for star 1 [M_sol]
1.55	mbar_voulue[1] : Baryon mass required for star 2 [M_sol]
###################################################################
400   	mermax : Maximum number of steps in the main iteration
0.5 	relax :	Relaxation factor in the main iteration
1	mermax_eqb : Maximum number of steps in Etoile_bin::equilibrium
1	prompt : 1 if no pause during the computation
0	graph : 1 if graphical outputs during the computation
1.E-6 	seuil : Threshold on the enthalpy relative change for ending the computation
2	fmer_stop : Step interval between pauses in the main iteration
5	fmer_save : Step interval between safeguards of the whole configuration
4	mermax_poisson : Maximum number of steps in Map_et::poisson
1.5	relax_poisson :  Relaxation factor in Map_et::poisson
4 	mermax_potvit : Maximum number of steps in Map_radial::poisson_compact
0.5	relax_potvit :  Relaxation factor in Map_radial::poisson_compact
3000	mer_masse : Step from which the baryon mass is forced to converge
0.5	aexp_masse : Exponent for the increase factor of the central enthalpy
8	fmer_udp_met : Step interval between metric updates
1	ind_rel_met : 1 if relaxation of the metric, 0 if not
0.65	relax_met : Relaxation factor of the metric (used only if ind_rel_met=1)
0.75	relax_omeg : Relaxation factor on Omega (orbital angular velocity)
0.7	fact_omeg_min : fact_omeg_min * omega = low bound in the omega search
1.3	fact_omeg_max : fact_omeg_max * omega = high bound in the omega search
0.	thres_adapt1 : threshold on dH/dr for the adaptation of the mapping in star 1
0.	thres_adapt2 : threshold on dH/dr for the adaptation of the mapping in star 2
0.6	reduce_shift : factor by which the initial analytical shift is reduced
