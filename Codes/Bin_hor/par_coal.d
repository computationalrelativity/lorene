# Parameters for the binary equilibrium computation by coal
###################################################################
static.d
0.015 	Initial omega
0 2.9	1 if search_masses and value of mass_area
1e-6	Precision for the virial
0.5	Relaxation	
10	Number of steps to go from omega_init to ``real'' omega.
10	Number of iteration when at the real omega
0	boundary condition for the shift
	0	 boundary_beta_cart()
	1	 vv_bound_cart()
