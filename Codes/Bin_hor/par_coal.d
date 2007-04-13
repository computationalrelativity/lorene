# Parameters for the binary equilibrium computation by coal
###################################################################
static.d
0.03 	Initial omega
1 3.	1 if search_masses and value of mass_area
1e-5	Precision for the virial
0.5	Relaxation	
10	Number of steps to go from omega_init to ``real'' omega.
10	Number of iteration when at the real omega
0	boundary condition for the shift
0	rotation parameter alpha (1 if corotational, 0 if irrotational :
	omega_loc = alpha * omega_orb)
