#################### PHYSICAL PARAMETERS ######################################
0.01     ampli_init_khi : initial amplitude of khi potential
0.     ampli_init_mu :  initial amplitude of mu potential
0.      ampli_tgam_dot : initial amplitude of u^{ij} = dtgam^{ij}/dt
#################### COMPUTATIONAL PARAMETERS #################################
0.001    pdt :   time step dt
10    nb_time_steps :  maximum number of time steps
2       niter_elliptic : number of iterations in the resolution of the ellip. eq.
0.8     relax_elliptic : relaxation factor for the elliptic equations
6       method_poisson_vect : method for solving vectorial Poisson equations
1.e-10  precis_init : precision in the resolution of initial data equations
1       nopause : 1 = no pause during output printing, 0 otherwise 
0       graph : 1 = graphical outputs during the computation, 0 = no graph
0       graph_init : 1 = graphical outputs for initial data
2       jmod_check_constraints : 1/frequency of constraint checking
5      jmod_save  : 1/frequency of saving monitoring quantities to file
0	verbose  : verbose output during evolution (1 = yes / 0 = no)
#################### MULTI-GRID PARAMETERS ###################################
1       symmetry_phi : 1 = symmetry phi --> phi + pi, 0 otherwise
6       nz : total number of domains
9       nt : number of collocation points in theta (the same in each domain)
8       np : number of collocation points in phi   (the same in each domain)
# Number of points and inner boundary of each domain:
17	0.      nr & min(r)  in domain 0  (nucleus)  	
17	1.      nr & min(r)  in domain 1
17	2.      nr & min(r)  in domain 2
17	4.      nr & min(r)  in domain 3
17	6.      nr & min(r)  in domain 4
9	8.      nr & min(r)  in domain 5
