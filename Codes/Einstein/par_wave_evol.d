#################### PHYSICAL PARAMETERS ######################################
0.1     ampli_init_khi : initial amplitude of khi potential
0.1     ampli_init_mu :  initial amplitude of mu potential
0.      ampli_tgam_dot : initial amplitude of u^{ij} = dtgam^{ij}/dt
#################### COMPUTATIONAL PARAMETERS #################################
0.01    pdt :   time step dt
1000    nb_time_steps :  maximum number of time steps
4       niter_elliptic : number of iterations in the resolution of the ellip. eq.
0.5     relax_elliptic : relaxation factor for the elliptic equations
1       method_poisson_vect : method for solving vectorial Poisson equations
1.e-10  precis_init : precision in the resolution of initial data equations
0       nopause : 1 = no pause during output printing, 0 otherwise 
1       graph : 1 = graphical outputs during the computation, 0 = no graph
1       graph_init : 1 = graphical outputs for initial data
#################### MULTI-GRID PARAMETERS ###################################
1       symmetry_phi : 1 = symmetry phi --> phi + pi, 0 otherwise
4       nz : total number of domains
17      nr : number of collocation points in r (the same in each domain)
9       nt : number of collocation points in theta (the same in each domain)
12      np : number of collocation points in phi   (the same in each domain)
# Inner boundary of each domain:
0.      min(r)  in domain 0  (nucleus)  	
1.      min(r)  in domain 1
2.      min(r)  in domain 2
4.      min(r)  in domain 3
8.      min(r)  in domain 4
