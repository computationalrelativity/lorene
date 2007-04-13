## Parameters for init_bh 
#######################################################
12	separation
6	nz : total number of domains
13	nt: number of points in theta (the same in each domain)
12	np: number of points in phi   (the same in each domain)
13       nr: number of points in r in the first domain
13       nr: number of points in r in all the other domains
0 1 2 4 8 16 boundaries of the domains
1e-6	Convergence treashold
0.5	Relaxation
0 0.3	boundary condition for the lapse / value of the coefficient
0	boundary condition for the psi

########################################################
For the lapse :
   0	boundary_nn_Dir(double)
   1	boundary_nn_Neu(double)
For Psi :
   0	boundary_psi_app_hor()

    
