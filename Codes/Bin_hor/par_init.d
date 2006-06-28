## Parameters for init_bh 
#######################################################
11.	separation
6	nz : total number of domains
9	nt: number of points in theta (the same in each domain)
8	np: number of points in phi   (the same in each domain)
9       nr: number of points in r in the first domain
9       nr: number of points in r in all the other domains
0 1 2 4 8 16 boundaries of the domains
1e-6	Convergence treashold
0.5	Relaxation
0 0.3	boundary condition for the lapse / value of the coefficient
1	boundary condition for the psi

########################################################
For the lapse :
   0	boundary_nn_Dir(double)
   1	boundary_nn_Neu_eff(double)
   2	boundary_nn_Dir_eff(double)   
   3	boundary_nn_Neu_kk()   
   4	boundary_nn_Dir_kk()
   
For Psi :
   0	boundary_psi_app_hor()
   1	boundary_psi_Neu_spat()
   2	boundary_psi_Dir_spat()
   3	boundary_psi_Neu_evol()
   4	boundary_psi_Dir_evol()

    
