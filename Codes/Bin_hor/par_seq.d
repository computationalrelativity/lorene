## Parameters for coal_seq
#######################################################
17       nt: number of points in theta (the same in each domain)
40       np: number of points in phi   (the same in each domain)
25       nr: number of points in r in the first domain
25       nr: number of points in r in all the other domains
1e-6    Convergence treashold
0.5     Relaxation
1 0.3   boundary condition for the lapse / value of the coefficient
1       boundary condition for the psi
0.03    Initial omega
1 2.6   1 if search_masses and value of irreductible mass
1e-6    Precision for the virial
10      Number of steps to go from omega_init to ``real'' omega.
10      Number of iteration when at the real omega
0       boundary condition for the shift
4	number of configurations to compute
8	closest congiguration (value of beta)
11	more distant configuration (value of beta)

########################################################
For the lapse :
   0    boundary_nn_Dir(double)
   1    boundary_nn_Neu_eff(double)
   2    boundary_nn_Dir_eff(double)
   3    boundary_nn_Neu_kk()
   4    boundary_nn_Dir_kk()

For Psi :
   0    boundary_psi_app_hor()
   1    boundary_psi_Neu_spat()
   2    boundary_psi_Dir_spat()
   3    boundary_psi_Neu_evol()
   4    boundary_psi_Dir_evol()

For the shift :
   0    boundary_beta_cart()
   1    vv_bound_cart()
