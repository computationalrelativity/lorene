# Parameters for the computation
###############################
4	nz: total number of domains
9	nt: number of points in theta (the same in each domain)
4	np: number of points in phi   (the same in each domain)
17	nr: number of points in r in the first domain
17	nr: number of points in r in all the other domains
1.	coordinate radius
0.3	relaxation lapse (1 -> no relaxation)
0.3	relaxation psi (1 -> no relaxation)
0.3	relaxation beta (1 -> no relaxation)
1e-8	threshold
20000	maximum of iterations
0.01	angular velocity
0. 0.	boost velocity in x direction / z direction
0 0.2	boundary condition for the lapse / value of the coefficient
1	boundary condition for psi
0	boundary condition for the shift
1	1 = solve for the lapse / 0 : not 
1	1 = solve for the psi / 0 : not 
1	1 = solve for the shift / 0 : not 

########################################################
For the lapse :
   0	boundary_nn_Dir(double)
   1	boundary_nn_Neu_eff(double)
   2	boundary_nn_Dir_eff(double)   
   3	boundary_nn_Neu_kk()   
   4	boundary_nn_Dir_kk()
   5    boundary_nn_Dir_lapl()
   6    boundary_nn_Neu_Cook()

For Psi :
   0	boundary_psi_app_hor()
   1	boundary_psi_Neu_spat()
   2	boundary_psi_Dir_spat()
   3	boundary_psi_Neu_evol()
   4	boundary_psi_Dir_evol()
   5	boundary_psi_Dir()

For the shift
   0	boundary_beta_x y,z (coordonnees horizon fixe)
   1	Berlin boundary condition
 

resu_lXX_P_sXX_F/NF_(bxXX/bzXX).d
