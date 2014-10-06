/*
 * Function et_rot_mag::equilibrium_mag_plus
 *
 * Computes rotating equilibirum with a magnetic field with extended features
 * (see file et_rot_mag.h for documentation)
 *
 */

/*
 *   Copyright (c) 2012 Pablo Cerda, Michael Gabler
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

char et_rot_mag_equil_plus_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2013/11/25 14:03:55  j_novak
 * Commented some variables to avoid warnings
 *
 * Revision 1.1  2012/08/12 17:48:35  p_cerda
 * Magnetstar: New classes for magnetstar. Allowing for non-equatorial symmetry in Etoile et al. Adding B_phi in Et_rot_mag.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "et_rot_mag.h"
#include "param.h"
#include "unites.h"
#include "graphique.h"

void Et_rot_mag::equilibrium_mag_plus(
     const Itbl& icontrol, const Tbl& control, Tbl& diff, 
     const int initial_j, 
     const Tbl an_j, 
     Cmp (*f_j)(const Cmp&, const Tbl), 
     Cmp (*)(const Cmp& x, const Tbl),
     const Tbl bn_j, 
     Cmp (*g_j)(const Cmp&, const Tbl), 
     Cmp (*N_j)(const Cmp& x, const Tbl),
     const double relax_mag) {
			     
    // Fundamental constants and units
    // -------------------------------
  using namespace Unites_mag ;
    
    // Grid parameters
    // ---------------
    
    //  const Mg3d* mg = mp.get_mg() ; 
    //  int nz = mg->get_nzone() ;	    // total number of domains
    //  int nzm1 = nz - 1 ; 
    
    // The following is required to initialize mp_prev as a Map_et:
    //Map_et& mp_et = dynamic_cast<Map_et&>(mp) ; 
        
    
    // Parameters to control the iteration
    // -----------------------------------
    
    int mer_max = icontrol(0) ; 
    //    int mer_rot = icontrol(1) ;
    //    int mer_change_omega = icontrol(2) ; 
    //    int mer_fix_omega = icontrol(3) ; 
    //    int mer_mass = icontrol(4) ; 
    int mermax_poisson = icontrol(5) ; 
    //    int delta_mer_kep = icontrol(6) ; 
    //    int mer_mag = icontrol(7) ;
    //    int mer_change_mag = icontrol(8) ;
    //    int mer_fix_mag = icontrol(9) ;

        
    double precis = control(0) ; 
    //    double omega_ini = control(1) ; 
    //    double relax = control(2) ;
    //    double relax_prev = double(1) - relax ;  
    double relax_poisson = control(3) ; 
    //    double thres_adapt = control(4) ; 
    //    double precis_adapt = control(5) ; 
    //    double Q_ini = control(6) ;
    //    double a_j_ini = control (7) ;

    // Error indicators
    // ----------------
    
    diff.set_etat_qcq() ; 
    double& diff_A_phi = diff.set(0) ; 

    // Parameters for the function Map_et::adapt
    // -----------------------------------------
    
    int niter ; 
    int adapt_flag = 1 ;    //  1 = performs the full computation, 
			    //  0 = performs only the rescaling by 
			    //      the factor alpha_r

   

 					   
    // Parameters for the Maxwell equations
    // -------------------------------------

    double precis_poisson = 1.e-16 ;     

    Param par_poisson_At ; // For scalar At Poisson equation
    Cmp ssjm1_At(mp) ;
    ssjm1_At.set_etat_zero() ;
    par_poisson_At.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_At.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_At.add_double(precis_poisson, 1) ; // required precision
    par_poisson_At.add_int_mod(niter, 0) ;  //  number of iterations actually used 
    par_poisson_At.add_cmp_mod( ssjm1_At ) ; 

    Param par_poisson_Avect ;  // For vector Aphi Poisson equation

    Cmp ssjm1_khi_mag(ssjm1_khi) ;
    Tenseur ssjm1_w_mag(ssjm1_wshift) ;

    par_poisson_Avect.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson_Avect.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson_Avect.add_double(precis_poisson, 1) ; // required precision
    par_poisson_Avect.add_cmp_mod( ssjm1_khi_mag ) ; 
    par_poisson_Avect.add_tenseur_mod( ssjm1_w_mag ) ; 
    par_poisson_Avect.add_int_mod(niter, 0) ;   

				   
    // Initializations
    // ---------------

    // Initial  magnetic quantities
    a_j = 0;
    
    update_metric() ;	// update of the metric coefficients

    equation_of_state() ;	// update of the density, pressure, etc...
    
    hydro_euler() ;	// update of the hydro quantities relative to the 
			//  Eulerian observer

    MHD_comput() ; // update of EM contributions to stress-energy tensor


    // Output files

    ofstream fichmulti("multipoles.d") ;    
    
		    
    ofstream fichconv("convergence.d") ;    // Output file for diff_A_phi
    fichconv << "#     diff_A_phi     GRV2    " << endl ; 
    
    ofstream fichfreq("frequency.d") ;    // Output file for  omega
    fichfreq << "#       f [Hz]" << endl ; 
    
    ofstream fichevol("evolution.d") ;    // Output file for various quantities
    fichevol << 
    "#       |dH/dr_eq/dH/dr_pole|      r_pole/r_eq	ent_c" 
    << endl ; 
    
    diff_A_phi = 1 ; 
    //    double err_grv2 = 1 ; 
    

    A_phi = 0. ;	
    A_phi.std_base_scal() ;
    A_t = 0.;	
    A_t.std_base_scal() ;
    j_phi = 0.;	
    j_phi.std_base_scal() ;

    Cmp A_phi_old = A_phi;
    Cmp A_phi_new = A_phi;
    Cmp A_t_old = A_t;
    Cmp A_t_new = A_t;
    Cmp j_phi_old = j_phi;
    Cmp j_phi_new = j_phi;

    //=========================================================================
    // 			Start of iteration
    //=========================================================================

    for(int mer=0 ; (diff_A_phi > precis) && (mer<mer_max) ; mer++ ) {

 	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
   
	fichconv << mer ;
	fichfreq << mer ;
	fichevol << mer ;
	
	A_t_old = A_t;
	A_phi_old = A_phi;
	j_phi_old = j_phi;

	//-----------------------------------------------
	// Computation of electromagnetic potentials :
	// -------------------------------------------

	magnet_comput_plus(adapt_flag, initial_j, 
			   an_j, f_j, bn_j, g_j, N_j, par_poisson_At, par_poisson_Avect) ;
	A_t_new = A_t;
	A_phi_new = A_phi;
	j_phi_new = j_phi;
	
	A_t = relax_mag*A_t_new + (1.-relax_mag)*A_t_old ;
	A_phi = relax_mag*A_phi_new + (1. - relax_mag)*A_phi_old ;
	

	double diff_A_phi_n = max(abs((A_phi_new - A_phi_old))).set(0);
	double max_Aphi = max(abs(A_phi)).set(0);
	double diff_j_phi_n = max(abs((j_phi_new - j_phi_old))).set(0);
	double max_jphi = max(abs(j_phi)).set(0);

	Tbl maphi = A_phi_new.multipole_spectrum();
	//	int nzmax = maphi.get_dim (1) -1;

	if (max_Aphi == 0) {
	  diff_A_phi = 100.;
	}else{
	  diff_A_phi = diff_A_phi_n / max_Aphi ;
	}
	cout << mer << " "<< diff_A_phi << " " << max(A_phi).set(0) << " " << min(A_phi).set(0) << endl;
	cout << mer << " "<< diff_j_phi_n << " " << max_jphi << endl;
	

	fichmulti << diff_A_phi<< " " 
		  <<maphi.set(0,0) <<  " " <<maphi.set(0,1) <<  " " 
		  <<maphi.set(0,2) <<  " " <<maphi.set(0,3) <<  " " 
		  <<maphi.set(0,4) << endl;
	

	//	des_coupe_y(A_phi, 0., nzet, "Magnetic field") ; 
	//		des_coupe_y(j_phi, 0., nzet, "Current") ; 

	

	fichconv << endl ;
	fichfreq << endl ;
	fichevol << endl ;
	fichconv.flush() ; 
	fichfreq.flush() ; 
	fichevol.flush() ; 

    } // End of main loop
    
    //=========================================================================
    // 			End of iteration
    //=========================================================================

    fichconv.close() ; 
    fichfreq.close() ; 
    fichevol.close() ; 
    fichmulti.close();

}
