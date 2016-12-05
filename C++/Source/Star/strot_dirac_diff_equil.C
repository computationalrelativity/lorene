/*
 *  Function Star_rot_Dirac_diff::equilibrium
 *
 *    (see file star_rot_dirac_diff.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Motoyuki Saijo
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

 

/*
 * $Header$
 *
 */


// C headers
#include <cmath>
#include <cassert>

// Lorene headers
#include "star_rot_dirac_diff.h"
#include "param.h"
#include "utilitaires.h"
#include "unites.h" 

namespace Lorene {
void Star_rot_Dirac_diff::equilibrium(double ent_c, double omega_c0, 
	         double fact_omega, int , const Tbl& ent_limit,
		 const Itbl& icontrol, const Tbl& control,
		 double, double, Tbl& diff){	
	

  // Fundamental constants and units
  // --------------------------------
  using namespace Unites ;

  // For the display 
  // ---------------
  char display_bold[]="x[1m" ; display_bold[0] = 27 ;
  char display_normal[] = "x[0m" ; display_normal[0] = 27 ;


  // Grid parameters
  // ----------------

  const Mg3d* mg = mp.get_mg() ;
  int nz = mg->get_nzone() ;    // total number of domains
  int nzm1 = nz - 1 ;

  // Index of the point at phi=0, theta=pi/2 at the surface of the star:
  assert(mg->get_type_t() == SYM) ; 
  int l_b = nzet - 1 ; 
  int i_b = mg->get_nr(l_b) - 1 ; 
  int j_b = mg->get_nt(l_b) - 1 ; 
  int k_b = 0 ;

  // Value of the enthalpy defining the surface of the star
  double ent_b = ent_limit(nzet-1) ;

  // Parameters to control the iteration
  // -----------------------------------

  int mer_max = icontrol(0) ;
  int mer_rot = icontrol(1) ;
  int mer_change_omega = icontrol(2) ; 
  int mer_fix_omega = icontrol(3) ; 
  //int mer_mass = icontrol(4) ; 
  int delta_mer_kep = icontrol(5) ; 
  
  // Protections:
  if (mer_change_omega < mer_rot) {
    cout << "Star_rot_Dirac_diff::equilibrium: mer_change_omega < mer_rot !" 
	 << '\n' ;
    cout << " mer_change_omega = " << mer_change_omega << '\n' ; 
    cout << " mer_rot = " << mer_rot << '\n' ; 
    abort() ; 
  }
  if (mer_fix_omega < mer_change_omega) {
    cout << "Star_rot_Dirac_diff::equilibrium: mer_fix_omega < mer_change_omega !" 
	 << '\n' ;
    cout << " mer_fix_omega = " << mer_fix_omega << '\n' ; 
    cout << " mer_change_omega = " << mer_change_omega << '\n' ; 
    abort() ; 
  }

  double precis = control(0) ; 
  double omega_ini = control(1) ; 
  double relax = control(2) ;
  double relax_prev = double(1) - relax ;  
      
  // Error indicators
  // ----------------

  diff.set_etat_qcq() ;
  double& diff_ent = diff.set(0) ; 

  double alpha_r = 1 ;

  // Initializations
  // ---------------

  // Initial angular velocities 
  double omega_c = 0 ;
  
  double accrois_omega = (omega_c0 - omega_ini) / 
                double(mer_fix_omega - mer_change_omega) ; 
 

  update_metric() ;    //update of the metric quantities 

  equation_of_state() ;  // update of the density, pressure,...etc

  hydro_euler() ; //update of the hydro quantities relative to the 
                  //  Eulerian observer

  // Quantities at the previous step :
  Scalar ent_prev = ent ;
  Scalar logn_prev = logn ;
  Scalar qqq_prev = qqq ;
  //  Vector beta_prev = beta ;
  // Sym_tensor_trans hh_prev = hh ;

  // Output files
  // -------------

  ofstream fichconv("convergence.d") ;    // Output file for diff_ent
  fichconv << "#     diff_ent    GRV2    max_triax   vit_triax" << '\n' ;   

  ofstream fichfreq("frequency.d") ;    // Output file for  omega
  fichfreq << "#       f [Hz]" << '\n' ; 

  ofstream fichevol("evolution.d") ;    // Output file for various quantities
  fichevol << 
    "#       |dH/dr_eq/dH/dr_pole|      r_pole/r_eq ent_c" 
	   << '\n' ; 

  diff_ent = 1 ;
  double err_grv2 = 1 ;



//=========================================================================
//                   Start of iteration
//=========================================================================

  for(int mer=0 ; (diff_ent > precis) && (mer<mer_max) ; mer++) {

 cout << "-----------------------------------------------" << '\n' ;
 cout << "step: " << mer << '\n' ;
 cout << "ent_c = " << display_bold << ent_c << display_normal
      << '\n' ;    
 cout << "diff_ent = " << display_bold << diff_ent << display_normal
      << '\n' ;    
 cout << "err_grv2 = " << err_grv2 << '\n' ;    
 fichconv << mer ;
 fichfreq << mer ;
 fichevol << mer ;


 // switch on rotation 
 if (mer >= mer_rot) {
          
   if (mer < mer_change_omega) {
     omega_c = omega_ini ; 
   }
   else {
     if (mer <= mer_fix_omega) {
       omega_c = omega_ini + accrois_omega * 
	 (mer - mer_change_omega) ;
     }
   }


 }


 //---------------------------------------------------//
 // Resolution of the Poisson equation for logn       //
 // Note: ln_f is due to the fluid part               //
 //       ln_q is due to the quadratic metric part    //
 //---------------------------------------------------//

 Scalar ln_f_new(mp) ;
 Scalar ln_q_new(mp) ;

 solve_logn_f( ln_f_new ) ;
 solve_logn_q( ln_q_new ) ;

 ln_f_new.std_spectral_base() ;
 ln_q_new.std_spectral_base() ;


 //--------------------------------------------------//
 // Resolution of the Poisson equation for shift     //
 //--------------------------------------------------//

 Vector beta_new(mp, CON, mp.get_bvect_spher()) ;

 solve_shift( beta_new ) ;

 //------------------------------------
 // Determination of the fluid velocity
 //------------------------------------

 if (mer > mer_fix_omega + delta_mer_kep) {
              
   omega_c *= fact_omega ;  // Increase of the angular velocity if 
 }              //  fact_omega != 1
 
 bool omega_trop_grand = false ;
 bool kepler = true ;

 while ( kepler ) {

   // Possible decrease of Omega to ensure a velocity < c

   bool superlum = true ;

   while ( superlum ){
   
     // Computation of Omega(r,theta)

	if (omega_c == 0.) {
	    omega_field = 0 ; 
	}
	else {
	    par_frot.set(0) = omega_c ; 
	    if (par_frot(2) != double(0)) {  // fixed a = R_eq / R_0
		par_frot.set(1) = ray_eq() / par_frot(2) ; 
	    }
	    double omeg_min = 0 ; 
	    double omeg_max = omega_c ; 
	    double precis1 = 1.e-14 ;
	    int nitermax1 = 200 ; 

	    fait_omega_field(omeg_min, omeg_max, precis1, nitermax1) ;
	}

     // New fluid velocity :
     //

     u_euler.set(1).set_etat_zero() ;
     u_euler.set(2).set_etat_zero() ;

     u_euler.set(3) = omega_field ;
     u_euler.set(3).annule(nzet,nzm1) ;   // nzet is defined in class Star
     u_euler.set(3).std_spectral_base() ;
     u_euler.set(3).mult_rsint() ;
     u_euler.set(3) += beta(3) ;
     u_euler.set(3).annule(nzet,nzm1) ;
     
     u_euler = u_euler / nn ; 


     // v2 (square of norm of u_euler)
     // -------------------------------

     v2 = contract(contract(gamma.cov(), 0, u_euler, 0), 0, u_euler, 0) ;

     // Is the new velocity larger than c in the equatorial plane ?

     superlum = false ;
     
     for (int l=0; l<nzet; l++) {
       for (int i=0; i<mg->get_nr(l); i++) {
	  
	 double u2 = v2.val_grid_point(l, 0, j_b, i) ; 
	 if (u2 >= 1.) {     // superluminal velocity
	   superlum = true ; 
	   cout << "U > c  for l, i : " << l << "  " << i 
		<< "   U = " << sqrt( u2 ) << '\n' ;  
	 }
       }
     }
     if ( superlum ) {
       cout << "**** VELOCITY OF LIGHT REACHED ****" << '\n' ; 
       omega /= fact_omega ;    // Decrease of Omega
       cout << "New rotation frequency : " 
	    << omega/(2.*M_PI) * f_unit <<  " Hz" << '\n' ; 
       omega_trop_grand = true ;  
     }
   }    // end of while ( superlum )


   // New computation of U (this time is not superluminal)
   // as well as of gam_euler, ener_euler,...etc

   hydro_euler() ;



 //--------------------------------//
 // First integral of motion       //
 //--------------------------------//

 Scalar mlngamma(mp) ;  //  -log( gam_euler )

 mlngamma = - log( gam_euler ) ;

 // Equatorial values of various potentials :
 double ln_f_b = ln_f_new.val_grid_point(l_b, k_b, j_b, i_b) ; 
 double ln_q_b = ln_q_new.val_grid_point(l_b, k_b, j_b, i_b) ;
 double mlngamma_b = mlngamma.val_grid_point(l_b, k_b, j_b, i_b) ;
 double primf_b = prim_field.val_grid_point(l_b, k_b, j_b, i_b) ;

 // Central values of various potentials :
 double ln_f_c = ln_f_new.val_grid_point(0,0,0,0) ;
 double ln_q_c = ln_q_new.val_grid_point(0,0,0,0) ;
 double mlngamma_c = 0 ;
 double primf_c = prim_field.val_grid_point(0,0,0,0) ; 

 // Scale factor to ensure that the (log of) enthalpy is equal to 
 // ent_b at the equator
 double alpha_r2 = ( ent_c - ent_b + mlngamma_c - mlngamma_b
		     + ln_q_c - ln_q_b + primf_c - primf_b) 
                     / ( ln_f_b - ln_f_c  ) ;

 alpha_r = sqrt(alpha_r2) ;

 cout << "alpha_r = " << alpha_r << '\n' ; 

 // Rescaling of the grid (no adaptation!)
 //---------------------------------------
 mp.homothetie(alpha_r) ;

 // Readjustment of logn :
 // -----------------------

 logn = alpha_r2 * ln_f_new + ln_q_new ;

 double logn_c = logn.val_grid_point(0,0,0,0) ;

 // First integral of motion -> (log of) enthalpy in all space
 // ----------------------------------------------------------

 ent = (ent_c + logn_c + mlngamma_c) - logn - mlngamma - prim_field;


 // --------------------------------------------------------------
 // Test: is the enthalpy negative somewhere in the equatorial plane
 //       inside the star?
 // --------------------------------------------------------

 kepler = false ; 
 for (int l=0; l<nzet; l++) {
   int imax = mg->get_nr(l) - 1 ;
   if (l == l_b) imax-- ;  // The surface point is skipped
   for (int i=0; i<imax; i++) { 
     if ( ent.val_grid_point(l, 0, j_b, i) < 0. ) {
       kepler = true ;
       cout << "ent < 0 for l, i : " << l << "  " << i 
	    << "   ent = " << ent.val_grid_point(l, 0, j_b, i) << '\n' ;  
     } 
   }
 }

 if ( kepler ) {
   cout << "**** KEPLERIAN VELOCITY REACHED ****" << '\n' ; 
   omega /= fact_omega ;    // Omega is decreased
   cout << "New central rotation frequency : " 
	<< omega_c/(2.*M_PI) * f_unit << " Hz" << '\n' ; 
   omega_trop_grand = true ;  
 }

 } // End of while ( kepler )

 if ( omega_trop_grand ) {   // fact_omega is decreased for the
                             //  next step 
   fact_omega = sqrt( fact_omega ) ; 
   cout << "**** New fact_omega : " << fact_omega << '\n' ; 
 }


//---------------------------------
 // Equation of state
 //---------------------------------

 equation_of_state() ;  // computes new values for nbar (n), ener (e),
                        // and press (p) from the new ent (H)

 hydro_euler() ;

 //---------------------------------------------//
 // Resolution of the Poisson equation for qqq  //
 //---------------------------------------------//

 Scalar q_new(mp) ;

 solve_qqq( q_new ) ;

 q_new.std_spectral_base() ;

 //----------------------------------------------//
 // Resolution of the Poisson equation for hh    //
 //----------------------------------------------//

 Sym_tensor_trans hij_new(mp, mp.get_bvect_spher(), flat) ;

 solve_hij( hij_new ) ;

 hh = hij_new ; 
 qqq = q_new ;
 beta = beta_new ; 

 //---------------------------------------
 // Calculate error of the GRV2 identity 
 //---------------------------------------

 err_grv2 = grv2() ;


 //--------------------------------------
 // Relaxations on some quantities....?
 //
 // ** On logn and qqq?
 //--------------------------------------

 if (mer >= 10) {
   logn = relax * logn + relax_prev * logn_prev ;

   qqq = relax * qqq + relax_prev * qqq_prev ;

 }

 // Update of the metric quantities :

   update_metric() ;

 //---------------------------
 //  Informations display
 // More to come later......
 //---------------------------

 // partial_display(cout) ;    // What is partial_display(cout) ? 
 fichfreq << "  " << omega_c / (2*M_PI) * f_unit ; 
 fichevol << "  " << ent_c ; 

 //-----------------------------------------------------------
 // Relative change in enthalpy with respect to previous step
 // ** Check:  Is diffrel(ent, ent_prev) ok?  
 //-----------------------------------------------------------

 Tbl diff_ent_tbl = diffrel( ent, ent_prev ) ;
 diff_ent = diff_ent_tbl(0) ;
 for (int l=1; l<nzet; l++) {
   diff_ent += diff_ent_tbl(l) ;
 }
 diff_ent /= nzet ;

 fichconv << " " << log10( fabs(diff_ent) + 1.e-16 ) ;
 fichconv << " " << log10( fabs(err_grv2) + 1.e-16 ) ;

 //------------------------------
 // Recycling for the next step
 //------------------------------

 ent_prev = ent ;
 logn_prev = logn ;
 qqq_prev = qqq ;

 fichconv << '\n' ;     
 fichfreq << '\n' ;     
 fichevol << '\n' ;
 fichconv.flush() ; 
 fichfreq.flush() ; 
 fichevol.flush() ; 
      
	 
  } // End of main loop				 
 
  //=================================================
  //             End of iteration
  //=================================================

  fichconv.close() ; 
  fichfreq.close() ; 
  fichevol.close() ; 


}
}
