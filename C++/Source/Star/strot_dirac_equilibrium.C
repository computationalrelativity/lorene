/*
 *  Function Star_rot_Dirac::equilibrium
 *
 *    (see file star_rot_dirac.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
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

char strot_dirac_equilibrium_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.11  2009/10/26 10:54:33  j_novak
 * Added the case of a NONSYM base in theta.
 *
 * Revision 1.10  2008/02/25 10:40:52  j_novak
 * Added the flag mer_hij to control the step from which the equation for h^{ij}
 * is being solved.
 *
 * Revision 1.9  2006/03/14 15:18:21  lm_lin
 *
 * Add convergence to a given baryon mass.
 *
 * Revision 1.8  2005/11/28 14:45:16  j_novak
 * Improved solution of the Poisson tensor equation in the case of a transverse
 * tensor.
 *
 * Revision 1.7  2005/09/16 14:04:49  j_novak
 * The equation for hij is now solved only for mer >  mer_fix_omega. It uses the
 * Poisson solver of the class Sym_tensor_trans.
 *
 * Revision 1.6  2005/04/20 14:26:29  j_novak
 * Removed some outputs.
 *
 * Revision 1.5  2005/04/05 09:24:05  j_novak
 * minor modifs
 *
 * Revision 1.4  2005/03/10 09:39:19  j_novak
 * The order of resolution has been changed in the iteration step.
 *
 * Revision 1.3  2005/02/17 17:30:09  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.2  2005/02/09 13:36:42  lm_lin
 *
 * Calculate GRV2 during iterations.
 *
 * Revision 1.1  2005/01/31 08:51:48  j_novak
 * New files for rotating stars in Dirac gauge (still under developement).
 *
 *
 * $Header$
 *
 */


// C headers
#include <math.h>
#include <assert.h>

// Lorene headers
#include "star_rot_dirac.h"

#include "utilitaires.h"
#include "unites.h" 

void Star_rot_Dirac::equilibrium(double ent_c, double omega0, 
	         double fact_omega, int , const Tbl& ent_limit,
		 const Itbl& icontrol, const Tbl& control,
		 double mbar_wanted, double aexp_mass, Tbl& diff){	
	

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
  int type_t = mg->get_type_t() ; 
  assert( ( type_t == SYM) || (type_t == NONSYM) ) ; 
  int l_b = nzet - 1 ; 
  int i_b = mg->get_nr(l_b) - 1 ; 
  int j_b = (type_t == SYM ? mg->get_nt(l_b) - 1 : mg->get_nt(l_b)/2 ) ; 
  int k_b = 0 ;

  // Value of the enthalpy defining the surface of the star
  double ent_b = ent_limit(nzet-1) ;

  // Parameters to control the iteration
  // -----------------------------------

  int mer_max = icontrol(0) ;
  int mer_rot = icontrol(1) ;
  int mer_change_omega = icontrol(2) ; 
  int mer_fix_omega = icontrol(3) ; 
  int mer_mass = icontrol(4) ; 
  int delta_mer_kep = icontrol(5) ; 
  int mer_hij = icontrol(6) ;
  
  // Protections:
  if (mer_change_omega < mer_rot) {
    cout << "Star_rot_Dirac::equilibrium: mer_change_omega < mer_rot !" 
	 << endl ;
    cout << " mer_change_omega = " << mer_change_omega << endl ; 
    cout << " mer_rot = " << mer_rot << endl ; 
    abort() ; 
  }
  if (mer_fix_omega < mer_change_omega) {
    cout << "Star_rot_Dirac::equilibrium: mer_fix_omega < mer_change_omega !" 
	 << endl ;
    cout << " mer_fix_omega = " << mer_fix_omega << endl ; 
    cout << " mer_change_omega = " << mer_change_omega << endl ; 
    abort() ; 
  }

  // In order to converge to a given baryon mass, shall the central
  // enthalpy be varied or Omega ?
  bool change_ent = true ; 
  if (mer_mass < 0) {
    change_ent = false ; 
    mer_mass = abs(mer_mass) ;
  }
  

  double precis = control(0) ; 
  double omega_ini = control(1) ; 
  double relax = control(2) ;
  double relax_prev = double(1) - relax ;  
      
  // Error indicators
  // ----------------

  diff.annule_hard() ;
  double& diff_ent = diff.set(0) ; 

  double alpha_r = 1 ;

  // Initializations
  // ---------------

  // Initial angular velocities 
  omega = 0 ;
  
  double accrois_omega = (omega0 - omega_ini) / 
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
  fichconv << "#     diff_ent    GRV2    max_triax   vit_triax" << endl ;   

  ofstream fichfreq("frequency.d") ;    // Output file for  omega
  fichfreq << "#       f [Hz]" << endl ; 

  ofstream fichevol("evolution.d") ;    // Output file for various quantities
  fichevol << 
    "#       |dH/dr_eq/dH/dr_pole|      r_pole/r_eq ent_c" 
	   << endl ; 

  diff_ent = 1 ;
  double err_grv2 = 1 ;



//=========================================================================
//                   Start of iteration
//=========================================================================

  for(int mer=0 ; (diff_ent > precis) && (mer<mer_max) ; mer++) {

 cout << "-----------------------------------------------" << endl ;
 cout << "step: " << mer << endl ;
 cout << "diff_ent = " << display_bold << diff_ent << display_normal
      << endl ;    
 cout << "err_grv2 = " << err_grv2 << endl ;    
 fichconv << mer ;
 fichfreq << mer ;
 fichevol << mer ;


 // switch on rotation 
 if (mer >= mer_rot) {
          
   if (mer < mer_change_omega) {
     omega = omega_ini ; 
   }
   else {
     if (mer <= mer_fix_omega) {
       omega = omega_ini + accrois_omega * 
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
              
   omega *= fact_omega ;  // Increase of the angular velocity if 
 }              //  fact_omega != 1
 
 bool omega_trop_grand = false ;
 bool kepler = true ;

 while ( kepler ) {

   // Possible decrease of Omega to ensure a velocity < c

   bool superlum = true ;

   while ( superlum ){
   
     // New fluid velocity :
     //

     u_euler.set(1).set_etat_zero() ;
     u_euler.set(2).set_etat_zero() ;

     u_euler.set(3) = omega ;
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
		<< "   U = " << sqrt( u2 ) << endl ;  
	 }
       }
     }
     if ( superlum ) {
       cout << "**** VELOCITY OF LIGHT REACHED ****" << endl ; 
       omega /= fact_omega ;    // Decrease of Omega
       cout << "New rotation frequency : " 
	    << omega/(2.*M_PI) * f_unit <<  " Hz" << endl ; 
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

 // Central values of various potentials :
 double ln_f_c = ln_f_new.val_grid_point(0,0,0,0) ;
 double ln_q_c = ln_q_new.val_grid_point(0,0,0,0) ;
 double mlngamma_c = 0 ;

 // Scale factor to ensure that the (log of) enthalpy is equal to 
 // ent_b at the equator
 double alpha_r2 = ( ent_c - ent_b + mlngamma_c - mlngamma_b
		     + ln_q_c - ln_q_b) / ( ln_f_b - ln_f_c  ) ;

 alpha_r = sqrt(alpha_r2) ;

 cout << "alpha_r = " << alpha_r << endl ; 

 // Rescaling of the grid (no adaptation!)
 //---------------------------------------
 mp.homothetie(alpha_r) ;

 // Readjustment of logn :
 // -----------------------

 logn = alpha_r2 * ln_f_new + ln_q_new ;

 double logn_c = logn.val_grid_point(0,0,0,0) ;

 // First integral of motion -> (log of) enthalpy in all space
 // ----------------------------------------------------------

 ent = (ent_c + logn_c + mlngamma_c) - logn - mlngamma ;


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
	    << "   ent = " << ent.val_grid_point(l, 0, j_b, i) << endl ;  
     } 
   }
 }

 if ( kepler ) {
   cout << "**** KEPLERIAN VELOCITY REACHED ****" << endl ; 
   omega /= fact_omega ;    // Omega is decreased
   cout << "New rotation frequency : " 
	<< omega/(2.*M_PI) * f_unit << " Hz" << endl ; 
   omega_trop_grand = true ;  
 }

 } // End of while ( kepler )

 if ( omega_trop_grand ) {   // fact_omega is decreased for the
                             //  next step 
   fact_omega = sqrt( fact_omega ) ; 
   cout << "**** New fact_omega : " << fact_omega << endl ; 
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

 if (mer > mer_hij )
   solve_hij( hij_new ) ;
 else
   hij_new.set_etat_zero() ;

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
 fichfreq << "  " << omega / (2*M_PI) * f_unit ; 
 fichevol << "  " << ent_c ; 

 
 //-----------------------------------------
 // Convergence towards a given baryon mass 
 //-----------------------------------------

 if (mer > mer_mass) {
	    
   double xx ; 
   if (mbar_wanted > 0.) {
     xx = mass_b() / mbar_wanted - 1. ;
     cout << "Discrep. baryon mass <-> wanted bar. mass : " << xx 
	  << endl ; 
   }
   else{
     xx = mass_g() / fabs(mbar_wanted) - 1. ;
     cout << "Discrep. grav. mass <-> wanted grav. mass : " << xx 
	  << endl ; 
   }
   double xprog = ( mer > 2*mer_mass) ? 1. : 
     double(mer-mer_mass)/double(mer_mass) ; 
   xx *= xprog ; 
   double ax = .5 * ( 2. + xx ) / (1. + xx ) ; 
   double fact = pow(ax, aexp_mass) ; 
   cout << "  xprog, xx, ax, fact : " << xprog << "  " <<
     xx << "  " << ax << "  " << fact << endl ; 
   
   if ( change_ent ) {
     ent_c *= fact ; 
   }
   else {
     if (mer%4 == 0) omega *= fact ; 
   }
 }
	

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

 fichconv << endl ;     
 fichfreq << endl ;     
 fichevol << endl ;
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
