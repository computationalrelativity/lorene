/*
 * Method of class Et_rot_bifluid to compute a static spherical configuration.
 *
 * (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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


char et_bfrot_equilibre_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/01/03 15:30:28  j_novak
 * Some comments modified.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2001/06/22  15:40:06  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Entetes du C++
#include <iostream.h>
#include <fstream.h>

// Headers C
#include <math.h>

// Headers Lorene
#include "et_rot_bifluid.h"
#include "param.h"

#include "graphique.h"
#include "utilitaires.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//                   Spherical Equilibrium (Not tested yet!!!)

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void Et_rot_bifluid::equilibrium_spher_bi(double ent_c, double ent2_c, 
				       double precis){
    
  // Fundamental constants and units
  // -------------------------------
#include "unites.h"	    
  // To avoid some compilation warnings
  if (ent_c < 0) {
    cout << f_unit << msol << mevpfm3 << endl ; 
  }    
    
  // Initializations
  // ---------------
    
  const Mg3d* mg = mp.get_mg() ; 
  int nz = mg->get_nzone() ;	    // total number of domains
  uuu = 0 ;
  uuu2 = 0 ; // needed for the EOS
  gam_euler = 1 ;
  gam_euler2 = 1 ;
    
  // Index of the point at phi=0, theta=pi/2 at the surface of the star:
  int l_b = nzet - 1 ; 
  int i_b = mg->get_nr(l_b) - 1 ; 
  int j_b = mg->get_nt(l_b) - 1 ; 
  int k_b = 0 ; 
    
  // Value of the enthalpy defining the surface of the star
  double ent_b = 0 ; 
    
  // Initialization of the enthalpy field to the constant value ent_c :
    
  ent = ent_c ; 
  ent.annule(nzet, nz-1) ; 

  ent2 = ent2_c ;
  ent2.annule(nzet, nz-1) ;

    
  // Corresponding profiles of baryon density, energy density and pressure
    
  equation_of_state() ; 
    
  // Initial metric 
  a_car = 1 ;	     // this value will remain unchanged in the Newtonian case
  beta_auto = 0 ;  // this value will remain unchanged in the Newtonian case
    

  // Auxiliary quantities
  // --------------------
    
  // Affine mapping for solving the Poisson equations
  Map_af mpaff(mp);	
        
  Param par_nul   ;	 // Param (null) for Map_af::poisson.

  Tenseur ent_jm1(mp) ;	// Enthalpy at previous step
  ent_jm1 = 0 ; 
    
  Tenseur ent2_jm1(mp) ;
  ent2_jm1 = 0 ;

  Tenseur source(mp) ; 
  Tenseur logn(mp) ; 
  Tenseur logn_quad(mp) ; 
  logn = 0 ; 
  logn_quad = 0 ; 

  Cmp dlogn(mp) ; 
  Cmp dbeta(mp) ; 
    
  double diff_ent = 1 ; 
  int mermax = 200 ;	    // Max number of iterations
    
  //=========================================================================
  // 			Start of iteration
  //=========================================================================

  for(int mer=0 ; (diff_ent > precis) && (mer<mermax) ; mer++ ) {

    double alpha_r = 1 ; 
	
    cout << "-----------------------------------------------" << endl ;
    cout << "step: " << mer << endl ;
    cout << "alpha_r: " << alpha_r << endl ;
    cout << "diff_ent = " << diff_ent << endl ;    

    //-----------------------------------------------------
    // Resolution of Poisson equation for ln(N)
    //-----------------------------------------------------

    // Matter part of ln(N)
    // --------------------
    if (relativistic) {
      source = a_car * (ener + 3*press) ;
    }
    else {
      source = nbar + nbar2 ; 
    }
	
    (source.set()).set_dzpuis(4) ; 
	
    source.set_std_base() ;	    // Sets the standard spectral bases. 
	
    logn_auto.set_etat_qcq() ; 

    mpaff.poisson(source(), par_nul, logn_auto.set()) ; 

    // NB: at this stage logn_auto is in machine units, not in c^2

    // Quadratic part of ln(N)
    // -----------------------

    if (relativistic) {
	    
      mpaff.dsdr(logn(), dlogn) ; 
      mpaff.dsdr(beta_auto(), dbeta) ; 
	    
      source = - dlogn * dbeta ; 
		      
      logn_quad.set_etat_qcq() ; 
	    
      mpaff.poisson(source(), par_nul, logn_quad.set()) ; 
	    	    
    }

    //-----------------------------------------------------
    // Computation of the new radial scale
    //-----------------------------------------------------

    // alpha_r (r = alpha_r r') is determined so that the enthalpy
    // takes the requested value ent_b at the stellar surface
	
    double nu_mat0_b  = logn_auto()(l_b, k_b, j_b, i_b) ; 
    double nu_mat0_c  = logn_auto()(0, 0, 0, 0) ; 

    double nu_quad0_b  = logn_quad()(l_b, k_b, j_b, i_b) ; 
    double nu_quad0_c  = logn_quad()(0, 0, 0, 0) ; 

    double alpha_r2 = ( ent_c - ent_b - nu_quad0_b + nu_quad0_c )
      / ( qpig*(nu_mat0_b - nu_mat0_c) ) ;
    double alpha2_r2 = ( ent2_c - ent_b - nu_quad0_b + nu_quad0_c )
      / ( qpig*(nu_mat0_b - nu_mat0_c) ) ;

    alpha_r = sqrt(alpha_r2) ;
    double alpha2_r = sqrt(alpha2_r2) ;
    alpha_r = (alpha_r > alpha2_r ? alpha_r : alpha2_r) ;
	
    // New radial scale
    mpaff.homothetie( alpha_r ) ; 

    //--------------------
    // First integral
    //--------------------

    // Gravitation potential in units c^2 :
    logn_auto = alpha_r2*qpig * logn_auto ;
    logn = logn_auto + logn_quad ;

    // Enthalpy in all space
    double logn_c = logn()(0, 0, 0, 0) ;
    ent = ent_c - logn() + logn_c ;
    //ent = 0.5*(ent + abs(ent)) let's keep negative values of ent

    ent2 = ent2_c - logn() + logn_c ;
    //ent2 = 0.5*(ent2 + abs(ent2)) ;

    //---------------------
    // Equation of state
    //---------------------
	
    equation_of_state() ; 
	
    if (relativistic) {
	    
      //----------------------------
      // Equation for beta = ln(AN)
      //----------------------------
	    
      mpaff.dsdr(logn(), dlogn) ; 
      mpaff.dsdr(beta_auto(), dbeta) ; 
	    
      source = 3 * qpig * a_car * press ;
	    
      source = source() 
	- 0.5 * (  dlogn * dlogn + dbeta * dbeta ) ;
	
      source.set_std_base() ;	    // Sets the standard spectral bases. 

      beta_auto.set_etat_qcq() ; 
	         
      mpaff.poisson(source(), par_nul, beta_auto.set()) ; 
	    

      // Metric coefficient A^2 update
	    
      a_car = exp(2*(beta_auto - logn)) ;
   	    
    }
	
    // Relative difference with enthalpy at the previous step
    // ------------------------------------------------------
	
    double xx = norme( diffrel(ent2(), ent2_jm1()) ) / nzet ;
    diff_ent = 0.5*(norme( diffrel(ent(), ent_jm1()) ) / nzet + xx) ; 
	
    // Next step
    // ---------
    
    ent_jm1 = ent ; 
    ent2_jm1 = ent2 ;


  }  // End of iteration loop 
    
  //=========================================================================
  // 			End of iteration
  //=========================================================================


  // The mapping is transfered to that of the star:
  // ----------------------------------------------
  mp = mpaff ; 
    
    
  // Sets value to all the Tenseur's of the star
  // -------------------------------------------
    
  // ... hydro 
  ent.annule(nzet, nz-1) ;	// enthalpy set to zero at the exterior of 
				// the star
  ener_euler = ener ; 
  s_euler = 3 * press ;
  u_euler = 0 ; 
    
  // ... metric
  nnn = exp( unsurc2 * logn ) ; 
  shift = 0 ; 

  // Info printing
  // -------------
    
  cout << endl 
       << "Characteristics of the star obtained by Et_rot_bifluid::equilibrium_spher : " 
       << endl 
       << "-----------------------------------------------------------------" 
       << endl ; 

  double ray = mp.val_r(l_b, 1., M_PI/2., 0) ; 
  cout << "Coordinate radius  :       " << ray / km << " km" << endl ; 

  double rcirc = ray * sqrt( a_car()(l_b, k_b, j_b, i_b) ) ; 
	
  double compact = qpig/(4.*M_PI) * mass_g() / rcirc ; 

  cout << "Circumferential radius R : " << rcirc/km  << " km" << endl ;
  cout << "Baryon mass :	     " << mass_b()/msol << " Mo" << endl ;
  cout << "Gravitational mass M :  " << mass_g()/msol << " Mo" << endl ;
  cout << "Compacity parameter GM/(c^2 R) : " << compact << endl ;


  //-----------------
  // Virial theorem
  //-----------------
    
  //... Pressure term

  source = qpig * a_car * sqrt(a_car) * s_euler ; 
  source.set_std_base() ;	    
  double vir_mat = source().integrale() ; 
    
  //... Gravitational term

  Cmp tmp = beta_auto() - logn() ; 

  source =  - ( logn().dsdr() * logn().dsdr() 
		- 0.5 * tmp.dsdr() * tmp.dsdr() ) 
    * sqrt(a_car())  ; 

  source.set_std_base() ;	    
  double vir_grav = source().integrale() ; 

  //... Relative error on the virial identity GRV3
    
  double grv3 = ( vir_mat + vir_grav ) / vir_mat ;

  cout << "Virial theorem GRV3 : " << endl ; 
  cout << "     3P term    : " << vir_mat << endl ; 
  cout << "     grav. term : " << vir_grav << endl ; 
  cout << "     relative error : " << grv3 << endl ; 
    
    
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//                        Axial Equilibrium

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void Et_rot_bifluid::equilibrium_bi
(double ent_c, double ent2_c, double omega0, double omega20, 
 const Tbl& ent_limit, const Tbl& ent2_limit, const Itbl& icontrol, 
 const Tbl& control, Tbl& diff) {
			     
  // Fundamental constants and units
  // -------------------------------
#include "unites.h"	    
  // To avoid some compilation warnings
  if (ent_c < 0) {
    cout << f_unit << msol << km << mevpfm3 << endl ; 
  }    
    
  // For the display 
  // ---------------
  char display_bold[]="x[1m" ; display_bold[0] = 27 ;
  char display_normal[] = "x[0m" ; display_normal[0] = 27 ;

  // Grid parameters
  // ---------------
    
  const Mg3d* mg = mp.get_mg() ; 
  int nz = mg->get_nzone() ;	    // total number of domains
  int nzm1 = nz - 1 ; 
    
    // Index of the point at phi=0, theta=pi/2 at the surface of the star:
  assert(mg->get_type_t() == SYM) ; 
  int l_b = nzet - 1 ; 
  int i_b = mg->get_nr(l_b) - 1 ; 
  int j_b = mg->get_nt(l_b) - 1 ; 
  int k_b = 0 ; 
    
  // Value of the enthalpies defining the surface of each fluid
  double ent_b = ent_limit(nzet-1) ;
  double ent2_b = ent2_limit(nzet-1) ;
  // This value is chosen so that the grid contain both fluids
  ent_b = (ent_b > ent2_b ? ent_b : ent2_b) ; 
    
  // Parameters to control the iteration
  // -----------------------------------
    
  int mer_max = icontrol(0) ; 
  int mer_rot = icontrol(1) ;
  int mer_change_omega = icontrol(2) ; 
  int mer_fix_omega = icontrol(3) ; 
  int mermax_poisson = icontrol(4) ; 
  int niter ;

  // Protections:
  if (mer_change_omega < mer_rot) {
    cout << "Et_rot_bifluid::equilibrium: mer_change_omega < mer_rot !" << endl ;
    cout << " mer_change_omega = " << mer_change_omega << endl ; 
    cout << " mer_rot = " << mer_rot << endl ; 
    abort() ; 
  }
  if (mer_fix_omega < mer_change_omega) {
    cout << "Et_rot_bifluid::equilibrium: mer_fix_omega < mer_change_omega !" 
	 << endl ;
    cout << " mer_fix_omega = " << mer_fix_omega << endl ; 
    cout << " mer_change_omega = " << mer_change_omega << endl ; 
    abort() ; 
  }

  double precis = control(0) ; 
  double omega_ini = control(1) ;
  double omega2_ini = control(2) ;
  double relax = control(3) ;
  double relax_prev = double(1) - relax ;  
  double relax_poisson = control(4) ; 

  // Error indicators
  // ----------------
    
  diff.set_etat_qcq() ; 
  double diff_ent ;
  double& diff_ent1 = diff.set(0) ; 
  double& diff_ent2 = diff.set(1) ; 
  double& diff_nuf = diff.set(2) ; 
  double& diff_nuq = diff.set(3) ; 
  //    double& diff_dzeta = diff.set(4) ; 
  //    double& diff_ggg = diff.set(5) ; 
  double& diff_shift_x = diff.set(6) ; 
  double& diff_shift_y = diff.set(7) ; 
    
  // Parameters for the function Map_et::poisson for nuf
  // ----------------------------------------------------

  double precis_poisson = 1.e-16 ;     

  Param par_poisson_nuf ; 
  par_poisson_nuf.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_nuf.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_nuf.add_double(precis_poisson, 1) ; // required precision
  par_poisson_nuf.add_int_mod(niter, 0) ;  //  number of iterations actually used 
  par_poisson_nuf.add_cmp_mod( ssjm1_nuf ) ; 
					   
  Param par_poisson_nuq ; 
  par_poisson_nuq.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_nuq.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_nuq.add_double(precis_poisson, 1) ; // required precision
  par_poisson_nuq.add_int_mod(niter, 0) ;  //  number of iterations actually used 
  par_poisson_nuq.add_cmp_mod( ssjm1_nuq ) ; 
					   
  Param par_poisson_tggg ; 
  par_poisson_tggg.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_tggg.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_tggg.add_double(precis_poisson, 1) ; // required precision
  par_poisson_tggg.add_int_mod(niter, 0) ;  //  number of iterations actually used 
  par_poisson_tggg.add_cmp_mod( ssjm1_tggg ) ; 
  double lambda_tggg ;
  par_poisson_tggg.add_double_mod( lambda_tggg ) ; 
    
  Param par_poisson_dzeta ; 
  double lambda_grv2 ;
  par_poisson_dzeta.add_double_mod( lambda_grv2 ) ; 
 					   
  // Parameters for the function Tenseur::poisson_vect
  // -------------------------------------------------

  Param par_poisson_vect ; 

  par_poisson_vect.add_int(mermax_poisson,  0) ;  // maximum number of iterations
  par_poisson_vect.add_double(relax_poisson,  0) ; // relaxation parameter
  par_poisson_vect.add_double(precis_poisson, 1) ; // required precision
  par_poisson_vect.add_cmp_mod( ssjm1_khi ) ; 
  par_poisson_vect.add_tenseur_mod( ssjm1_wshift ) ; 
  par_poisson_vect.add_int_mod(niter, 0) ;   

 					   
    // Initializations
    // ---------------

    // Initial angular velocities
  omega = 0 ; 
  omega2 = 0 ;
  
  double accrois_omega = (omega0 - omega_ini) / 
    double(mer_fix_omega - mer_change_omega) ; 
  double accrois_omega2 = (omega20 - omega2_ini) / 
    double(mer_fix_omega - mer_change_omega) ; 


  update_metric() ;	// update of the metric coefficients

  equation_of_state() ;	// update of the densities, pressure, etc...
    
  hydro_euler() ;	// update of the hydro quantities relative to the 
  //  Eulerian observer

    // Quantities at the previous step : 	
  Tenseur ent_prev = ent ;	
  Tenseur ent2_prev = ent2 ;
  Tenseur logn_prev = logn ;	    
  Tenseur dzeta_prev = dzeta ;	    
    
    // Creation of uninitialized tensors:
  Tenseur source_nuf(mp) ;    // source term in the equation for nuf
  Tenseur source_nuq(mp) ;    // source term in the equation for nuq
  Tenseur source_dzf(mp) ;	// matter source term in the eq. for dzeta
  Tenseur source_dzq(mp) ;	// quadratic source term in the eq. for dzeta
  Tenseur source_tggg(mp) ;	// source term in the eq. for tggg
  Tenseur source_shift(mp, 1, CON, mp.get_bvect_cart()) ;  
  // source term for shift
  Tenseur mlngamma(mp) ;	// centrifugal potential
  Tenseur mlngamma2(mp) ;	// centrifugal potential
    
  // Preparations for the Poisson equations:
  // --------------------------------------
  if (nuf.get_etat() == ETATZERO) {
    nuf.set_etat_qcq() ; 
    nuf.set() = 0 ; 
  }
    
  if (relativistic) {
    if (nuq.get_etat() == ETATZERO) {
      nuq.set_etat_qcq() ; 
      nuq.set() = 0 ; 
    }

    if (tggg.get_etat() == ETATZERO) {
      tggg.set_etat_qcq() ; 
      tggg.set() = 0 ; 
    }
	
    if (dzeta.get_etat() == ETATZERO) {
      dzeta.set_etat_qcq() ; 
      dzeta.set() = 0 ; 
    }
  }
		    
  ofstream fichconv("convergence.d") ;    // Output file for diff_ent
  fichconv << "#     diff_ent     GRV2      " << endl ; 
    
  ofstream fichfreq("frequency.d") ;    // Output file for  omega
  fichfreq << "#       f1 [Hz]     f2 [Hz]" << endl ; 
    
  ofstream fichevol("evolution.d") ;    // Output file for various quantities
  fichevol << 
    "#           r_pole/r_eq	     ent_c    ent2_c" 
	   << endl ; 
    
  diff_ent = 1 ; 
  double err_grv2 = 1 ; 
    
    //=========================================================================
    // 			Start of iteration
    //=========================================================================

  for(int mer=0 ; (diff_ent > precis) && (mer<mer_max) ; mer++ ) {

    cout << "-----------------------------------------------" << endl ;
    cout << "step: " << mer << endl ;
    cout << "diff_ent = " << display_bold << diff_ent << display_normal
	 << endl ;    
    cout << "err_grv2 = " << err_grv2 << endl ;    
    fichconv << mer ;
    fichfreq << mer ;
    fichevol << mer ;
      
    if (mer >= mer_rot) {
	
      if (mer < mer_change_omega) {
	omega = omega_ini ; 
	omega2 = omega2_ini ;
      }
      else {
	if (mer <= mer_fix_omega) {
	  omega = omega_ini + accrois_omega * 
	    (mer - mer_change_omega) ;
	  omega2 = omega2_ini + accrois_omega2 * 
	    (mer - mer_change_omega) ;
	}
      }
	
    }

    //-----------------------------------------------
    //  Sources of the Poisson equations
    //-----------------------------------------------
	
    // Source for nu
    // -------------
    Tenseur beta = log(bbb) ; 
    beta.set_std_base() ; 

    if (relativistic) {
      source_nuf =  qpig * a_car *( ener_euler + s_euler ) ; 

      source_nuq = ak_car - flat_scalar_prod(logn.gradient_spher(), 
					     logn.gradient_spher() + beta.gradient_spher()) ; 
    }
    else {
      source_nuf = qpig * (nbar + nbar2); 

      source_nuq = 0 ; 
    }
    source_nuf.set_std_base() ; 	
    source_nuq.set_std_base() ; 	

    // Source for dzeta
    // ----------------
    source_dzf = 2 * qpig * a_car * (s_euler - 2*press ) ;
    source_dzf.set_std_base() ; 
  
    source_dzq = 1.5 * ak_car - flat_scalar_prod(logn.gradient_spher(),
						 logn.gradient_spher() ) ;	    
    source_dzq.set_std_base() ; 	
	
    // Source for tggg
    // ---------------
	
    source_tggg = 4 * qpig * nnn * a_car * bbb * press ;
    source_tggg.set_std_base() ; 
	
    (source_tggg.set()).mult_rsint() ; 
	

    // Source for shift
    // ----------------
	    
    // Matter term (u_euler is NOT the same as in Etoile_rot): 
    source_shift = (-4*qpig) * nnn * a_car * u_euler ;

    // Quadratic terms:
    Tenseur vtmp =  6 * beta.gradient_spher() - 2 * logn.gradient_spher() ;
    vtmp.change_triad(mp.get_bvect_cart()) ; 

    Tenseur squad  = nnn * flat_scalar_prod(tkij, vtmp) ;     

    // The addition of matter terms and quadratic terms is performed
    //  component by component because u_euler is contravariant,
    //  while squad is covariant. 
		
    if (squad.get_etat() == ETATQCQ) {
      for (int i=0; i<3; i++) {
	source_shift.set(i) += squad(i) ; 
      }
    }

    source_shift.set_std_base() ; 	

    //----------------------------------------------
    // Resolution of the Poisson equation for nuf 
    //----------------------------------------------

    source_nuf().poisson(par_poisson_nuf, nuf.set()) ; 
	
    cout << "Test of the Poisson equation for nuf :" << endl ; 
    Tbl err = source_nuf().test_poisson(nuf(), cout, true) ; 
    diff_nuf = err(0, 0) ; 
	
    if (relativistic) {
	    
      //----------------------------------------------
      // Resolution of the Poisson equation for nuq 
      //----------------------------------------------

      source_nuq().poisson(par_poisson_nuq, nuq.set()) ; 
	    
      cout << "Test of the Poisson equation for nuq :" << endl ; 
      err = source_nuq().test_poisson(nuq(), cout, true) ;
      diff_nuq = err(0, 0) ; 
	
      //---------------------------------------------------------
      // Resolution of the vector Poisson equation for the shift
      //---------------------------------------------------------


      if (source_shift.get_etat() != ETATZERO) {

	for (int i=0; i<3; i++) {
	  if(source_shift(i).dz_nonzero()) {
	    assert( source_shift(i).get_dzpuis() == 4 ) ; 
	  }
	  else{
	    (source_shift.set(i)).set_dzpuis(4) ; 
	  }
	}

      }

      double lambda_shift = double(1) / double(3) ; 

      if ( mg->get_np(0) == 1 ) {
	lambda_shift = 0 ; 
      }
	
      source_shift.poisson_vect(lambda_shift, par_poisson_vect, 
				shift, w_shift, khi_shift) ;      
	    
      cout << "Test of the Poisson equation for shift_x :" << endl ; 
      err = source_shift(0).test_poisson(shift(0), cout, true) ;
      diff_shift_x = err(0, 0) ; 
	
      cout << "Test of the Poisson equation for shift_y :" << endl ; 
      err = source_shift(1).test_poisson(shift(1), cout, true) ;
      diff_shift_y = err(0, 0) ; 
	    
      // Computation of tnphi and nphi from the Cartesian components
      //  of the shift
      // -----------------------------------------------------------
	    
      fait_nphi() ; 
	
    }

    //------------------------------------------------
    // Determination of the fluid velocities U1 and U2
    //------------------------------------------------
	
    Cmp tmp = omega - nphi() ; 
    tmp.annule(nzm1) ; 
    tmp.std_base_scal() ;
	
    tmp.mult_rsint() ;	    //  Multiplication by r sin(theta)
	
    uuu = bbb() / nnn() * tmp ; 
	
    if (uuu.get_etat() == ETATQCQ) {
      // Same basis as (Omega -N^phi) r sin(theta) :
      ((uuu.set()).va).set_base( (tmp.va).base ) ;   
    }
	
    tmp = omega2 - nphi() ; 
    tmp.annule(nzm1) ; 
    tmp.std_base_scal() ;
	
    tmp.mult_rsint() ;	    //  Multiplication by r sin(theta)
	
    uuu2 = bbb() / nnn() * tmp ; 
	
    if (uuu2.get_etat() == ETATQCQ) {
      // Same basis as (Omega -N^phi) r sin(theta) :
      ((uuu2.set()).va).set_base( (tmp.va).base ) ;   
    }
	
    // Is one of the new velocities larger than c in the equatorial plane ?
	
    bool superlum = false ; 
	
    for (int l=0; l<nzet; l++) {
      for (int i=0; i<mg->get_nr(l); i++) {
	    
	double u1 = uuu()(l, 0, j_b, i) ; 
	double u2 = uuu2()(l, 0, j_b, i) ;
	if ((u1 >= 1.) || (u2>=1.)) {	    // superluminal velocity
	  superlum = true ; 
	  cout << "U > c  for l, i : " << l << "  " << i 
	       << "   U1 = " << u1 << endl ;
	  cout << "   U2 = " << u2 << endl ;
	}
      }
    }
    if ( superlum ) {
      cout << "**** VELOCITY OF LIGHT REACHED ****" << endl ; 
      abort() ;
    }
	
    // New computation of gam_euler, ener_euler, etc...
    // ------------------------------------------------
	
    hydro_euler() ; 
	

    //------------------------------------------------------
    //	First integral of motion 
    //------------------------------------------------------
	
    // Centrifugal potential : 
    if (relativistic) {
      mlngamma = - log( gam_euler ) ;
      mlngamma2 = - log( gam_euler2) ;
    }
    else {
      mlngamma = - 0.5 * uuu*uuu ;
      mlngamma2 = -0.5 * uuu2*uuu2 ;
    }
	
    // Equatorial values of various potentials :
    double nuf_b  = nuf()(l_b, k_b, j_b, i_b) ; 
    double nuq_b  = nuq()(l_b, k_b, j_b, i_b) ; 
    double mlngamma_b  = mlngamma()(l_b, k_b, j_b, i_b) ; 
    double mlngamma2_b  = mlngamma2()(l_b, k_b, j_b, i_b) ; 
	
    // Central values of various potentials :
    double nuf_c = nuf()(0,0,0,0) ; 
    double nuq_c = nuq()(0,0,0,0) ; 
    double mlngamma_c = 0 ;
    double mlngamma2_c = 0 ;
	
    // Scale factor to ensure that the enthalpy is equal to ent_b at 
    //  the equator for the "outer" fluid
    double alpha_r2 = ( ent_c - ent_b + mlngamma_c - mlngamma_b
			+ nuq_c - nuq_b) / ( nuf_b - nuf_c  ) ;
    double alpha2_r2 = ( ent2_c - ent_b + mlngamma2_c - mlngamma2_b
			 + nuq_c - nuq_b) / ( nuf_b - nuf_c  ) ;
    double alpha_r = sqrt(alpha_r2) ;
    alpha_r = (alpha_r > sqrt(alpha2_r2) ? alpha_r : sqrt(alpha2_r2)) ;
    alpha_r2 = alpha_r * alpha_r ;
    cout << "alpha_r = " << alpha_r << endl ; 

    // Rescaling of the grid (no adaptation!)
    //---------------------------------------
    mp.homothetie(alpha_r) ;
	
    // Readjustment of nu :
    // -------------------
	
    logn = alpha_r2 * nuf + nuq ;
    double nu_c =  logn()(0,0,0,0) ;
	
    // First integral	--> enthalpy in all space
    //-----------------
	
    ent = (ent_c + nu_c + mlngamma_c) - logn - mlngamma ;
    ent2 = (ent2_c + nu_c + mlngamma2_c) - logn - mlngamma2 ;
	
    //----------------------------------------------------
    // Equation of state  
    //----------------------------------------------------
	
    equation_of_state() ; 	// computes new values for nbar1,2 , ener (e) 
				// and press (p) from the new ent,ent2 
	
    //---------------------------------------------------------
    // Matter source terms in the gravitational field equations	
    //---------------------------------------------------------

    //## Computation of tnphi and nphi from the Cartesian components
    //  of the shift for the test in hydro_euler():
	    
    fait_nphi() ; 

    hydro_euler() ;		// computes new values for ener_euler (E), 
				// s_euler (S) and u_euler (U^i)

    if (relativistic) {

      //-------------------------------------------------------
      //	2-D Poisson equation for tggg
      //-------------------------------------------------------

      mp.poisson2d(source_tggg(), mp.cmp_zero(), par_poisson_tggg,
		   tggg.set()) ; 
	    
      //-------------------------------------------------------
      //	2-D Poisson equation for dzeta
      //-------------------------------------------------------

      mp.poisson2d(source_dzf(), source_dzq(), par_poisson_dzeta,
		   dzeta.set()) ; 
	    
      err_grv2 = lambda_grv2 - 1; 
      cout << "GRV2: " << err_grv2 << endl ; 
	    
    }


    //---------------------------------------
    // Computation of the metric coefficients (except for N^phi)
    //---------------------------------------

    // Relaxations on nu and dzeta :  

    logn = relax * logn + relax_prev * logn_prev ;

    dzeta = relax * dzeta + relax_prev * dzeta_prev ; 

    // Update of the metric coefficients N, A, B and computation of K_ij :

    update_metric() ; 
	
    //-----------------------
    //  Informations display
    //-----------------------

    partial_display(cout) ; 
    fichfreq << "  " << omega / (2*M_PI) * f_unit ; 
    fichfreq << "  " << omega2 / (2*M_PI) * f_unit ; 
    fichevol << "  " << ray_pole() / ray_eq() ; 
    fichevol << "  " << ent_c ; 
    fichevol << "  " << ent2_c ; 

    //-------------------------------------------------------------
    //  Relative change in enthalpies with respect to previous step 
    //-------------------------------------------------------------

    Tbl diff_ent_tbl = diffrel( ent(), ent_prev() ) ; 
    diff_ent1 = diff_ent_tbl(0) ; 
    for (int l=1; l<nzet; l++) {
      diff_ent1 += diff_ent_tbl(l) ; 
    }
    diff_ent1 /= nzet ; 
    diff_ent_tbl = diffrel( ent2(), ent2_prev() ) ;
    diff_ent2 = diff_ent_tbl(0) ; 
    for (int l=1; l<nzet; l++) {
      diff_ent2 += diff_ent_tbl(l) ; 
    }
    diff_ent2 /= nzet ;
    diff_ent = 0.5*(diff_ent1 + diff_ent2) ; 
	
    fichconv << "  " << log10( fabs(diff_ent) + 1.e-16 ) ;
    fichconv << "  " << log10( fabs(err_grv2) + 1.e-16 ) ;

    //------------------------------
    //  Recycling for the next step
    //------------------------------
	
    ent_prev = ent ; 
    ent2_prev = ent2 ; 
    logn_prev = logn ; 
    dzeta_prev = dzeta ; 
	
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
    

}

