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
 * Revision 1.16  2004/09/01 10:56:05  r_prix
 * added option of converging baryon-mass to equilibrium_bi()
 *
 * Revision 1.15  2004/08/30 09:54:20  r_prix
 * experimental version of Kepler-limit finder for 2-fluid stars
 *
 * Revision 1.14  2004/03/25 10:29:03  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.13  2003/12/11 12:43:35  r_prix
 * activated adaptive grid for 2-fluid star (taken from Etoile_rot)
 *
 * Revision 1.12  2003/12/04 14:28:26  r_prix
 * allow for the case of "slow-rot-style" EOS inversion, in which we need to adapt
 * the inner domain to n_outer=0 instead of mu_outer=0 ...
 * (this should only be used for comparison to analytic slow-rot solution!)
 *
 * Revision 1.11  2003/11/25 12:49:44  j_novak
 * Modified headers to compile on IRIX. Changed Mapping to be Map_af (speed
 * enhancement).
 *
 * Revision 1.10  2003/11/20 14:01:25  r_prix
 * changed member names to better conform to Lorene coding standards:
 * J_euler -> j_euler, EpS_euler -> enerps_euler, Delta_car -> delta_car
 *
 * Revision 1.9  2003/11/19 22:01:57  e_gourgoulhon
 * -- Relaxation on logn and dzeta performed only if mer >= 10.
 * -- err_grv2 is now evaluated also in the Newtonian case.
 *
 * Revision 1.8  2003/11/18 18:38:11  r_prix
 * use of new member EpS_euler: matter sources in equilibrium() and global quantities
 * no longer distinguish Newtonian/relativistic, as all terms should have the right limit...
 *
 * Revision 1.7  2003/11/17 13:49:43  r_prix
 * - moved superluminal check into hydro_euler()
 * - removed some warnings
 *
 * Revision 1.6  2003/11/13 12:14:35  r_prix
 * *) removed all use of etoile-type specific u_euler, press
 *   and use 3+1 components of Tmunu instead
 *
 * Revision 1.5  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.4  2002/04/05 09:09:36  j_novak
 * The inversion of the EOS for 2-fluids polytrope has been modified.
 * Some errors in the determination of the surface were corrected.
 *
 * Revision 1.3  2002/01/16 15:03:28  j_novak
 * *** empty log message ***
 *
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

// Headers C
#include <math.h>

// Headers Lorene
#include "et_rot_bifluid.h"
#include "param.h"
#include "unites.h"	    

#include "graphique.h"
#include "utilitaires.h"

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//                        Axial Equilibrium

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

void Et_rot_bifluid::equilibrium_bi
(double ent_c, double ent2_c, double omega0, double omega20, 
 const Tbl& ent_limit, const Tbl& ent2_limit, const Itbl& icontrol, 
 const Tbl& control, Tbl& diff,
 int mer_mass, double mbar1_wanted, double mbar2_wanted, double aexp_mass) 
{
			     
  // Fundamental constants and units
  // -------------------------------
  using namespace Unites ;

  // For the display 
  // ---------------
  char display_bold[]="x[1m" ; display_bold[0] = 27 ;
  char display_normal[] = "x[0m" ; display_normal[0] = 27 ;

  // Grid parameters
  // ---------------
    
  const Mg3d* mg = mp.get_mg() ; 

  // The following is required to initialize mp_prev as a Map_af
  Map_et& mp_et = dynamic_cast<Map_et&>(mp) ;  // reference
    
  // Index of the point at phi=0, theta=pi/2 at the surface of the star:
  assert(mg->get_type_t() == SYM) ; 
  int l_b = nzet - 1 ; 
  int i_b = mg->get_nr(l_b) - 1 ; 
  int j_b = mg->get_nt(l_b) - 1 ; 
  int k_b = 0 ; 
    
  // Value of the enthalpies defining the surface of each fluid
  double ent1_b = ent_limit(nzet-1) ;
  double ent2_b = ent2_limit(nzet-1) ;

  // This value is chosen so that the grid contain both fluids
  //  double ent_b = (ent1_b > ent2_b ? ent1_b : ent2_b) ; 
    
  // Parameters to control the iteration
  // -----------------------------------
    
  int mer_max = icontrol(0) ; 
  int mer_rot = icontrol(1) ;
  int mer_change_omega = icontrol(2) ; 
  int mer_fix_omega = icontrol(3) ; 
  int mermax_poisson = icontrol(4) ; 
  int nzadapt = icontrol(5);		// number of domains for adaptive grid
  int kepler_fluid = icontrol(6); 	// fluid-index for kepler-limit (0=none,3=both)
  int kepler_wait_steps = icontrol(7);
  
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
  // some additional stuff for adaptive grid:
  double thres_adapt = control(5) ; 
  double precis_adapt = control(6) ; 
  double kepler_factor = control(7);
  if (kepler_factor <= 1.0)
    {
      cout << "ERROR: Kepler factor has to be greater than 1!!\n";
      abort();
    }


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
    

  // Parameters for the function Map_et::adapt
  // -----------------------------------------
    
  Param par_adapt ; 
  int nitermax = 100 ;  
  int adapt_flag = 1 ;    //  1 = performs the full computation, 
			    //  0 = performs only the rescaling by 
			    //      the factor alpha_r
  int nz_search = nzet + 1 ;  // Number of domains for searching the enthalpy
				//  isosurfaces
  double alpha_r ; 
  double reg_map = 1. ; // 1 = regular mapping, 0 = contracting mapping

  par_adapt.add_int(nitermax, 0) ; // maximum number of iterations to
				     // locate zeros by the secant method
  par_adapt.add_int(nzadapt, 1) ; // number of domains where the adjustment 
				    // to the isosurfaces of ent is to be 
				    // performed
  par_adapt.add_int(nz_search, 2) ;	// number of domains to search for
					// the enthalpy isosurface
  par_adapt.add_int(adapt_flag, 3) ; //  1 = performs the full computation, 
				       //  0 = performs only the rescaling by 
				       //      the factor alpha_r
  par_adapt.add_int(j_b, 4) ; //  theta index of the collocation point 
			        //  (theta_*, phi_*)
  par_adapt.add_int(k_b, 5) ; //  theta index of the collocation point 
			        //  (theta_*, phi_*)

  par_adapt.add_int_mod(niter, 0) ;  //  number of iterations actually used in 
				       //  the secant method
    
  par_adapt.add_double(precis_adapt, 0) ; // required absolute precision in 
					     // the determination of zeros by 
					     // the secant method
  par_adapt.add_double(reg_map, 1)	;  // 1. = regular mapping, 0 = contracting mapping
    
  par_adapt.add_double(alpha_r, 2) ;	// factor by which all the radial 
					   // distances will be multiplied 
    	   
  par_adapt.add_tbl(ent_limit, 0) ;	// array of values of the field ent 
				        // to define the isosurfaces. 			   


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
  double lbda_grv2 ;
  par_poisson_dzeta.add_double_mod( lbda_grv2 ) ; 
 					   
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

  Map_et mp_prev = mp_et;
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

  Tenseur *outer_ent_p;		// pointer to the enthalpy field of the outer fluid
    
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
  fichevol << "#           r_pole/r_eq	     ent_c    ent2_c" << endl ; 
    
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

    // common source term for relativistic and Newtonian ! (enerps_euler has the right limit)
    source_nuf =  qpig * a_car * enerps_euler;

    if (relativistic) 
      source_nuq = ak_car - flat_scalar_prod(logn.gradient_spher(),logn.gradient_spher() + beta.gradient_spher()) ; 
    else 
      source_nuq = 0 ; 

    source_nuf.set_std_base() ; 	
    source_nuq.set_std_base() ; 	

    if (relativistic) {
      // Source for dzeta
      // ----------------
      source_dzf = 2 * qpig * a_car * sphph_euler;
      source_dzf.set_std_base() ; 
      
      source_dzq = 1.5 * ak_car - flat_scalar_prod(logn.gradient_spher(),logn.gradient_spher() ) ;	    
      source_dzq.set_std_base() ; 	
      
      // Source for tggg
      // ---------------
      
      source_tggg = 2 * qpig * nnn * a_car * bbb * (s_euler - sphph_euler);
      source_tggg.set_std_base() ; 
      
      (source_tggg.set()).mult_rsint() ; 
      
      
      // Source for shift
      // ----------------
      
      // Matter term 
      source_shift = (-4*qpig) * nnn * a_car * j_euler;
      
      // Quadratic terms:
      Tenseur vtmp =  3 * beta.gradient_spher() - logn.gradient_spher() ;
      vtmp.change_triad(mp.get_bvect_cart()) ; 
      
      Tenseur squad  = 2 * nnn * flat_scalar_prod(tkij, vtmp) ;     
      
      // The addition of matter terms and quadratic terms is performed
      //  component by component because u_euler is contravariant,
      //  while squad is covariant. 
      
      if (squad.get_etat() == ETATQCQ) {
	for (int i=0; i<3; i++) {
	  source_shift.set(i) += squad(i) ; 
	}
      }
      
      source_shift.set_std_base() ; 	
    }
    //----------------------------------------------
    // Resolution of the Poisson equation for nuf 
    //----------------------------------------------
    
    source_nuf().poisson(par_poisson_nuf, nuf.set()) ; 
    
//     cout << "Test of the Poisson equation for nuf :" << endl ; 
//     Tbl err = source_nuf().test_poisson(nuf(), cout, true) ; 
//     diff_nuf = err(0, 0) ; 
    diff_nuf = 0 ;
    
    if (relativistic) {
      
      //----------------------------------------------
      // Resolution of the Poisson equation for nuq 
      //----------------------------------------------

      source_nuq().poisson(par_poisson_nuq, nuq.set()) ; 
	    
//        cout << "Test of the Poisson equation for nuq :" << endl ; 
//        err = source_nuq().test_poisson(nuq(), cout, true) ;
//        diff_nuq = err(0, 0) ; 
      diff_nuq = 0 ;
	
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
	    
//        cout << "Test of the Poisson equation for shift_x :" << endl ; 
//        err = source_shift(0).test_poisson(shift(0), cout, true) ;
//        diff_shift_x = err(0, 0) ; 
      diff_shift_x = 0 ; 
	
//        cout << "Test of the Poisson equation for shift_y :" << endl ; 
//        err = source_shift(1).test_poisson(shift(1), cout, true) ;
//        diff_shift_y = err(0, 0) ; 
      diff_shift_y = 0 ;
	    
      // Computation of tnphi and nphi from the Cartesian components
      //  of the shift
      // -----------------------------------------------------------
	    
      fait_nphi() ; 
	
    }


    //----------------------------------------
    // Shall we search for the Kepler limit?
    //----------------------------------------
    bool kepler = false;
    bool too_fast = false;

    if ( (kepler_fluid > 0) && (mer > mer_fix_omega + kepler_wait_steps) )
      {
	if (kepler_fluid & 0x01)
	  omega *= kepler_factor;
	if (kepler_fluid & 0x02)
	  omega2 *= kepler_factor;
      }


    // ============================================================
    kepler = true;
    while (kepler)	
      {	

	// New computation of delta_car, gam_euler, enerps_euler etc...
	// ------------------------------------------------------
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
	
	// Central values of various potentials :
	double nuf_c = nuf()(0,0,0,0) ; 
	double nuq_c = nuq()(0,0,0,0) ; 

	// Scale factor to ensure that the enthalpy is equal to ent_b at 
	//  the equator for the "outer" fluid
	double alpha_r2 = 0;
	
	int j=j_b;
    
	// Boundary values of various potentials :
	double nuf_b  = nuf()(l_b, k_b, j, i_b) ; 
	double nuq_b  = nuq()(l_b, k_b, j, i_b) ; 
	double mlngamma_b  = mlngamma()(l_b, k_b, j, i_b) ; 
	double mlngamma2_b  = mlngamma2()(l_b, k_b, j, i_b) ; 


	// RP: "hack": adapt the radius correctly if using "slow-rot-style" EOS inversion
	// 
	if ( eos.identify() == 2 ) // only applies to Eos_bf_poly_newt
	  {
	    const Eos_bf_poly_newt &eos0 = dynamic_cast<const Eos_bf_poly_newt&>(eos);
	    if (eos0.get_typeos() == 5)
	      {
		double vn_b = uuu()(l_b, k_b, j, i_b);
		double vp_b = uuu2()(l_b, k_b, j, i_b);
		double D2_b = (vp_b - vn_b)*(vp_b - vn_b);
		double kdet = eos0.get_kap3() + eos0.get_beta()*D2_b;
		double kaps1 = kdet / ( eos0.get_kap2() - kdet );
		double kaps2 = kdet / ( eos0.get_kap1() - kdet );
		
		ent1_b = kaps1 * ( ent2_c - ent_c - mlngamma2_b + mlngamma_b );
		ent2_b = kaps2 * ( ent_c - ent2_c - mlngamma_b + mlngamma2_b );
		
		cout << "**********************************************************************\n";
		cout << "DEBUG: Rescaling domain for slow-rot-style EOS inversion \n";
		cout << "DEBUG: ent1_b = " << ent1_b << "; ent2_b = " << ent2_b << endl;
		cout << "**********************************************************************\n";
		
		adapt_flag = 0;	// don't do adaptive-grid if using slow-rot-style inversion!
	      }
	  }
    
	double alpha1_r2 = ( ent_c - ent1_b - mlngamma_b + nuq_c - nuq_b) / ( nuf_b - nuf_c  ) ;
	double alpha2_r2 = ( ent2_c - ent2_b  - mlngamma2_b + nuq_c - nuq_b) / ( nuf_b - nuf_c  ) ;
	
	cout << "DEBUG: j= "<< j<<" ; alpha1 = " << alpha1_r2 <<" ; alpha2 = " << alpha2_r2 << endl;

	int outer_fluid = (alpha1_r2 > alpha2_r2) ? 1 : 2;  // index of 'outer' fluid (at equator!)
	
	outer_ent_p = (outer_fluid == 1) ? (&ent) : (&ent2);
	
	alpha_r2 = (outer_fluid == 1) ? alpha1_r2 : alpha2_r2 ;
	
	alpha_r = sqrt(alpha_r2);
	
	cout << "alpha_r = " << alpha_r << endl ; 
	
	// Readjustment of nu :
	// -------------------
	
	logn = alpha_r2 * nuf + nuq ;
	double nu_c =  logn()(0,0,0,0) ;
	
	
	// First integral	--> enthalpy in all space
	//-----------------
	
	ent = (ent_c + nu_c) - logn - mlngamma ;
	ent2 = (ent2_c + nu_c) - logn - mlngamma2 ;


	// now let's try to figure out if we have overstepped the Kepler-limit
	// (FIXME) we assume that the enthalpy of the _outer_ fluid being negative
	// inside the star 
	kepler = false;
	for (int l=0; l<nzet; l++) {
	  int imax = mg->get_nr(l) - 1 ;
	  if (l == l_b) imax-- ;	// The surface point is skipped
	  for (int i=0; i<imax; i++) { 
	    if ( (*outer_ent_p)()(l, 0, j_b, i) < 0. ) {
	      kepler = true;
	      cout << "(outer) ent < 0 for l, i : " << l << "  " << i 
		   << "   ent = " << (*outer_ent_p)()(l, 0, j_b, i) << endl ;  
	    } 
	  }
	}
	
	if ( kepler ) 
	  {
	    cout << "**** KEPLERIAN VELOCITY REACHED ****" << endl ; 
	    if (kepler_fluid & 0x01)
	      omega /= kepler_factor ;    // Omega is decreased
	    if (kepler_fluid & 0x02)
	      omega2 /= kepler_factor;

	    cout << "New rotation frequencies : " 
		 << "Omega = " << omega/(2.*M_PI) * f_unit << " Hz; " 
		 << "Omega2 = " << omega2/(2.*M_PI) * f_unit << endl ; 

	    too_fast = true;
	  }

      } /* while kepler */

	
    if ( too_fast ) 
      {	// fact_omega is decreased for the next step 
	kepler_factor = sqrt( kepler_factor ) ; 
	cout << "**** New fact_omega : " << kepler_factor << endl ; 
      }
    // ============================================================


    // Cusp-check: shall the adaptation still be performed?
    // ------------------------------------------
    double dent_eq = (*outer_ent_p)().dsdr()(l_b, k_b, j_b, i_b) ; 
    double dent_pole = (*outer_ent_p)().dsdr()(l_b, k_b, 0, i_b) ;
    double rap_dent = fabs( dent_eq / dent_pole ) ; 
    cout << "| dH/dr_eq / dH/dr_pole | = " << rap_dent << endl ; 
    
    if ( rap_dent < thres_adapt ) {
      adapt_flag = 0 ;	// No adaptation of the mapping 
      cout << "******* FROZEN MAPPING  *********" << endl ; 
    }
    
    // Rescaling of the grid and adaption to (outer) enthalpy surface
    //---------------------------------------
    if (adapt_flag  && (nzadapt > 0) ) 
      {
	mp_prev = mp_et ; 
	
	mp.adapt( (*outer_ent_p)(), par_adapt) ; 
	
	mp_prev.homothetie(alpha_r) ;
	
	mp.reevaluate(&mp_prev, nzet+1, ent.set()) ; 
	mp.reevaluate(&mp_prev, nzet+1, ent2.set()) ; 
      }
    else
      mp.homothetie (alpha_r);


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

      mp.poisson2d(source_tggg(), mp.cmp_zero(), par_poisson_tggg, tggg.set()) ; 
	    
      //-------------------------------------------------------
      //	2-D Poisson equation for dzeta
      //-------------------------------------------------------

      mp.poisson2d(source_dzf(), source_dzq(), par_poisson_dzeta, dzeta.set()) ; 
	    
      err_grv2 = lbda_grv2 - 1; 
      cout << "GRV2: " << err_grv2 << endl ; 
	    
    }
    else {
      err_grv2 = grv2() ; 
    }


    //---------------------------------------
    // Computation of the metric coefficients (except for N^phi)
    //---------------------------------------

    // Relaxations on nu and dzeta :  

	if (mer >= 10) {
    	logn = relax * logn + relax_prev * logn_prev ;

    	dzeta = relax * dzeta + relax_prev * dzeta_prev ; 
	}

    // Update of the metric coefficients N, A, B and computation of K_ij :

    update_metric() ; 
	
    //-----------------------
    //  Informations display
    //-----------------------

    //    partial_display(cout) ; 
    fichfreq << "  " << omega / (2*M_PI) * f_unit ; 
    fichfreq << "  " << omega2 / (2*M_PI) * f_unit ; 
    fichevol << "  " << ray_pole() / ray_eq() ; 
    fichevol << "  " << ent_c ; 
    fichevol << "  " << ent2_c ; 


    //-----------------------------------------
    // Convergence towards given baryon masses  (if mer_mass > 0)
    //-----------------------------------------
    
    if ((mer_mass>0) && (mer > mer_mass)) {
      
      double xx, xprog, ax, fact; 

      // fluid 1
      xx = mass_b1() / mbar1_wanted - 1. ;
      cout << "Discrep. baryon mass1 <-> wanted bar. mass1 : " << xx << endl ; 

      xprog = ( mer > 2*mer_mass) ? 1. : double(mer - mer_mass)/double(mer_mass) ; 
      xx *= xprog ; 
      ax = 0.5 * ( 2. + xx ) / (1. + xx ) ; 
      fact = pow(ax, aexp_mass) ; 
      cout << "Fluid1:  xprog, xx, ax, fact : " << xprog << "  " << xx << "  " << ax << "  " << fact << endl ; 
      ent_c *= fact ; 

      // fluid 2
      xx = mass_b2() / mbar2_wanted - 1. ;
      cout << "Discrep. baryon mass2 <-> wanted bar. mass2 : " << xx << endl ; 

      xprog = ( mer > 2*mer_mass) ? 1. : double(mer - mer_mass)/double(mer_mass) ; 
      xx *= xprog ; 
      ax = 0.5 * ( 2. + xx ) / (1. + xx ) ; 
      fact = pow(ax, aexp_mass) ; 
      cout << "Fluid2: xprog, xx, ax, fact : " << xprog << "  " << xx << "  " << ax << "  " << fact << endl ; 
      ent2_c *= fact ; 

    } /* if mer > mer_mass */
	




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

