/*
 * Function et_rot_mag::equilibrium_mag
 *
 * Computes rotating equilibirum with a magnetic field
 * (see file et_rot_mag.h for documentation)
 *
 */

/*
 *   Copyright (c) 2002 Eric Gourgoulhon
 *   Copyright (c) 2002 Emmanuel Marcq
 *   Copyright (c) 2002 Jerome Novak
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

char et_rot_mag_equil_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.14  2003/10/03 15:58:47  j_novak
 * Cleaning of some headers
 *
 * Revision 1.13  2002/10/16 14:36:36  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.12  2002/10/11 11:47:35  j_novak
 * Et_rot_mag::MHD_comput is now virtual.
 * Use of standard constructor for Tenseur mtmp in Et_rot_mag::equilibrium_mag
 *
 * Revision 1.11  2002/06/05 15:15:59  j_novak
 * The case of non-adapted mapping is treated.
 * parmag.d and parrot.d have been merged.
 *
 * Revision 1.10  2002/06/03 13:23:16  j_novak
 * The case when the mapping is not adapted is now treated
 *
 * Revision 1.9  2002/06/03 13:00:45  e_marcq
 *
 * conduc parameter read in parmag.d
 *
 * Revision 1.6  2002/05/17 15:08:01  e_marcq
 *
 * Rotation progressive plug-in, units corrected, Q and a_j new member data
 *
 * Revision 1.5  2002/05/16 10:02:09  j_novak
 * Errors in stress energy tensor corrected
 *
 * Revision 1.4  2002/05/15 09:53:59  j_novak
 * First operational version
 *
 * Revision 1.3  2002/05/14 13:38:36  e_marcq
 *
 *
 * Unit update, new outputs
 *
 * Revision 1.1  2002/05/10 09:26:52  j_novak
 * Added new class Et_rot_mag for magnetized rotating neutron stars (under development)
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "et_rot_mag.h"
#include "param.h"

void Et_rot_mag::equilibrium_mag(double ent_c, double omega0, 
     double fact_omega, int nzadapt, const Tbl& ent_limit, 
     const Itbl& icontrol, const Tbl& control, double mbar_wanted, 
     double aexp_mass, Tbl& diff, const double Q0, const double a_j0, 
     Cmp (*f_j)(const Cmp&, const double), 
     Cmp (*M_j)(const Cmp& x, const double)) {
			     
    // Fundamental constants and units
    // -------------------------------
    #include "unites_mag.h"	    
    // To avoid some compilation warnings
    if (ent_c < 0) {
	cout << f_unit << msol << km << mevpfm3 << endl ; 
	cout << mag_unit << elec_unit << endl ;
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
    
    // The following is required to initialize mp_prev as a Map_et:
    Map_et& mp_et = dynamic_cast<Map_et&>(mp) ; 
        
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
    int mer_mass = icontrol(4) ; 
    int mermax_poisson = icontrol(5) ; 
    int delta_mer_kep = icontrol(6) ; 
    int mer_mag = icontrol(7) ;
    int mer_change_mag = icontrol(8) ;
    int mer_fix_mag = icontrol(9) ;

    // Protections:
    if (mer_change_omega < mer_rot) {
	cout << "Etoile_rot::equilibrium: mer_change_omega < mer_rot !" << endl ;
	cout << " mer_change_omega = " << mer_change_omega << endl ; 
	cout << " mer_rot = " << mer_rot << endl ; 
	abort() ; 
    }
    if (mer_fix_omega < mer_change_omega) {
	cout << "Etoile_rot::equilibrium: mer_fix_omega < mer_change_omega !" 
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
    double relax_poisson = control(3) ; 
    double thres_adapt = control(4) ; 
    double precis_adapt = control(5) ; 
    double Q_ini = control(6) ;
    double a_j_ini = control (7) ;

    // Error indicators
    // ----------------
    
    diff.set_etat_qcq() ; 
    double& diff_ent = diff.set(0) ; 

    // Parameters for the function Map_et::adapt
    // -----------------------------------------
    
    Param par_adapt ; 
    int nitermax = 100 ;  
    int niter ; 
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

 					   
    // Parameters for the Maxwell equations
    // -------------------------------------

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

    // Initial angular velocity / magnetic quantities
    omega = 0 ; 
    Q = 0 ;
    a_j = 0 ;

    double accrois_omega = (omega0 - omega_ini) / 
			    double(mer_fix_omega - mer_change_omega) ; 
    double accrois_Q = (Q0 - Q_ini) /
                            double(mer_fix_mag - mer_change_mag);
    double accrois_a_j = (a_j0 - a_j_ini) / 
                            double(mer_fix_mag - mer_change_mag); 

    update_metric() ;	// update of the metric coefficients

    equation_of_state() ;	// update of the density, pressure, etc...
    
    hydro_euler() ;	// update of the hydro quantities relative to the 
			//  Eulerian observer

    MHD_comput() ; // update of EM contributions to stress-energy tensor


    // Quantities at the previous step : 	
    Map_et mp_prev = mp_et ; 
    Tenseur ent_prev = ent ;	    
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
    fichconv << "#     diff_ent     GRV2    " << endl ; 
    
    ofstream fichfreq("frequency.d") ;    // Output file for  omega
    fichfreq << "#       f [Hz]" << endl ; 
    
    ofstream fichevol("evolution.d") ;    // Output file for various quantities
    fichevol << 
    "#       |dH/dr_eq/dH/dr_pole|      r_pole/r_eq	ent_c" 
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
	    }
	    else {
		if (mer <= mer_fix_omega) {
		    omega = omega_ini + accrois_omega * 
			      (mer - mer_change_omega) ;
		}
	    }

	}

	if (mer >= mer_mag) {
	  if (mer < mer_change_mag) {
	    Q   = Q_ini ;
	    a_j = a_j_ini ;
	  }
	  else {
	    if (mer <= mer_fix_mag) {
	      Q = Q_ini + accrois_Q * (mer - mer_change_mag) ;
	      a_j = a_j_ini + accrois_a_j * (mer - mer_change_mag) ;
	    }
	  }
	}


	//-----------------------------------------------
	// Computation of electromagnetic potentials :
	// -------------------------------------------
	magnet_comput(adapt_flag, 
		      f_j, par_poisson_At, par_poisson_Avect) ;

	MHD_comput() ; // computes EM contributions to T_{mu,nu}

	//-----------------------------------------------
	//  Sources of the Poisson equations
	//-----------------------------------------------
	
	// Source for nu
	// -------------
	Tenseur beta = log(bbb) ; 
	beta.set_std_base() ; 

	if (relativistic) {
	  source_nuf =  qpig * a_car *( ener_euler + s_euler) ; 
	  
	  source_nuq = ak_car - flat_scalar_prod(logn.gradient_spher(), 
		       logn.gradient_spher() + beta.gradient_spher()) 
	               + qpig * a_car * 2*E_em ;
	}
	else {
	    source_nuf = qpig * nbar ; 

	    source_nuq = 0 ; 
	}
	source_nuf.set_std_base() ; 	
	source_nuq.set_std_base() ; 	

	// Source for dzeta
	// ----------------
	source_dzf = 2 * qpig * a_car * (press + (ener_euler+press) * uuu*uuu ) ;
	source_dzf.set_std_base() ; 
  
	source_dzq = 2 * qpig * a_car * E_em + 1.5 * ak_car - 
	  flat_scalar_prod(logn.gradient_spher(), logn.gradient_spher() ) ;  
	source_dzq.set_std_base() ; 	
	
	// Source for tggg
	// ---------------
	
	source_tggg = 4 * qpig * nnn * a_car * bbb * press ;
	source_tggg.set_std_base() ; 
	
	(source_tggg.set()).mult_rsint() ; 
	

	// Source for shift
	// ----------------
	    
	// Matter term: 
	
	Cmp tjpem(Jp_em()) ;
	tjpem.div_rsint() ;

	source_shift = (-4*qpig) * nnn * a_car  * (ener_euler + press)
				* u_euler ;

	// Quadratic terms:
	Tenseur vtmp =  6 * beta.gradient_spher() - 2 * logn.gradient_spher() ;
	Tenseur mtmp(mp, 1, COV, mp.get_bvect_spher()) ;
	if (tjpem.get_etat() == ETATZERO) mtmp.set_etat_zero() ;
	else {
	  mtmp.set_etat_qcq() ;
	  mtmp.set(0) = 0 ;
	  mtmp.set(1) = 0 ;
	  mtmp.set(2) = (-4*qpig)*tjpem*nnn()*a_car()/b_car() ;
	}
	mtmp.change_triad(mp.get_bvect_cart()) ; 

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
	if (mtmp.get_etat() == ETATQCQ) {
	  if (source_shift.get_etat() == ETATZERO) {
	    source_shift.set_etat_qcq() ;
	    for (int i=0; i<3; i++) {
	      source_shift.set(i) = mtmp(i) ;
	      source_shift.set(i).va.coef_i() ;
	    }
	  }
	  else
	    for (int i=0; i<3; i++) 
	      source_shift.set(i) += mtmp(i) ; 
	}

      	source_shift.set_std_base() ; 	

	//----------------------------------------------
	// Resolution of the Poisson equation for nuf 
	//----------------------------------------------

	source_nuf().poisson(par_poisson_nuf, nuf.set()) ; 
		
	if (relativistic) {
	    
	    //----------------------------------------------
	    // Resolution of the Poisson equation for nuq 
	    //----------------------------------------------

	    source_nuq().poisson(par_poisson_nuq, nuq.set()) ; 
	    
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
	    //##
	    // source_shift.dec2_dzpuis() ;    // dzpuis 4 -> 2

	    double lambda_shift = double(1) / double(3) ; 

	    if ( mg->get_np(0) == 1 ) {
		lambda_shift = 0 ; 
	    }
	
	    source_shift.poisson_vect(lambda_shift, par_poisson_vect, 
				      shift, w_shift, khi_shift) ;      
	    
	    // Computation of tnphi and nphi from the Cartesian components
	    //  of the shift
	    // -----------------------------------------------------------
	    
	    fait_nphi() ; 
	
		//##	cout.precision(10) ;
		//		cout << "nphi : " << nphi()(0, 0, 0, 0)
	    //			 << "  " << nphi()(l_b, k_b, j_b, i_b/2)
	    //			 << "  " << nphi()(l_b, k_b, j_b, i_b) << endl ;

	}

	//-----------------------------------------
	// Determination of the fluid velociy U
	//-----------------------------------------
	
	if (mer > mer_fix_omega + delta_mer_kep) {
		
    	    omega *= fact_omega ;  // Increase of the angular velocity if 
	}			   //  fact_omega != 1
	
	bool omega_trop_grand = false ; 
	bool kepler = true ; 

	while ( kepler ) {
	
	    // Possible decrease of Omega to ensure a velocity < c 
	
	    bool superlum = true ; 
	
	    while ( superlum ) {
	
		// New fluid velocity U :
	
		Cmp tmp = omega - nphi() ; 
		tmp.annule(nzm1) ; 
		tmp.std_base_scal() ;
    
		tmp.mult_rsint() ;	    //  Multiplication by r sin(theta)

		uuu = bbb() / nnn() * tmp ; 
    
		if (uuu.get_etat() == ETATQCQ) {
		    // Same basis as (Omega -N^phi) r sin(theta) :
		    ((uuu.set()).va).set_base( (tmp.va).base ) ;   
		}


		// Is the new velocity larger than c in the equatorial plane ?
	
		superlum = false ; 
	
		for (int l=0; l<nzet; l++) {
		    for (int i=0; i<mg->get_nr(l); i++) {
		
			double u1 = uuu()(l, 0, j_b, i) ; 
			if (u1 >= 1.) {	    // superluminal velocity
			    superlum = true ; 
			    cout << "U > c  for l, i : " << l << "  " << i 
				 << "   U = " << u1 << endl ;  
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
	    }	// end of while ( superlum )

	
	    // New computation of U (which this time is not superluminal)
	    //  as well as of gam_euler, ener_euler, etc...
	    // -----------------------------------

	    hydro_euler() ; 
	

	    //------------------------------------------------------
	    //	First integral of motion 
	    //------------------------------------------------------
	
	    // Centrifugal potential : 
	    if (relativistic) {
		mlngamma = - log( gam_euler ) ;
	    }
	    else {
		mlngamma = - 0.5 * uuu*uuu ; 
	    }

	    Tenseur mag(mp) ;
	    if (is_conduct()) {
	      mag = mu0*M_j(A_phi, a_j) ;}
	    else{
	      mag = mu0*M_j(omega*A_phi-A_t, a_j) ;}

	    // Equatorial values of various potentials :
	    double nuf_b  = nuf()(l_b, k_b, j_b, i_b) ; 
	    double nuq_b  = nuq()(l_b, k_b, j_b, i_b) ; 
	    double mlngamma_b  = mlngamma()(l_b, k_b, j_b, i_b) ; 
	    double mag_b = mag()(l_b, k_b, j_b, i_b) ; 

	    // Central values of various potentials :
	    double nuf_c = nuf()(0,0,0,0) ; 
	    double nuq_c = nuq()(0,0,0,0) ; 
	    double mlngamma_c = 0 ;
	    double mag_c = mag()(0,0,0,0) ;
	
	    // Scale factor to ensure that the enthalpy is equal to ent_b at 
	    //  the equator
	    double alpha_r2 = ( ent_c - ent_b + mlngamma_c - mlngamma_b
				+ nuq_c - nuq_b + mag_c - mag_b) 
	      / ( nuf_b - nuf_c  ) ;
	    alpha_r = sqrt(alpha_r2) ;
	    cout << "alpha_r = " << alpha_r << endl ; 

	    // Readjustment of nu :
	    // -------------------

	    logn = alpha_r2 * nuf + nuq ;
	    double nu_c =  logn()(0,0,0,0) ;


	    // First integral	--> enthalpy in all space
	    //-----------------
	    ent = (ent_c + nu_c + mlngamma_c + mag_c) - logn - mlngamma
	      - mag ;

	    // Test: is the enthalpy negative somewhere in the equatorial plane
	    //  inside the star ? If yes, this means that the Keplerian velocity
	    //  has been overstep.

	    kepler = false ; 
	    for (int l=0; l<nzet; l++) {
		int imax = mg->get_nr(l) - 1 ;
		if (l == l_b) imax-- ;	// The surface point is skipped
		for (int i=0; i<imax; i++) { 
		    if ( ent()(l, 0, j_b, i) < 0. ) {
			kepler = true ;
			cout << "ent < 0 for l, i : " << l << "  " << i 
			     << "   ent = " << ent()(l, 0, j_b, i) << endl ;  
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

	}   // End of while ( kepler )
	
	if ( omega_trop_grand ) {	// fact_omega is decreased for the
					//  next step 
	    fact_omega = sqrt( fact_omega ) ; 
	    cout << "**** New fact_omega : " << fact_omega << endl ; 
	}

	//----------------------------------------------------
	// Adaptation of the mapping to the new enthalpy field
	//----------------------------------------------------
    
	// Shall the adaptation be performed (cusp) ?
	// ------------------------------------------
	
	double dent_eq = ent().dsdr()(l_b, k_b, j_b, i_b) ; 
	double dent_pole = ent().dsdr()(l_b, k_b, 0, i_b) ;
	double rap_dent = fabs( dent_eq / dent_pole ) ; 
	cout << "| dH/dr_eq / dH/dr_pole | = " << rap_dent << endl ; 
	
	if ( rap_dent < thres_adapt ) {
	    adapt_flag = 0 ;	// No adaptation of the mapping 
	    cout << "******* FROZEN MAPPING  *********" << endl ; 
	}
	else{
	    adapt_flag = 1 ;	// The adaptation of the mapping is to be
				//  performed
	}

	mp_prev = mp_et ; 

	mp.adapt(ent(), par_adapt) ; 

 	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------
	
	mp_prev.homothetie(alpha_r) ; 
	
	mp.reevaluate(&mp_prev, nzet+1, ent.set()) ; 

	//----------------------------------------------------
	// Equation of state  
	//----------------------------------------------------
	
	equation_of_state() ; 	// computes new values for nbar (n), ener (e) 
				// and press (p) from the new ent (H)
	
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

	//	partial_display(cout) ; 
	fichfreq << "  " << omega / (2*M_PI) * f_unit ; 
	fichevol << "  " << rap_dent ; 
	fichevol << "  " << ray_pole() / ray_eq() ; 
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
	

	//------------------------------------------------------------
	//  Relative change in enthalpy with respect to previous step 
	//------------------------------------------------------------

	Tbl diff_ent_tbl = diffrel( ent(), ent_prev() ) ; 
	diff_ent = diff_ent_tbl(0) ; 
	for (int l=1; l<nzet; l++) {
	    diff_ent += diff_ent_tbl(l) ; 
	}
	diff_ent /= nzet ; 
	
	fichconv << "  " << log10( fabs(diff_ent) + 1.e-16 ) ;
	fichconv << "  " << log10( fabs(err_grv2) + 1.e-16 ) ;

	//------------------------------
	//  Recycling for the next step
	//------------------------------
	
	ent_prev = ent ; 
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
