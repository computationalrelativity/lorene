/*
 * Function Etoile_rot::equilibrium
 *
 * (see file etoile.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.9  2016/12/05 16:17:54  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.8  2014/10/13 08:52:57  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:13:09  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2005/10/05 15:15:31  j_novak
 * Added a Param* as parameter of Etoile_rot::equilibrium
 *
 * Revision 1.5  2004/03/25 10:29:06  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2003/11/19 21:30:57  e_gourgoulhon
 * -- Relaxation on logn and dzeta performed only if mer >= 10.
 * -- err_grv2 is now evaluated also in the Newtonian case
 *
 * Revision 1.3  2003/10/27 10:54:43  e_gourgoulhon
 * Changed local variable name lambda_grv2 to lbda_grv2 in order not
 * to shadow method name.
 *
 * Revision 1.2  2002/10/16 14:36:36  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.19  2000/11/23  15:44:10  eric
 * Ajout de l'argument ent_limit.
 *
 * Revision 2.18  2000/11/19  22:34:15  eric
 * Correction erreur ds convergence vers une masse baryonique fixee.
 *
 * Revision 2.17  2000/11/18  17:13:18  eric
 * Modifs pour permettre np=1 (axisymetrie). En particulier,
 *  lambda_shift est mis a zero si np=1.
 *
 * Revision 2.16  2000/11/10  15:17:50  eric
 * Ajout des arguments icontrol(7) (delta_mer_kep) et control(6) (precis_adapt)
 * Creation des fichiers freqeuncy.d et evolution.d
 *
 * Revision 2.15  2000/11/08  15:21:16  eric
 * Appel de fait_nphi() avant hydro_euler() pour le test sur u_euler.
 *
 * Revision 2.14  2000/10/25  15:13:36  eric
 * omega est initialise a zero.
 *
 * Revision 2.13  2000/10/23  14:02:55  eric
 * Modif de Map_et::adapt: on y rentre desormais avec nz_search
 *   dans le cas present nz_search = nzet + 1).
 *
 * Revision 2.12  2000/10/23  13:47:47  dorota
 * Ajout en sortie (dans diff(7)) de vit_triax.
 * Suppression de l'initialisation a zero de omega_ini (!)
 *
 * Revision 2.11  2000/10/20  13:56:43  eric
 * Ecriture dans le fichier convergence.d de la vitesse de developpement
 *  de la perturbation triaxiale (vit_triax).
 *
 * Revision 2.10  2000/10/20  13:11:23  eric
 * Ajout de l'argument nzadapt.
 *
 * Revision 2.9  2000/10/17  16:00:24  eric
 * Ajout de la perturbation triaxiale.
 *
 * Revision 2.8  2000/10/12  15:33:22  eric
 * Ajout de l'appel a fait_nphi() pour le calcul de tnphi et nphi.
 * Emploi de la nouvelle version de Tenseur::set_std_base() : on n'a plus
 *  besoin de tester l'etat du tenseur avant.
 *
 * Revision 2.7  2000/10/11  15:14:09  eric
 * Ajout des equations pour tggg et dzeta --> 1ere version complete !
 *
 * Revision 2.6  2000/10/06  15:06:14  eric
 * Version relativiste avec le lapse et le shift uniquement.
 * Ca converge.
 *
 * Revision 2.5  2000/09/18  16:15:26  eric
 * Premiers termes relativistes.
 *
 * Revision 2.4  2000/08/31  15:39:02  eric
 * Appel du nouvel operateur Cmp::mult_rsint pour le calcul de uuu.
 *
 * Revision 2.3  2000/08/25  12:27:55  eric
 * modifs mineures (fichconv).
 *
 * Revision 2.2  2000/08/18  14:01:24  eric
 * Premiere version operationnelle (testee en spherique Newtonien)
 *
 * Revision 2.1  2000/08/17  12:39:47  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/07/21  16:30:32  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "etoile.h"
#include "param.h"

#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
void Etoile_rot::equilibrium(double ent_c, double omega0, double fact_omega, 
			     int nzadapt, const Tbl& ent_limit, const Itbl& icontrol,
			     const Tbl& control, double mbar_wanted, 
			     double aexp_mass, Tbl& diff, Param*) {
			     
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
    int mer_triax = icontrol(6) ; 
    int delta_mer_kep = icontrol(7) ; 

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
    double ampli_triax = control(5) ; 
    double precis_adapt = control(6) ; 


    // Error indicators
    // ----------------
    
    diff.set_etat_qcq() ; 
    double& diff_ent = diff.set(0) ; 
    double& diff_nuf = diff.set(1) ; 
    double& diff_nuq = diff.set(2) ; 
//    double& diff_dzeta = diff.set(3) ; 
//    double& diff_ggg = diff.set(4) ; 
    double& diff_shift_x = diff.set(5) ; 
    double& diff_shift_y = diff.set(6) ; 
    double& vit_triax = diff.set(7) ; 
    
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

    // Initial angular velocity
    omega = 0 ; 
  
    double accrois_omega = (omega0 - omega_ini) / 
			    double(mer_fix_omega - mer_change_omega) ; 


    update_metric() ;	// update of the metric coefficients

    equation_of_state() ;	// update of the density, pressure, etc...
    
    hydro_euler() ;	// update of the hydro quantities relative to the 
			//  Eulerian observer

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
    fichconv << "#     diff_ent     GRV2       max_triax      vit_triax" << endl ; 
    
    ofstream fichfreq("frequency.d") ;    // Output file for  omega
    fichfreq << "#       f [Hz]" << endl ; 
    
    ofstream fichevol("evolution.d") ;    // Output file for various quantities
    fichevol << 
    "#       |dH/dr_eq/dH/dr_pole|      r_pole/r_eq	ent_c" 
    << endl ; 
    
    diff_ent = 1 ; 
    double err_grv2 = 1 ; 
    double max_triax_prev = 0 ;	    // Triaxial amplitude at previous step
    
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
	    source_nuf = qpig * nbar ; 

	    source_nuq = 0 ; 
	}
	source_nuf.set_std_base() ; 	
	source_nuq.set_std_base() ; 	

	// Source for dzeta
	// ----------------
	source_dzf = 2 * qpig * a_car * (press + (ener_euler+press) * uuu*uuu ) ;
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
	    
	// Matter term: 
	source_shift = (-4*qpig) * nnn * a_car  * (ener_euler + press)
				* u_euler ;

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
	
	//---------------------------------------
	// Triaxial perturbation of nuf
	//---------------------------------------
	
	if (mer == mer_triax) {
	
	    if ( mg->get_np(0) == 1 ) {
		cout << 
		"Etoile_rot::equilibrium: np must be stricly greater than 1"
		<< endl << " to set a triaxial perturbation !" << endl ; 
		abort() ; 
	    }
	
	    const Coord& phi = mp.phi ; 
	    const Coord& sint = mp.sint ; 
	    Cmp perturb(mp) ; 
	    perturb = 1 + ampli_triax * sint*sint * cos(2*phi) ;
	    nuf.set() = nuf() * perturb ;  
	    
	    nuf.set_std_base() ;    // set the bases for spectral expansions
				    //  to be the standard ones for a 
				    //  scalar field

	}
	
	// Monitoring of the triaxial perturbation
	// ---------------------------------------
	
	Valeur& va_nuf = nuf.set().va ; 
	va_nuf.coef() ;		// Computes the spectral coefficients
	double max_triax = 0 ; 

	if ( mg->get_np(0) > 1 ) {

	    for (int l=0; l<nz; l++) {	// loop on the domains
		for (int j=0; j<mg->get_nt(l); j++) {
		    for (int i=0; i<mg->get_nr(l); i++) {
		
			// Coefficient of cos(2 phi) : 
			double xcos2p = (*(va_nuf.c_cf))(l, 2, j, i) ; 

			// Coefficient of sin(2 phi) : 
			double xsin2p = (*(va_nuf.c_cf))(l, 3, j, i) ; 
		    
			double xx = sqrt( xcos2p*xcos2p + xsin2p*xsin2p ) ; 
		    
			max_triax = ( xx > max_triax ) ? xx : max_triax ; 		    
		    }
		} 
	    }
	    
	}
	
	cout << "Triaxial part of nuf : " << max_triax << endl ; 

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
	    //##
	    // source_shift.dec2_dzpuis() ;    // dzpuis 4 -> 2

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
	
	    // Equatorial values of various potentials :
	    double nuf_b  = nuf()(l_b, k_b, j_b, i_b) ; 
	    double nuq_b  = nuq()(l_b, k_b, j_b, i_b) ; 
	    double mlngamma_b  = mlngamma()(l_b, k_b, j_b, i_b) ; 

	    // Central values of various potentials :
	    double nuf_c = nuf()(0,0,0,0) ; 
	    double nuq_c = nuq()(0,0,0,0) ; 
	    double mlngamma_c = 0 ;
	
	    // Scale factor to ensure that the enthalpy is equal to ent_b at 
	    //  the equator
	    double alpha_r2 = ( ent_c - ent_b + mlngamma_c - mlngamma_b
				+ nuq_c - nuq_b) / ( nuf_b - nuf_c  ) ;
	    alpha_r = sqrt(alpha_r2) ;
	    cout << "alpha_r = " << alpha_r << endl ; 

	    // Readjustment of nu :
	    // -------------------

	    logn = alpha_r2 * nuf + nuq ;
	    double nu_c =  logn()(0,0,0,0) ;

	    // First integral	--> enthalpy in all space
	    //-----------------

	    ent = (ent_c + nu_c + mlngamma_c) - logn - mlngamma ;

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

//##	if (mer >= mer_triax) {
//	    des_coupe_y(ent(), 0., 1, "ent before adapt", &(ent()) ) ; 
//	    des_coupe_z(ent(), 0., 1, "ent before adapt (EQUAT)", &(ent()) ) ; 
//##	}

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

//##	if (mer >= mer_triax) {
//	    des_coupe_y(ent(), 0., 1, "ent after adapt", &(ent()) ) ; 
//	    des_coupe_z(ent(), 0., 1, "ent after adapt (EQUAT)", &(ent()) ) ; 
//##	}

 	//----------------------------------------------------
	// Computation of the enthalpy at the new grid points
	//----------------------------------------------------
	
	mp_prev.homothetie(alpha_r) ; 
	
	mp.reevaluate(&mp_prev, nzet+1, ent.set()) ; 

//##	if (mer >= mer_triax) {
//	    des_coupe_y(ent(), 0., 1, "ent after reevaluate", &(ent()) ) ; 
//	    des_coupe_z(ent(), 0., 1, "ent after reevaluate (EQUAT)", &(ent()) ) ; 
//##	}


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

	partial_display(cout) ; 
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
	fichconv << "  " << log10( fabs(max_triax) + 1.e-16 ) ;

	vit_triax = 0 ;
	if ( (mer > mer_triax+1) && (max_triax_prev > 1e-13) ) {
		vit_triax = (max_triax - max_triax_prev) / max_triax_prev ; 
	}

	fichconv << "  " << vit_triax ;
	
	//------------------------------
	//  Recycling for the next step
	//------------------------------
	
	ent_prev = ent ; 
	logn_prev = logn ; 
	dzeta_prev = dzeta ; 
	max_triax_prev = max_triax ; 
	
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
}
