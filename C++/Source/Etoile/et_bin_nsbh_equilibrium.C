/*
 *  Method of class Etoile to compute a static spherical configuration
 *   of a neutron star in a NS-BH binary system.
 *
 *    (see file etoile.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Keisuke Taniguchi
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

char et_bin_nsbh_equilibrium_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.13  2014/10/13 08:52:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.12  2014/10/06 15:13:08  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.11  2008/09/26 08:38:45  p_grandclement
 * get rid of desaliasing
 *
 * Revision 1.10  2006/09/05 13:39:45  p_grandclement
 * update of the bin_ns_bh project
 *
 * Revision 1.9  2006/06/01 12:47:53  p_grandclement
 * update of the Bin_ns_bh project
 *
 * Revision 1.8  2006/04/25 07:21:58  p_grandclement
 * Various changes for the NS_BH project
 *
 * Revision 1.7  2006/03/30 07:33:47  p_grandclement
 * *** empty log message ***
 *
 * Revision 1.6  2005/10/18 13:12:33  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.5  2005/08/29 15:10:17  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.4  2004/03/25 10:29:04  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.3  2003/10/24 12:34:06  k_taniguchi
 * Change the notation as it should be
 *
 * Revision 1.2  2003/10/21 11:49:33  k_taniguchi
 * Change the class from Etoile_bin to sub-class Et_bin_nsbh.
 *
 * Revision 1.1  2003/10/20 15:01:55  k_taniguchi
 * Computation of an equilibrium configuration of a neutron star
 * in a NS-BH binary system.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "etoile.h"
#include "map.h"
#include "nbr_spx.h"
#include "et_bin_nsbh.h"
#include "param.h"

#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {
void Et_bin_nsbh::equilibrium_nsbh(bool adapt, double ent_c, int& niter, int mermax,
				   int mermax_poisson, double relax_poisson,
				   int mermax_potvit, double relax_potvit,
				   Tbl& diff) {
				   
    // Fundamental constants and units
    // -------------------------------
  using namespace Unites ;

    // Initializations
    // --------------
   
    const Mg3d* mg = mp.get_mg() ;
    int nz = mg->get_nzone() ;	    // total number of domains

    // The following is required to initialize mp_prev as a Map_et:
    Map_et& mp_et = dynamic_cast<Map_et&>(mp) ;

    // Error indicators
    // ----------------
    double& diff_ent = diff.set(0) ;
    double& diff_vel_pot = diff.set(1) ;
    double& diff_lapse = diff.set(2) ;
    double& diff_confpsi = diff.set(3) ;
    double& diff_shift_x = diff.set(4) ;
    double& diff_shift_y = diff.set(5) ;
    double& diff_shift_z = diff.set(6) ;

   
    // Parameters for the function Map_et::poisson for n_auto
    // ------------------------------------------------------
    double precis_poisson = 1.e-16 ;

    Param par_poisson1 ;
    par_poisson1.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson1.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson1.add_double(precis_poisson, 1) ; // required precision
    par_poisson1.add_int_mod(niter, 0) ; // number of iterations actually used
    par_poisson1.add_cmp_mod( ssjm1_lapse ) ;

    // Parameters for the function Map_et::poisson for confpsi_auto
    // ------------------------------------------------------------

    Param par_poisson2 ;
    par_poisson2.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par_poisson2.add_double(relax_poisson,  0) ; // relaxation parameter
    par_poisson2.add_double(precis_poisson, 1) ; // required precision
    par_poisson2.add_int_mod(niter, 0) ; // number of iterations actually used
    par_poisson2.add_cmp_mod( ssjm1_confpsi ) ;

    // Parameters for the function Tenseur::poisson_vect
    // -------------------------------------------------

    Param par_poisson_vect ;
    par_poisson_vect.add_int(mermax_poisson, 0) ; 
                                            // maximum number of iterations
    par_poisson_vect.add_double(relax_poisson, 0) ; // relaxation parameter
    par_poisson_vect.add_double(precis_poisson, 1) ; // required precision
    par_poisson_vect.add_cmp_mod( ssjm1_khi ) ;
    par_poisson_vect.add_tenseur_mod( ssjm1_wshift ) ;
    par_poisson_vect.add_int_mod(niter, 0) ;
    
    // Parameters for the adaptation
    Param par_adapt ; 
    int nitermax = 100 ;  
    int niter_adapt ; 
    int adapt_flag = (adapt) ? 1 : 0 ;  
    int nz_search = nzet + 1 ;
    double precis_secant = 1.e-14 ; 
    double alpha_r ; 
    double reg_map = 1. ;   
    int k_b ;
    int j_b ; 
    Tbl ent_limit(nzet) ; 
    
    par_adapt.add_int(nitermax, 0) ; 
    par_adapt.add_int(nzet, 1) ;
    par_adapt.add_int(nz_search, 2) ;	
    par_adapt.add_int(adapt_flag, 3) ;
    par_adapt.add_int(j_b, 4) ;
    par_adapt.add_int(k_b, 5) ; 
    par_adapt.add_int_mod(niter_adapt, 0) ; 
    par_adapt.add_double(precis_secant, 0) ; 
    par_adapt.add_double(reg_map, 1)	; 
    par_adapt.add_double(alpha_r, 2) ;
    par_adapt.add_tbl(ent_limit, 0) ;
    
    // External potential
    // See Eq (99) from Gourgoulhon et al. (2001)
    // -----------------------------------------
   

    Tenseur ent_jm1 = ent ;	// Enthalpy at previous step
    Tenseur source(mp) ;    // source term in the equation for logn_auto
			    // and beta_auto
    Tenseur source_shift(mp, 1, CON, ref_triad) ;  // source term in the
						   //  equation for shift_auto

    //=========================================================================
    // 			Start of iteration
    //=========================================================================

    for(int mer=0 ; mer<mermax ; mer++ ) {

	cout << "-----------------------------------------------" << endl ;
	cout << "step: " << mer << endl ;
	cout << "diff_ent = " << diff_ent << endl ;
	//-----------------------------------------------------
	// Resolution of the elliptic equation for the velocity
	// scalar potential
	//-----------------------------------------------------

	if (irrotational) {
	    diff_vel_pot = velocity_potential(mermax_potvit, precis_poisson, 
					      relax_potvit) ; 
	}

	// Equation de la surface
	//if (adapt) {
	
	   // Rescaling of the radius : (Be carefull !)
	   int nt = mg->get_nt(nzet-1) ;
	   int np = mg->get_np(nzet-1) ;
	   int nr = mg->get_nr(nzet-1) ;
	  
	   // valeurs au centre
	   double hc = exp(ent_c) ;
	   double gamma_c = exp(loggam())(0,0,0,0) ;
	   double gamma_0_c = exp(-pot_centri())(0,0,0,0) ;
	   double n_auto_c = n_auto()(0,0,0,0) ;
	   double n_comp_c = n_comp()(0,0,0,0) ;   
	 
	   double alpha_square = 0 ;
	   double constante = 0;
	   for (int k=0; k<np; k++) {
	    for (int j=0; j<nt; j++) {
		
		  // valeurs au bord
	          double gamma_b = exp(loggam())(nzet-1,k,j,nr-1) ;
	          double gamma_0_b = exp(-pot_centri())(nzet-1,k,j,nr-1) ;
	          double n_auto_b = n_auto()(nzet-1,k,j,nr-1) ;
	          double n_comp_b = n_comp()(nzet-1,k,j,nr-1) ;
  
	   // Les solutions :
	   double alpha_square_courant = (gamma_0_c*gamma_b*n_comp_b - hc*gamma_c*gamma_0_b*n_comp_c) /
	                         (hc*gamma_c*gamma_0_b*n_auto_c-gamma_0_c*gamma_b*n_auto_b) ;
	   double constante_courant = gamma_b*(n_comp_b+alpha_square_courant*n_auto_b)/gamma_0_b ;

		if (alpha_square_courant > alpha_square) {
		    alpha_square = alpha_square_courant ; 
		    k_b = k ; 
		    j_b = j ;
		    constante = constante_courant ; 
		}
	    }
	}

	   alpha_r = sqrt(alpha_square) ;
	   cout << "Adaptation : " << k_b << " " << j_b << " " << alpha_r << endl ;
	    
	   // Le potentiel : 
	   Tenseur potentiel (constante*exp(-loggam-pot_centri)/(n_auto*alpha_square+n_comp)) ;
	   potentiel.set_std_base() ;
	   for (int l=nzet+1 ; l<nz ; l++)
	    	potentiel.set().va.set(l) = 1 ;

	   Map_et mp_prev = mp_et ; 
	   ent = log(potentiel) ;
	   ent.set_std_base() ;
           ent().va.smooth(nzet, (ent.set().va)) ;
	
	   ent_limit.set_etat_qcq() ; 
	   for (int l=0; l<nzet; l++) {	// loop on domains inside the star
	       ent_limit.set(l) = ent()(l, k_b, j_b, nr-1) ; 
	   } 
	   
	   // On adapte :  
	   mp.adapt(ent(), par_adapt, 4) ;
	   mp_prev.homothetie(alpha_r) ;

	   for (int l=nzet ; l<nz-1 ; l++)
	     mp.resize(l, 1./alpha_r) ;
	   mp.reevaluate_symy (&mp_prev, nzet, ent.set()) ;
	   


	// Equation of state
	//----------------------------------------------------
	equation_of_state() ; 	// computes new values for nbar (n), ener (e)
				// and press (p) from the new ent (H)

	//---------------------------------------------------------
	// Matter source terms in the gravitational field equations
	//---------------------------------------------------------
	hydro_euler() ;	// computes new values for ener_euler (E),
				// s_euler (S) and u_euler (U^i)
				
				
	//-------------------------------------------------
	//  Relative change in enthalpy
	//-------------------------------------------------

	Tbl diff_ent_tbl = diffrel( ent(), ent_jm1() ) ;
	diff_ent = diff_ent_tbl(0) ;
	for (int l=1; l<nzet; l++) {
	    diff_ent += diff_ent_tbl(l) ;
	}
	diff_ent /= nzet ;

	ent_jm1 = ent ;
	
	
	//--------------------------------------------------------
	// Poisson equation for n_auto
	//--------------------------------------------------------

	// Source 
	// See Eq (50) from Gourgoulhon et al. (2001)
	// ------------------------------------------

	Tenseur confpsi_q = pow(confpsi, 4.) ;
	Tenseur confpsi_c = pow(confpsi, 5.) ;
	
	if (relativistic) {
	    Tenseur tmp = flat_scalar_prod(tkij_tot, tkij_auto) ;
	    Tenseur kk (mp) ;
	    kk = 0 ;
	    Tenseur tmp2(mp) ;
	    tmp2.set_etat_qcq() ;
	    for (int i=0 ; i<3 ; i++) {
		tmp2.set() = tmp(i, i) ;
		kk = kk + tmp2 ;
	    }
	    
	    source = qpig * nnn * confpsi_q * (ener_euler + s_euler)
		+ nnn * confpsi_q * kk
		- 2.*flat_scalar_prod(d_confpsi_auto+d_confpsi_comp, d_n_auto) /
		confpsi ;
	}
	else {
	    cout <<
		"WARNING : Et_bin_nsbh is for the relativistic calculation"
		 << endl ;
	    abort() ;
	}

	source.set_std_base() ;	

	// Resolution of the Poisson equation
	// ----------------------------------
	Cmp n_auto_old (n_auto()) ;
	source().poisson(par_poisson1, n_auto.set()) ;
	n_auto.set() = n_auto() + 0.5 ; 
	
	// Difference pas précédent
	// -----------------------------------------------------
	
	Tbl tdiff_lapse = diffrel(n_auto(), n_auto_old) ;
	cout <<
	"Relative difference on n_auto : "
	<< endl ;
	for (int l=0; l<nz; l++) {
	    cout << tdiff_lapse(l) << "  " ;
	}
	cout << endl ;
	diff_lapse = max(abs(tdiff_lapse)) ; 
	    
	if (relativistic) {
         
	
	    //--------------------------------------------------------
	    // Poisson equation for confpsi_auto 
	    //--------------------------------------------------------

	    // Source
	    // See Eq (51) from Gourgoulhon et al. (2001)
	    // ------------------------------------------

	    Tenseur tmp = flat_scalar_prod(tkij_tot, tkij_auto) ;
	    Tenseur kk (mp) ;
	    kk = 0 ;
	    Tenseur tmp2(mp) ;
	    tmp2.set_etat_qcq() ;
	    for (int i=0 ; i<3 ; i++) {
		tmp2.set() = tmp(i, i) ;
		kk = kk + tmp2 ;
	    }

	    source = -0.5 * qpig * confpsi_c * ener_euler
		- 0.125 * confpsi_c * kk ;

	    source.set_std_base() ; 	
	    
	    // Resolution of the Poisson equation 
	    // ----------------------------------
	    Cmp psi_old (confpsi_auto()) ;
	    source().poisson(par_poisson2, confpsi_auto.set()) ;
	    confpsi_auto.set() = confpsi_auto() + 0.5 ; 
	    
	    
	    // Check: has the Poisson equation been correctly solved ?
	    // -----------------------------------------------------

	    Tbl tdiff_confpsi = diffrel(confpsi_auto(), psi_old) ;
	    cout << 
		"Relative difference on confpsi_auto : "
		 << endl ;
	    for (int l=0; l<nz; l++) {
		cout << tdiff_confpsi(l) << "  " ;
	    }
	    cout << endl ;
	    diff_confpsi = max(abs(tdiff_confpsi)) ;
	    
	    //--------------------------------------------------------
	    // Vector Poisson equation for shift_auto 
	    //--------------------------------------------------------

	    // Source
	    // See Eq (52) from Gourgoulhon et al. (2001)
	    // ------
	 Tenseur vtmp = d_n_auto -6. * nnn * d_confpsi_auto / confpsi ;
	    source_shift = 4.*qpig * nnn *confpsi_q *(ener_euler + press)
	       * u_euler ;
        if (tkij_tot.get_etat() != ETATZERO)	    
	source_shift = source_shift + 2.* flat_scalar_prod(tkij_tot, vtmp) ;
	    source_shift.set_std_base() ;
	    // Resolution of the Poisson equation 
	    // ----------------------------------
	    // Filter for the source of shift vector 
	for (int i=0 ; i<3 ; i++)
         if (source_shift(i).get_etat() !=  ETATZERO)
           source_shift.set(i).va.coef_i() ;	   

for (int i=0; i<3; i++)
	      if ((source_shift(i).get_etat() != ETATZERO) && (source_shift(i).va.c->t[nz-1]->get_etat() != ETATZERO)) 
		source_shift.set(i).filtre(4) ;
	    for (int i=0; i<3; i++) {
		if(source_shift(i).dz_nonzero()) {
		    assert( source_shift(i).get_dzpuis() == 4 ) ;
		}
		else{
		    (source_shift.set(i)).set_dzpuis(4) ;
		}
	    }
	    //##
	    // source_shift.dec2_dzpuis() ;    // dzpuis 4 -> 2

	    double lambda_shift = double(1) / double(3) ;
	    // ON DOIT CHANGER DE TRIADE
	    source_shift.change_triad(mp.get_bvect_cart()) ;
	    Tenseur shift_old (shift_auto) ;  
	    source_shift.poisson_vect_oohara(lambda_shift, par_poisson_vect,
				      shift_auto, khi_shift) ;
	   shift_auto.change_triad(ref_triad) ;

	    // Check: has the equation for shift_auto been correctly solved ?
	    // --------------------------------------------------------------

	   

	    Tbl tdiff_shift_x = diffrel(shift_auto(0), shift_old(0)) ;
	    Tbl tdiff_shift_y = diffrel(shift_auto(1), shift_old(1)) ;
	    Tbl tdiff_shift_z = diffrel(shift_auto(2), shift_old(2)) ;

	    cout <<
		"Relative difference on shift_auto : "
		 << endl ; 
	    cout << "x component : " ;
	    for (int l=0; l<nz; l++) {
		cout << tdiff_shift_x(l) << "  " ;
	    }
	    cout << endl ;
	    cout << "y component : " ;
	    for (int l=0; l<nz; l++) {
		cout << tdiff_shift_y(l) << "  " ;
	    }
	    cout << endl ;
	    cout << "z component : " ;
	    for (int l=0; l<nz; l++) {
		cout << tdiff_shift_z(l) << "  " ;
	    }
	    cout << endl ;

	    diff_shift_x = max(abs(tdiff_shift_x)) ;
	    diff_shift_y = max(abs(tdiff_shift_y)) ;
	    diff_shift_z = max(abs(tdiff_shift_z)) ;
        } // End of relativistic equations

	
    } // End of main loop

    //=========================================================================
    // 			End of iteration
    //=========================================================================
}

// Truc pourri
void Et_bin_nsbh::equilibrium_nsbh (double, int, int, double,
				  int, double, double, const Tbl&, Tbl&) {
				  
		cout << "Not implemented !" << endl ;
		abort() ;
}
}
