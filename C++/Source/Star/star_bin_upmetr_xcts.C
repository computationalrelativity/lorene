/*
 * Methods of Star_bin_xcts::update_metric
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

char star_bin_upmetr_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2010/06/15 08:13:01  m_bejger
 * Some more corrections: Psi, chi
 *
 * Revision 1.2  2010/06/04 20:01:59  m_bejger
 * Corrected definitions of lapse, Psi4; added definition of gamma
 *
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// Headers Lorene
#include "cmp.h"
#include "star.h"
#include "graphique.h"
#include "utilitaires.h"

//----------------------------------//
//	 Version without relaxation     //
//----------------------------------//

void Star_bin_xcts::update_metric(const Star_bin_xcts& comp) {

    // Computation of quantities coming from the companion
    // ---------------------------------------------------
    
    const Map& mp_comp (comp.get_mp()) ;
     
    if ( (comp.Psi_auto).get_etat() == ETATZERO ) {	    
		Psi_comp.set_etat_zero() ;
    
    } else {	
		Psi_comp.set_etat_qcq() ;
		Psi_comp.import(comp.Psi_auto) ;
		Psi_comp.std_spectral_base() ;
    } 
  
    beta_comp.set_etat_qcq() ; 
    beta_comp.set_triad(mp.get_bvect_cart()) ;
    
    Vector comp_beta(comp.beta_auto) ;
    comp_beta.change_triad(mp_comp.get_bvect_cart()) ;
    comp_beta.change_triad(mp.get_bvect_cart()) ;

    assert ( *(beta_comp.get_triad()) == *(comp_beta.get_triad())) ;

    (beta_comp.set(1)).import( comp_beta(1) ) ;  
    (beta_comp.set(2)).import( comp_beta(2) ) ;  
    (beta_comp.set(3)).import( comp_beta(3) ) ;  

    beta_comp.std_spectral_base() ;   

    if ( (comp.chi_auto).get_etat()  == ETATZERO ) {
		chi_comp.set_etat_zero() ;

    } else	{
		chi_comp.set_etat_qcq() ;  
		chi_comp.import( comp.chi_auto ) ;		
		chi_comp.std_spectral_base() ;
    } 

// Conformal factor Psi
// --------------------

    Psi = Psi_auto * Psi_comp ; 
    Psi.std_spectral_base() ; 

    Scalar psi4 = pow(Psi_auto*Psi_comp, 4.) ; 
    psi4.std_spectral_base() ; 

// Function chi = NPsi
// --------------------

    chi = chi_auto * chi_comp ; 
    chi.std_spectral_base() ; 
   
// Lapse function N
// ----------------

    nn = chi_auto*chi_comp / ( Psi_auto*Psi_comp ) ; 
    nn.std_spectral_base() ; 

// Shift vector 
// -------------

    beta = beta_auto + beta_comp ;
    
    Sym_tensor tens_gamma(( flat.con() ) / psi4) ;
    gamma = tens_gamma ;

// Extrinsic curvature (haij_auto and hacar_auto)
//-----------------------------------------------

    extrinsic_curvature() ;

   
// The derived quantities are obsolete
// -----------------------------------

    del_deriv() ;
    
}

//----------------------------------//
//	  Version with relaxation       //
//----------------------------------//

void Star_bin_xcts::update_metric(const Star_bin_xcts& comp,
			     				  const Star_bin_xcts& star_jm1, 
			     				  double relax) {


    // Computation of quantities coming from the companion
    // ---------------------------------------------------
  
    const Map& mp_comp (comp.get_mp()) ;

    if ( (comp.Psi_auto).get_etat() == ETATZERO ) {
		Psi_comp.set_etat_zero() ;
    
    } else {

		Psi_comp.set_etat_qcq() ;
		Psi_comp.import( comp.Psi_auto ) ;
		Psi_comp.std_spectral_base() ;

    }

    beta_comp.set_etat_qcq() ; 
    beta_comp.set_triad(mp.get_bvect_cart()) ;

    Vector comp_beta(comp.beta_auto) ;
    comp_beta.change_triad(mp_comp.get_bvect_cart()) ;
    comp_beta.change_triad(mp.get_bvect_cart()) ;

    assert ( *(beta_comp.get_triad()) == *(comp_beta.get_triad())) ;

    (beta_comp.set(1)).import( comp_beta(1) ) ;  
    (beta_comp.set(2)).import( comp_beta(2) ) ;  
    (beta_comp.set(3)).import( comp_beta(3) ) ;  

    beta_comp.std_spectral_base() ;   
 
    if ( (comp.chi_auto).get_etat()  == ETATZERO ) {
    	chi_comp.set_etat_zero() ;

    } else {

		chi_comp.set_etat_qcq() ;
		chi_comp.import( comp.chi_auto ) ;
 		chi_comp.std_spectral_base() ;
 		
   }	
  
// Relaxation on Psi_comp, beta_comp, chi_comp
// -------------------------------------------
    double relaxjm1 = 1. - relax ; 
    
    Psi_comp = relax * Psi_comp 
    			+ relaxjm1 * (star_jm1.Psi_comp) ; 
    
    beta_comp = relax * beta_comp 
    			+ relaxjm1 * (star_jm1.beta_comp) ; 

    chi_comp = relax * chi_comp 
    			+ relaxjm1 * (star_jm1.chi_comp) ;

// Conformal factor Psi
// --------------------

    Psi = Psi_auto * Psi_comp ; 
    Psi.std_spectral_base() ; 

    Scalar psi4 = pow(Psi_auto*Psi_comp, 4.) ; 
    psi4.std_spectral_base() ; 

// Function chi = NPsi
// --------------------

    chi = chi_auto * chi_comp ; 
    chi.std_spectral_base() ; 
   
// Lapse function N
// ----------------

    nn = chi_auto*chi_comp / ( Psi_auto*Psi_comp ) ; 
    nn.std_spectral_base() ; 

// Shift vector
// ------------
	    
    beta = beta_auto + beta_comp ;
        
    Sym_tensor tens_gamma((flat.con()) / psi4) ;
    gamma = tens_gamma ;

    // Extrinsic curvature (haij_auto and hacar_auto)
    extrinsic_curvature() ;

   
// The derived quantities are obsolete
// -----------------------------------

    del_deriv() ;

}

