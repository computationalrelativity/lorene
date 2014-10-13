/*
 * Methods of class Binaire for estimating the error in the Hamiltionian
 *  and momentum constraints
 *
 * (see file binaire.h for documentation).
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


char binaire_constr_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:52:44  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2004/03/25 10:28:59  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.1.1.1  2001/11/20 15:19:30  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  2000/03/13  17:05:34  eric
 * *** empty log message ***
 *
 * Revision 2.0  2000/03/13  14:26:08  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "binaire.h"
#include "unites.h"



		//----------------------------------------------//
		//	    Hamiltonian constraint		//
		//----------------------------------------------//

namespace Lorene {
double Binaire::ham_constr() const {
    
  using namespace Unites ;

    if (p_ham_constr == 0x0) {	    // A new computation is required
	

	Tenseur lap_alpha1( star1.get_mp() ) ; 
	Tenseur lap_alpha2( star2.get_mp() ) ; 

	Tenseur source1( star1.get_mp() ) ; 
	Tenseur source2( star2.get_mp() ) ; 

	Tenseur* p_lap_alpha[2] ; 
	Tenseur* p_source[2] ; 
	p_lap_alpha[0] = &lap_alpha1 ; 
	p_lap_alpha[1] = &lap_alpha2 ; 
	p_source[0] = &source1 ; 
	p_source[1] = &source2 ; 


	// Computation of the l.h.s. and r.h.s. of the Hamiltonian
	// constraint in each star. 
	// -------------------------------------------------------
	
	double som = 0 ; 
	
	for (int i=0; i<2; i++) {
	    
	    // Laplacian of alpha = ln(A) 
	    // --------------------------

	    Tenseur alpha_auto = et[i]->get_beta_auto()	
				 - et[i]->get_logn_auto() ; 
				 
	    *(p_lap_alpha[i]) = alpha_auto().laplacien() ; 
	    
	    // Right-hand-side of the Hamiltonian constraint
	    // ---------------------------------------------
	    
	    const Tenseur& a_car = et[i]->get_a_car() ; 
	    const Tenseur& ener_euler = et[i]->get_ener_euler() ; 

	    Tenseur d_alpha_auto = et[i]->get_d_beta_auto()	
				 - et[i]->get_d_logn_auto() ; 
	    
	    Tenseur d_alpha_comp = et[i]->get_d_beta_comp()	
				 - et[i]->get_d_logn_comp() ; 
	    
	    const Tenseur& akcar_auto = et[i]->get_akcar_auto() ; 
	    const Tenseur& akcar_comp = et[i]->get_akcar_comp() ; 
	    
	    *(p_source[i]) = - qpig * a_car * ener_euler 
			     - 0.25 * ( akcar_auto + akcar_comp ) 
			     - 0.5 * 
				(   flat_scalar_prod(d_alpha_auto, d_alpha_auto)
			          + flat_scalar_prod(d_alpha_auto, d_alpha_comp)
				) ;

	    // Relative difference
	    // -------------------
	    Tbl diff = diffrel( (*(p_lap_alpha[i]))(), (*(p_source[i]))() ) ; 
	    
	    cout << 
	    "Binaire::ham_constr : relative difference Lap(alpha) <-> source : "
	    << endl << diff << endl ; 
	    
	    som += max( abs(diff) ) ; 

	}

	
	// Total error
	// -----------
	p_ham_constr = new double ; 

	*p_ham_constr = 0.5 * som  ; 
	    
    }
        
    return *p_ham_constr ; 
    
}


		//----------------------------------------------//
		//		Momentum constraint		//
		//----------------------------------------------//

const Tbl& Binaire::mom_constr() const {

  using namespace Unites ;

    if (p_mom_constr == 0x0) {	    // A new computation is required
	
	Tenseur divk1( star1.get_mp(), 1, CON, ref_triad ) ; 
	Tenseur divk2( star2.get_mp(), 1, CON, ref_triad ) ; 

	Tenseur source1( star1.get_mp(), 1, CON, ref_triad ) ; 
	Tenseur source2( star2.get_mp(), 1, CON, ref_triad ) ; 

	Tenseur* p_divk[2] ; 
	Tenseur* p_source[2] ; 
	p_divk[0] = &divk1 ; 
	p_divk[1] = &divk2 ; 
	p_source[0] = &source1 ; 
	p_source[1] = &source2 ; 


	// Computation of the l.h.s. and r.h.s. of the momentum
	// constraint in each star. 
	// -------------------------------------------------------
	
	double somx = 0 ; 
	double somy = 0 ; 
	double somz = 0 ; 
	
	for (int i=0; i<2; i++) {
	
	    // (flat space) divergence of K^{ij}
	    // ---------------------------------
	    
	    const Tenseur& a_car = et[i]->get_a_car() ; 
	    Tenseur kij_auto = et[i]->get_tkij_auto() / a_car ; 
	    
	    kij_auto.dec2_dzpuis() ; // dzpuis : 2 --> 0 
				     // so that in the external domain, kij_auto
				     // contains now exactly K^{ij}
    
	    // The gradient of K^{ij} is computed on the local triad:
	    kij_auto.change_triad( (et[i]->get_mp()).get_bvect_cart() ) ; 
    
	    *(p_divk[i]) = contract( kij_auto.gradient(), 0, 1) ; 
	    
	    // Back to the Reference triad : 
	    p_divk[i]->change_triad( ref_triad ) ; 
	    kij_auto.change_triad( ref_triad ) ; 
	
	    // Right-hand-side of the momentum constraint
	    // ------------------------------------------
	    
	    const Tenseur& u_euler = et[i]->get_u_euler() ; 
	    const Tenseur& ener_euler = et[i]->get_ener_euler() ; 
	    const Tenseur& press = et[i]->get_press() ; 
	
	    
	    Tenseur d_alpha =   et[i]->get_d_beta_auto()	
			      - et[i]->get_d_logn_auto() 
			      + et[i]->get_d_beta_comp()	
			      - et[i]->get_d_logn_comp() ; 
			      
	    *(p_source[i]) = 2 * qpig * (ener_euler + press) * u_euler
			     - 5 * contract(kij_auto, 1, d_alpha, 0) ; 

	    // Relative differences
	    // --------------------
	    Tbl diffx = diffrel( (*(p_divk[i]))(0), (*(p_source[i]))(0)) ;
	    Tbl diffy = diffrel( (*(p_divk[i]))(1), (*(p_source[i]))(1)) ;
	    Tbl diffz = diffrel( (*(p_divk[i]))(2), (*(p_source[i]))(2)) ;

	    cout << "Binaire::mom_constr : norme div(K) : " << endl ;
	    cout << "X component : " << norme( (*(p_divk[i]))(0) ) << endl ;  
	    cout << "Y component : " << norme( (*(p_divk[i]))(1) ) << endl ;  
	    cout << "Z component : " << norme( (*(p_divk[i]))(2) ) << endl ;  
	
	    cout << "Binaire::mom_constr : norme source : " << endl ;
	    cout << "X component : " << norme( (*(p_source[i]))(0) ) << endl ;  
	    cout << "Y component : " << norme( (*(p_source[i]))(1) ) << endl ;  
	    cout << "Z component : " << norme( (*(p_source[i]))(2) ) << endl ;  
	
	
	    cout << 
	    "Binaire::mom_constr : rel. diff. div(K) <-> source : "
	    << endl ;
	    cout << "X component : " << diffx  << endl ;  
	    cout << "Y component : " << diffy  << endl ;  
	    cout << "Z component : " << diffz  << endl ;  
	
	
	    somx += max( abs(diffx) ) ;    
	    somy += max( abs(diffy) ) ;    
	    somz += max( abs(diffz) ) ;    
	}   
	
	// Total error
	// -----------

	p_mom_constr = new Tbl(3) ; 
	p_mom_constr->set_etat_qcq() ;
	
	p_mom_constr->set(0) = 0.5 * somx ;  
	p_mom_constr->set(1) = 0.5 * somy ;  
	p_mom_constr->set(2) = 0.5 * somz ; 
	
	
    }    
    
    return *p_mom_constr ;     

}
}
