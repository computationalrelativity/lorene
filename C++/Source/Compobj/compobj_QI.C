/*
 *  Methods of the class Compobj_QI
 *
 *    (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2012 Claire Some, Eric Gourgoulhon
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

char compobj_QI_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2012/11/20 16:28:48  c_some
 * -- tkij is created on the Cartesian triad.
 * -- implemented method extrinsic_curvature()
 *
 * Revision 1.1  2012/11/16 16:14:11  c_some
 * New class Compobj_QI
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>

// Lorene headers
#include "compobj.h"

                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
Compobj_QI::Compobj_QI(Map& map_i) :
		Compobj(map_i) ,
		a_car(map_i) ,
		bbb(map_i) ,
		b_car(map_i) ,
		nphi(map_i) ,
 		tkij(map_i, COV, map_i.get_bvect_cart()),		
		ak_car(map_i) 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

    // Initialization to a flat metric : 
    a_car = 1 ;
    a_car.std_spectral_base() ; 
    bbb = 1 ;
    bbb.std_spectral_base() ; 
    b_car = bbb*bbb ;
    nphi = 0 ;   
    tkij.set_etat_zero() ; 
    ak_car = 0 ; 

}

// Copy constructor
// --------------------
Compobj_QI::Compobj_QI(const Compobj_QI& co) :
		Compobj(co), 
		a_car(co.a_car) , 
		bbb(co.bbb) , 
		b_car(co.b_car) ,
 		nphi(co.nphi) ,
		tkij(co.tkij) ,
		ak_car(co.ak_car) 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}


// Constructor from a file
// -----------------------
Compobj_QI::Compobj_QI(Map& map_i, FILE* fich) :
		Compobj(map_i) , 
		a_car(map_i, *(map_i.get_mg()), fich) , 
		bbb(map_i, *(map_i.get_mg()), fich) , 
		b_car(map_i) , 
		nphi(map_i, *(map_i.get_mg()), fich) , 
 		tkij(map_i, COV, map_i.get_bvect_cart()),		
		ak_car(map_i) 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    Scalar nn_file(mp, *(mp.get_mg()), fich) ; 
    nn = nn_file ;
    
    b_car = bbb*bbb ;
     
    // Initialization of gamma_ij:
    update_metric() ; 
    
    // Computation of tkij and ak_car:
    extrinsic_curvature() ; 
}

			    //------------//
			    // Destructor //
			    //------------//

Compobj_QI::~Compobj_QI(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Compobj_QI::del_deriv() const {

    Compobj::del_deriv() ; 

   	if (p_angu_mom != 0x0) delete p_angu_mom ; 
    if (p_r_isco != 0x0) delete p_r_isco ;
    if (p_f_isco != 0x0) delete p_f_isco ;
    if (p_lspec_isco != 0x0) delete p_lspec_isco ;
    if (p_espec_isco != 0x0) delete p_espec_isco ;

    Compobj_QI::set_der_0x0() ; 
}			    


void Compobj_QI::set_der_0x0() const {

    p_angu_mom = 0x0 ; 
    p_r_isco = 0x0 ;
    p_f_isco = 0x0 ;
    p_lspec_isco = 0x0 ;
    p_espec_isco = 0x0 ;
 	 
}			    

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Compobj_QI
// --------------------------------
void Compobj_QI::operator=(const Compobj_QI& co) {

    // Assignment of quantities common to all the derived classes of Compobj
    Compobj::operator=(co) ;	    
    
    a_car = co.a_car ; 
    bbb = co.bbb ; 
    b_car = co.b_car ; 
    nphi = co.nphi ; 
    tkij = co.tkij ; 
    ak_car = co.ak_car ;

    del_deriv() ;  // Deletes all derived quantities
}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Compobj_QI::sauve(FILE* fich) const {

	a_car.sauve(fich) ; 
	bbb.sauve(fich) ;
	nphi.sauve(fich) ;

	nn.sauve(fich) ;     

}

// Printing
// --------

ostream& Compobj_QI::operator>>(ostream& ost) const {
   
    Compobj::operator>>(ost) ; 
    
    ost << endl << "Axisymmetric stationary compact object in quasi-isotropic coordinates (class Compobj_QI) " << endl ; 
    ost << "A^2 : " << a_car << endl ; 
    ost << "B^2 : " << b_car << endl ; 
    ost << "nphi : " << nphi << endl ; 
	
    return ost ; 
      
}

// Updates the 3-metric and the shift

void Compobj_QI::update_metric() {

    Sym_tensor gam(mp, COV, mp.get_bvect_spher()) ; 
    gam.set(1,1) = a_car ; 
    gam.set(1,2) = 0 ; 
    gam.set(1,3) = 0 ; 
    gam.set(2,2) = a_car ; 
    gam.set(2,3) = 0 ; 
    gam.set(3,3) = b_car ;
    
    gamma = gam ;

	assert(*(beta.get_triad()) == mp.get_bvect_spher()) ; 
	
	beta.set(1) = 0 ;
	beta.set(2) = 0 ;
	Scalar nphi_ortho(nphi) ; 
	nphi_ortho.mult_rsint() ;
	beta.set(3) = - nphi_ortho ; 
	
    // Tensor B^{-2} K_{ij} and Scalar A^2 K_{ij} K^{ij}
    // -------------------------------------------------
    
    extrinsic_curvature() ; 
    
  
    // The derived quantities are no longer up to date : 
    // -----------------------------------------------

    del_deriv() ;  
	
}


// Updates the extrinsic curvature

void Compobj_QI::extrinsic_curvature() {

	// ---------------------------------------
	// Special treatment for axisymmetric case
	// ---------------------------------------
	
 	if ( (mp.get_mg())->get_np(0) == 1) {
 	
 		tkij.set_etat_zero() ;		// initialisation
				
		// Computation of K_xy
		// -------------------
		
		Scalar dnpdr = nphi.dsdr() ; 		// d/dr (N^phi)
 		Scalar dnpdt = nphi.srdsdt() ; 		// 1/r d/dtheta (N^phi)
 		
 		// What follows is valid only for a mapping of class Map_radial :	
		assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ;
		
		if (dnpdr.get_etat() == ETATQCQ) {
		    // multiplication by sin(theta)
		    dnpdr.set_spectral_va() = (dnpdr.get_spectral_va()).mult_st() ;	
 		}
		
		if (dnpdt.get_etat() == ETATQCQ) {
		    // multiplication by cos(theta)
		    dnpdt.set_spectral_va() = (dnpdt.get_spectral_va()).mult_ct() ;	
 		}
 	
		Scalar tmp = dnpdr + dnpdt ;
 	
		tmp.mult_rsint() ;	// multiplication by r sin(theta)
 	
		tkij.set(1,2) = - 0.5 * tmp / nn ; 	// component (x,y)
 	
 	
		// Computation of K_yz
		// -------------------
 	
		dnpdr = nphi.dsdr() ; 		// d/dr (N^phi)
 		dnpdt = nphi.srdsdt() ; 		// 1/r d/dtheta (N^phi)
 		
		if (dnpdr.get_etat() == ETATQCQ) {
		    // multiplication by cos(theta)
		    dnpdr.set_spectral_va() = (dnpdr.get_spectral_va()).mult_ct() ;	
 		}
		
		if (dnpdt.get_etat() == ETATQCQ) {
		    // multiplication by sin(theta)
		    dnpdt.set_spectral_va() = (dnpdt.get_spectral_va()).mult_st() ;	
 		}
 	
		tmp = dnpdr - dnpdt ;
		
		tmp.mult_rsint() ;	// multiplication by r sin(theta)
 		
		tkij.set(2,3) = - 0.5 * tmp / nn ; 	// component (y,z)
 	
		// The other components are set to zero
		// ------------------------------------
		tkij.set(1,1) = 0 ;	// component (x,x)
		tkij.set(1,3) = 0 ;     // component (x,z)
		tkij.set(2,2) = 0 ;    	// component (y,y)
		tkij.set(3,3) = 0 ;     // component (z,z)
 	
 	}
    else {

    // ------------
    // General case
    // ------------

    	// Gradient (Cartesian components) of the shift
    	// D_j N^i
    
    	Tensor dn = - beta.derive_cov( mp.flat_met_cart() ) ;
    
    	// Trace of D_j N^i = divergence of N^i :
    	Scalar divn = contract(dn, 0, 1) ;
    
    	if (divn.get_etat() == ETATQCQ) {
    
		// Computation of B^{-2} K_{ij}
		// ----------------------------
		tkij.set_etat_qcq() ;
		for (int i=1; i<=3; i++) {
		    for (int j=i; j<=3; j++) {
			  tkij.set(i, j) = dn(i, j) + dn(j, i)  ;
		    }
		    tkij.set(i, i) -= double(2) /double(3) * divn ;
		}
    
		tkij = - 0.5 * tkij / nn ;
	
    	}
    	else{
		assert( divn.get_etat() == ETATZERO ) ;
		tkij.set_etat_zero() ;
    	}
   }
    
    // Computation of A^2 K_{ij} K^{ij}
    // --------------------------------
        
    ak_car = 0 ;
    
    for (int i=1; i<=3; i++) {
		for (int j=1; j<=3; j++) {
	
	    ak_car += tkij(i, j) * tkij(i, j) ;
	
		}
    }
    
    ak_car = b_car * ak_car ;
    
	del_deriv() ; 

}
