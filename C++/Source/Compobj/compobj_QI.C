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

 

/*
 * $Id$
 * $Log$
 * Revision 1.10  2016/12/05 16:17:49  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2014/10/13 08:52:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/05/16 11:55:18  o_straub
 * fixed: GYOTO output from compobj & compobj_QI
 *
 * Revision 1.7  2013/07/25 19:44:11  o_straub
 * calculation of the marginally bound radius
 *
 * Revision 1.6  2013/04/04 15:32:32  e_gourgoulhon
 * Better computation of the ISCO
 *
 * Revision 1.5  2013/04/04 08:53:47  e_gourgoulhon
 * Minor improvements
 *
 * Revision 1.4  2013/04/03 12:10:13  e_gourgoulhon
 * Added member kk to Compobj; suppressed tkij
 *
 * Revision 1.3  2012/11/22 16:04:51  c_some
 * Minor modifications
 *
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
#include <cmath>
#include <cstdio>

// Lorene headers
#include "compobj.h"
#include "nbr_spx.h"
#include "utilitaires.h"



                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
namespace Lorene {
Compobj_QI::Compobj_QI(Map& map_i) :
		Compobj(map_i) ,
		a_car(map_i) ,
		bbb(map_i) ,
		b_car(map_i) ,
		nphi(map_i) ,
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
		ak_car(map_i) 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    Scalar nn_file(mp, *(mp.get_mg()), fich) ; 
    nn = nn_file ;
    
    b_car = bbb*bbb ;
     
    // Initialization of gamma_ij:
    update_metric() ; 
    
    // Computation of K_ij and ak_car:
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
    if (p_r_mb != 0x0) delete p_r_mb ;
    if (p_r_isco != 0x0) delete p_r_isco ;
    if (p_f_isco != 0x0) delete p_f_isco ;
    if (p_lspec_isco != 0x0) delete p_lspec_isco ;
    if (p_espec_isco != 0x0) delete p_espec_isco ;

    Compobj_QI::set_der_0x0() ; 
}			    


void Compobj_QI::set_der_0x0() const {

    p_angu_mom = 0x0 ; 
    p_r_mb = 0x0 ;
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

            
// Save in a file for GYOTO input  
// -------------------------------

// Redefinition of Compobj::gyoto_data

void Compobj_QI::gyoto_data(const char* file_name) const {
    
    FILE* file_out = fopen(file_name, "w") ;
    double total_time = 0. ; // for compatibility
    double RISCO=r_isco(0) ;
    double RMB=r_mb(0);

    fwrite_be(&total_time, sizeof(double), 1, file_out) ;
    mp.get_mg()->sauve(file_out) ;
    mp.sauve(file_out) ;
    nn.sauve(file_out) ;
    beta.sauve(file_out) ;
    gamma.cov().sauve(file_out) ;
    gamma.con().sauve(file_out) ;
    kk.sauve(file_out) ;
    fwrite_be(&RISCO, sizeof(double), 1, file_out) ;
    fwrite_be(&RMB, sizeof(double), 1, file_out) ;
    fclose(file_out) ;    

    cout << "WRITING RISCO TO GYOTO FILE : " << RISCO << endl ; 
    cout << "WRITING RMB TO GYOTO FILE   : " << RMB << endl ; 
    cout << "WRITING TO GYOTO FILE - end of part " << endl ; 
}






// Printing
// --------

ostream& Compobj_QI::operator>>(ostream& ost) const {
   
    Compobj::operator>>(ost) ; 
    
    ost << endl << "Axisymmetric stationary compact object in quasi-isotropic coordinates (class Compobj_QI) " << endl ; 

    ost << "Central values of various fields : " << endl ; 
    ost << "-------------------------------- " << endl ; 
    ost << "   metric coefficient A^2   : " << a_car.val_grid_point(0,0,0,0) << endl ; 
    ost << "   metric coefficient B^2   : " << b_car.val_grid_point(0,0,0,0) << endl ; 
    ost << "   metric coefficient N^phi : " << nphi.val_grid_point(0,0,0,0) << endl ; 
    ost << "   A^2 K_{ij} K^{ij}        = " << ak_car.val_grid_point(0,0,0,0) << endl << endl ; 

    
    double RISCO = r_isco(0, &ost) ; 
    ost << "Coordinate r at the innermost stable circular orbit (ISCO) : " << 
        RISCO  << endl ;  
     // ost << "Circumferential radius of the innermost stable circular orbit (ISCO) : " << 
     // bbb.val_point(risco, M_PI/2, 0)*risco << endl ;  
     //	ost << "Orbital frequency at the ISCO : " << f_isco(0) << endl ; 
     // ost << "Specific energy of a particle on the ISCO : " << espec_isco(0) << endl ;	
     // ost << "Specific angular momentum of a particle on the ISCO : " << lspec_isco(0) << endl ;	
	

    double RMB = r_mb(0, &ost) ; 
    ost << "Coordinate r at the marginally bound circular orbit (R_mb) : " << RMB  << endl ; 


  //  ost << "A^2 : " << a_car << endl ; 
  //  ost << "B^2 : " << b_car << endl ; 
  //  ost << "nphi : " << nphi << endl ; 
	
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
// -------------------------------

void Compobj_QI::extrinsic_curvature() {

	// Special treatment for axisymmetric case:
	
 	if ( (mp.get_mg())->get_np(0) == 1) {
 	
		Scalar dnpdr = nphi.dsdr() ; 		// d/dr (N^phi)
 		Scalar dnpdt = nphi.dsdt() ; 		// d/dtheta (N^phi)
 		
 		// What follows is valid only for a mapping of class Map_radial :	
		assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ;
		
        dnpdr.mult_rsint() ;    // multiplication by r sin(theta)
        kk.set(1,3) = - b_car * dnpdr / (2*nn) ; 
        
        dnpdt.mult_sint() ; // multiplication by sin(theta)
        kk.set(2,3) = - b_car * dnpdt / (2*nn) ; 
        kk.set(2,3).inc_dzpuis(2) ;  // to have the same dzpuis as kk(1,3)
        
        kk.set(1,1) = 0 ; 
        kk.set(1,2) = 0 ; 
        kk.set(2,2) = 0 ; 
        kk.set(3,3) = 0 ; 
	}
    else {

    // General case:

        Compobj::extrinsic_curvature() ; 
   }
    
    // Computation of A^2 K_{ij} K^{ij}
    // --------------------------------
        
    ak_car = 2 * ( kk(1,3)*kk(1,3) +  kk(2,3)*kk(2,3) ) / b_car ;
    
	del_deriv() ; 

}
}
