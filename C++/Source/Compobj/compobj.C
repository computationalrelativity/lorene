/*
 *  Methods of the class Compobj
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

char compobj_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2012/11/22 16:04:51  c_some
 * Minor modifications
 *
 * Revision 1.2  2012/11/20 16:24:09  c_some
 * Added computation of ADM mass (method mass_q())
 *
 * Revision 1.1  2012/11/15 16:20:51  c_some
 * New class Compobj
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>
#include <cmath>

// Lorene headers
#include "compobj.h"
#include "nbr_spx.h"

                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
Compobj::Compobj(Map& map_i) :
		mp(map_i) , 
		nn(map_i) , 
		beta(map_i, CON, map_i.get_bvect_spher()) ,
 		gamma(map_i.flat_met_spher()) ,
		ener_euler(map_i) ,
		mom_euler(map_i, CON, map_i.get_bvect_spher()) ,
		stress_euler(map_i, COV, map_i.get_bvect_spher()) 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

	// Some initialisations:
	nn = 1 ; 
    	nn.std_spectral_base() ; 

	beta.set_etat_zero() ;
	ener_euler = 0 ; 
	mom_euler.set_etat_zero() ;
	stress_euler.set_etat_zero() ;
	
}

// Copy constructor
// --------------------
Compobj::Compobj(const Compobj& co) :
		mp(co.mp) , 
		nn(co.nn) , 
		beta(co.beta) ,
 		gamma(co.gamma) ,
		ener_euler(co.ener_euler) ,
		mom_euler(co.mom_euler) ,
		stress_euler(co.stress_euler) 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}


// Constructor from a file
// -----------------------
Compobj::Compobj(Map& map_i, FILE* fich) :
		mp(map_i) , 
		nn(map_i, *(map_i.get_mg()), fich) , 
		beta(map_i, map_i.get_bvect_spher(), fich) ,
 		gamma(map_i, fich) ,
		ener_euler(map_i, *(map_i.get_mg()), fich) ,
		mom_euler(map_i,  map_i.get_bvect_spher(), fich) ,
		stress_euler(map_i, map_i.get_bvect_spher(), fich)  
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}

			    //------------//
			    // Destructor //
			    //------------//

Compobj::~Compobj(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Compobj::del_deriv() const {

    if (p_mass_g != 0x0) delete p_mass_g ; 

    Compobj::set_der_0x0() ; 
}			    


void Compobj::set_der_0x0() const {

    p_mass_g = 0x0 ; 

}			    

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Compobj
// ----------------------------
void Compobj::operator=(const Compobj& co) {

	assert( &(co.mp) == &mp ) ;		    // Same mapping
    
	nn = co.nn ;
	beta = co.beta ;
 	gamma = co.gamma ; 
	ener_euler = co.ener_euler ;
	mom_euler = co.mom_euler ;
	stress_euler = co.stress_euler ; 

    	del_deriv() ;  // Deletes all derived quantities
}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Compobj::sauve(FILE* fich) const {

		nn.sauve(fich) ; 
		beta.sauve(fich) ;
 		gamma.sauve(fich) ;
		ener_euler.sauve(fich) ;
		mom_euler.sauve(fich) ;
		stress_euler.sauve(fich) ;

}

// Printing
// --------

ostream& operator<<(ostream& ost, const Compobj& co)  {
    co >> ost ;
    return ost ;
}


ostream& Compobj::operator>>(ostream& ost) const {
    
    ost << endl << "Compact object (class Compobj) " << endl ; 
    ost << "Mapping : " << mp << endl ; 
    ost << "Central values of various fields : " << endl ; 
    ost << "-------------------------------- " << endl ; 
    ost << "   lapse function : N_c = " << nn.val_grid_point(0,0,0,0) << endl ; 
    ost << "   metric components gamma_{ij} : " << endl
    << "    ( " << gamma.cov()(1,1).val_grid_point(0,0,0,0) << "  " 
    		<< gamma.cov()(1,2).val_grid_point(0,0,0,0) << "  " 
     		<< gamma.cov()(1,3).val_grid_point(0,0,0,0) << " )" << endl  
    << "    ( " << gamma.cov()(2,1).val_grid_point(0,0,0,0) << "  " 
    		<< gamma.cov()(2,2).val_grid_point(0,0,0,0) << "  " 
     		<< gamma.cov()(2,3).val_grid_point(0,0,0,0) << " )" << endl  
    << "    ( " << gamma.cov()(3,1).val_grid_point(0,0,0,0) << "  " 
    		<< gamma.cov()(3,2).val_grid_point(0,0,0,0) << "  " 
     		<< gamma.cov()(3,3).val_grid_point(0,0,0,0) << " )" << endl ; 
    ost << "   energy density / Eulerian observer : E_c = " << ener_euler.val_grid_point(0,0,0,0) << endl ; 
    ost << "   components of the stress tensor S_{ij} / Eulerian observer : " << endl
    << "    ( " << stress_euler(1,1).val_grid_point(0,0,0,0) << "  " 
    		<< stress_euler(1,2).val_grid_point(0,0,0,0) << "  " 
     		<< stress_euler(1,3).val_grid_point(0,0,0,0) << " )" << endl  
    << "    ( " << stress_euler(2,1).val_grid_point(0,0,0,0) << "  " 
    		<< stress_euler(2,2).val_grid_point(0,0,0,0) << "  " 
     		<< stress_euler(2,3).val_grid_point(0,0,0,0) << " )" << endl  
    << "    ( " << stress_euler(3,1).val_grid_point(0,0,0,0) << "  " 
    		<< stress_euler(3,2).val_grid_point(0,0,0,0) << "  " 
     		<< stress_euler(3,3).val_grid_point(0,0,0,0) << " )" << endl ; 

	ost << endl << "ADM mass : " << mass_g() << endl ; 
	 	
    return ost ; 
      
}


			    //-------------------------//
			    //	Computational methods  //
			    //-------------------------//
			    
/// Gravitational mass
double Compobj::mass_g() const {

    if (p_mass_g == 0x0) {    // a new computation is required
		
        const Sym_tensor& gam_dd = gamma.cov() ;  // components \gamma_{ij} of the 3-metric
        Metric_flat ff(mp, *(gam_dd.get_triad())) ;
    
        Vector ww = gam_dd.derive_con(ff).trace(1,2).up(0,ff) 
                    - gam_dd.trace(ff).derive_con(ff) ; 

        p_mass_g = new double( ww.flux(__infinity, ff) / (16.* M_PI) ) ; 
	}
            
    return *p_mass_g ; 

} 
	

