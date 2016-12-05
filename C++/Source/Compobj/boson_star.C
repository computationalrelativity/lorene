/*
 *  Methods of the class Boson_star
 *
 *    (see file boson_star.h for documentation).
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
 * Revision 1.5  2016/12/05 16:17:49  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2012/12/03 15:27:30  c_some
 * Small changes
 *
 * Revision 1.2  2012/11/23 15:44:10  c_some
 * Small changes + method update_ener_mom
 *
 * Revision 1.1  2012/11/22 16:04:12  c_some
 * New class Boson_star
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>

// Lorene headers
#include "boson_star.h"
#include "utilitaires.h"

                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
namespace Lorene {
Boson_star::Boson_star(Map& mpi, double m, int k) :
			 Star_QI(mpi) ,
			 rphi(mpi), 
			 iphi(mpi),
			 omega(0),
			 kkk(k),
			 mmm(m),
			 m2(m*m) 	 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;

    // Initialization of the scalar field to zero
    rphi = 0 ;   
    iphi = 0 ;   
    
}

// Copy constructor
// --------------------
Boson_star::Boson_star(const Boson_star& st) :
			 Star_QI(st), 
			 rphi(st.rphi), 
			 iphi(st.iphi),
			 omega(st.omega),
			 kkk(st.kkk),		 
			 mmm(st.mmm),
			 m2(st.m2)		 
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}


// Constructor from a file
// -----------------------
Boson_star::Boson_star(Map& mpi, FILE* fich) :
			 Star_QI(mpi, fich) , 
			 rphi(mpi, *(mpi.get_mg()), fich) ,
			 iphi(mpi, *(mpi.get_mg()), fich)
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    fread_be(&omega, sizeof(double), 1, fich) ;		
    fread_be(&kkk, sizeof(int), 1, fich) ;		
    fread_be(&mmm, sizeof(double), 1, fich) ;	
    m2 = mmm*mmm ;	
        
}

			    //------------//
			    // Destructor //
			    //------------//

Boson_star::~Boson_star(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void Boson_star::del_deriv() const {

    Star_QI::del_deriv() ; 

    Boson_star::set_der_0x0() ; 
}			    


void Boson_star::set_der_0x0() const {
 	 
}			    

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Boson_star
// --------------------------------
void Boson_star::operator=(const Boson_star& st) {

    // Assignment of quantities common to all the derived classes of Star_QI
    Star_QI::operator=(st) ;	    
    
    rphi = st.rphi ; 
    iphi = st.iphi ; 
    omega = st.omega ; 
    kkk = st.kkk ; 
   	mmm = st.mmm ; 
   	m2 = st.m2 ; 

    del_deriv() ;  // Deletes all derived quantities
}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Boson_star::sauve(FILE* fich) const {

	Star_QI::sauve(fich) ; 
    rphi.sauve(fich) ; 
    iphi.sauve(fich) ; 
    fwrite_be(&omega, sizeof(double), 1, fich) ;		
    fwrite_be(&kkk, sizeof(int), 1, fich) ;		
    fwrite_be(&mmm, sizeof(double), 1, fich) ;		

}

// Printing
// --------

ostream& Boson_star::operator>>(ostream& ost) const {
   
    Star_QI::operator>>(ost) ; 
    
    ost << endl << "Axisymmetric stationary boson star in quasi-isotropic coordinates (class Boson_star) " << endl ;
    
    ost << "Boson mass : " << mmm << endl ; //## which units ? 
    
    ost << "omega = " << omega << " ,   k = " << kkk << endl ; 
    
    ost << "Central value of the scalar field : " << rphi.val_grid_point(0,0,0,0) << " + i " << iphi.val_grid_point(0,0,0,0)<< endl ;  

    // ost << "Real part of the scalar field : " << rphi << endl ; 
    // ost << "Imaginary part of the scalar field : " << iphi << endl ; 
    
    //## ost << "ADM mass : " << adm_mass() << endl ; 	
    ost << "Gravitational mass : " << mass_g() << endl ; 	
    return ost ; 
      
}

//------------------------
// Computational routines
//------------------------

	/* Computes the 3+1 components of the energy-momentum tensor (E, P_i and S_{ij})
	 *  from the values of the scalar field and the metric 
	 */
void Boson_star::update_ener_mom() {

	const Metric_flat& ff = mp.flat_met_spher() ;
	
	Scalar mod_phi2 = rphi*rphi + iphi*iphi ; 
	mod_phi2.inc_dzpuis(4) ; // mod_phi2 multiplied by r^4 in the last domain
	
	Sym_tensor dphidphi = rphi.derive_cov(ff) * rphi.derive_cov(ff) 
							+ iphi.derive_cov(ff) * iphi.derive_cov(ff) ;
	
	Scalar tmp = (omega - kkk*nphi)/nn  ;
	Scalar tmp2 = tmp*tmp ; 
	Scalar dphi2 = dphidphi.trace(gamma) ;

    ener_euler = 0.5*( (tmp2 + m2)*mod_phi2 + dphi2 ) ;
    
	mom_euler.set(1) = 0 ;   // p^(r) = 0
	mom_euler.set(2) = 0 ;   // p^(theta) = 0
	mom_euler.set(3) = -kkk*tmp * mod_phi2 ;

	stress_euler = dphidphi - 0.5*( (m2-tmp2)*mod_phi2 + dphi2 )*gamma.cov() ;
}


}
