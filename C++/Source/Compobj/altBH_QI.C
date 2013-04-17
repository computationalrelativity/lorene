/*
 *  Methods of the class AltBH_QI
 *
 *    (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2013 Odele Straub, Eric Gourgoulhon
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

char altBH_QI_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2013/04/17 13:02:47  e_gourgoulhon
 * Added member krphi + method extrinsic_curvature
 *
 * Revision 1.1  2013/04/16 15:27:27  e_gourgoulhon
 * New class AltBH_QI
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>

// Lorene headers
#include "compobj.h"
#include "unites.h"
#include "nbr_spx.h"

                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
AltBH_QI::AltBH_QI(Map& mpi, const char* file_name, double a_spin_i) :
			 Compobj_QI(mpi) , 
			 a_spin(a_spin_i),
			 krphi(mpi)
{

    ifstream file(file_name) ; 
    if ( !file.good() ) {
        cerr << "Problem in opening the file " << file_name << endl ;
        abort() ;
    }
    
    file.getline(description1, 256) ;
    file.getline(description2, 256) ;
    cout << "description1 : " << description1 << endl ; 
    cout << "description2 : " << description2 << endl ; 
//     description1[0] = " " ; 
//     description2[0] = " " ; 
//     cout << "description1 : " << description1 << endl ; 
//     cout << "description2 : " << description2 << endl ; 

    const Mg3d* mg = mp.get_mg() ; 
    double tmp ; 
    file.ignore(1000,'\n') ;
    int nz = mg->get_nzone() ; 
    cout << "nz : " << nz << endl ; 
    a_car.set_etat_qcq() ; // Memory allocation for A^2
    nn.set_etat_qcq() ; // Memory allocation for N
    nphi.allocate_all() ; // Memory allocation for N^phi
    krphi.allocate_all() ; // Memory allocation for K_rphi
    double r_inner = 0 ;
    for (int l=1; l<nz; l++) {
        cout << "l = " << l << endl ; 
        int nr = mg->get_nr(l) ; 
        int nt = mg->get_nt(l) ; 
        int np = mg->get_np(l) ;
        double* r_iso = new double[nr] ;
        double* r_areal = new double[nr] ;
        double* psi4 = new double[nr] ;
        double* alpha = new double[nr] ; 
        double* Krphi = new double[nr] ; 
        double* beta_phi = new double[nr] ; 
        int nr_max = nr ; 
        if (l==nz-1) {
            nr_max = nr - 1 ;
            r_areal[nr-1] = __infinity ; 
            r_iso[nr-1] = __infinity ;
            psi4[nr-1] = 1 ; 
            alpha[nr-1] = 1 ; 
            Krphi[nr-1] = 0 ; 
            beta_phi[nr-1] = 0 ; 
        }
        for (int i=0; i<nr_max; i++) {
            file >> r_areal[i] ; 
            file >> r_iso[i] ; 
            file >> tmp ; 
            file >> psi4[i] ; 
            file >> alpha[i] ; 
            file >> tmp ; 
            file >> Krphi[i] ; 
            file >> beta_phi[i] ; 
            cout << "r_iso, psi4, beta_phi : " << r_iso[i] << "  " << psi4[i] << "  " << beta_phi[i] << endl ; 
        }
        
        if (l==1) r_inner = r_iso[0] ; 
        
        for (int k=0; k<np; k++) {
            for (int j=0; j<nt; j++) {
                for (int i=0; i<nr; i++) {
                    a_car.set_grid_point(l,k,j,i) = psi4[i] ;
                    nn.set_grid_point(l,k,j,i) = alpha[i] ;
                    nphi.set_grid_point(l,k,j,i) = - a_spin * beta_phi[i] / psi4[i] / (r_iso[i]*r_iso[i]) ;
                    krphi.set_grid_point(l,k,j,i) = a_spin * Krphi[i] / r_iso[i] ;
                }
            }
        }
    }
    file.close() ; 
    
    mp.homothetie(r_inner / mp.val_r(1,-1.,0.,0.)) ; 
    
    // Set the shift to zero in domain l=0:
    nphi.annule(0,0) ;
    
    a_car.std_spectral_base() ;
    nn.std_spectral_base() ;
    nphi.std_spectral_base() ;
    krphi.std_spectral_base() ;
    
    b_car = a_car ; // slow rotation limit: B^2 = A^2 
    
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}

// Copy constructor
// --------------------
AltBH_QI::AltBH_QI(const AltBH_QI& other) :
			 Compobj_QI(other),
			 krphi(other.krphi)
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
}


// Constructor from a file
// -----------------------
AltBH_QI::AltBH_QI(Map& mpi, FILE* ) :
			 Compobj_QI(mpi),
			 krphi(mpi)
{
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    // Read of the saved fields:
    // ------------------------

}

			    //------------//
			    // Destructor //
			    //------------//

AltBH_QI::~AltBH_QI(){

    del_deriv() ; 

}


			//----------------------------------//
			// Management of derived quantities //
			//----------------------------------//

void AltBH_QI::del_deriv() const {

    Compobj_QI::del_deriv() ; 


    AltBH_QI::set_der_0x0() ; 
}			    


void AltBH_QI::set_der_0x0() const {
 	 
}			    

			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another AltBH_QI
// --------------------------------
void AltBH_QI::operator=(const AltBH_QI& other) {

    // Assignment of quantities common to all the derived classes of Compobj_QI
    Compobj_QI::operator=(other) ;	    
    
    del_deriv() ;  // Deletes all derived quantities
}	

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void AltBH_QI::sauve(FILE* ) const {

    
}

// Printing
// --------

ostream& AltBH_QI::operator>>(ostream& ost) const {

 	using namespace Unites ;
	
	Compobj_QI::operator>>(ost) ; 
    
    ost << endl << "Alternative black hole  spacetime in quasi-isotropic coordinates (class AltBH_QI) " << endl ; 
    ost << description1 << endl ; 
    ost << description2 << endl ; 
   
    return ost ; 
      
}

			    //-------------------------//
			    //	Computational methods  //
			    //-------------------------//
			    
// Updates the extrinsic curvature
// -------------------------------

void AltBH_QI::extrinsic_curvature() {

    // Special treatment for axisymmetric case:
    
    if ( (mp.get_mg())->get_np(0) == 1) {
    
        // What follows is valid only for a mapping of class Map_radial :   
        assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ;
        
        Scalar tmp = krphi ; 
        tmp.mult_sint() ;  // multiplication by sin(theta)
        kk.set(1,3) = tmp ; 
        
        kk.set(2,3) = 0 ; 

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
