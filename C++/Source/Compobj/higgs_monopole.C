/*
 *  Methods of the class HiggsMonopole
 *
 *    (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2014 Marie Leroy,  Eric Gourgoulhon
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

char HiggsMonopole_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/13 08:52:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/01/31 15:35:27  e_gourgoulhon
 * Reading file in constructor
 *
 * Revision 1.1  2014/01/29 16:31:42  e_gourgoulhon
 * New class HiggsMonopole
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

                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
// --------------------
namespace Lorene {
HiggsMonopole::HiggsMonopole(Map& mpi, const char* file_name) :
			 Compobj(mpi) , 
			 hh(mpi) ,
             grr(mpi) , 
             press(mpi)
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

    file.ignore(1000,'\n') ;

    const Mg3d* mg = mp.get_mg() ; 
    double tmp ; 
    int nz = mg->get_nzone() ; 
    cout << "nz : " << nz << endl ; 
    nn.set_etat_qcq() ; // Memory allocation for N
    grr.allocate_all() ; // Memory allocation for g_rr
    hh.allocate_all() ; // Memory allocation for h
    press.allocate_all() ; // Memory allocation for press
    double explamb_last, nn_last, hh_last, press_last ; 
    for (int l=0; l<nz; l++) {
        cout << "l = " << l << endl ; 
        int nr = mg->get_nr(l) ; 
        int nt = mg->get_nt(l) ; 
        int np = mg->get_np(l) ;
        double u_r ; // u = r m_p
        double* explamb_tab = new double[nr] ;
        double* nn_tab = new double[nr] ;
        double* hh_tab = new double[nr] ;
        double* press_tab = new double[nr] ; 
        int i_min = 0 ; 
        if (l>0) {
            i_min = 1 ;
            explamb_tab[0] = explamb_last ; 
            nn_tab[0] = nn_last ; 
            hh_tab[0] = hh_last ;
            press_tab[0] = press_last ;
        }
        for (int i=i_min; i<nr; i++) {
            file >> u_r ; 
            file >> explamb_tab[i] ; 
            file >> nn_tab[i] ; 
            file >> hh_tab[i] ; 
            file >> press_tab[i] ; 
        }
        explamb_last = explamb_tab[nr-1] ; 
        nn_last = nn_tab[nr-1] ; 
        hh_last = hh_tab[nr-1] ; 
        press_last = press_tab[nr-1] ; 
 
        for (int i=0; i<nr; i++) {
            cout << " explamb, nn, hh  : " << explamb_tab[i] << "  " 
                    << nn_tab[i] << "  " << hh_tab[i] << endl ; 
        }
        
        for (int k=0; k<np; k++) {
            for (int j=0; j<nt; j++) {
                for (int i=0; i<nr; i++) {
                    grr.set_grid_point(l,k,j,i) = explamb_tab[i]*explamb_tab[i] ;
                    nn.set_grid_point(l,k,j,i) = nn_tab[i] ;
                    hh.set_grid_point(l,k,j,i) = hh_tab[i] ;
                    press.set_grid_point(l,k,j,i) = press_tab[i] ;
                }
            }
        }
        
        delete[] explamb_tab ;
        delete[] nn_tab ;
        delete[] hh_tab ;
        delete[] press_tab ;

    }
    file.close() ; 
    
    nn.std_spectral_base() ;  // Sets bases for spectral expansions
    grr.std_spectral_base() ;  
    hh.std_spectral_base() ;  
    press.std_spectral_base() ;  

   // Pointers of derived quantities initialized to zero : 
    // set_der_0x0() ;
}

// Copy constructor
// --------------------
HiggsMonopole::HiggsMonopole(const HiggsMonopole& other) :
			 Compobj(other),
			 hh(other.hh) , 
             grr(other.grr) , 
             press(other.press)
{
    // Pointers of derived quantities initialized to zero : 
    // set_der_0x0() ;
}

			    //------------//
			    // Destructor //
			    //------------//

HiggsMonopole::~HiggsMonopole(){

    // del_deriv() ; 

}

// Printing
// --------

ostream& HiggsMonopole::operator>>(ostream& ost) const {

 	using namespace Unites ;
	
    ost << endl << "Higgs monopole" << endl ; 

    ost << description1 << endl ; 
    ost << description2 << endl ; 

	Compobj::operator>>(ost) ; 

    return ost ; 
      
}



}
