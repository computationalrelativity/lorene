/*
 *  Function Scalar::asymptot
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2002 Eric Gourgoulhon (for previous class Cmp)
 *   Copyright (c) 1999-2001 Philippe Grandclement (for previous class Cmp)
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


char scalar_asymptot_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:16:15  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2004/02/21 17:05:14  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 * Method Scalar::set_point renamed Scalar::set_grid_point.
 *
 * Revision 1.2  2003/10/08 14:24:09  j_novak
 * replaced mult_r_zec with mult_r_ced
 *
 * Revision 1.1  2003/09/25 07:18:00  j_novak
 * Method asymptot implemented.
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "tensor.h"

Valeur** Scalar::asymptot(int n0, const int flag) const {
    
    assert(n0 >= 0) ; 
    const Mg3d& mg = *(mp->get_mg()) ; 
    int nz = mg.get_nzone() ; 
    int nzm1 = nz-1 ; 
    assert(mg.get_type_r(nzm1) == UNSURR) ; 
    int np = mg.get_np(nzm1) ; 
    int nt = mg.get_nt(nzm1) ; 
    int nr = mg.get_nr(nzm1) ; 
    int nrm1 = nr-1 ; 
        
    Valeur** vu = new Valeur*[n0+1] ;
    for (int h=0; h<=n0; h++) {
	vu[h] = new Valeur(mg.get_angu()) ; 
    }

    Scalar uu = *this ; 

    int precis = 2 ; 

    // The terms 1/r^h with h < dzpuis are null :
    // -----------------------------------------
    for (int h=0; h<dzpuis; h++) {

	vu[h]->set_etat_zero() ; 	
				
    }

    // Terms 1/r^h with h >= dzpuis :
    // -----------------------------
    for (int h=dzpuis; h<=n0; h++) {
	
	// Value on the sphere S^2 at infinity
	// -----------------------------------	
	vu[h]->set_etat_c_qcq() ; 
	vu[h]->c->set_etat_qcq() ; 
	for (int l=0; l<nzm1; l++) {
	    vu[h]->c->t[l]->set_etat_zero() ; 
	}
	vu[h]->c->t[nzm1]->set_etat_qcq() ; 

	for (int k=0; k<np; k++) {
	    for (int j=0; j<nt; j++) {
		vu[h]->set(nzm1, k, j, 0) = uu.val_grid_point(nzm1, k, j, nrm1) ; 
	    }
	}
	
	vu[h]->set_base( uu.va.base ) ; 
	
	// Printing
	// --------
	if (flag != 0) {
	    cout << "Term in 1/r^" << h << endl ; 
	    cout << "-------------" << endl ; 

	    double ndec = 0 ; 
	    double vmin = (*vu[h])(nzm1, 0, 0, 0) ; 
	    double vmax = vmin ; 

	    cout << "            Values at the point (phi_k, theta_j) : " << endl ; 
	    cout.precision(precis) ;
	    cout.setf(ios::showpoint) ;
	    for (int k=0; k<np; k++) {
		cout <<  " k=" << k << " : "  ;
	    for (int j=0; j<nt; j++) {
		double xx = (*vu[h])(nzm1, k, j, 0) ;
		cout <<  " " << setw(precis) << xx ;	
		ndec += fabs(xx) ; 	
		vmin = ( xx < vmin ) ? xx : vmin ; 
		vmax = ( xx > vmax ) ? xx : vmax ; 
	    }
	    cout << endl;
	}
	ndec /= np*nt ; 
	cout << "Minimum value on S^2 : " << vmin << endl ; 
	cout << "Maximum value on S^2 : " << vmax << endl ; 
	cout << "L^1 norm on S^2     : " << ndec << endl ;  
	}		
	// The value at infinity is substracted
	// ------------------------------------
	for (int k=0; k<np; k++) {
	    for (int j=0; j<nt; j++) {
		double v_inf = (*vu[h])(nzm1, k, j, 0) ; 
		for (int i=0; i<nr; i++) {
		   uu.set_grid_point(nzm1, k, j, i) -= v_inf ;  
		}
	    }
	}
	
	// Mutliplication by r
	// -------------------
	
	uu.mult_r_ced() ; 
	
    } // End of loop on h  (development order)	   
	
    return vu ; 
    
}
