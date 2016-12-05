/*
 *  Function Cmp::asymptot
 *
 */

/*
 *   Copyright (c) 1999-2002 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:17:48  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:52:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:03  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2002/10/16 14:36:34  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/05/22 13:57:25  f_limousin
 * Corrected error in determination of min and max values on the sphere.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2000/11/15  13:24:43  phil
 * modification output
 *
 * Revision 2.1  2000/11/15  13:16:01  phil
 * ajout gestion affichage
 *
 * Revision 2.0  2000/03/25  12:53:03  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cmath>

// Headers Lorene
#include "cmp.h"

namespace Lorene {
Valeur** Cmp::asymptot(int n0, const int flag) const {
    
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

    Cmp uu = *this ; 

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
		vu[h]->set(nzm1, k, j, 0) = uu(nzm1, k, j, nrm1) ; 
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
		   uu.set(nzm1, k, j, i) -= v_inf ;  
		}
	    }
	}
	
	// Mutliplication by r
	// -------------------
	
	uu.mult_r_zec() ; 
	
    } // End of loop on h  (development order)	   
	
    return vu ; 
    
}
}
