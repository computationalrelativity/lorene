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


 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:11  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2010/01/26 16:47:40  e_gourgoulhon
 * Added the Scalar version.
 *
 * Revision 1.2  2003/12/19 16:21:49  j_novak
 * Shadow hunt
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/03/22  10:08:55  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "cmp.h"
#include "scalar.h"

namespace Lorene {
Cmp raccord_c1(const Cmp& uu, int l1) {
    
    const Map_radial* mpi = dynamic_cast<const Map_radial*>( uu.get_mp() ) ; 

    if (mpi == 0x0) {
	cout << 
	"raccord_c1 : The mapping does not belong to the class Map_radial !"
	    << endl ; 
	abort() ;
    }

    const Map_radial& mp = *mpi ; 

    const Mg3d& mg = *(mp.get_mg()) ; 
    
#ifndef NDEBUG
    int nz = mg.get_nzone() ; 
#endif
    assert(l1 > 0) ; 
    assert(l1 < nz-1) ; 
    
    int l0 = l1 - 1 ;	// index of left domain
    int l2 = l1 + 1 ;	// index of right domain
    
    
    Mtbl dxdr = mp.dxdr ; 
    Mtbl r2 = mp.r * mp.r ; 
     
    Cmp resu = uu ; 
    
    Valeur& va = resu.va ; 
    
    va.coef_i() ; 
    va.set_etat_c_qcq() ; 
    
    va.c->set_etat_qcq() ; 
    va.c->t[l1]->set_etat_qcq() ; 
    
    int np = mg.get_np(l1) ; 
    int nt = mg.get_nt(l1) ; 
    assert( mg.get_np(l0) == np ) ;    
    assert( mg.get_nt(l0) == nt ) ;    
    assert( mg.get_np(l2) == np ) ;    
    assert( mg.get_nt(l2) == nt ) ;    
    
    int nr0 = mg.get_nr(l0) ; 
    int nr1 = mg.get_nr(l1) ; 
    
    for (int k=0; k<np; k++) {
	for (int j=0; j<nt; j++) {
	    double lam0 = uu(l0, k, j, nr0-1) ;    
	    double lam1 = uu.dsdr()(l0, k, j, nr0-1) / dxdr(l1, k, j, 0) ;    
	    double mu0 = uu(l2, k, j, 0) ;    
	    double mu1 = uu.dsdr()(l2, k, j, 0) / dxdr(l1, k, j, nr1-1) ;  
	   
	    if ( mg.get_type_r(l2) == UNSURR ) {
		mu1 /= r2(l2, k, j, 0) ; 
	    }
	   
	    double s0 = 0.25 * (mu0 + lam0) ;  
	    double s1 = 0.25 * (mu1 + lam1) ;  
	    double d0 = 0.25 * (mu0 - lam0) ;  
	    double d1 = 0.25 * (mu1 - lam1) ;  
	    
	    double a0 = 2. * s0 - d1 ; 
	    double a1 = 3. * d0 - s1 ; 
	    double a2 = d1 ; 
	    double a3 = s1 - d0 ; 
	    
	    for (int i=0; i<nr1; i++) {
		double x = mg.get_grille3d(l1)->x[i] ; 
    		va.set(l1, k, j, i) = a0 + x * ( a1 + x * ( a2 + x * a3 ) ) ; 
	    }
	     
	}
    }
    
    return resu ; 
    
}

/*
 * Scalar version 
 */
Scalar raccord_c1(const Scalar& uu, int l1) {

    Cmp cuu(uu) ; 
    return Scalar( raccord_c1(cuu, l1) ) ; 
   
}
}
