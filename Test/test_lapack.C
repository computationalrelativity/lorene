/*
 * Test code for LORENE Matrice class and LAPACK
 */
 
/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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
 * Revision 1.4  2014/10/13 08:54:08  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:07:36  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2002/10/16 14:37:19  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:31  e_gourgoulhon
 * LORENE
 *
 *
 * $Header$
 *
 */

// C headers:
#include <cmath>

// Lorene headers
#include "matrice.h"

using namespace Lorene ;

int main() {
    
    Matrice a(3, 3) ; 
    a.set_etat_qcq() ; 
    
    for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	    a.set(i, j) = 0 ; 
	}
    }
    a.set(0, 0) = 1 ;     
    a.set(1, 1) = 1 ;     
    a.set(2, 2) = 1 ;     
    a.set(1, 2) = 1 ; 
    a.set(2, 1) = -1 ; 
    a.set(0, 2) = 2 ; 
    a.set(1, 0) = -3 ; 

    Tbl x0(3) ; 
    x0.set_etat_qcq() ; 
    
    for (int i=0; i<3; i++) {
	x0.set(i) = 1+i ; 
    }
	
    Tbl y(3) ; 
    y.set_etat_qcq() ; 

    for (int i=0; i<3; i++) {
	y.set(i) = 0 ; 
	for (int j=0; j<3; j++) {
	    y.set(i) += a(i, j) * x0(j) ; 
	}
    }
    
    a.set_band(2, 2) ; 
    a.set_lu() ; 
    
    Tbl x = a.inverse(y) ; 
    
    cout << "Matrix a : " << endl ; 
    cout << a << endl ; 
    
    cout << "Det(a) = " << a.determinant() << endl ; 
    
    cout << "R.H.S. y : " << endl ; 
    cout << y << endl ; 
    
    
    cout << "Test of resolution of the system A.x = y : " << endl ; 
    for (int i=0; i<3; i++) {
	cout << x(i) - x0(i) << endl ; 
    }
    
  return EXIT_SUCCESS ; 
 
}
