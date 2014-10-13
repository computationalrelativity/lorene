/*
 * Test code for Lorene class Tbl and PGPLOT
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
 * Revision 1.5  2014/10/13 08:54:08  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:07:36  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2008/09/22 19:08:56  j_novak
 * Minor modifs.
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
#include "tbl.h"
#include "graphique.h"

using namespace Lorene ;

int main() {
    
    int np = 100 ;
     
    Tbl a(np) ; 
    a.set_etat_qcq() ; 
    
    double h = 2*M_PI / double(np-1) ; 
    
    for (int i=0; i<np; i++) {
	a.set(i) = sin(h*i) ; 
    }
    
    float* uutab = new float[np] ; 
    
    for (int i=0; i<np; i++) {
	uutab[i] = float(a(i)) ; 
    }
       
    
    float xmin = 0. ; 
    float xmax = float(2*M_PI) ; 
    
    des_profile(uutab, np, xmin, xmax, "x", "sin(x)", "Test") ; 
    
    delete [] uutab ; 

    return EXIT_SUCCESS ; 
 
}
