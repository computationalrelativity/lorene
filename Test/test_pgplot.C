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
 * Revision 1.1  2001/11/20 15:19:31  e_gourgoulhon
 * Initial revision
 *
 *
 * $Header$
 *
 */

// C++ headers:
#include <iostream.h>

// C headers:
#include <math.h>

// Lorene headers
#include "tbl.h"
#include "graphique.h"

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
	uutab[i] = a(i) ; 
    }
       
    
    float xmin = 0. ; 
    float xmax = 2*M_PI ; 
    
    des_profile(uutab, np, xmin, xmax, "x", "sin(x)", "Test") ; 
    
    delete [] uutab ; 

    return EXIT_SUCCESS ; 
 
}
