/*
 *  Test code for Mg3d class.
 *
 */

/*
 *   Copyright (c) 2001 Eric Gourgoulhon.
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

char test_mg3d_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/12/11 06:50:14  e_gourgoulhon
 * test code for Mg3d class
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <iostream.h>

// C headers
#include <stdlib.h>

// Lorene headers
#include "grilles.h"
#include "type_parite.h"


int main() {

        int nz = 4 ;
        int nr = 17 ;
        int nt = 9 ;
        int np = 4 ;

        Mg3d mg(nz, nr, nt, np, SYM, NONSYM) ;

        cout << "Multi-grid:" << endl ;
        cout << mg << endl ;
	
	return EXIT_SUCCESS ;

}

