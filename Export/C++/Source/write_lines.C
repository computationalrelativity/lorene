/*
 * Writes a C array in an ostream with a fixed number of items per line
 *
 */

/*
 *   Copyright (c) 2002  Eric Gourgoulhon
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


char write_lines_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/01/09 11:08:00  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.1  2002/01/11 17:03:02  e_gourgoulhon
 * Exportation of binary neutron stars configuration to a Cartesian grid
 *
 *
 * $Header$
 *
 */

#ifdef OBSOLETE_HEADERS

#include <iostream.h>

#else

#include <iostream>
using namespace std ;
#endif

void write_lines(ostream& fich, int dpl, const double* pdata, int np) {

        int nlines = np / dpl ;   // number of filled lines
        int reste = np - nlines * dpl ;	// number of remaining data

	for (int line = 0; line < nlines; line++) {
	    for (int i=0; i<dpl; i++) {
		fich << *pdata << "  " ;
		pdata++ ;
	    }
	    fich << endl ;
	}
	for (int i=0; i<reste; i++) {
	    fich << *pdata << "  " ;
		pdata++ ;
	}
	fich << endl ;

}

