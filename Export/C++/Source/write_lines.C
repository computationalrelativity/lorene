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


 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:31  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:54:06  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2006/09/12 08:04:07  j_novak
 * Removal of the include path Export/C++/Include, updating of the relevant
 * source files in Export/C++/Source.
 *
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

#include "headcpp.h"

namespace Lorene {
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

}
