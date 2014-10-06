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
 * Revision 1.5  2014/10/06 15:12:54  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2003/01/09 11:07:54  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.3  2001/12/12 12:33:44  e_gourgoulhon
 * Installation instructions now in files INSTALL and INSTALL_linux
 *
 * Revision 1.2  2001/12/12 09:23:46  e_gourgoulhon
 * Parameter compact added to the simplified constructor of class Mg3d
 *
 * Revision 1.1  2001/12/11 06:50:14  e_gourgoulhon
 * test code for Mg3d class
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>

// Lorene headers
#include "headcpp.h"
#include "grilles.h"
#include "type_parite.h"


int main() {

        int nz = 4 ;
        int nr = 17 ;
        int nt = 9 ;
        int np = 4 ;

	{
	  Mg3d mg(nz, nr, nt, np, SYM, NONSYM, true) ;
	  cout << "Multi-grid:" << endl ;
	  cout << mg << endl ;
	}
	
	return EXIT_SUCCESS ;

}

