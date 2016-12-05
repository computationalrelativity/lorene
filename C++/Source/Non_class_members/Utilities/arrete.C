/*
 *  Stops the execution of a code, until the 'Enter' case is hit.
 *
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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
 * Revision 1.4  2016/12/05 16:18:10  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:31  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.1  1999/12/15  15:40:50  eric
 * L'argument par defaut a=0 est desormais precise dans le prototypage
 * declare dans utilitaires.h
 *
 * Revision 2.0  1999/03/02  14:04:18  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// headers du C++
#include"headcpp.h"

namespace Lorene {
void arrete(int a) {
    char c ;
    
    if (a == 0) {
    
	cout.flush() ;
	cout << "Continue = 'return'" << endl ;
	cin.get(c) ;

    }
    
}
}
