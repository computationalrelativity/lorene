/*
 * Stops the execution. In debug mode: ask for continuation. In normal mode:
 * abort the main program.
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
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


char c_est_pas_fait_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2001/11/29 16:17:54  e_gourgoulhon
 * minor modifs
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/01/15  09:10:39  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/01/15  09:10:39  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers du C++
#include"headcpp.h"
#include <stdlib.h>

void c_est_pas_fait(char * fichier) {

#ifdef NDEBUG
    cout.flush() ;
    cout << "Routine not ready in " << fichier << " !" << endl ;
    abort() ;
#else
    cout.flush() ;
    cout << "Routine not ready in " << fichier << " !" << endl ;
    cout << "Next = 'return',  abort = '0'" << endl ;
    char c ;
    cin.get(c) ;
    if (c == '0') {
	abort() ;
    }
#endif    
}
