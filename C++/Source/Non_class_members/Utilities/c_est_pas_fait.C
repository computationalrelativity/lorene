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


 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:18:11  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:31  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:16:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2008/08/19 06:42:01  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
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
#include <cstdlib>

namespace Lorene {
void c_est_pas_fait(const char * fichier) {

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
}
