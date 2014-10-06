/*
 *  Simple code to search for the enthalpy corresponding
 *  to a given value of a global parameter 
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon
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

char chent_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/06 15:09:45  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.1  2003/09/08 15:12:06  e_gourgoulhon
 * Initial version.
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdio>

int main(int argc, char** argv){

    if (argc < 2) {
        cout <<
        "chent : the value of the global quantity Q for which the enthalpy is searched "
            << endl << " must be given in argument !" << endl ;
	return -1 ; 
    }

    char* qchar = argv[1] ;
    
    double qq ;
    sscanf(qchar, "%lf", &qq) ;
    
    cout << "Search of ent for Q = " << qq << " : " << endl ; 

    double ent1, ent2, qq1,  qq2 ;
    cout << "  ent1, Q1 ? "  ;
    cin >> ent1 ;
    cin >> qq1 ;
    cout << endl << "  ent2, Q2 ? "  ;
    cin >> ent2 ;
    cin >> qq2 ;

    double y1 = qq1 - qq ;
    double y2 = qq2 - qq ;

    double ent = ent1 - y1 /(y2 - y1) * (ent2 - ent1) ;
    cout.precision(16) ;
    cout << endl << "ent for Q = " << qq << " : "
         << ent << endl ;
	
    return 0 ; 

}
