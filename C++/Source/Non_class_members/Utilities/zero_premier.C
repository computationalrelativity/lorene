/*
 *  Locates the first zero of a function in a given interval.
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
 * Revision 1.4  2016/12/05 16:18:11  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:53:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/01/04  10:57:51  eric
 * Le test f1*f2 < 0. est remplace par f1*f2 <= double(0).
 *
 * Revision 1.1  1999/12/24  13:00:10  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers Lorene 
#include "headcpp.h"
#include "param.h"
//****************************************************************************

namespace Lorene {

bool zero_premier(double (*f)(double, const Param&), const Param& par,
		  double a, double b, int n, double& a0, double& b0) {

    double dx = (b-a)/double(n) ;     

    a0 = a ; 
    b0 = a0 + dx ; 

    double f1 = f(a0, par) ;
    bool trouve = false ; 
    
    for (int i=0; i<n; i++) {
	double f2 = f(b0, par) ;
	if (f1*f2 <= double(0)) {	    // on a encadre le zero
	    trouve = true ; 
	    break ; 
	} 

	// On passe au sous-intervalle suivant :
	a0 = b0 ; 
	f1 = f2 ;  
	b0 += dx ;
    } 
    
    return trouve ; 

}
}
