/*
 * Search for a zero of a function in a given interval, by means of a
 *  secant method.
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
 * Revision 1.7  2016/12/05 16:18:11  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:53:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/07/04 12:09:06  j_novak
 * New argument in zerosec(): a boolean (false by default) for aborting if the number of iteration is greater than the max.
 *
 * Revision 1.4  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/04/11 09:19:46  j_novak
 * Back to old version of zerosec
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.6  2001/10/17  08:16:47  eric
 * In case there is not a single zero in the interval, the found
 * zero is displayed in the warning message.
 *
 * Revision 1.5  2000/01/04  13:20:34  eric
 * Test final f0 != double(0) remplace par fabs(f0) > 1.e-15 .
 *
 * Revision 1.4  1999/12/20  09:46:08  eric
 * Anglicisation des messages.
 *
 * Revision 1.3  1999/12/17  10:08:46  eric
 * Le test final fabs(f0) > 1.e-14 est remplace par f0 != 0.
 *
 * Revision 1.2  1999/12/17  09:37:40  eric
 * Ajout de assert(df != 0).
 *
 * Revision 1.1  1999/12/15  09:41:34  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cmath>
#include <cassert>

// Headers C++
#include <exception>

// Headers Lorene 
#include "headcpp.h"
#include "param.h"
//****************************************************************************

namespace Lorene {

double zerosec(double (*f)(double, const Param&), const Param& parf, 
	       double x1, double x2, double precis, int nitermax, int& niter, 
	       bool abor) {
    
    double f0_prec, f0, x0, x0_prec, dx, df ;

// Teste si un zero unique existe dans l'intervalle [x_1,x_2]

    bool warning = false ; 
    
    f0_prec = f(x1, parf) ;
    f0 = f(x2, parf) ;
    if ( f0*f0_prec > 0.) {
	warning = true ; 
	cout << 
      "WARNING: zerosec: there does not exist a unique zero of the function" 
	<< endl ;
	cout << "  between x1 = " << x1 << " ( f(x1)=" << f0_prec << " )" << endl ; 
	cout << "      and x2 = " << x2 << " ( f(x2)=" << f0 << " )" << endl ;
    }

// Choisit la borne avec la plus petite valeur de |f(x)| comme la valeur la
//  "plus recente" de x0

    if ( fabs(f0) < fabs(f0_prec) ) {  // On a bien choisi f0_prec et f0
	x0_prec = x1 ;
	x0 = x2 ;
    }
    else {  // il faut interchanger f0_prec et f0
	x0_prec = x2 ;
	x0 = x1 ;
	double swap = f0_prec ;
	f0_prec = f0 ;
	f0 = swap ;	
    }

// Debut des iterations de la methode de la secante
    
    niter = 0 ;
    do {
	df = f0 - f0_prec ;
	assert(df != double(0)) ; 
	dx = (x0_prec - x0) * f0 / df ;
	x0_prec = x0 ;
	f0_prec = f0 ;
	x0 += dx ;
	f0 = f(x0, parf) ;
	niter++ ;
	if (niter > nitermax) {
	    cout << "zerosec: Maximum number of iterations has been reached ! " 
	    << endl ;
	    if (abor)
	      abort () ;
	    else {
	      warning = true ;
	      f0 = 0. ;
	    }
	}
    }
    while ( ( fabs(dx) > precis ) && ( fabs(f0) > 1.e-15 ) ) ;

    if (warning) {
	cout << "      A zero may have been found at x0 = " << x0 << endl ; 
    }

    return x0 ;
}  



}
