/*
 * Search for a zero of a function in a given interval, by means of a
 *  secant method. The input parameters x1 and x2 must be such that 
 *  f(x1)*f(x2) < 0 . The obtained root is then necessarily in the 
 *  interval [x1,x2].
 *
 */

/*
 *   Copyright (c) 2002 Jerome Novak
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


char zerosec_borne_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/06 15:16:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.4  2002/05/02 15:16:22  j_novak
 * Added functions for more general bi-fluid EOS
 *
 * Revision 1.3  2002/04/18 09:17:17  j_novak
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */
 
// Headers C
#include <cstdlib>
#include <cmath>
#include <cassert>

// Headers Lorene 
#include "headcpp.h"
#include "param.h"
//****************************************************************************

double zerosec_b(double (*f)(double, const Param&), const Param& parf, 
    double x1, double x2, double precis, int nitermax, int& niter) {
    
    double f0_moins, f0, x0, x0_moins, dx, df , fnew, xnew;

// Teste si un zero unique existe dans l'intervalle [x_1,x_2]

    f0_moins = f(x1, parf) ;
    f0 = f(x2, parf) ;
    if ( f0*f0_moins > 0.) {
	cout << 
      "WARNING: zerosec: there may not exist a zero of the function" 
	<< endl ;
	cout << "  between x1 = " << x1 << " ( f(x1)=" << f0_moins << " )" << endl ; 
	cout << "      and x2 = " << x2 << " ( f(x2)=" << f0 << " )" << endl ;
    }

// Choisit la borne avec la plus grande valeur de f(x) (borne positive) 
//  comme la valeur la de x0

    if ( f0_moins < f0) {  // On a bien choisi f0_moins et f0
	x0_moins = x1 ;
	x0 = x2 ;
    }
    else {  // il faut interchanger f0_moins et f0
	x0_moins = x2 ;
	x0 = x1 ;
	double swap = f0_moins ;
	f0_moins = f0 ;
	f0 = swap ;	
    }

// Debut des iterations de la methode de la secante
    
    niter = 0 ;
    do {
	df = f0 - f0_moins ;
	assert(df != double(0)) ; 
	xnew = (x0_moins*f0 - x0*f0_moins)/df ; ;
	fnew = f(xnew, parf) ;
	if (fnew < 0.) {
	  dx = x0_moins - xnew ;
	  x0_moins = xnew ;
	  f0_moins = fnew ;
	}
	else {
	  dx = x0 - xnew ; 
	  x0 = xnew ;
	  f0 = fnew ;
	}
	niter++ ;
	if (niter > nitermax) {
//cout << "zerosec: Maximum number of iterations has been reached ! " 
	  //	    << endl ;
	  //cout << x0_moins << ", " << xnew << ", " << x0 << endl ;
	  //cout << f0_moins << ", " << fnew << ", " << f0 << endl ;
	
	    return xnew ;
	}
    }
    while ( ( fabs(dx) > precis ) && ( fabs(fnew) > 1.e-15 ) ) ;

    return xnew ;
}  



