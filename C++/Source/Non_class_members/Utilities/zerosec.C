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


char zerosec_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/04/08 07:37:29  j_novak
 * zerosec method changed!
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

 
// Headers C++
#include <iostream.h>

// Headers C
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Headers Lorene 
#include "param.h"
//****************************************************************************

double zerosec(double (*f)(double, const Param&), const Param& parf, 
    double x1, double x2, double precis, int nitermax, int& niter) {
    
    double f0_moins, f0, x0, x0_moins, dx, df , fnew, xnew;

// Teste si un zero unique existe dans l'intervalle [x_1,x_2]

    f0_moins = f(x1, parf) ;
    f0 = f(x2, parf) ;
    if ( f0*f0_moins > 0.) {
	cout << 
      "WARNING: zerosec: there does not exist a unique zero of the function" 
	<< endl ;
	cout << "  between x1 = " << x1 << " ( f(x1)=" << f0_moins << " )" << endl ; 
	cout << "      and x2 = " << x2 << " ( f(x2)=" << f0 << " )" << endl ;
	abort() ;
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
	    cout << "zerosec: Maximum number of iterations has been reached ! " 
	    << endl ;
	cout << x0_moins << ", " << xnew << ", " << x0 << endl ;
	cout << f0_moins << ", " << fnew << ", " << f0 << endl ;
	
	    abort () ;
	}
    }
    while ( ( fabs(dx) > precis ) && ( fabs(fnew) > 1.e-15 ) ) ;

    return xnew ;
}  



