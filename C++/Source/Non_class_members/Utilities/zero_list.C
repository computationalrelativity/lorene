/*
 * Function zero_list 
 * Locates approximatively all the zeros of a function in a given interval
 * 
 * (see file utilitaires.h for documentation)
 *
 */


/*
 *   Copyright (c) 2003 Eric Gourgoulhon
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
 * Revision 1.3  2016/12/05 16:18:11  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:53:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2003/09/08 20:22:02  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */


// Headers Lorene 
#include "param.h"
#include "tbl.h"

//****************************************************************************

namespace Lorene {

void zero_list( double (*f)(double, const Param&), const Param& par,
		double xmin, double xmax, int nsub, 
		Tbl*& az, Tbl*& bz ) {
		
    int nzero = 0 ; 
    
    double dx = (xmax-xmin) / double(nsub) ;     
    double f1 = f(xmin, par) ;
    double x1 = xmin ; 
    double x2 = xmin + dx ; 
    
    double* borne_inf = new double[nsub] ;   // At maximum nsub zeros
    double* borne_sup = new double[nsub] ;	    
    
    for (int i=0; i<nsub; i++) {
	double f2 = f(x2, par) ;
	if (f1*f2 < 0.) {	    // A zero has been found
	    borne_inf[nzero] = x1 ; 
	    borne_sup[nzero] = x2 ; 
	    nzero += 1 ;	
	} 
	// Next sub-interval :
	x1 = x2 ; 
	f1 = f2 ;  
	x2 += dx ;
    } 

    // Result:

    az = new Tbl(nzero) ; 
    bz = new Tbl(nzero) ; 
    
    if (nzero > 0) {

	az->set_etat_qcq() ; 
	bz->set_etat_qcq() ; 

	for (int i=0; i<nzero; i++) {
	    az->set(i) = borne_inf[i] ;
	    bz->set(i) = borne_sup[i] ;
	}
    }
    
    delete [] borne_inf ; 
    delete [] borne_sup ; 
    
}  
}
