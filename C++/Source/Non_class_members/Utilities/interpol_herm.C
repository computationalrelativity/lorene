/*
 * Hermite interpolation of degree 3
 *
 */

/*
 *   Copyright (c) 2000-2002 Eric Gourgoulhon
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


char interpol_herm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/11/21 16:14:51  m_bejger
 * Added the linear interpolation
 *
 * Revision 1.3  2003/05/15 09:42:12  e_gourgoulhon
 * Added the new function interpol_herm_der
 *
 * Revision 1.2  2002/09/09 13:00:40  e_gourgoulhon
 * Modification of declaration of Fortran 77 prototypes for
 * a better portability (in particular on IBM AIX systems):
 * All Fortran subroutine names are now written F77_* and are
 * defined in the new file C++/Include/proto_f77.h.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.0  2000/11/22  19:31:42  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tbl.h"

// Prototypes of F77 subroutines
#include "proto_f77.h"

// Linear interpolation 
void interpol_linear(const Tbl& xtab, const Tbl& ytab, 
                     double x, int& i, double& y) {

	assert(ytab.dim == xtab.dim) ;
	//assert(dytab.dim == xtab.dim) ;	
	
	int np = xtab.get_dim(0) ;
	
	F77_huntm(xtab.t, &np, &x, &i) ;
	
	i-- ; 	// Fortran --> C
	
	int i1 = i + 1 ;
	
	// double dx  = xtab(i1) - xtab(i) ;
        double y1  = ytab(i) ;
        double y2  = ytab(i1) ;

        double x1  = xtab(i) ;
        double x2  = xtab(i1) ;
        double x12 = x1-x2 ;
  

        double a  = (y1-y2)/x12 ;
        double b  = (x1*y2-y1*x2)/x12 ;
	
        y  = x*a+b ; 

}

// Version returning the function and its first derivative 
void interpol_herm(const Tbl& xtab, const Tbl& ytab, const Tbl& dytab,
		   double x, int& i, double& y, double& dy) {

	assert(ytab.dim == xtab.dim) ;
	assert(dytab.dim == xtab.dim) ;	
	
	int np = xtab.get_dim(0) ;
	
	F77_huntm(xtab.t, &np, &x, &i) ;
	
	i-- ; 	// Fortran --> C
	
	int i1 = i + 1 ;
	
	double dx = xtab(i1) - xtab(i) ;

	double u = (x - xtab(i)) / dx ;
	double u2 = u*u ;
	double u3 = u2*u ;
	
	y =   ytab(i) * ( 2.*u3 - 3.*u2 + 1.)
	    + ytab(i1) * ( 3.*u2 - 2.*u3)
     	    + dytab(i) * dx * ( u3 - 2.*u2 + u )
     	    - dytab(i1) * dx * ( u2 - u3 ) ;
     	
 	dy =   6. * ( ytab(i) / dx * ( u2 - u )
 	     - ytab(i1) / dx * ( u2 - u ) )
     	     + dytab(i) * ( 3.*u2 - 4.*u + 1. )
     	     + dytab(i1) * ( 3.*u2 - 2.*u ) ;
	
		   		
}


// Version returning the second derivative 
void interpol_herm_der(const Tbl& xtab, const Tbl& ytab, const Tbl& dytab,
		   double x, int& i, double& y, double& dy, double& ddy) {

	assert(ytab.dim == xtab.dim) ;
	assert(dytab.dim == xtab.dim) ;	
	
	int np = xtab.get_dim(0) ;
	
	F77_huntm(xtab.t, &np, &x, &i) ;
	
	i-- ; 	// Fortran --> C
	
	int i1 = i + 1 ;
	
	double dx = xtab(i1) - xtab(i) ;

	double u = (x - xtab(i)) / dx ;
	double u2 = u*u ;
	double u3 = u2*u ;
	
	y =   ytab(i) * ( 2.*u3 - 3.*u2 + 1.)
	    + ytab(i1) * ( 3.*u2 - 2.*u3)
     	    + dytab(i) * dx * ( u3 - 2.*u2 + u )
     	    - dytab(i1) * dx * ( u2 - u3 ) ;
     	
 	dy =   6. * ( ytab(i) - ytab(i1) ) * ( u2 - u ) / dx 
     	     + dytab(i) * ( 3.*u2 - 4.*u + 1. )
     	     + dytab(i1) * ( 3.*u2 - 2.*u ) ;
	     
	ddy = 6 * ( ( ytab(i) - ytab(i1) ) * ( 2.*u - 1. ) / dx
		+  dytab(i) * (6.*u - 4.)
		+  dytab(i1) * (6.*u - 2.) ) / dx ; 
		   		
}


