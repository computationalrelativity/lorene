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
 * Revision 1.5  2011/09/27 15:38:11  j_novak
 * New function for 2D interpolation added. The computation of 1st derivative is
 * still missing.
 *
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


void interpol_herm_2d(const Tbl& xtab, const Tbl& ytab, const Tbl& ftab, 
		      const Tbl& dfdxtab, const Tbl& dfdytab, const Tbl& d2fdxdytab,
		      double x, double y, double& f, double& dfdy) {

  assert(ytab.dim == xtab.dim) ;
  assert(ftab.dim == xtab.dim) ;
  assert(dfdxtab.dim == xtab.dim) ;
  assert(dfdytab.dim == xtab.dim) ;
  assert(d2fdxdytab.dim == xtab.dim) ;
  
  int nbp1, nbp2;
  nbp1 = xtab.get_dim(1);
  nbp2 = xtab.get_dim(0);
  
  int i,j;
  i = 0;
  while ((xtab(i,0) <= x) && (nbp2 > i)) {
    i = i + 1;
  }
  if (i != 0) {
    i = i - 1;
  }
  j = 0;
  while ((ytab(i,j) < y) && (nbp1 > j)) {
    j = j + 1;
  }
  if (j != 0) {
    j = j - 1;
  }
  
  int i1 = i+1 ; int j1 = j+1 ;

  double dx = xtab(i1, j) - xtab(i, j) ;
  double dy = ytab(i, j1) - ytab(i, j) ;

  double u = (x - xtab(i, j)) / dx ;
  double v = (y - ytab(i, j)) / dy ;

  double u2 = u*u ; double v2 = v*v ;
  double u3 = u2*u ; double v3 = v2*v ;

  double psi0_u = 2.*u3 - 3.*u2 + 1. ;
  double psi0_1mu = -2.*u3 + 3.*u2 ;
  double psi1_u = u3 - 2.*u2 + u ;
  double psi1_1mu = -u3 + u2 ;

  double psi0_v = 2.*v3 - 3.*v2 + 1. ;
  double psi0_1mv = -2.*v3 + 3.*v2 ;
  double psi1_v = v3 - 2.*v2 + v ;
  double psi1_1mv = -v3 + v2 ;

  f = ftab(i, j) * psi0_u * psi0_v
    + ftab(i1, j) * psi0_1mu * psi0_v 
    + ftab(i, j1) * psi0_u * psi0_1mv
    + ftab(i1, j1)  * psi0_1mu * psi0_1mv ;

  f += (dfdxtab(i, j) * psi1_u * psi0_v
	- dfdxtab(i1, j) * psi1_1mu * psi0_v
	+ dfdxtab(i, j1) * psi1_u * psi0_1mv
	- dfdxtab(i1, j1) * psi1_1mu * psi0_1mv) * dx ;

  f += (dfdytab(i, j) * psi0_u * psi1_v
	+ dfdytab(i1, j) * psi0_1mu * psi1_v
	- dfdytab(i, j1) * psi0_u * psi1_1mv
	- dfdytab(i1, j1) * psi0_1mu * psi1_1mv) * dy ;
  
  f += (d2fdxdytab(i, j) * psi1_u * psi1_v
	- d2fdxdytab(i1, j) * psi1_1mu * psi1_v
	- d2fdxdytab(i, j1) * psi1_u * psi1_1mv 
	+ d2fdxdytab(i1, j1) * psi1_1mu * psi1_1mv) * dx * dy ;


  return ;

}
