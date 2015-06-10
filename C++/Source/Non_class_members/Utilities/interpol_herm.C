/*
 * Hermite interpolation functions.
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
 * Revision 1.12  2015/06/10 14:39:18  a_sourie
 * New class Eos_bf_tabul for tabulated 2-fluid EoSs and associated functions for the computation of rotating stars with such EoSs.
 *
 * Revision 1.11  2015/01/27 14:22:38  j_novak
 * New methods in Eos_tabul to correct for EoS themro consistency (optional).
 *
 * Revision 1.10  2014/10/13 08:53:32  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2013/12/12 16:07:30  j_novak
 * interpol_herm_2d outputs df/dx, used to get the magnetization.
 *
 * Revision 1.8  2012/09/04 14:53:28  j_novak
 * Replacement of the FORTRAN version of huntm by a C one.
 *
 * Revision 1.7  2011/10/04 16:05:19  j_novak
 * Update of Eos_mag class. Suppression of loge, re-definition of the derivatives
 * and use of interpol_herm_2d.
 *
 * Revision 1.6  2011/10/03 13:44:45  j_novak
 * Updated the y-derivative for the 2D version
 *
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
#include "math.h"

namespace Lorene {

  //---------------------------------------------------------------
  // Value bracketting in an ordered table (from Numerical Recipes)
  //---------------------------------------------------------------
void huntm(const Tbl& xx, double& x, int& i_low) {

  assert (xx.get_etat() == ETATQCQ) ;
  int nx = xx.get_taille() ;
  bool ascend = ( xx(nx-1) > xx(0) ) ;
  int i_hi ;
  if ( (i_low < 0)||(i_low>=nx) ) {
    i_low = -1 ;
    i_hi = nx ;
  }
  else {
    int inc = 1 ;
    if ( (x >= xx(i_low)) == ascend ) {
	if (i_low == nx -1) return ;
	i_hi = i_low + 1 ;
	while ( (x >= xx(i_hi)) == ascend ) {
	  i_low = i_hi ;
	  inc += inc ;
	  i_hi = i_low + inc ;
	  if (i_hi >= nx) {
	    i_hi = nx ;
	    break ;
	  }
	}
    } else {
      if (i_low == 0) {
	i_low = -1 ;
	return ;
      }
      i_hi = i_low-- ;
      while ( (x < xx(i_low)) == ascend ) {
	i_hi = i_low ;
	inc += inc ;
	if ( inc >= i_hi ) {
	  i_low = 0 ;
	  break ;
	}
	else i_low = i_hi - inc ;
      }
    }
  }
  while ( (i_hi - i_low) > 1) {
    int i_med = (i_hi + i_low) / 2 ;
    if ( (x>=xx(i_med)) == ascend ) i_low = i_med ;
    else i_hi = i_med ;
  }
  if (x == xx(nx-1)) i_low = nx-2 ;
  if (x == xx(0)) i_low = 0 ;
  return ;
}

  //---------------------
  // Linear interpolation 
  //---------------------
void interpol_linear(const Tbl& xtab, const Tbl& ytab, 
                     double x, int& i, double& y) {

	assert(ytab.dim == xtab.dim) ;
	//assert(dytab.dim == xtab.dim) ;	
	
	huntm(xtab, x, i) ;
	
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

  //------------------------------------------------------------
  // Cubic Hermite interpolation, returning the first derivative
  //------------------------------------------------------------
void interpol_herm(const Tbl& xtab, const Tbl& ytab, const Tbl& dytab,
		   double x, int& i, double& y, double& dy) {

	assert(ytab.dim == xtab.dim) ;
	assert(dytab.dim == xtab.dim) ;	

	huntm(xtab, x, i) ;
	
	int i1 = i + 1 ;
	
	double dx = xtab(i1) - xtab(i) ;

	double u = (x - xtab(i)) / dx ;
	double u2 = u*u ;
	double u3 = u2*u ;
	
	y =   ytab(i) * ( 2.*u3 - 3.*u2 + 1.)
	    + ytab(i1) * ( 3.*u2 - 2.*u3)
     	    + dytab(i) * dx * ( u3 - 2.*u2 + u )
     	    - dytab(i1) * dx * ( u2 - u3 ) ;
   	//cout << "X = " << xtab(i) << "    " << x << "    " << xtab(i+1) << endl;
   	//cout << "Y = " << ytab(i) << "    " << y << "    " << ytab(i+1) << endl;
 	dy =   6. * ( ytab(i) / dx * ( u2 - u )
 	     - ytab(i1) / dx * ( u2 - u ) )
     	     + dytab(i) * ( 3.*u2 - 4.*u + 1. )
     	     + dytab(i1) * ( 3.*u2 - 2.*u ) ;
  	//cout << "DY = " <<  dytab(i) << "    " << dy << "    " << dytab(i+1) << endl;
		//cout << xtab(i) << "   " << x << "   " << xtab(i1) << endl ;   		
		//cout << ytab(i) << "   " << y << "   " << ytab(i1) << endl ;
}


  //-------------------------------------------------------------
  // Cubic Hermite interpolation, returning the second derivative
  //-------------------------------------------------------------
void interpol_herm_der(const Tbl& xtab, const Tbl& ytab, const Tbl& dytab,
		   double x, int& i, double& y, double& dy, double& ddy) {

	assert(ytab.dim == xtab.dim) ;
	assert(dytab.dim == xtab.dim) ;	
		
	huntm(xtab, x, i) ;
	
	//	i-- ; 	// Fortran --> C
	
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

  //----------------------------------------------
  // Bi-cubic Hermite interpolation, for 2D arrays
  //----------------------------------------------
void interpol_herm_2d(const Tbl& xtab, const Tbl& ytab, const Tbl& ftab, 
		      const Tbl& dfdxtab, const Tbl& dfdytab, const Tbl& d2fdxdytab,
		      double x, double y, double& f, double& dfdx, double& dfdy) {

  assert(ytab.dim == xtab.dim) ;
  assert(ftab.dim == xtab.dim) ;
  assert(dfdxtab.dim == xtab.dim) ;
  assert(dfdytab.dim == xtab.dim) ;
  assert(d2fdxdytab.dim == xtab.dim) ;
  
  int nbp1, nbp2;
  nbp1 = xtab.get_dim(0);
  nbp2 = xtab.get_dim(1);
  
  int i_near = 0 ;
  int j_near = 0 ;

  while ((xtab(i_near,0) <= x) && (nbp2 > i_near)) {
    i_near++; 
  }
  if (i_near != 0) {
    i_near-- ; 
  }
  j_near = 0;
  while ((ytab(i_near,j_near) < y) && (nbp1 > j_near)) {
    j_near++ ;
  }
  if (j_near != 0) {
    j_near-- ; 
  }
  
  int i1 = i_near+1 ; int j1 = j_near+1 ;

  double dx = xtab(i1, j_near) - xtab(i_near, j_near) ;
  double dy = ytab(i_near, j1) - ytab(i_near, j_near) ;

  double u = (x - xtab(i_near, j_near)) / dx ;
  double v = (y - ytab(i_near, j_near)) / dy ;

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

  f = ftab(i_near, j_near) * psi0_u * psi0_v
    + ftab(i1, j_near) * psi0_1mu * psi0_v 
    + ftab(i_near, j1) * psi0_u * psi0_1mv
    + ftab(i1, j1)  * psi0_1mu * psi0_1mv ;


  f += (dfdxtab(i_near, j_near) * psi1_u * psi0_v
	- dfdxtab(i1, j_near) * psi1_1mu * psi0_v
	+ dfdxtab(i_near, j1) * psi1_u * psi0_1mv
	- dfdxtab(i1, j1) * psi1_1mu * psi0_1mv) * dx ;

  f += (dfdytab(i_near, j_near) * psi0_u * psi1_v
	+ dfdytab(i1, j_near) * psi0_1mu * psi1_v
	- dfdytab(i_near, j1) * psi0_u * psi1_1mv
	- dfdytab(i1, j1) * psi0_1mu * psi1_1mv) * dy ;
  
  f += (d2fdxdytab(i_near, j_near) * psi1_u * psi1_v
	- d2fdxdytab(i1, j_near) * psi1_1mu * psi1_v
	- d2fdxdytab(i_near, j1) * psi1_u * psi1_1mv 
	+ d2fdxdytab(i1, j1) * psi1_1mu * psi1_1mv) * dx * dy ;

	
  double dpsi0_u = 6.*(u2 - u) ;
  double dpsi0_1mu = 6.*(u2 - u) ;
  double dpsi1_u = 3.*u2 - 4.*u + 1. ;
  double dpsi1_1mu = 3.*u2 - 2.*u ;

  dfdx = (ftab(i_near, j_near) * dpsi0_u * psi0_v
    - ftab(i1, j_near) * dpsi0_1mu * psi0_v 
    + ftab(i_near, j1) * dpsi0_u * psi0_1mv
	  - ftab(i1, j1)  * dpsi0_1mu * psi0_1mv ) / dx;

  dfdx += (dfdxtab(i_near, j_near) * dpsi1_u * psi0_v
	+ dfdxtab(i1, j_near) * dpsi1_1mu * psi0_v
	+ dfdxtab(i_near, j1) * dpsi1_u * psi0_1mv
	+ dfdxtab(i1, j1) * dpsi1_1mu * psi0_1mv) ;

  dfdx += (dfdytab(i_near, j_near) * dpsi0_u * psi1_v
	- dfdytab(i1, j_near) * dpsi0_1mu * psi1_v
	- dfdytab(i_near, j1) * dpsi0_u * psi1_1mv
	+ dfdytab(i1, j1) * dpsi0_1mu * psi1_1mv) * dy /dx ;
  
  dfdx += (d2fdxdytab(i_near, j_near) * dpsi1_u * psi1_v
	+ d2fdxdytab(i1, j_near) * dpsi1_1mu * psi1_v
	- d2fdxdytab(i_near, j1) * dpsi1_u * psi1_1mv 
	- d2fdxdytab(i1, j1) * dpsi1_1mu * psi1_1mv) * dy ;


  double dpsi0_v = 6.*(v2 - v) ;
  double dpsi0_1mv = 6.*(v2 - v) ;
  double dpsi1_v = 3.*v2 - 4.*v + 1. ;
  double dpsi1_1mv = 3.*v2 - 2.*v ;


  dfdy = (ftab(i_near, j_near) * psi0_u * dpsi0_v
    + ftab(i1, j_near) * psi0_1mu * dpsi0_v 
    - ftab(i_near, j1) * psi0_u * dpsi0_1mv
	  - ftab(i1, j1)  * psi0_1mu * dpsi0_1mv) / dy ;

  dfdy += (dfdxtab(i_near, j_near) * psi1_u * dpsi0_v
	- dfdxtab(i1, j_near) * psi1_1mu * dpsi0_v
	- dfdxtab(i_near, j1) * psi1_u * dpsi0_1mv
	+ dfdxtab(i1, j1) * psi1_1mu * dpsi0_1mv) * dx / dy ;

  dfdy += (dfdytab(i_near, j_near) * psi0_u * dpsi1_v
	+ dfdytab(i1, j_near) * psi0_1mu * dpsi1_v
	+ dfdytab(i_near, j1) * psi0_u * dpsi1_1mv
	+ dfdytab(i1, j1) * psi0_1mu * dpsi1_1mv) ;
  
  dfdy += (d2fdxdytab(i_near, j_near) * psi1_u * dpsi1_v
	- d2fdxdytab(i1, j_near) * psi1_1mu * dpsi1_v
	+ d2fdxdytab(i_near, j1) * psi1_u * dpsi1_1mv 
	- d2fdxdytab(i1, j1) * psi1_1mu * dpsi1_1mv) * dx ;

  return ;

}






void interpol_herm_2d_sans(const Tbl& xtab, const Tbl& ytab, const Tbl& ftab, 
		      const Tbl& dfdxtab, const Tbl& dfdytab, 
		      double x, double y, double& f, double& dfdx, double& dfdy) {

  assert(ytab.dim == xtab.dim) ;
  assert(ftab.dim == xtab.dim) ;
  assert(dfdxtab.dim == xtab.dim) ;
  assert(dfdytab.dim == xtab.dim) ;
  
  int nbp1, nbp2;
  nbp1 = xtab.get_dim(0);
  nbp2 = xtab.get_dim(1);
  
  int i_near = 0 ;
  int j_near = 0 ;

  while ((xtab(i_near,0) <= x) && (nbp2 > i_near)) {
    i_near++; 
  }
  if (i_near != 0) {
    i_near-- ; 
  }
  j_near = 0;
  while ((ytab(i_near,j_near) < y) && (nbp1 > j_near)) {
    j_near++ ;
  }
  if (j_near != 0) {
    j_near-- ; 
  }
  
  int i1 = i_near+1 ; int j1 = j_near+1 ;

  double dx = xtab(i1, j_near) - xtab(i_near, j_near) ;
  double dy = ytab(i_near, j1) - ytab(i_near, j_near) ;

  double u = (x - xtab(i_near, j_near)) / dx ;
  double v = (y - ytab(i_near, j_near)) / dy ;

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

  f = ftab(i_near, j_near) * psi0_u * psi0_v
    + ftab(i1, j_near) * psi0_1mu * psi0_v 
    + ftab(i_near, j1) * psi0_u * psi0_1mv
    + ftab(i1, j1)  * psi0_1mu * psi0_1mv ;


  f += (dfdxtab(i_near, j_near) * psi1_u * psi0_v
	- dfdxtab(i1, j_near) * psi1_1mu * psi0_v
	+ dfdxtab(i_near, j1) * psi1_u * psi0_1mv
	- dfdxtab(i1, j1) * psi1_1mu * psi0_1mv) * dx ;

  f += (dfdytab(i_near, j_near) * psi0_u * psi1_v
	+ dfdytab(i1, j_near) * psi0_1mu * psi1_v
	- dfdytab(i_near, j1) * psi0_u * psi1_1mv
	- dfdytab(i1, j1) * psi0_1mu * psi1_1mv) * dy ;

  double dpsi0_u = 6.*(u2 - u) ;
  double dpsi0_1mu = 6.*(u2 - u) ;
  double dpsi1_u = 3.*u2 - 4.*u + 1. ;
  double dpsi1_1mu = 3.*u2 - 2.*u ;

  dfdx = (ftab(i_near, j_near) * dpsi0_u * psi0_v
    - ftab(i1, j_near) * dpsi0_1mu * psi0_v 
    + ftab(i_near, j1) * dpsi0_u * psi0_1mv
	  - ftab(i1, j1)  * dpsi0_1mu * psi0_1mv ) / dx;

  dfdx += (dfdxtab(i_near, j_near) * dpsi1_u * psi0_v
	+ dfdxtab(i1, j_near) * dpsi1_1mu * psi0_v
	+ dfdxtab(i_near, j1) * dpsi1_u * psi0_1mv
	+ dfdxtab(i1, j1) * dpsi1_1mu * psi0_1mv) ;

  dfdx += (dfdytab(i_near, j_near) * dpsi0_u * psi1_v
	- dfdytab(i1, j_near) * dpsi0_1mu * psi1_v
	- dfdytab(i_near, j1) * dpsi0_u * psi1_1mv
	+ dfdytab(i1, j1) * dpsi0_1mu * psi1_1mv) * dy /dx ;


  double dpsi0_v = 6.*(v2 - v) ;
  double dpsi0_1mv = 6.*(v2 - v) ;
  double dpsi1_v = 3.*v2 - 4.*v + 1. ;
  double dpsi1_1mv = 3.*v2 - 2.*v ;


  dfdy = (ftab(i_near, j_near) * psi0_u * dpsi0_v
    + ftab(i1, j_near) * psi0_1mu * dpsi0_v 
    - ftab(i_near, j1) * psi0_u * dpsi0_1mv
	  - ftab(i1, j1)  * psi0_1mu * dpsi0_1mv) / dy ;

  dfdy += (dfdxtab(i_near, j_near) * psi1_u * dpsi0_v
	- dfdxtab(i1, j_near) * psi1_1mu * dpsi0_v
	- dfdxtab(i_near, j1) * psi1_u * dpsi0_1mv
	+ dfdxtab(i1, j1) * psi1_1mu * dpsi0_1mv) * dx / dy ;

  dfdy += (dfdytab(i_near, j_near) * psi0_u * dpsi1_v
	+ dfdytab(i1, j_near) * psi0_1mu * dpsi1_v
	+ dfdytab(i_near, j1) * psi0_u * dpsi1_1mv
	+ dfdytab(i1, j1) * psi0_1mu * dpsi1_1mv) ;
  

  return ;

}



/*interpolation for functions of 3 variables : hermite interpolation on the 2 last variables and linear interpolation on the first one*/

/*
xtab/x refers to delta_car (logdelta_car)
ytab/y refers to mu_n (logent1)
ztab/z refers to mu_p (logent2)
ftab/f refers to psi (logp) or alpha (dlpsddelta_car or logalpha)
dfdytab/dfdy refers to dpsi/dmu_n (dlpsdlent1) or dalpha/dmu_n (d2lpsdlent1ddelta_car or dlalphadlent1)
dfdztab/dfdz refers to dpsi/dmu_p (dlpsdlent2) or dalpha/dmu_p (dl2psdlent2ddelta_car or dlalphadlent2)
d2fdytabdztab refers to d2psi/dmu_ndmu_p (d2lpsdlent1dlent2) or d2alpha/dmu_ndmu_p (d3lpsdlent1dlent2ddelta_car or d2lalphadlent1dlent2)
*/

/*this routine provides the interpolated values of f, dfdy and dfdz at point (x,y,z) via the use of adapted tables*/

void interpol_mixed_3d(const Tbl& xtab, const Tbl& ytab, const Tbl& ztab, const Tbl& ftab, 
		      const Tbl& dfdytab, const Tbl& dfdztab, const Tbl& d2fdydztab,
		      double x, double y, double z, double& f, double& dfdy, double& dfdz) 


{
  assert(ytab.dim == xtab.dim) ; 
  assert(ztab.dim == xtab.dim) ;
  assert(ftab.dim == xtab.dim) ;
  assert(dfdytab.dim == xtab.dim) ;
  assert(dfdztab.dim == xtab.dim) ;
  assert(d2fdydztab.dim == xtab.dim) ;
  
  int nbp1, nbp2, nbp3;
  nbp1 = xtab.get_dim(2) ; // \Delta^{2}
  nbp2 = xtab.get_dim(1) ; // \mu_n
  nbp3 = xtab.get_dim(0) ; // \mu_p


  /* we can put these tests directly in the code (before calling interpol_mixed_3d) 
  assert(x >= xtab(0,0,0)) ;
  assert(x <= xtab(nbp1,0,0)) ;
  assert(y >= ytab(0,0,0)) ;
  assert(y <= ytab(0,nbp2,0)) ;
  assert(z >= ztab(0,0,0)) ;
  assert(z <= ztab(0,0,nbp3)) ;
  */

  int i_near = 0 ; 
  int j_near = 0 ;
  int k_near = 0 ;


  // look for the positions of (x,y,z) in the tables
  while ( ( xtab(i_near,0,0) <= x ) && ( ( nbp1-1 ) > i_near ) ) {
      i_near++ ;
  }
  if (i_near != 0) { 
      i_near -- ; 
  }

  while ( ( ytab(i_near,j_near, 0) <= y ) && ( ( nbp2-1 ) > j_near ) ) {
      j_near++ ;
  }
  if (j_near != 0) {
      j_near -- ; 
  }

  while ( ( ztab( i_near, j_near, k_near) <= z) && ( ( nbp3-1 ) > k_near ) ) {
      k_near++ ;
  }
  if (k_near != 0) {
      k_near-- ; 
  }

  int i1 = i_near + 1 ;
  int j1 = j_near + 1 ;
  int k1 = k_near + 1 ;
//cout << j_near<< "   " << k_near<< endl;
  // Hermite interpolation on the slice i_near 
/*cout << "y_dessous= " << ytab(i_near, j_near, k_near) << endl ;
 cout << "y= " << y << endl ;
cout << "y_dessus= " << ytab(i_near, j1, k_near) << endl ;
 cout << "z_dessus= " << ztab(i_near, j_near, k1) << endl ;
 cout << "z= " << z << endl ;
 cout << "z_dessous= " << ztab(i_near, j_near, k_near) << endl ;
*/
//cout << j_near << "   " << k_near << endl ;
 

/*
bool hum1 =   (( ytab(i_near, j1, k_near) - ytab(i_near, j_near, k_near) ) == (ytab(i_near, j1, k1) - ytab(i_near, j_near, k1) ) ) ;
bool hum2 = ((ztab(i_near, j_near, k1) - ztab(i_near, j_near, k_near) ) == ( ztab(i_near, j1, k1) - ztab(i_near, j1, k_near) ))  ;
if (hum1 == false ) { 
cout << setprecision(16);
cout <<  ytab(i_near, j1, k_near) - ytab(i_near, j_near, k_near) <<  "   " << ytab(i_near, j1, k1) - ytab(i_near, j_near, k1)<< endl ;
cout << j_near << "    " << k_near << endl ;
}

if (hum2 == false ) { 
cout << setprecision(16);
cout <<  ztab(i_near, j_near, k1) - ztab(i_near, j_near, k_near) <<  "   " << ztab(i_near, j1, k1) - ztab(i_near, j1, k_near)<< endl ;
cout << j_near << "    " << k_near << endl ;
}*/

  double dy_i_near = ytab(i_near, j1, k_near) - ytab(i_near, j_near, k_near) ;
  double dz_i_near = ztab(i_near, j_near, k1) - ztab(i_near, j_near, k_near) ;
  double u_i_near = (y - ytab(i_near, j_near, k_near)) / dy_i_near ;
  double v_i_near = (z - ztab(i_near, j_near, k_near)) / dz_i_near ;
/*
//lineaire
  double dy_anal = 1. ; // log(2.) ;
  double dz_anal = 1. ; //log(2.) ;
//en log
cout << y << "   " << z << endl ;
  double u_i_anal = y /log(2.) ; //(y - 1.) / 1.
  double v_i_anal = z /log(2.) ; //(z - 1.) / 1. 

cout << setprecision(16);
cout <<"interpol  " <<  u_i_near <<"    " << u_i_anal<< "   "<<  v_i_near <<  "   "  << v_i_anal << endl ;
*/
  double u2_i_near = u_i_near * u_i_near ; 
  double v2_i_near = v_i_near * v_i_near ;
  double u3_i_near = u2_i_near * u_i_near ; 
  double v3_i_near = v2_i_near * v_i_near ;

 /* 
//lineaire
double u2_i_anal =  (y - 1. ) * (y - 1. ) ; //
  double v2_i_anal =  (z - 1.) * (z - 1.) ; // 
  double u3_i_anal = (y - 1. ) * (y - 1. )* (y - 1. )  ; //
  double v3_i_anal = (z - 1.) * (z - 1.)*  (z - 1.)  ; // 

cout << " interpol2" << u2_i_near << "    "  << u2_i_anal << "    " << u3_i_near  << "       "  << u3_i_anal << endl ;
*/

  double psi0_u_i_near = double(2)*u3_i_near - double(3)*u2_i_near + double(1) ;
  double psi0_1mu_i_near = -double(2)*u3_i_near + double(3)*u2_i_near ;
  double psi1_u_i_near = u3_i_near - double(2)*u2_i_near + u_i_near ;
  double psi1_1mu_i_near = -u3_i_near + u2_i_near ;
  double psi0_v_i_near = double(2)*v3_i_near - double(3)*v2_i_near + double(1) ;
  double psi0_1mv_i_near = -double(2)*v3_i_near + double(3)*v2_i_near ;
  double psi1_v_i_near = v3_i_near - double(2)*v2_i_near + v_i_near ;
  double psi1_1mv_i_near = -v3_i_near + v2_i_near ;

  double f_i_near ;
/*// lineaire
 double psi0_u_i_anal = double(2)*(y - 1. ) * (y - 1. )* (y - 1. ) - double(3)*(y - 1. ) *(y - 1. )  + double(1) ;
  double psi0_1mu_i_anal = -double(2)* (y - 1. ) * (y - 1. )* (y - 1. ) + double(3)*(y - 1. ) * (y - 1. );
  double psi1_u_i_anal = (y - 1. ) * (y - 1. )* (y - 1. ) - double(2)*(y - 1. ) * (y - 1.) + (y - 1. ) ;
  double psi1_1mu_i_anal = -(y - 1. ) * (y - 1. )* (y - 1. ) +(y - 1. ) * (y - 1. ) ;
//en log
 double psi0_u_i_anal = double(2)*y /log(2.) * y /log(2.) * y /log(2.) - double(3)*  y /log(2.) * y /log(2.)  + double(1) ;
  double psi0_1mu_i_anal = -double(2)*y /log(2.)*y /log(2.)*y /log(2.) + double(3)*y /log(2.)*y /log(2.);
  double psi1_u_i_anal = y /log(2.)*y /log(2.)*y /log(2.) - double(2)* y /log(2.) * y /log(2.) + y /log(2.);
  double psi1_1mu_i_anal = -y /log(2.)*y /log(2.)*y /log(2.) +y /log(2.)*y /log(2.) ;
// lineaire
 double psi0_v_i_anal = double(2)*(z - 1. ) * (z - 1. )* (z - 1. ) - double(3)*(z - 1. ) *(z - 1. )  + double(1) ;
  double psi0_1mv_i_anal = -double(2)* (z - 1. ) * (z - 1. )* (z - 1. ) + double(3)*(z - 1. ) * (z - 1. );
  double psi1_v_i_anal = (z - 1. ) * (z - 1. )* (z - 1. ) - double(2)*(z - 1. ) * (z - 1.) + (z - 1. ) ;
  double psi1_1mv_i_anal = -(z - 1. ) * (z - 1. )* (z - 1. ) +(z - 1. ) * (z - 1. ) ;

cout << " Interpol3   " << psi0_u_i_near << "    " <<  psi0_u_i_anal << "    " << psi0_1mu_i_near << "  " << psi0_1mu_i_anal << "   " << 
psi1_u_i_near << "  " << psi1_u_i_anal  << "  " <<  psi1_1mu_i_near << "  " <<  psi1_1mu_i_anal << endl;

cout << " Interpol4   " << psi0_v_i_near << "    " <<  psi0_v_i_anal << "    " << psi0_1mv_i_near << "  " << psi0_1mv_i_anal << "   " << 
psi1_v_i_near << "  " << psi1_v_i_anal  << "  " <<  psi1_1mv_i_near << "  " <<  psi1_1mv_i_anal << endl;

cout<< setprecision(16);
cout << "u_i_near " << u_i_near << "  psi1_u_i_near  " << psi1_u_i_near << endl;
cout<< " ftab(i_near, j_near, k_near) " << ftab(i_near, j_near, k_near)  << "  log ftab(i_near, j_near, k_near)  " << log(ftab(i_near, j_near, k_near) ) << endl ;
cout<< " dfdytab(i_near, j1, k_near) " << dfdytab(i_near, j1, k_near)   << "  log dfdytab(i_near, j1, k_near)   " << dfdytab(i_near, j1, k_near) / ftab(i_near, j1, k_near)  << endl ;
cout<< " dfdztab(i_near, j1, k_near) " << dfdztab(i_near, j1, k_near)   << "  log dfdztab(i_near, j1, k_near)   " << dfdztab(i_near, j1, k_near) / ftab(i_near, j1, k_near)  << endl ;

*/

  f_i_near = ftab(i_near, j_near, k_near) * psi0_u_i_near * psi0_v_i_near
    + ftab(i_near, j1, k_near) * psi0_1mu_i_near * psi0_v_i_near
    + ftab(i_near, j_near, k1) * psi0_u_i_near * psi0_1mv_i_near
    + ftab(i_near, j1, k1)  * psi0_1mu_i_near * psi0_1mv_i_near ;


  f_i_near += (dfdytab(i_near, j_near, k_near) * psi1_u_i_near * psi0_v_i_near
	- dfdytab(i_near, j1, k_near) * psi1_1mu_i_near * psi0_v_i_near
	+ dfdytab(i_near, j_near, k1) * psi1_u_i_near * psi0_1mv_i_near
	- dfdytab(i_near, j1, k1) * psi1_1mu_i_near * psi0_1mv_i_near) * dy_i_near ;

  f_i_near += (dfdztab(i_near, j_near, k_near) * psi0_u_i_near * psi1_v_i_near
	+ dfdztab(i_near, j1, k_near) * psi0_1mu_i_near * psi1_v_i_near
	- dfdztab(i_near, j_near, k1) * psi0_u_i_near * psi1_1mv_i_near
	- dfdztab(i_near, j1, k1) * psi0_1mu_i_near * psi1_1mv_i_near) * dz_i_near ;

  f_i_near += (d2fdydztab(i_near, j_near, k_near) * psi1_u_i_near * psi1_v_i_near
	- d2fdydztab(i_near, j1, k_near) * psi1_1mu_i_near * psi1_v_i_near
	- d2fdydztab(i_near, j_near, k1) * psi1_u_i_near * psi1_1mv_i_near 
	+ d2fdydztab(i_near, j1, k1) * psi1_1mu_i_near * psi1_1mv_i_near) * dy_i_near * dz_i_near ;


//cout << f_i_near << endl;


/*
cout << " Y " << y << endl;
cout << " Z " << z << endl;
cout << "dy " << dy_i_near << "  " <<  ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))<< endl ;
cout << "dz " << dz_i_near << "  " <<  log(ztab(i_near, j_near, k1))-log(ztab(i_near, j_near, k_near))<< endl ;

cout << " u  " << u_i_near << "  " <<  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))<< endl ;
cout << " v  " << v_i_near << "  " <<  (log(z)-log(ztab(i_near, j_near, k_near)))/(log(ztab(i_near, j_near, k1))-log(ztab(i_near, j_near, k_near)))<< endl ;

cout << " u2  " << u2_i_near << "  " <<  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))  * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) << endl ;
cout << " psi0_u_i_near " << psi0_u_i_near << "  " <<  double(2)* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) - (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) *  double(3) +double(1) << endl ;

cout << "psi0_1mu_i_near" << psi0_1mu_i_near << "  " << -double(2) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) + double(3)* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))*  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))<< endl ;

cout <<"psi1_u_i_near " << psi1_u_i_near  << "   " <<  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))*  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) - (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) *  double(2) + (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) << endl ;


cout << "psi1_1mu_i_near " <<  psi1_1mu_i_near << "   " << - (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))*  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) + (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))*  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) << endl ;


cout << " ftab(i_near, j_near, k_near) " << ftab(i_near, j_near, k_near) << "   " << log(ftab(i_near, j_near, k_near)) << endl;
cout << " dfdytab(i_near, j_near, k_near) " << dfdytab(i_near, j_near, k_near) <<  "   " << dfdytab(i_near, j_near, k_near)  * ytab(i_near, j_near, k_near) / ftab(i_near, j_near, k_near) << endl ;
cout << "(d2fdydztab(i_near, j_near, k_near)  " << d2fdydztab(i_near, j_near, k_near) << "    " << 
	ytab(i_near, j_near, k_near)  * ztab(i_near, j_near, k_near) / ftab(i_near, j_near, k_near) * ( d2fdydztab(i_near, j_near, k_near) - dfdytab(i_near, j_near, k_near) * dfdztab(i_near, j_near, k_near) / 
 ftab(i_near, j_near, k_near) ) << endl ;
cout << "  f = " << f_i_near << "   " << log(f_i_near) << endl;




*/

  double dpsi0_u_i_near = 6.*(u2_i_near - u_i_near) ;
  double dpsi0_1mu_i_near = 6.*(u2_i_near - u_i_near) ;
  double dpsi1_u_i_near = 3.*u2_i_near - 4.*u_i_near + 1. ;
  double dpsi1_1mu_i_near = 3.*u2_i_near - 2.*u_i_near ;

  double dfdy_i_near;
 
  dfdy_i_near = (ftab(i_near, j_near, k_near) * dpsi0_u_i_near * psi0_v_i_near
	  - ftab(i_near, j1, k_near) * dpsi0_1mu_i_near * psi0_v_i_near
	  + ftab(i_near, j_near, k1) * dpsi0_u_i_near * psi0_1mv_i_near
	  - ftab(i_near, j1, k1)  * dpsi0_1mu_i_near * psi0_1mv_i_near ) / dy_i_near;

  dfdy_i_near += (dfdytab(i_near, j_near, k_near) * dpsi1_u_i_near * psi0_v_i_near
	+ dfdytab(i_near, j1, k_near) * dpsi1_1mu_i_near * psi0_v_i_near
	+ dfdytab(i_near, j_near, k1) * dpsi1_u_i_near * psi0_1mv_i_near
	+ dfdytab(i_near, j1, k1) * dpsi1_1mu_i_near * psi0_1mv_i_near) ;

  dfdy_i_near += (dfdztab(i_near, j_near, k_near) * dpsi0_u_i_near * psi1_v_i_near
	- dfdztab(i_near, j1, k_near) * dpsi0_1mu_i_near * psi1_v_i_near
	- dfdztab(i_near, j_near, k1) * dpsi0_u_i_near * psi1_1mv_i_near
	+ dfdztab(i_near, j1, k1) * dpsi0_1mu_i_near * psi1_1mv_i_near) * dz_i_near /dy_i_near ;
  
 dfdy_i_near += (d2fdydztab(i_near, j_near, k_near) * dpsi1_u_i_near * psi1_v_i_near
	+ d2fdydztab(i_near, j1, k_near) * dpsi1_1mu_i_near * psi1_v_i_near
	- d2fdydztab(i_near, j_near, k1) * dpsi1_u_i_near * psi1_1mv_i_near
	- d2fdydztab(i_near, j1, k1) * dpsi1_1mu_i_near * psi1_1mv_i_near) * dz_i_near ;


  double dpsi0_v_i_near = 6.*(v2_i_near - v_i_near) ;
  double dpsi0_1mv_i_near = 6.*(v2_i_near - v_i_near) ;
  double dpsi1_v_i_near= 3.*v2_i_near - 4.*v_i_near + 1. ;
  double dpsi1_1mv_i_near= 3.*v2_i_near - 2.*v_i_near ;

  double dfdz_i_near;

  dfdz_i_near = (ftab(i_near, j_near, k_near) * psi0_u_i_near * dpsi0_v_i_near
    + ftab(i_near, j1, k_near) * psi0_1mu_i_near * dpsi0_v_i_near
    - ftab(i_near, j_near, k1) * psi0_u_i_near * dpsi0_1mv_i_near
    - ftab(i_near, j1, k1)  * psi0_1mu_i_near * dpsi0_1mv_i_near) / dz_i_near ;

  dfdz_i_near += (dfdytab(i_near, j_near, k_near) * psi1_u_i_near * dpsi0_v_i_near
	- dfdytab(i_near, j1, k_near) * psi1_1mu_i_near * dpsi0_v_i_near
	- dfdytab(i_near, j_near, k1) * psi1_u_i_near * dpsi0_1mv_i_near
	+ dfdytab(i_near, j1, k1) * psi1_1mu_i_near * dpsi0_1mv_i_near) * dy_i_near / dz_i_near ;

  dfdz_i_near += (dfdztab(i_near, j_near, k_near) * psi0_u_i_near * dpsi1_v_i_near
	+ dfdztab(i_near, j1, k_near) * psi0_1mu_i_near * dpsi1_v_i_near
	+ dfdztab(i_near, j_near, k1) * psi0_u_i_near * dpsi1_1mv_i_near
	+ dfdztab(i_near, j1, k1) * psi0_1mu_i_near * dpsi1_1mv_i_near) ;
  
  dfdz_i_near += (d2fdydztab(i_near, j_near, k_near) * psi1_u_i_near* dpsi1_v_i_near
	- d2fdydztab(i_near, j1, k_near) * psi1_1mu_i_near * dpsi1_v_i_near
	+ d2fdydztab(i_near, j_near, k1) * psi1_u_i_near* dpsi1_1mv_i_near
	- d2fdydztab(i_near, j1, k1) * psi1_1mu_i_near * dpsi1_1mv_i_near) * dy_i_near;


  // Hermite interpolation on the slice i1
// cette recherche est inutile car meme couverture du plan mu1, mu2 pour des delta differents
j_near = 0; 
k_near = 0; 
 while ( ( ytab(i1,j_near, 0) <= y ) && ( ( nbp2-1 ) > j_near ) ) {
      j_near++ ;
  }
  if (j_near != 0) {
      j_near -- ; 
  }

  while ( ( ztab( i1, j_near, k_near) <= z) && ( ( nbp3-1 ) > k_near ) ) {
      k_near++ ;
  }
  if (k_near != 0) {
      k_near-- ; 
  }


j1 = j_near + 1 ;
k1 = k_near + 1 ;


//cout << "i1 " << j_near << "   " << k_near << endl ;


/*
 cout << "y_dessous= " << ytab(i1, j_near, k_near) << endl ;
 cout << "y= " << y << endl ;
cout << "y_dessus= " << ytab(i1, j1, k_near) << endl ;
 cout << "z_dessus= " << ztab(i1, j_near, k1) << endl ;
 cout << "z= " << z << endl ;
 cout << "z_dessous= " << ztab(i1, j_near, k_near) << endl ;
*/
  double dy_i1 = ytab(i1, j1, k_near) - ytab(i1, j_near, k_near) ;
  double dz_i1 = ztab(i1, j_near, k1) - ztab(i1, j_near, k_near) ;

  double u_i1 = (y - ytab(i1, j_near, k_near)) / dy_i1 ;
  double v_i1 = (z - ztab(i1, j_near, k_near)) / dz_i1 ;

  double u2_i1 = u_i1 * u_i1 ; 
  double v2_i1 = v_i1 * v_i1 ;
  double u3_i1 = u2_i1 * u_i1 ; 
  double v3_i1 = v2_i1 * v_i1 ;

  double psi0_u_i1 = 2.*u3_i1 - 3.*u2_i1 + 1. ;
  double psi0_1mu_i1 = -2.*u3_i1 + 3.*u2_i1 ;
  double psi1_u_i1 = u3_i1 - 2.*u2_i1 + u_i1 ;
  double psi1_1mu_i1 = -u3_i1 + u2_i1 ;

  double psi0_v_i1 = 2.*v3_i1 - 3.*v2_i1 + 1. ;
  double psi0_1mv_i1 = -2.*v3_i1 + 3.*v2_i1 ;
  double psi1_v_i1 = v3_i1 - 2.*v2_i1 + v_i1 ;
  double psi1_1mv_i1 = -v3_i1 + v2_i1 ;

  double f_i1;
 
  f_i1 = ftab(i1, j_near, k_near) * psi0_u_i1 * psi0_v_i1
    + ftab(i1, j1, k_near) * psi0_1mu_i1 * psi0_v_i1
    + ftab(i1, j_near, k1) * psi0_u_i1 * psi0_1mv_i1
    + ftab(i1, j1, k1)  * psi0_1mu_i1 * psi0_1mv_i1 ;

  f_i1 += (dfdytab(i1, j_near, k_near) * psi1_u_i1 * psi0_v_i1
	- dfdytab(i1, j1, k_near) * psi1_1mu_i1 * psi0_v_i1
	+ dfdytab(i1, j_near, k1) * psi1_u_i1 * psi0_1mv_i1
	- dfdytab(i1, j1, k1) * psi1_1mu_i1 * psi0_1mv_i1) * dy_i1 ;

  f_i1 += (dfdztab(i1, j_near, k_near) * psi0_u_i1 * psi1_v_i1
	+ dfdztab(i1, j1, k_near) * psi0_1mu_i1 * psi1_v_i1
	- dfdztab(i1, j_near, k1) * psi0_u_i1 * psi1_1mv_i1
	- dfdztab(i1, j1, k1) * psi0_1mu_i1 * psi1_1mv_i1) * dz_i1 ;
  
  f_i1 += (d2fdydztab(i1, j_near, k_near) * psi1_u_i1 * psi1_v_i1
	- d2fdydztab(i1, j1, k_near) * psi1_1mu_i1 * psi1_v_i1
	- d2fdydztab(i1, j_near, k1) * psi1_u_i1 * psi1_1mv_i1 
	+ d2fdydztab(i1, j1, k1) * psi1_1mu_i1 * psi1_1mv_i1) * dy_i1 * dz_i1 ;

  double dpsi0_u_i1 = 6.*(u2_i1 - u_i1) ;
  double dpsi0_1mu_i1 = 6.*(u2_i1 - u_i1) ;
  double dpsi1_u_i1 = 3.*u2_i1 - 4.*u_i1 + 1. ;
  double dpsi1_1mu_i1 = 3.*u2_i1 - 2.*u_i1 ;

  double dfdy_i1;
 
  dfdy_i1 = (ftab(i1, j_near, k_near) * dpsi0_u_i1 * psi0_v_i1
	  - ftab(i1, j1, k_near) * dpsi0_1mu_i1 * psi0_v_i1
	  + ftab(i1, j_near, k1) * dpsi0_u_i1 * psi0_1mv_i1
	  - ftab(i1, j1, k1)  * dpsi0_1mu_i1 * psi0_1mv_i1 ) / dy_i1;

  dfdy_i1 += (dfdytab(i1, j_near, k_near) * dpsi1_u_i1 * psi0_v_i1
	+ dfdytab(i1, j1, k_near) * dpsi1_1mu_i1 * psi0_v_i1
	+ dfdytab(i1, j_near, k1) * dpsi1_u_i1 * psi0_1mv_i1
	+ dfdytab(i1, j1, k1) * dpsi1_1mu_i1 * psi0_1mv_i1) ;

  dfdy_i1 += (dfdztab(i1, j_near, k_near) * dpsi0_u_i1 * psi1_v_i1
	- dfdztab(i1, j1, k_near) * dpsi0_1mu_i1 * psi1_v_i1
	- dfdztab(i1, j_near, k1) * dpsi0_u_i1 * psi1_1mv_i1
	+ dfdztab(i1, j1, k1) * dpsi0_1mu_i1 * psi1_1mv_i1) * dz_i1 /dy_i1 ;
  
 dfdy_i1 += (d2fdydztab(i1, j_near, k_near) * dpsi1_u_i1 * psi1_v_i1
	+ d2fdydztab(i1, j1, k_near) * dpsi1_1mu_i1 * psi1_v_i1
	- d2fdydztab(i1, j_near, k1) * dpsi1_u_i1 * psi1_1mv_i1
	- d2fdydztab(i1, j1, k1) * dpsi1_1mu_i1 * psi1_1mv_i1) * dz_i1 ;


  double dpsi0_v_i1 = 6.*(v2_i1 - v_i1) ;
  double dpsi0_1mv_i1 = 6.*(v2_i1 - v_i1) ;
  double dpsi1_v_i1= 3.*v2_i1 - 4.*v_i1 + 1. ;
  double dpsi1_1mv_i1= 3.*v2_i1 - 2.*v_i1 ;

  double dfdz_i1;

  dfdz_i1 = (ftab(i1, j_near, k_near) * psi0_u_i1 * dpsi0_v_i1
    + ftab(i1, j1, k_near) * psi0_1mu_i1 * dpsi0_v_i1
    - ftab(i1, j_near, k1) * psi0_u_i1 * dpsi0_1mv_i1
    - ftab(i1, j1, k1)  * psi0_1mu_i1 * dpsi0_1mv_i1) / dz_i1 ;

  dfdz_i1 += (dfdytab(i1, j_near, k_near) * psi1_u_i1 * dpsi0_v_i1
	- dfdytab(i1, j1, k_near) * psi1_1mu_i1 * dpsi0_v_i1
	- dfdytab(i1, j_near, k1) * psi1_u_i1 * dpsi0_1mv_i1
	+ dfdytab(i1, j1, k1) * psi1_1mu_i1 * dpsi0_1mv_i1) * dy_i1 / dz_i1 ;

  dfdz_i1 += (dfdztab(i1, j_near, k_near) * psi0_u_i1 * dpsi1_v_i1
	+ dfdztab(i1, j1, k_near) * psi0_1mu_i1 * dpsi1_v_i1
	+ dfdztab(i1, j_near, k1) * psi0_u_i1 * dpsi1_1mv_i1
	+ dfdztab(i1, j1, k1) * psi0_1mu_i1 * dpsi1_1mv_i1) ;
  
  dfdz_i1 += (d2fdydztab(i1, j_near, k_near) * psi1_u_i1* dpsi1_v_i1
	- d2fdydztab(i1, j1, k_near) * psi1_1mu_i1 * dpsi1_v_i1
	+ d2fdydztab(i1, j_near, k1) * psi1_u_i1* dpsi1_1mv_i1
	- d2fdydztab(i1, j1, k1) * psi1_1mu_i1 * dpsi1_1mv_i1) * dy_i1;


  /* linear interpolation on the first variable */

       double x1  = xtab(i_near, 0, 0) ;
       double x2  = xtab(i1, 0, 0) ;
   
       double x12 = x1-x2 ;
   
       //for f
       double y1 = f_i_near;
       double y2 = f_i1;
       double a  = (y1-y2)/x12 ;
       double b  = (x1*y2-y1*x2)/x12 ;
	
       f  = x*a+b ; 
/*cout << " x1 " << x1 << endl ;
cout << " x2 " << x2 << endl ;
cout << " x " << x << endl ;
cout << " y1 " << y1 << endl ;
cout << " y2 " << y2 << endl ;
cout << " y " <<f << endl ;*/
       //for df/dy
       double y1_y = dfdy_i_near;
       double y2_y = dfdy_i1;
       double a_y  = (y1_y-y2_y)/x12 ;
       double b_y  = (x1*y2_y-y1_y*x2)/x12 ;
	
       dfdy  = x*a_y+b_y ; 

       //for df/dz
       double y1_z = dfdz_i_near;
       double y2_z = dfdz_i1;
       double a_z  = (y1_z-y2_z)/x12 ;
       double b_z  = (x1*y2_z-y1_z*x2)/x12 ;
	
       dfdz  = x*a_z+b_z ; 
/*cout << "i_near  " << i_near << endl;
cout << "j_near  " << j_near << endl;
cout << "k_near  " << k_near << endl;
abort () ;*/
//cout << dfdy_i_near  << "  " <<dfdy_i1 << "  " << dfdz_i_near << "  " << dfdz_i1<< endl;
       return;

}



void interpol_mixed_3d_mod(const Tbl& xtab, const Tbl& ytab, const Tbl& ztab, const Tbl& ftab, 
		      const Tbl& dfdytab, const Tbl& dfdztab, 
		      double x, double y, double z, double& f, double& dfdy, double& dfdz) 


{
  assert(ytab.dim == xtab.dim) ; 
  assert(ztab.dim == xtab.dim) ;
  assert(ftab.dim == xtab.dim) ;
  assert(dfdytab.dim == xtab.dim) ;
  assert(dfdztab.dim == xtab.dim) ;
   
  int nbp1, nbp2, nbp3;
  nbp1 = xtab.get_dim(2) ; // \Delta^{2}
  nbp2 = xtab.get_dim(1) ; // \mu_n
  nbp3 = xtab.get_dim(0) ; // \mu_p


  /* we can put these tests directly in the code (before calling interpol_mixed_3d) 
  assert(x >= xtab(0,0,0)) ;
  assert(x <= xtab(nbp1,0,0)) ;
  assert(y >= ytab(0,0,0)) ;
  assert(y <= ytab(0,nbp2,0)) ;
  assert(z >= ztab(0,0,0)) ;
  assert(z <= ztab(0,0,nbp3)) ;
  */

  int i_near = 0 ; 
  int j_near = 0 ;
  int k_near = 0 ;


  // look for the positions of (x,y,z) in the tables
  while ( ( xtab(i_near,0,0) <= x ) && ( ( nbp1-1 ) > i_near ) ) {
      i_near++ ;
  }
  if (i_near != 0) { 
      i_near -- ; 
  }

  while ( ( ytab(i_near,j_near, 0) <= y ) && ( ( nbp2-1 ) > j_near ) ) {
      j_near++ ;
  }
  if (j_near != 0) {
      j_near -- ; 
  }

  while ( ( ztab( i_near, j_near, k_near) <= z) && ( ( nbp3-1 ) > k_near ) ) {
      k_near++ ;
  }
  if (k_near != 0) {
      k_near-- ; 
  }

  int i1 = i_near + 1 ;
  int j1 = j_near + 1 ;
  int k1 = k_near + 1 ;
//cout << j_near<< "   " << k_near<< endl;
  // Hermite interpolation on the slice i_near 
/*cout << "y_dessous= " << ytab(i_near, j_near, k_near) << endl ;
 cout << "y= " << y << endl ;
cout << "y_dessus= " << ytab(i_near, j1, k_near) << endl ;
 cout << "z_dessus= " << ztab(i_near, j_near, k1) << endl ;
 cout << "z= " << z << endl ;
 cout << "z_dessous= " << ztab(i_near, j_near, k_near) << endl ;
*/
//cout << j_near << "   " << k_near << endl ;
 

/*
bool hum1 =   (( ytab(i_near, j1, k_near) - ytab(i_near, j_near, k_near) ) == (ytab(i_near, j1, k1) - ytab(i_near, j_near, k1) ) ) ;
bool hum2 = ((ztab(i_near, j_near, k1) - ztab(i_near, j_near, k_near) ) == ( ztab(i_near, j1, k1) - ztab(i_near, j1, k_near) ))  ;
if (hum1 == false ) { 
cout << setprecision(16);
cout <<  ytab(i_near, j1, k_near) - ytab(i_near, j_near, k_near) <<  "   " << ytab(i_near, j1, k1) - ytab(i_near, j_near, k1)<< endl ;
cout << j_near << "    " << k_near << endl ;
}

if (hum2 == false ) { 
cout << setprecision(16);
cout <<  ztab(i_near, j_near, k1) - ztab(i_near, j_near, k_near) <<  "   " << ztab(i_near, j1, k1) - ztab(i_near, j1, k_near)<< endl ;
cout << j_near << "    " << k_near << endl ;
}*/

  double dy_i_near = ytab(i_near, j1, k_near) - ytab(i_near, j_near, k_near) ;
  double dz_i_near = ztab(i_near, j_near, k1) - ztab(i_near, j_near, k_near) ;
  double u_i_near = (y - ytab(i_near, j_near, k_near)) / dy_i_near ;
  double v_i_near = (z - ztab(i_near, j_near, k_near)) / dz_i_near ;
/*
//lineaire
  double dy_anal = 1. ; // log(2.) ;
  double dz_anal = 1. ; //log(2.) ;
//en log
cout << y << "   " << z << endl ;
  double u_i_anal = y /log(2.) ; //(y - 1.) / 1.
  double v_i_anal = z /log(2.) ; //(z - 1.) / 1. 

cout << setprecision(16);
cout <<"interpol  " <<  u_i_near <<"    " << u_i_anal<< "   "<<  v_i_near <<  "   "  << v_i_anal << endl ;
*/


  double u2_i_near = u_i_near * u_i_near ; 
  double v2_i_near = v_i_near * v_i_near ;
  double u3_i_near = u2_i_near * u_i_near ; 
  double v3_i_near = v2_i_near * v_i_near ;

 /* 
//lineaire
double u2_i_anal =  (y - 1. ) * (y - 1. ) ; //
  double v2_i_anal =  (z - 1.) * (z - 1.) ; // 
  double u3_i_anal = (y - 1. ) * (y - 1. )* (y - 1. )  ; //
  double v3_i_anal = (z - 1.) * (z - 1.)*  (z - 1.)  ; // 

cout << " interpol2" << u2_i_near << "    "  << u2_i_anal << "    " << u3_i_near  << "       "  << u3_i_anal << endl ;
*/

  double psi0_u_i_near = double(2)*u3_i_near - double(3)*u2_i_near + double(1) ;
  double psi0_1mu_i_near = -double(2)*u3_i_near + double(3)*u2_i_near ;
  double psi1_u_i_near = u3_i_near - double(2)*u2_i_near + u_i_near ;
  double psi1_1mu_i_near = -u3_i_near + u2_i_near ;
  double psi0_v_i_near = double(2)*v3_i_near - double(3)*v2_i_near + double(1) ;
  double psi0_1mv_i_near = -double(2)*v3_i_near + double(3)*v2_i_near ;
  double psi1_v_i_near = v3_i_near - double(2)*v2_i_near + v_i_near ;
  double psi1_1mv_i_near = -v3_i_near + v2_i_near ;

  double f_i_near ;
/*// lineaire
 double psi0_u_i_anal = double(2)*(y - 1. ) * (y - 1. )* (y - 1. ) - double(3)*(y - 1. ) *(y - 1. )  + double(1) ;
  double psi0_1mu_i_anal = -double(2)* (y - 1. ) * (y - 1. )* (y - 1. ) + double(3)*(y - 1. ) * (y - 1. );
  double psi1_u_i_anal = (y - 1. ) * (y - 1. )* (y - 1. ) - double(2)*(y - 1. ) * (y - 1.) + (y - 1. ) ;
  double psi1_1mu_i_anal = -(y - 1. ) * (y - 1. )* (y - 1. ) +(y - 1. ) * (y - 1. ) ;
//en log
 double psi0_u_i_anal = double(2)*y /log(2.) * y /log(2.) * y /log(2.) - double(3)*  y /log(2.) * y /log(2.)  + double(1) ;
  double psi0_1mu_i_anal = -double(2)*y /log(2.)*y /log(2.)*y /log(2.) + double(3)*y /log(2.)*y /log(2.);
  double psi1_u_i_anal = y /log(2.)*y /log(2.)*y /log(2.) - double(2)* y /log(2.) * y /log(2.) + y /log(2.);
  double psi1_1mu_i_anal = -y /log(2.)*y /log(2.)*y /log(2.) +y /log(2.)*y /log(2.) ;
// lineaire
 double psi0_v_i_anal = double(2)*(z - 1. ) * (z - 1. )* (z - 1. ) - double(3)*(z - 1. ) *(z - 1. )  + double(1) ;
  double psi0_1mv_i_anal = -double(2)* (z - 1. ) * (z - 1. )* (z - 1. ) + double(3)*(z - 1. ) * (z - 1. );
  double psi1_v_i_anal = (z - 1. ) * (z - 1. )* (z - 1. ) - double(2)*(z - 1. ) * (z - 1.) + (z - 1. ) ;
  double psi1_1mv_i_anal = -(z - 1. ) * (z - 1. )* (z - 1. ) +(z - 1. ) * (z - 1. ) ;

cout << " Interpol3   " << psi0_u_i_near << "    " <<  psi0_u_i_anal << "    " << psi0_1mu_i_near << "  " << psi0_1mu_i_anal << "   " << 
psi1_u_i_near << "  " << psi1_u_i_anal  << "  " <<  psi1_1mu_i_near << "  " <<  psi1_1mu_i_anal << endl;

cout << " Interpol4   " << psi0_v_i_near << "    " <<  psi0_v_i_anal << "    " << psi0_1mv_i_near << "  " << psi0_1mv_i_anal << "   " << 
psi1_v_i_near << "  " << psi1_v_i_anal  << "  " <<  psi1_1mv_i_near << "  " <<  psi1_1mv_i_anal << endl;
*/

  f_i_near = ftab(i_near, j_near, k_near) * psi0_u_i_near * psi0_v_i_near
    + ftab(i_near, j1, k_near) * psi0_1mu_i_near * psi0_v_i_near
    + ftab(i_near, j_near, k1) * psi0_u_i_near * psi0_1mv_i_near
    + ftab(i_near, j1, k1)  * psi0_1mu_i_near * psi0_1mv_i_near ;


  f_i_near += (dfdytab(i_near, j_near, k_near) * psi1_u_i_near * psi0_v_i_near
	- dfdytab(i_near, j1, k_near) * psi1_1mu_i_near * psi0_v_i_near
	+ dfdytab(i_near, j_near, k1) * psi1_u_i_near * psi0_1mv_i_near
	- dfdytab(i_near, j1, k1) * psi1_1mu_i_near * psi0_1mv_i_near) * dy_i_near ;

  f_i_near += (dfdztab(i_near, j_near, k_near) * psi0_u_i_near * psi1_v_i_near
	+ dfdztab(i_near, j1, k_near) * psi0_1mu_i_near * psi1_v_i_near
	- dfdztab(i_near, j_near, k1) * psi0_u_i_near * psi1_1mv_i_near
	- dfdztab(i_near, j1, k1) * psi0_1mu_i_near * psi1_1mv_i_near) * dz_i_near ;

  


/*

cout<< setprecision(16);

cout << " Y " << y << endl;
cout << " Z " << z << endl;
cout << "dy " << dy_i_near << "  " <<  ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))<< endl ;
cout << "dz " << dz_i_near << "  " <<  log(ztab(i_near, j_near, k1))-log(ztab(i_near, j_near, k_near))<< endl ;

cout << " u  " << u_i_near << "  " <<  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))<< endl ;
cout << " v  " << v_i_near << "  " <<  (log(z)-log(ztab(i_near, j_near, k_near)))/(log(ztab(i_near, j_near, k1))-log(ztab(i_near, j_near, k_near)))<< endl ;

cout << " u2  " << u2_i_near << "  " <<  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))  * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) << endl ;
cout << " psi0_u_i_near " << psi0_u_i_near << "  " <<  double(2)* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) - (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) *  double(3) +double(1) << endl ;

cout << "psi0_1mu_i_near" << psi0_1mu_i_near << "  " << -double(2) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) + double(3)* (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))*  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))<< endl ;

cout <<"psi1_u_i_near " << psi1_u_i_near  << "   " <<  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))*  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) - (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) *  double(2) + (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) << endl ;


cout << "psi1_1mu_i_near " <<  psi1_1mu_i_near << "   " << - (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))*  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) * (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) + (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near)))*  (log(y) - log(ytab(i_near, j_near, k_near)))/(log(ytab(i_near, j1, k_near))-log(ytab(i_near, j_near, k_near))) << endl ;


cout << " ftab(i_near, j_near, k_near) " << ftab(i_near, j_near, k_near) << "   " << log(ftab(i_near, j_near, k_near)) << endl;
cout << " dfdytab(i_near, j_near, k_near) " << dfdytab(i_near, j_near, k_near) <<  "   " << dfdytab(i_near, j_near, k_near)  * ytab(i_near, j_near, k_near) / ftab(i_near, j_near, k_near) << endl ;
cout << "(d2fdydztab(i_near, j_near, k_near)  " << d2fdydztab(i_near, j_near, k_near) << "    " << 
	ytab(i_near, j_near, k_near)  * ztab(i_near, j_near, k_near) / ftab(i_near, j_near, k_near) * ( d2fdydztab(i_near, j_near, k_near) - dfdytab(i_near, j_near, k_near) * dfdztab(i_near, j_near, k_near) / 
 ftab(i_near, j_near, k_near) ) << endl ;
cout << "  f = " << f_i_near << "   " << log(f_i_near) << endl;




*/

  double dpsi0_u_i_near = 6.*(u2_i_near - u_i_near) ;
  double dpsi0_1mu_i_near = 6.*(u2_i_near - u_i_near) ;
  double dpsi1_u_i_near = 3.*u2_i_near - 4.*u_i_near + 1. ;
  double dpsi1_1mu_i_near = 3.*u2_i_near - 2.*u_i_near ;

  double dfdy_i_near;
 
  dfdy_i_near = (ftab(i_near, j_near, k_near) * dpsi0_u_i_near * psi0_v_i_near
	  - ftab(i_near, j1, k_near) * dpsi0_1mu_i_near * psi0_v_i_near
	  + ftab(i_near, j_near, k1) * dpsi0_u_i_near * psi0_1mv_i_near
	  - ftab(i_near, j1, k1)  * dpsi0_1mu_i_near * psi0_1mv_i_near ) / dy_i_near;

  dfdy_i_near += (dfdytab(i_near, j_near, k_near) * dpsi1_u_i_near * psi0_v_i_near
	+ dfdytab(i_near, j1, k_near) * dpsi1_1mu_i_near * psi0_v_i_near
	+ dfdytab(i_near, j_near, k1) * dpsi1_u_i_near * psi0_1mv_i_near
	+ dfdytab(i_near, j1, k1) * dpsi1_1mu_i_near * psi0_1mv_i_near) ;

  dfdy_i_near += (dfdztab(i_near, j_near, k_near) * dpsi0_u_i_near * psi1_v_i_near
	- dfdztab(i_near, j1, k_near) * dpsi0_1mu_i_near * psi1_v_i_near
	- dfdztab(i_near, j_near, k1) * dpsi0_u_i_near * psi1_1mv_i_near
	+ dfdztab(i_near, j1, k1) * dpsi0_1mu_i_near * psi1_1mv_i_near) * dz_i_near /dy_i_near ;
  


  double dpsi0_v_i_near = 6.*(v2_i_near - v_i_near) ;
  double dpsi0_1mv_i_near = 6.*(v2_i_near - v_i_near) ;
  double dpsi1_v_i_near= 3.*v2_i_near - 4.*v_i_near + 1. ;
  double dpsi1_1mv_i_near= 3.*v2_i_near - 2.*v_i_near ;

  double dfdz_i_near;

  dfdz_i_near = (ftab(i_near, j_near, k_near) * psi0_u_i_near * dpsi0_v_i_near
    + ftab(i_near, j1, k_near) * psi0_1mu_i_near * dpsi0_v_i_near
    - ftab(i_near, j_near, k1) * psi0_u_i_near * dpsi0_1mv_i_near
    - ftab(i_near, j1, k1)  * psi0_1mu_i_near * dpsi0_1mv_i_near) / dz_i_near ;

  dfdz_i_near += (dfdytab(i_near, j_near, k_near) * psi1_u_i_near * dpsi0_v_i_near
	- dfdytab(i_near, j1, k_near) * psi1_1mu_i_near * dpsi0_v_i_near
	- dfdytab(i_near, j_near, k1) * psi1_u_i_near * dpsi0_1mv_i_near
	+ dfdytab(i_near, j1, k1) * psi1_1mu_i_near * dpsi0_1mv_i_near) * dy_i_near / dz_i_near ;

  dfdz_i_near += (dfdztab(i_near, j_near, k_near) * psi0_u_i_near * dpsi1_v_i_near
	+ dfdztab(i_near, j1, k_near) * psi0_1mu_i_near * dpsi1_v_i_near
	+ dfdztab(i_near, j_near, k1) * psi0_u_i_near * dpsi1_1mv_i_near
	+ dfdztab(i_near, j1, k1) * psi0_1mu_i_near * dpsi1_1mv_i_near) ;
  
  




  // Hermite interpolation on the slice i1
// cette recherche est inutile car meme couverture du plan mu1, mu2 pour des delta differents
j_near = 0; 
k_near = 0; 
 while ( ( ytab(i1,j_near, 0) <= y ) && ( ( nbp2-1 ) > j_near ) ) {
      j_near++ ;
  }
  if (j_near != 0) {
      j_near -- ; 
  }

  while ( ( ztab( i1, j_near, k_near) <= z) && ( ( nbp3-1 ) > k_near ) ) {
      k_near++ ;
  }
  if (k_near != 0) {
      k_near-- ; 
  }


j1 = j_near + 1 ;
k1 = k_near + 1 ;


//cout << "i1 " << j_near << "   " << k_near << endl ;


/*
 cout << "y_dessous= " << ytab(i1, j_near, k_near) << endl ;
 cout << "y= " << y << endl ;
cout << "y_dessus= " << ytab(i1, j1, k_near) << endl ;
 cout << "z_dessus= " << ztab(i1, j_near, k1) << endl ;
 cout << "z= " << z << endl ;
 cout << "z_dessous= " << ztab(i1, j_near, k_near) << endl ;
*/
  double dy_i1 = ytab(i1, j1, k_near) - ytab(i1, j_near, k_near) ;
  double dz_i1 = ztab(i1, j_near, k1) - ztab(i1, j_near, k_near) ;

  double u_i1 = (y - ytab(i1, j_near, k_near)) / dy_i1 ;
  double v_i1 = (z - ztab(i1, j_near, k_near)) / dz_i1 ;

  double u2_i1 = u_i1 * u_i1 ; 
  double v2_i1 = v_i1 * v_i1 ;
  double u3_i1 = u2_i1 * u_i1 ; 
  double v3_i1 = v2_i1 * v_i1 ;

  double psi0_u_i1 = 2.*u3_i1 - 3.*u2_i1 + 1. ;
  double psi0_1mu_i1 = -2.*u3_i1 + 3.*u2_i1 ;
  double psi1_u_i1 = u3_i1 - 2.*u2_i1 + u_i1 ;
  double psi1_1mu_i1 = -u3_i1 + u2_i1 ;

  double psi0_v_i1 = 2.*v3_i1 - 3.*v2_i1 + 1. ;
  double psi0_1mv_i1 = -2.*v3_i1 + 3.*v2_i1 ;
  double psi1_v_i1 = v3_i1 - 2.*v2_i1 + v_i1 ;
  double psi1_1mv_i1 = -v3_i1 + v2_i1 ;

  double f_i1;
 
  f_i1 = ftab(i1, j_near, k_near) * psi0_u_i1 * psi0_v_i1
    + ftab(i1, j1, k_near) * psi0_1mu_i1 * psi0_v_i1
    + ftab(i1, j_near, k1) * psi0_u_i1 * psi0_1mv_i1
    + ftab(i1, j1, k1)  * psi0_1mu_i1 * psi0_1mv_i1 ;

  f_i1 += (dfdytab(i1, j_near, k_near) * psi1_u_i1 * psi0_v_i1
	- dfdytab(i1, j1, k_near) * psi1_1mu_i1 * psi0_v_i1
	+ dfdytab(i1, j_near, k1) * psi1_u_i1 * psi0_1mv_i1
	- dfdytab(i1, j1, k1) * psi1_1mu_i1 * psi0_1mv_i1) * dy_i1 ;

  f_i1 += (dfdztab(i1, j_near, k_near) * psi0_u_i1 * psi1_v_i1
	+ dfdztab(i1, j1, k_near) * psi0_1mu_i1 * psi1_v_i1
	- dfdztab(i1, j_near, k1) * psi0_u_i1 * psi1_1mv_i1
	- dfdztab(i1, j1, k1) * psi0_1mu_i1 * psi1_1mv_i1) * dz_i1 ;
  
  
  double dpsi0_u_i1 = 6.*(u2_i1 - u_i1) ;
  double dpsi0_1mu_i1 = 6.*(u2_i1 - u_i1) ;
  double dpsi1_u_i1 = 3.*u2_i1 - 4.*u_i1 + 1. ;
  double dpsi1_1mu_i1 = 3.*u2_i1 - 2.*u_i1 ;

  double dfdy_i1;
 
  dfdy_i1 = (ftab(i1, j_near, k_near) * dpsi0_u_i1 * psi0_v_i1
	  - ftab(i1, j1, k_near) * dpsi0_1mu_i1 * psi0_v_i1
	  + ftab(i1, j_near, k1) * dpsi0_u_i1 * psi0_1mv_i1
	  - ftab(i1, j1, k1)  * dpsi0_1mu_i1 * psi0_1mv_i1 ) / dy_i1;

  dfdy_i1 += (dfdytab(i1, j_near, k_near) * dpsi1_u_i1 * psi0_v_i1
	+ dfdytab(i1, j1, k_near) * dpsi1_1mu_i1 * psi0_v_i1
	+ dfdytab(i1, j_near, k1) * dpsi1_u_i1 * psi0_1mv_i1
	+ dfdytab(i1, j1, k1) * dpsi1_1mu_i1 * psi0_1mv_i1) ;

  dfdy_i1 += (dfdztab(i1, j_near, k_near) * dpsi0_u_i1 * psi1_v_i1
	- dfdztab(i1, j1, k_near) * dpsi0_1mu_i1 * psi1_v_i1
	- dfdztab(i1, j_near, k1) * dpsi0_u_i1 * psi1_1mv_i1
	+ dfdztab(i1, j1, k1) * dpsi0_1mu_i1 * psi1_1mv_i1) * dz_i1 /dy_i1 ;
  
 


  double dpsi0_v_i1 = 6.*(v2_i1 - v_i1) ;
  double dpsi0_1mv_i1 = 6.*(v2_i1 - v_i1) ;
  double dpsi1_v_i1= 3.*v2_i1 - 4.*v_i1 + 1. ;
  double dpsi1_1mv_i1= 3.*v2_i1 - 2.*v_i1 ;

  double dfdz_i1;

  dfdz_i1 = (ftab(i1, j_near, k_near) * psi0_u_i1 * dpsi0_v_i1
    + ftab(i1, j1, k_near) * psi0_1mu_i1 * dpsi0_v_i1
    - ftab(i1, j_near, k1) * psi0_u_i1 * dpsi0_1mv_i1
    - ftab(i1, j1, k1)  * psi0_1mu_i1 * dpsi0_1mv_i1) / dz_i1 ;

  dfdz_i1 += (dfdytab(i1, j_near, k_near) * psi1_u_i1 * dpsi0_v_i1
	- dfdytab(i1, j1, k_near) * psi1_1mu_i1 * dpsi0_v_i1
	- dfdytab(i1, j_near, k1) * psi1_u_i1 * dpsi0_1mv_i1
	+ dfdytab(i1, j1, k1) * psi1_1mu_i1 * dpsi0_1mv_i1) * dy_i1 / dz_i1 ;

  dfdz_i1 += (dfdztab(i1, j_near, k_near) * psi0_u_i1 * dpsi1_v_i1
	+ dfdztab(i1, j1, k_near) * psi0_1mu_i1 * dpsi1_v_i1
	+ dfdztab(i1, j_near, k1) * psi0_u_i1 * dpsi1_1mv_i1
	+ dfdztab(i1, j1, k1) * psi0_1mu_i1 * dpsi1_1mv_i1) ;
  



  /* linear interpolation on the first variable */

       double x1  = xtab(i_near, 0, 0) ;
       double x2  = xtab(i1, 0, 0) ;
   
       double x12 = x1-x2 ;
   
       //for f
       double y1 = f_i_near;
       double y2 = f_i1;
       double a  = (y1-y2)/x12 ;
       double b  = (x1*y2-y1*x2)/x12 ;
	
       f  = x*a+b ; 
// cout << " x1 " << x1 << " x2 " << x2 << " x " << x <<  " y1 " << y1 << " y2 " << y2 << " y " <<f << endl ;
       //for df/dy
       double y1_y = dfdy_i_near;
       double y2_y = dfdy_i1;
       double a_y  = (y1_y-y2_y)/x12 ;
       double b_y  = (x1*y2_y-y1_y*x2)/x12 ;
	
       dfdy  = x*a_y+b_y ; 

       //for df/dz
       double y1_z = dfdz_i_near;
       double y2_z = dfdz_i1;
       double a_z  = (y1_z-y2_z)/x12 ;
       double b_z  = (x1*y2_z-y1_z*x2)/x12 ;
	
       dfdz  = x*a_z+b_z ; 
/*cout << "i_near  " << i_near << endl;
cout << "j_near  " << j_near << endl;
cout << "k_near  " << k_near << endl;
abort () ;*/
       return;

}



































void interpol_herm_2d_new_avec( double y, double z, 
				double mu1_11, double mu1_21, double mu2_11, double mu2_12, 
				double p_11, double p_21,double  p_12, double p_22,
				double n1_11, double n1_21, double n1_12,double  n1_22, 
				double n2_11, double n2_21,double  n2_12, double n2_22, 
				double cross_11, double cross_21, double cross_12, double cross_22,
				double& f, double& dfdy, double& dfdz) {


  double dy = mu1_21 - mu1_11 ;
  double dz = mu2_12 - mu2_11;

  double u = (y - mu1_11) / dy ;
  double v = (z - mu2_11) / dz ;

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

  f = p_11 * psi0_u * psi0_v
    + p_21 * psi0_1mu * psi0_v 
    + p_12 * psi0_u * psi0_1mv
    + p_22  * psi0_1mu * psi0_1mv ;


  f += (n1_11 * psi1_u * psi0_v
	- n1_21 * psi1_1mu * psi0_v
	+ n1_12* psi1_u * psi0_1mv
	- n1_22 * psi1_1mu * psi0_1mv) * dy ;

  f += (n2_11 * psi0_u * psi1_v
	+ n2_21 * psi0_1mu * psi1_v
	-  n2_12 * psi0_u * psi1_1mv
	- n2_22* psi0_1mu * psi1_1mv) * dz ;
  
  f += (cross_11 * psi1_u * psi1_v
	- cross_21 * psi1_1mu * psi1_v
	- cross_12 * psi1_u * psi1_1mv 
	+ cross_22 * psi1_1mu * psi1_1mv) * dy * dz ;

  double dpsi0_u = 6.*(u2 - u) ;
  double dpsi0_1mu = 6.*(u2 - u) ;
  double dpsi1_u = 3.*u2 - 4.*u + 1. ;
  double dpsi1_1mu = 3.*u2 - 2.*u ;

  dfdy = (p_11 * dpsi0_u * psi0_v
    - p_21 * dpsi0_1mu * psi0_v 
    + p_12 * dpsi0_u * psi0_1mv
	  - p_22 * dpsi0_1mu * psi0_1mv ) / dy;

  dfdy += (n1_11* dpsi1_u * psi0_v
	+ n1_21 * dpsi1_1mu * psi0_v
	+ n1_12 * dpsi1_u * psi0_1mv
	+ n1_22 * dpsi1_1mu * psi0_1mv) ;

  dfdy += (n2_11 * dpsi0_u * psi1_v
	- n2_21 * dpsi0_1mu * psi1_v
	-  n2_12 * dpsi0_u * psi1_1mv
	+ n2_22 * dpsi0_1mu * psi1_1mv) * dz /dy ;
  
  dfdy += (cross_11 * dpsi1_u * psi1_v
	+ cross_21* dpsi1_1mu * psi1_v
	- cross_12 * dpsi1_u * psi1_1mv 
	- cross_22 * dpsi1_1mu * psi1_1mv) * dz ;


  double dpsi0_v = 6.*(v2 - v) ;
  double dpsi0_1mv = 6.*(v2 - v) ;
  double dpsi1_v = 3.*v2 - 4.*v + 1. ;
  double dpsi1_1mv = 3.*v2 - 2.*v ;


  dfdz = (p_11* psi0_u * dpsi0_v
    + p_21 * psi0_1mu * dpsi0_v 
    - p_12 * psi0_u * dpsi0_1mv
	  - p_22  * psi0_1mu * dpsi0_1mv) / dz ;

  dfdz += (n1_11 * psi1_u * dpsi0_v
	- n1_21 * psi1_1mu * dpsi0_v
	- n1_12 * psi1_u * dpsi0_1mv
	+ n1_22 * psi1_1mu * dpsi0_1mv) * dy / dz ;

  dfdz += (n2_11 * psi0_u * dpsi1_v
	+ n2_21 * psi0_1mu * dpsi1_v
	+  n2_12 * psi0_u * dpsi1_1mv
	+ n2_22 * psi0_1mu * dpsi1_1mv) ;
  
  dfdz += (cross_11 * psi1_u * dpsi1_v
	- cross_21 * psi1_1mu * dpsi1_v
	+ cross_12 * psi1_u * dpsi1_1mv 
	- cross_22 * psi1_1mu * dpsi1_1mv) * dy ;

  return ;

}







void interpol_herm_2d_new_sans( double y, double z, 
				double mu1_11, double mu1_21, double mu2_11, double mu2_12, 
				double p_11, double p_21,double  p_12, double p_22,
				double n1_11, double n1_21, double n1_12,double  n1_22, 
				double n2_11, double n2_21,double  n2_12, double n2_22, 
				double& f, double& dfdy, double& dfdz) {


  double dy = mu1_21 - mu1_11 ;
  double dz = mu2_12 - mu2_11;

  double u = (y - mu1_11) / dy ;
  double v = (z - mu2_11) / dz ;

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

  f = p_11 * psi0_u * psi0_v
    + p_21 * psi0_1mu * psi0_v 
    + p_12 * psi0_u * psi0_1mv
    + p_22  * psi0_1mu * psi0_1mv ;


  f += (n1_11 * psi1_u * psi0_v
	- n1_21 * psi1_1mu * psi0_v
	+ n1_12* psi1_u * psi0_1mv
	- n1_22 * psi1_1mu * psi0_1mv) * dy ;

  f += (n2_11 * psi0_u * psi1_v
	+ n2_21 * psi0_1mu * psi1_v
	-  n2_12 * psi0_u * psi1_1mv
	- n2_22* psi0_1mu * psi1_1mv) * dz ;

  double dpsi0_u = 6.*(u2 - u) ;
  double dpsi0_1mu = 6.*(u2 - u) ;
  double dpsi1_u = 3.*u2 - 4.*u + 1. ;
  double dpsi1_1mu = 3.*u2 - 2.*u ;

  dfdy = (p_11 * dpsi0_u * psi0_v
    - p_21 * dpsi0_1mu * psi0_v 
    + p_12 * dpsi0_u * psi0_1mv
	  - p_22 * dpsi0_1mu * psi0_1mv ) / dy;

  dfdy += (n1_11* dpsi1_u * psi0_v
	+ n1_21 * dpsi1_1mu * psi0_v
	+ n1_12 * dpsi1_u * psi0_1mv
	+ n1_22 * dpsi1_1mu * psi0_1mv) ;

  dfdy += (n2_11 * dpsi0_u * psi1_v
	- n2_21 * dpsi0_1mu * psi1_v
	-  n2_12 * dpsi0_u * psi1_1mv
	+ n2_22 * dpsi0_1mu * psi1_1mv) * dz /dy ;



  double dpsi0_v = 6.*(v2 - v) ;
  double dpsi0_1mv = 6.*(v2 - v) ;
  double dpsi1_v = 3.*v2 - 4.*v + 1. ;
  double dpsi1_1mv = 3.*v2 - 2.*v ;


  dfdz = (p_11* psi0_u * dpsi0_v
    + p_21 * psi0_1mu * dpsi0_v 
    - p_12 * psi0_u * dpsi0_1mv
	  - p_22  * psi0_1mu * dpsi0_1mv) / dz ;

  dfdz += (n1_11 * psi1_u * dpsi0_v
	- n1_21 * psi1_1mu * dpsi0_v
	- n1_12 * psi1_u * dpsi0_1mv
	+ n1_22 * psi1_1mu * dpsi0_1mv) * dy / dz ;

  dfdz += (n2_11 * psi0_u * dpsi1_v
	+ n2_21 * psi0_1mu * dpsi1_v
	+  n2_12 * psi0_u * dpsi1_1mv
	+ n2_22 * psi0_1mu * dpsi1_1mv) ;

   //cout << "-alpha , devrait etre negatif !!" << f << endl;
  return ;

}























/*interpolation for functions of 3 variables : hermite interpolation on the 2 last variables and linear interpolation on the first one*/

/*
xtab/x refers to delta_car (logdelta_car)
ytab/y refers to mu_n (logent1)
ztab/z refers to mu_p (logent2)
ftab/f refers to psi (logp) or alpha (dlpsddelta_car or logalpha)
dfdytab/dfdy refers to dpsi/dmu_n (dlpsdlent1) or dalpha/dmu_n (d2lpsdlent1ddelta_car or dlalphadlent1)
dfdztab/dfdz refers to dpsi/dmu_p (dlpsdlent2) or dalpha/dmu_p (dl2psdlent2ddelta_car or dlalphadlent2)
d2fdytabdztab refers to d2psi/dmu_ndmu_p (d2lpsdlent1dlent2) or d2alpha/dmu_ndmu_p (d3lpsdlent1dlent2ddelta_car or d2lalphadlent1dlent2)
*/

/*this routine provides the interpolated values of f, dfdy and dfdz at point (x,y,z) via the use of adapted tables*/

void interpol_mixed_3d_new(double m_1, double m_2, const Tbl& xtab, const Tbl& ytab, const Tbl& ztab, 
			const Tbl& ftab, const Tbl& dfdytab, const Tbl& dfdztab, const Tbl& d2fdydztab,
			const Tbl& dlpsddelta_car, const Tbl&  d2lpsdlent1ddelta_car, const Tbl& d2lpsdlent2ddelta_car, 
			const Tbl&  mu2_P, const Tbl&  n_p_P,   const Tbl& press_P,
			const Tbl& mu1_N, const Tbl& n_n_N, const Tbl& press_N,
			const Tbl& delta_car_n0, const Tbl& mu1_n0, const Tbl& mu2_n0,
			const Tbl& delta_car_p0, const Tbl& mu1_p0, const Tbl& mu2_p0, 
		        double x, double y, double z, double& f, double& dfdy, double& dfdz, double& alpha)

{

    assert(ytab.dim == xtab.dim) ; 
    assert(ztab.dim == xtab.dim) ;
    assert(ftab.dim == xtab.dim) ;
    assert(dfdytab.dim == xtab.dim) ;
    assert(dfdztab.dim == xtab.dim) ;
    assert(d2fdydztab.dim == xtab.dim) ;
  
    int nbp1, nbp2, nbp3;
    nbp1 = xtab.get_dim(2) ; // \Delta^{2}
    nbp2 = xtab.get_dim(1) ; // \mu_n
    nbp3 = xtab.get_dim(0) ; // \mu_p

    int i_near = 0 ; 
    int j_near = 0 ;
    int k_near = 0 ;

    // look for the positions of (x,y,z) in the tables
    while ( ( xtab(i_near,0,0) <= x ) && ( ( nbp1-1 ) > i_near ) ) {
	i_near++ ;
    }
    if (i_near != 0) { 
	i_near -- ; 
    }

    while ( ( ytab(i_near,j_near, 0) <= y ) && ( ( nbp2-1 ) > j_near ) ) {
	j_near++ ;
    }
    if (j_near != 0) {
	j_near -- ; 
    }

    while ( ( ztab( i_near, j_near, k_near) <= z) && ( ( nbp3-1 ) > k_near ) ) {
	k_near++ ;
    }
    if (k_near != 0) {
	k_near-- ; 
    }
  
    int i1 = i_near + 1 ;
    int j1 = j_near + 1 ;
    int k1 = k_near + 1 ;

    
    
    // Remarque : cette recherche implique qu'on est au moins supérieure a la premiere valeur (il faudrait mettre un abort sinon!)
    // Pour vérifier que la recherche se déroule comme prévue : 
  //  cout << " DEBUG mode : " << endl ;
    if ( ( xtab( i_near, j_near, k_near) > x ) || (x > xtab( i1, j_near, k_near) ) ) {
      cout << "mauvais positionnement de x dans xtab " << endl ;
      cout << xtab( i_near, j_near, k_near) << "  " << x <<  "  " << xtab( i1, j_near, k_near) << endl;
      abort();
    }
 
    if ( ( ytab( i_near, j_near, k_near) > y ) || (y > ytab( i1, j1, k_near) ) ) {
      cout << "mauvais positionnement de y dans ytab " << endl ;
      cout << ytab( i_near, j_near, k_near) * 1.009000285 * 1.66e-27 * 3e8 * 3e8 /(1.6e-19 * 1e6) << "  " << y <<   "  " << ytab( i1, j1, k_near) * 1.009000285 * 1.66e-27 * 3e8 * 3e8 /(1.6e-19 * 1e6)  << endl;
      abort();
    }
  
    if ( ( ztab( i_near, j_near, k_near) > z ) || ( z > ztab( i1, j_near, k1) ) ){
      cout << "mauvais positionnement de z dans ztab " << endl ;
      cout << ztab( i_near, j_near, k_near) << "  " << z <<  "  " << ztab( i1, j_near, k1) << endl;
      abort();
    }
    
    
  double f_i_near =0. ;
  double dfdy_i_near =0. ;
  double dfdz_i_near =0. ;
  double alpha_i_near = 0.;
  double f_i1 =0. ;
  double dfdy_i1 =0. ;
  double dfdz_i1 =0. ;
  double alpha_i1 = 0.;
   
 
  int n_deltaN = delta_car_n0.get_dim(1) ;
  int n_mu1N = delta_car_n0.get_dim(0) ;
  int n_deltaP = delta_car_p0.get_dim(1) ;
  int n_mu2P = delta_car_p0.get_dim(0) ;

  //-----------------------------------------
  //----------------TRANCHE  i --------------
  //-----------------------------------------
 
  //REPERAGE DES ZONES:
  int Placei = 0 ; // 0 = 2fluides, 1 = fluideN seul, 2 = fluideP seul

  int i_nearN_i = 0;
  int j_nearN_i = 0;
  int i_nearP_i = 0;
  int j_nearP_i = 0;

  
  if ( y > m_1 ) // Dans cette zone : soit deux fluides, soit fluideN seul (ie np=0)
  {
	// on enregistre l'indice correpondant au delta_car le plus proche de xtab(i_near, j_near, k_near) 
	while ( ( (delta_car_n0)(i_nearN_i,0) <= xtab(i_near, j_near, k_near) ) && ( ( n_deltaN-1 ) > i_nearN_i ) ) {
  	  i_nearN_i++ ;
  	}
  	if (i_nearN_i != 0) { 
	  i_nearN_i -- ; 
 	}
	// on localise maintenant ou se situe mu_n dans la table de mun 
	while ( ( (mu1_n0)(i_nearN_i,j_nearN_i) <= y ) && ( ( n_mu1N-1 ) > j_nearN_i ) ) {
	  j_nearN_i++ ;
  	}
  	if (j_nearN_i != 0) { 
	  j_nearN_i -- ; 
  	}

 
  	// Pour vérifier le positionnement : 
  	//cout << " DEBUG mode = " << (delta_car_n0)(i_nearN_i,0) << "  " <<  xtab( i_near, j_near, k_near)  << endl ;
// 	if ( (delta_car_n0)(i_nearN_i,0) != xtab( i_near, j_near, k_near) ) {
// 	  cout << "c'est un test : " << endl;
// 	  abort();
// 	}
// 	if ( ( (delta_car_n0)(i_nearN_i,0) > xtab(i_near, j_near, k_near)  ) || (xtab(i_near, j_near, k_near)  > (delta_car_n0)(i_nearN_i+1,0) ) ) {
// 	  cout << "mauvais positionnement de delta_car_i dans delta_car_n0 (courbe limite np = 0) " << endl ;
// 	  cout << (delta_car_n0)(i_nearN_i,0) << "  " << xtab(i_near, j_near, k_near) <<  "  " <<  (delta_car_n0)(i_nearN_i+1,0) << endl;
// 	abort();
// 	}
// 	if ( ( (mu1_n0)(i_nearN_i,j_nearN_i) > y  ) || (y > (mu1_n0)(i_nearN_i,j_nearN_i+1) ) ) {
// 	  cout << "mauvais positionnement demu_n dans mu1_n0 (courbe limite np = 0) " << endl ;
// 	  cout << (mu1_n0)(i_nearN_i,j_nearN_i) << "  " << y <<  "  " <<  (mu1_n0)(i_nearN_i,j_nearN_i+1) << endl;
// 	abort();
// 	}
    
    
	// on regarde alors ou on est par rapport a la courbe np = 0.
	double aN_i, bN_i;
	aN_i = ((mu2_n0)(i_nearN_i,j_nearN_i+1) -  (mu2_n0)(i_nearN_i,j_nearN_i) ) / ((mu1_n0)(i_nearN_i,j_nearN_i+1) -  (mu1_n0)(i_nearN_i,j_nearN_i) ) ;
	bN_i = (mu2_n0)(i_nearN_i,j_nearN_i)  - aN_i * (mu1_n0)(i_nearN_i,j_nearN_i) ;
	double zN_i = aN_i * y + bN_i ;

	  if (zN_i <  z) { // Zone ou deux fluides
		Placei = 0;
	  }
	  else {	// Zone ou fluide N seul
		Placei = 1 ;
	  }	
// 	  cout << y * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "  " << 
// 	z* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << "   "  << m_1 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13)<< 
// 	m_2 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13)<<endl;
// 	cout << " Placei = " << Placei<< endl ;

  }
  
  else //if ( y <= m_1) // Dans cette zone : soit deux fluides, soit fluideP seul (ie nn=0), soit zéro fluide !!!!
  {
      
    //cout << y << "    " << m_1 << endl ;
    
	if ( z <= m_2) {
	   Placei = 3 ; // NO fluid at all !
	}
	else {
	
	  // on enregistre l'indice correpondant au delta_car le plus proche de xx
	  while ( ( (delta_car_p0)(i_nearP_i,0) <= xtab(i_near, j_near, k_near)) && ( ( n_deltaP-1 ) > i_nearP_i ) ) {
      		i_nearP_i++ ;
	  }
	  if (i_nearP_i != 0) { 
      		i_nearP_i -- ; 
	  }
	
	  // on localise maintenant ou se situe z dans la table de mup 
	  while ( ( (mu2_p0)(i_nearP_i,j_nearP_i) <= z ) && ( ( n_mu2P-1 ) > j_nearP_i ) ) {
		 j_nearP_i++ ;
	  }
	  if (j_nearP_i != 0) { 
		j_nearP_i -- ; 
	  }

	  // Pour vérifier le positionnement : 
	  // cout << " DEBUG mode = " << (delta_car_p0)(i_nearP_i,0) << "  " <<  xtab( i_near, j_near, k_near)  << endl ;
// 	  if ( (delta_car_p0)(i_nearP_i,0) != xtab( i_near, j_near, k_near) ) {
// 	      cout << "c'est un test : " << endl;
// 	      abort();
// 	  }
// 	    if ( ( (delta_car_p0)(i_nearP_i,0) > xtab(i_near, j_near, k_near)  ) || (xtab(i_near, j_near, k_near)  > (delta_car_p0)(i_nearP_i+1,0) ) ) {
// 	      cout << "mauvais positionnement de delta_car_i dans delta_car_p0 (courbe limite nn = 0) " << endl ;
// 	      cout << (delta_car_p0)(i_nearP_i,0) << "  " << xtab(i_near, j_near, k_near) <<  "  " <<  (delta_car_p0)(i_nearP_i+1,0) << endl;
// 	      abort();
// 	  }
// 	  if ( ( (mu2_p0)(i_nearP_i,j_nearP_i) > y  ) || (y > (mu2_p0)(i_nearP_i,j_nearP_i+1) ) ) {
// 	      cout << "mauvais positionnement demu_p dans mu2_p0 (courbe limite nn = 0) " << endl ;
// 	      cout << (mu2_p0)(i_nearP_i,j_nearP_i) << "  " << y <<  "  " <<  (mu2_p0)(i_nearP_i,j_nearP_i+1) << endl;
// 	      abort();
// 	  }
  		
	
	// on regarde alors ou on est par rapport a la droite nn = 0.
	  double aP_i, bP_i;
	  aP_i = ( (mu2_p0)(i_nearP_i,j_nearP_i+1) -  (mu2_p0)(i_nearP_i,j_nearP_i) ) / ( (mu1_p0)(i_nearP_i,j_nearP_i+1) -  (mu1_p0)(i_nearP_i,j_nearP_i) ) ;
	  bP_i = (mu2_p0)(i_nearP_i,j_nearP_i)  - aP_i * (mu1_p0)(i_nearP_i,j_nearP_i) ;
	  double yP_i = (z- bP_i) /aP_i ;
		

	  if (yP_i <  y) { // Zone ou deux fluides
	    Placei = 0;
	  }
	  else {	// Zone ou fluide 2 seul
	    Placei = 2 ;
	  }
	  
	 // cout << yP_i * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "   " << y * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "  " << 
	 // z* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) <<  endl;
		  
	}
// 	cout << y * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "  " << 
// 	z* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << "   "  << m_1 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13)<< 
// 	m_2 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13)<<endl;
// 	cout << " Placei = " << Placei<< endl ;
	
}	
//FIN DU REPERAGE -----------------------------------------------

  // traitement au cas par cas

 if (Placei == 3 ) { // NO FLUID
  		f_i_near = 0. ;
		dfdy_i_near = 0. ;
		dfdz_i_near = 0. ;
		alpha_i_near = 0.;
  }
 else if (Placei == 1 ) {	
		//cout << "fluide N seul" << endl;
		alpha_i_near = 0.;
		dfdz_i_near = 0. ;
		int i = 0;
		interpol_herm(mu1_N, press_N, n_n_N,
		  y,  i, f_i_near , dfdy_i_near ) ;
 		if (f_i_near < 0.) { 
		  cout << " INTERPOLATION FLUID N --> negative pressure " << endl ;
		  abort();
		  // f_i_near = 0.;
		}
 		if (dfdy_i_near < 0.) { 
		  cout << " INTERPOLATION FLUID N --> negative density " << endl ;
		  abort();
		  // dfdy_i_near = 0.;
		}
  }
  else if (Placei == 2 ) {	
		//cout << "fluide P seul" << endl;
		alpha_i_near = 0.;
		dfdy_i_near = 0. ;
		int i =0;
		interpol_herm( mu2_P, press_P, n_p_P,
		   z,  i, f_i_near,  dfdz_i_near) ;
 		if (f_i_near < 0.) { 
		  cout << " INTERPOLATION FLUID P --> negative pressure " << endl ;
		  abort();
		  // f_i_near = 0.;
		}
 		if (dfdz_i_near < 0.) { 
		  cout << " INTERPOLATION FLUID P --> negative density " << endl ;
		  abort();
		  // dfdz_i_near = 0.;
		}
  }
  else if (Placei == 0 ) {	
    /*
	// les deux fluides sont presents
	// on regarde alors si on a des termes nuls montrant qu'on est a la surface
					
	if (    (dfdytab( i_near, j_near, k_near) > 0. ) && (dfdztab( i_near, j_near, k_near) > 0. ) &&
  		(dfdytab( i_near, j1, k_near) > 0. ) && (dfdztab( i_near, j1, k_near) > 0. ) && 
		(dfdytab( i_near, j_near, k1) > 0. ) && (dfdztab( i_near, j_near, k1) > 0. ) && 
		(dfdytab( i_near, j1, k1) > 0. ) && (dfdztab( i_near, j1, k1) > 0. ) ) {

		//cout << "DEBUG : LOIN DES COURBES LIMITES" << endl;
		interpol_mixed_3d(xtab, ytab, ztab, ftab, 
		      dfdytab, dfdztab, d2fdydztab,
		      xtab(i_near, j_near, k_near), y, z, f_i_near, dfdy_i_near , dfdz_i_near ) ; 
		
		//cout << "f_i_near " << f_i_near << " dfdy_i_near  " << dfdy_i_near << " dfdz_i_near " << dfdz_i_near << endl ;      
		      
		if (f_i_near < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS --> negative pressure " << endl ;
//		  abort();
		}
 		if (dfdy_i_near < 0.) { 
//		   cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS --> negative density for fluid N " << endl ;
//		  abort();
		}
		if (dfdz_i_near < 0.) { 

//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS --> negative density for fluid P " << endl ;
//		  abort();
		}
		    
		double der1 = 0., der2 = 0.;
		interpol_mixed_3d_mod(xtab, ytab, ztab, dlpsddelta_car, 
		     d2lpsdlent1ddelta_car, d2lpsdlent2ddelta_car, 
		      xtab(i_near, j_near, k_near), y, z, alpha_i_near,  der1 , der2 ) ;
		
		alpha_i_near = - alpha_i_near ;
		//cout << " alpha_i_near" <<   alpha_i_near << endl ;
		
		}

 	else if ( (dfdytab( i_near, j_near, k_near) == 0. ) && (dfdztab( i_near, j_near, k_near) == 0. ) &&
		(dfdytab( i_near, j1, k_near) > 0. ) && (dfdztab( i_near, j1, k_near) == 0. ) && 
		(dfdytab( i_near, j_near, k1) == 0. ) && (dfdztab( i_near, j_near, k1) > 0. ) && 
 		(dfdytab( i_near, j1, k1) > 0. ) && (dfdztab( i_near, j1, k1) > 0. ) ) {
	

		//cout << "Croisement des courbes limites - Tranche i" << endl;
//		//********** point pris en haut a gauche
		double aP_inf, bP_inf;
		aP_inf = ( mu2_p0(i_nearP_i,j_nearP_i+1) - mu2_p0(i_nearP_i,j_nearP_i) ) / 
			 ( mu1_p0(i_nearP_i,j_nearP_i+1) - mu1_p0(i_nearP_i,j_nearP_i) ) ;
		bP_inf = mu2_p0(i_nearP_i,j_nearP_i)  - aP_inf * mu1_p0(i_nearP_i,j_nearP_i) ;
		
 		double mu2_nul_inf_Left = ztab(i_near, j1, k1) ;
		double mu1_nul_inf_Left = (mu2_nul_inf_Left - bP_inf) / aP_inf ;
//
//		//********** point pris en bas a droite
//		
		double aN_inf, bN_inf;
 		aN_inf = (mu2_n0(i_nearN_i,j_nearN_i +1) -  mu2_n0(i_nearN_i,j_nearN_i ) ) / 
			 (mu1_n0(i_nearN_i,j_nearN_i +1) -  mu1_n0(i_nearN_i,j_nearN_i ) ) ;
		bN_inf = mu2_n0(i_nearN_i,j_nearN_i )  - aN_inf * mu1_n0(i_nearN_i,j_nearN_i ) ;
 

		double mu1_nul_inf_Right= ytab(i_near, j1, k1) ;
      		double mu2_nul_inf_Right = aN_inf * mu1_nul_inf_Right + bN_inf ;
//
		double mu1_11, mu1_21, mu2_11, mu2_12 ; 
		double p_11,p_21,p_12,p_22 ; 
		double n1_11,n1_21,n1_12, n1_22, n2_11 ,n2_21, n2_12, n2_22 ; 
		double cross_11 , cross_21, cross_12, cross_22 ;
		double dp2sddelta_car_11,dp2sddelta_car_21,dp2sddelta_car_12, dp2sddelta_car_22 ; 
		double d2psdent1ddelta_car_11 ,d2psdent1ddelta_car_21, d2psdent1ddelta_car_12, d2psdent1ddelta_car_22 ;
 		double	d2psdent2ddelta_car_11,d2psdent2ddelta_car_21, d2psdent2ddelta_car_12, d2psdent2ddelta_car_22 ; 
		

		mu1_11 = mu1_nul_inf_Left ;
		mu1_21 = mu1_nul_inf_Right ;
		mu2_11 = mu2_nul_inf_Right ;
		mu2_12 = mu2_nul_inf_Left ; 
		p_21 =  ftab(i_near,j1, k_near) ;
		p_12 =  ftab(i_near,j_near, k1) ;
		p_22 =  ftab(i_near,j1, k1) ; 
 		n1_21 =  dfdytab(i_near,j1, k_near) ;
		n1_12 =  dfdytab(i_near,j_near, k1) ;
		n1_22 =  dfdytab(i_near,j1, k1) ; 
 		n2_21 =  dfdztab(i_near,j1, k_near) ;
		n2_12 =  dfdztab(i_near,j_near, k1) ;
		n2_22 =  dfdztab(i_near,j1, k1) ; 
		cross_21 =  d2fdydztab(i_near,j1, k_near) ; 
		cross_12 =  d2fdydztab(i_near,j_near, k1) ; 
		cross_22 =  d2fdydztab(i_near,j1, k1) ; 
		dp2sddelta_car_21 = dlpsddelta_car(i_near,j1, k_near) ;
		dp2sddelta_car_12 = dlpsddelta_car(i_near,j_near, k1) ;
		dp2sddelta_car_22 = dlpsddelta_car(i_near,j1, k1) ; 
		d2psdent1ddelta_car_21 = d2lpsdlent1ddelta_car(i_near,j1, k_near) ;
		d2psdent1ddelta_car_12 = d2lpsdlent1ddelta_car(i_near,j_near, k1) ;
 		d2psdent1ddelta_car_22 = d2lpsdlent1ddelta_car(i_near,j1, k1) ; 	
		d2psdent2ddelta_car_21 = d2lpsdlent2ddelta_car(i_near,j1, k_near) ;
		d2psdent2ddelta_car_12 = d2lpsdlent2ddelta_car(i_near,j_near, k1) ;
		d2psdent2ddelta_car_22 = d2lpsdlent2ddelta_car(i_near,j1, k1) ; 
		p_11 = 0.;
		n1_11 = 0. ;
		n2_11 =  0.;
		cross_11 = 0.;
		dp2sddelta_car_11 = 0. ;
		d2psdent1ddelta_car_11 = 0. ;	
 		d2psdent2ddelta_car_11 = 0. ;
		
//		// changmeent
 		interpol_herm_2d_new_avec( y, z, 
				mu1_11, mu1_21,  mu2_11,mu2_12, 
				p_11,p_21,  p_12,  p_22,
				n1_11, n1_21,  n1_12,  n1_22, 
				n2_11, n2_21,  n2_12,  n2_22, 
				cross_11, cross_21, cross_12, cross_22,
				f_i_near, dfdy_i_near, dfdz_i_near) ;
		
				
			//	cout << "f_i_near " << f_i_near << " dfdy_i_near  " << dfdy_i_near << " dfdz_i_near " << dfdz_i_near << endl ;
				
		double der1= 0., der2=0.;
		interpol_herm_2d_new_sans( y, z, 
				mu1_11, mu1_21,  mu2_11, mu2_12, 
 				dp2sddelta_car_11, dp2sddelta_car_21, dp2sddelta_car_12, dp2sddelta_car_22,
				d2psdent1ddelta_car_11, d2psdent1ddelta_car_21, d2psdent1ddelta_car_12,	d2psdent1ddelta_car_22,	
				d2psdent2ddelta_car_11, d2psdent2ddelta_car_21, d2psdent2ddelta_car_12, d2psdent2ddelta_car_22,
				alpha_i_near, der1, der2) ;
		alpha_i_near = - alpha_i_near ;
		
		if (f_i_near < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS CAS 2 --> negative pressure " << endl ;
//		  abort();
		}
 		if (dfdy_i_near < 0.) { 
//		   cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS 2 --> negative density for fluid N " << endl ;
//		  abort();
		}
		if (dfdz_i_near < 0.) { 

//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS 2 --> negative density for fluid P " << endl ;
//		  abort();
		}
		
		
		//cout << " alpha_i_near" <<   alpha_i_near << endl ;
	}
	else if ( (dfdytab( i_near, j_near, k_near) == 0. ) && (dfdztab( i_near, j_near, k_near) > 0. ) &&
		(dfdytab( i_near, j_near, k1) == 0. ) && (dfdztab( i_near, j_near, k1) > 0. )  && 
		(dfdytab( i_near, j1, k_near) > 0. ) && (dfdztab( i_near, j1, k_near) > 0. ) && 
		(dfdytab( i_near, j1, k1) > 0. ) && (dfdztab( i_near, j1, k1) > 0. ) ) {
//
		//cout << "cas Fluide de Neutrons seul / Courbe limite - Tranche i" << endl;
////			
		double aP_inf, bP_inf;
 		aP_inf = ( mu2_p0(i_nearP_i,j_nearP_i+1) - mu2_p0(i_nearP_i,j_nearP_i) ) / 
			 ( mu1_p0(i_nearP_i,j_nearP_i+1) - mu1_p0(i_nearP_i,j_nearP_i) ) ;
		bP_inf = mu2_p0(i_nearP_i,j_nearP_i)  - aP_inf * mu1_p0(i_nearP_i,j_nearP_i) ;
		
		double mu2_nul_inf = ztab(i_near, j1, k1) ;				
		double mu1_nul_inf = (mu2_nul_inf - bP_inf)/aP_inf;
		//cout << mu1_nul_inf * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "  " << mu2_nul_inf* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << endl;
		double mu1_11, mu1_21, mu2_11, mu2_12 ; 
		double p_11,p_21,p_12,p_22 ; 
		double n1_11,n1_21,n1_12, n1_22, n2_11 ,n2_21, n2_12, n2_22 ; 
		double cross_11 , cross_21, cross_12, cross_22 ;
 		double dp2sddelta_car_11,dp2sddelta_car_21,dp2sddelta_car_12, dp2sddelta_car_22 ; 
		double d2psdent1ddelta_car_11 ,d2psdent1ddelta_car_21, d2psdent1ddelta_car_12, d2psdent1ddelta_car_22 ;
		double	d2psdent2ddelta_car_11,d2psdent2ddelta_car_21, d2psdent2ddelta_car_12, d2psdent2ddelta_car_22 ; 
//////		

		//cout <<  " aP_inf = " << aP_inf << endl ;
		//cout <<  " bP_inf = " <<  bP_inf << endl ;
		//cout <<  " mu2_nul_inf " << mu2_nul_inf  << endl ;
		//cout <<  " mu1_nul_inf " << mu1_nul_inf << endl ;

		mu1_11 = mu1_nul_inf ;
		mu1_21 = ytab(i_near,j1, k_near) ;
		mu2_11 = ztab(i_near,j1, k_near) ;
		mu2_12 = mu2_nul_inf ; 
		p_21 =  ftab(i_near,j1, k_near) ;
		p_12 =  ftab(i_near,j_near, k1) ;
		p_22 =  ftab(i_near,j1, k1) ; 
		n1_21 =  dfdytab(i_near,j1, k_near) ;
 		n1_12 =  dfdytab(i_near,j_near, k1) ;
		n1_22 =  dfdytab(i_near,j1, k1) ; 
		n2_21 =  dfdztab(i_near,j1, k_near) ;
		n2_12 =  dfdztab(i_near,j_near, k1) ;
		n2_22 =  dfdztab(i_near,j1, k1) ; 
 		cross_21 =  d2fdydztab(i_near,j1, k_near) ;
		cross_12 =  d2fdydztab(i_near,j_near, k1) ;
		cross_22 =  d2fdydztab(i_near,j1, k1) ; 
		p_11 =  ftab(i_near,j_near, k_near) ;	
		n1_11 = 0. ; //dfdytab(i_near,j_near, k_near) ;
		//cout << "n1_11 " << n1_11 << endl ;
		n2_11 = dfdztab(i_near,j_near, k_near) ;
 		cross_11 = 0.; //d2fdydztab(i_near,j_near, k_near) ;
		//cout << "cross_11 " << cross_11 << endl ;
		dp2sddelta_car_21 = dlpsddelta_car(i_near,j1, k_near) ;
		dp2sddelta_car_12 = dlpsddelta_car(i_near,j_near, k1) ;
		dp2sddelta_car_22 = dlpsddelta_car(i_near,j1, k1) ; 
		d2psdent1ddelta_car_21 = d2lpsdlent1ddelta_car(i_near,j1, k_near) ;
		d2psdent1ddelta_car_12 = d2lpsdlent1ddelta_car(i_near,j_near, k1) ;
		d2psdent1ddelta_car_22 = d2lpsdlent1ddelta_car(i_near,j1, k1) ; 
		d2psdent2ddelta_car_21 = d2lpsdlent2ddelta_car(i_near,j1, k_near) ;
		d2psdent2ddelta_car_12 = d2lpsdlent2ddelta_car(i_near,j_near, k1) ;
		d2psdent2ddelta_car_22 = d2lpsdlent2ddelta_car(i_near,j1, k1) ; 
		dp2sddelta_car_11 = 0. ;//dlpsddelta_car(i_near,j_near, k_near)  ;
		//cout << dp2sddelta_car_11 << endl ; 
		d2psdent1ddelta_car_11 =0. ;// d2lpsdlent1ddelta_car(i_near,j_near, k_near) ;	
		//cout << d2psdent1ddelta_car_11 << endl ;
 		d2psdent2ddelta_car_11 = 0. ; //d2lpsdlent2ddelta_car(i_near,j_near, k_near)  ;
		//cout << d2psdent2ddelta_car_11 << endl ;
		
		interpol_herm_2d_new_avec( y, z, 
				mu1_11, mu1_21,  mu2_11,mu2_12, 
 				p_11,p_21,  p_12,  p_22,
				n1_11, n1_21,  n1_12,  n1_22, 
				n2_11, n2_21,  n2_12,  n2_22, 
				cross_11, cross_21, cross_12, cross_22,
				f_i_near, dfdy_i_near, dfdz_i_near) ;
		
		//cout << "f_i_near " << f_i_near << " dfdy_i_near  " << dfdy_i_near << " dfdz_i_near " << dfdz_i_near << endl ;
//		
		double der1= 0., der2=0.;
		
		//cout << "y   " << y * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< " z  " << z * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< endl ;
		interpol_herm_2d_new_sans( y, z, 
				mu1_11, mu1_21,  mu2_11, mu2_12, 
				dp2sddelta_car_11, dp2sddelta_car_21, dp2sddelta_car_12, dp2sddelta_car_22,
				d2psdent1ddelta_car_11, d2psdent1ddelta_car_21, d2psdent1ddelta_car_12,	d2psdent1ddelta_car_22,	
 				d2psdent2ddelta_car_11, d2psdent2ddelta_car_21, d2psdent2ddelta_car_12, d2psdent2ddelta_car_22,
				alpha_i_near, der1, der2) ;
		alpha_i_near = - alpha_i_near ;
		//cout << " alpha_i_near " << alpha_i_near << endl ;
	
				if (f_i_near < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS CAS 3 --> negative pressure " << endl ;
//		  abort();
		}
 		if (dfdy_i_near < 0.) { 
//		   cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS 3 --> negative density for fluid N " << endl ;
//		  abort();
		}
		if (dfdz_i_near < 0.) { 

//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS 3 --> negative density for fluid P " << endl ;
//		  abort();
		}
		
//		
	}
	else {
	  
*/

//cout << "TRAITEMENT DE LA SURFACE PAS ENCORE REALISE - Tranche i" << endl;
		interpol_mixed_3d(xtab, ytab, ztab, ftab, 
		      dfdytab, dfdztab, d2fdydztab,
		      xtab(i_near, j_near, k_near), y, z, f_i_near, dfdy_i_near , dfdz_i_near ) ; 
		      
		      //cout << "f_i_near " << f_i_near << " dfdy_i_near  " << dfdy_i_near << " dfdz_i_near " << dfdz_i_near << endl ;
		      
		double der1 = 0., der2 = 0.;
		interpol_mixed_3d_mod(xtab, ytab, ztab, dlpsddelta_car, 
		     d2lpsdlent1ddelta_car, d2lpsdlent2ddelta_car, 
		      xtab(i_near, j_near, k_near), y, z, alpha_i_near,  der1 , der2 ) ;
		 //     cout << "alpha " << alpha  << endl; 
		alpha_i_near = - alpha_i_near ;
		
/*		if (f_i_near < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS CAS NON TRAITE --> negative pressure " << endl ;
//		  abort();
		}
 		if (dfdy_i_near < 0.) { 
//		   cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS NON TRAITE --> negative density for fluid N " << endl ;
//		  abort();
		}
		if (dfdz_i_near < 0.) { 

//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS NON TRAITE --> negative density for fluid P " << endl ;
//		  abort();
		}
*/		
		
		//cout << " alpha_i_near" <<   alpha_i_near << endl ;
//			//cout << "surface non 6" << endl;
//		// Zone en plein milieu de la table a deux fluides
//	


/*}*/


}



  //-----------------------------------------
  //--------------TRANCHE  i+1 --------------
  //-----------------------------------------
  

  
  
 //REPERAGE DES ZONES:
  int Placei1 = 0 ; // 0 = 2fluides, 1 = fluideN seul, 2 = fluideP seul, 3 = 0 fluide !

  int i_nearN_i1 = 0;
  int j_nearN_i1 = 0;
  int i_nearP_i1 = 0;
  int j_nearP_i1 = 0; 
  
  if ( y > m_1 ) // Dans cette zone : soit deux fluides, soit fluideN seul (ie np=0)
  {
	
	// on enregistre l'indice correpondant au delta_car le plus proche de xtab(i_near, j_near, k_near) 
	while ( ( (delta_car_n0)(i_nearN_i1,0) <= xtab(i_near+1, j_near, k_near) ) && ( ( n_deltaN-1 ) > i_nearN_i1 ) ) {
  	  i_nearN_i1++ ;
  	}
  	if (i_nearN_i1 != 0) { 
	  i_nearN_i1 -- ; 
	}
	// on localise maintenant ou se situe mu_n dans la table de mun 	
	while ( ( (mu1_n0)(i_nearN_i1,j_nearN_i1) <= y ) && ( ( n_mu1N-1 ) > j_nearN_i1 ) ) {
	  j_nearN_i1++ ;
  	}
  	if (j_nearN_i1 != 0) { 
	  j_nearN_i1 -- ; 
  	}

  	// Pour vérifier le positionnement : 
//   	cout << " DEBUG mode = " << (delta_car_n0)(i_nearN_i1,0) << "  " <<  xtab( i_near+1, j_near, k_near)  << endl ;
// 	if ( (delta_car_n0)(i_nearN_i1,0) != xtab( i_near+1, j_near, k_near) ) {
// 	  cout << "c'est un test : " << endl;
// 	  abort();
// 	}
// 	if ( ( (delta_car_n0)(i_nearN_i1,0) > xtab(i_near+1, j_near, k_near)  ) || (xtab(i_near+1, j_near, k_near)  > (delta_car_n0)(i_nearN_i1+1,0) ) ) {
// 	  cout << "mauvais positionnement de delta_car_i+1 dans delta_car_n0 (courbe limite np = 0) " << endl ;
// 	  cout << (delta_car_n0)(i_nearN_i1,0) << "  " << xtab(i_near+1, j_near, k_near) <<  "  " <<  (delta_car_n0)(i_nearN_i1+1,0) << endl;
// 	abort();
// 	}
// 	if ( ( (mu1_n0)(i_nearN_i1,j_nearN_i1) > y  ) || (y > (mu1_n0)(i_nearN_i1,j_nearN_i1+1) ) ) {
// 	  cout << "mauvais positionnement demu_n dans mu1_n0 (courbe limite np = 0) " << endl ;
// 	  cout << (mu1_n0)(i_nearN_i1,j_nearN_i1) << "  " << y <<  "  " <<  (mu1_n0)(i_nearN_i1,j_nearN_i1+1) << endl;
// 	abort();
// 	}
    
	// on regarde alors ou on est par rapport a la courbe np = 0.
	double aN_i1, bN_i1;
	aN_i1 = ((mu2_n0)(i_nearN_i1,j_nearN_i1+1) -  (mu2_n0)(i_nearN_i1,j_nearN_i1) ) / ((mu1_n0)(i_nearN_i1,j_nearN_i1+1) -  (mu1_n0)(i_nearN_i1,j_nearN_i1) ) ;
	bN_i1 = (mu2_n0)(i_nearN_i1,j_nearN_i1)  - aN_i1 * (mu1_n0)(i_nearN_i1,j_nearN_i1) ;
	double zN_i1 = aN_i1 * y + bN_i1 ;

	if (zN_i1 <  z) {	// Zone ou deux fluides
	  Placei1 = 0;
	}
	else {			// Zone ou fluide N seul
	  Placei1 = 1 ;
	}	
		
// 	cout << y * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "  " << 
// 	z* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << "   "  << m_1 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13)<< 
// 	m_2 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13)<<endl;
// 	cout << " Placei1 = " << Placei1<< endl ;
  }
 
  else //if ( y <= m_1) // Dans cette zone : soit deux fluides, soit fluideP seul (ie nn=0), soit zéro fluide !!!!
  {
      	if ( z <= m_2) {
	   Placei1 = 3 ; // NO fluid at all !
	}
	else {
		
	  
	  // on enregistre l'indice correpondant au delta_car le plus proche de xx
	  while ( ( (delta_car_p0)(i_nearP_i1,0) <= xtab(i_near+1, j_near, k_near)) && ( ( n_deltaP-1 ) > i_nearP_i1) ) {
	    i_nearP_i1++ ;
	  }
	  if (i_nearP_i1 != 0) { 
	    i_nearP_i1 -- ; 
	  }
	  // on localise maintenant ou se situe z dans la table de mup 		
	  while ( ( (mu2_p0)(i_nearP_i1,j_nearP_i1) <= z ) && ( ( n_mu2P-1 ) > j_nearP_i1 ) ) {
	    j_nearP_i1++ ;
	  }
	  if (j_nearP_i1 != 0) { 
	    j_nearP_i1 -- ; 
	  }

	  // Pour vérifier le positionnement : 
// 	  cout << " DEBUG mode = " << (delta_car_p0)(i_nearP_i1,0) << "  " <<  xtab( i_near+1, j_near, k_near)  << endl ;
// 	  if ( (delta_car_p0)(i_nearP_i1,0) != xtab( i_near+1, j_near, k_near) ) {
// 	      cout << "c'est un test : " << endl;
// 	      abort();
// 	  }
// 	    if ( ( (delta_car_p0)(i_nearP_i1,0) > xtab(i_near+1, j_near, k_near)  ) || (xtab(i_near+1, j_near, k_near)  > (delta_car_p0)(i_nearP_i1+1,0) ) ) {
// 	      cout << "mauvais positionnement de delta_car_i+1 dans delta_car_p0 (courbe limite nn = 0) " << endl ;
// 	      cout << (delta_car_p0)(i_nearP_i1,0) << "  " << xtab(i_near+1, j_near, k_near) <<  "  " <<  (delta_car_p0)(i_nearP_i1+1,0) << endl;
// 	      abort();
// 	  }
// 	  if ( ( (mu2_p0)(i_nearP_i1,j_nearP_i1) > y  ) || (y > (mu2_p0)(i_nearP_i1,j_nearP_i1+1) ) ) {
// 	      cout << "mauvais positionnement demu_p dans mu2_p0 (courbe limite nn = 0) " << endl ;
// 	      cout << (mu2_p0)(i_nearP_i1,j_nearP_i1) << "  " << y <<  "  " <<  (mu2_p0)(i_nearP_i1,j_nearP_i1+1) << endl;
// 	      abort();
// 	  }
	
	  // on regarde alors ou on est par rapport a la droite nn = 0.
	  double aP_i1, bP_i1;
	  aP_i1 = ( (mu2_p0)(i_nearP_i1,j_nearP_i1+1) -  (mu2_p0)(i_nearP_i1,j_nearP_i1) ) / ( (mu1_p0)(i_nearP_i1,j_nearP_i1+1) -  (mu1_p0)(i_nearP_i1,j_nearP_i1) ) ;
	  bP_i1 = (mu2_p0)(i_nearP_i1,j_nearP_i1)  - aP_i1 * (mu1_p0)(i_nearP_i1,j_nearP_i1) ;
	  double yP_i1 = (z- bP_i1) /aP_i1 ;
		
	  if (yP_i1 <  y) {   // Zone ou deux fluides
	    Placei1 = 0;
	  }
	  else {  // Zone ou fluide 2 seul
	    Placei1 = 2 ;
	  }
		
	 //cout << yP_i1 * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "   " << y * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "  " << 
	 //z* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) <<  endl;
		  
	}
	
// 	cout << y * 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13)<< "  " << 
// 	z* 2.99792458E+8 * 2.99792458E+8 * 1.66e-27 /(1.6021892E-13) << "   "  << m_1 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13)<< 
// 	m_2 *2.99792458E+8 * 2.99792458E+8 * 1.66e-27  / (1.6021892E-13)<<endl;
// 	cout << " Placei1 = " << Placei1<< endl ;
	
}	




//FIN DU REPERAGE -----------------------------------------------


  // traitement au cas par cas
  if (Placei1 == 3 ) { // NO FLUID
  		f_i1 = 0. ;
		dfdy_i1 = 0. ;
		dfdz_i1 = 0. ;
		alpha_i1 = 0.;
  }
  else if (Placei1 == 1 ) {	
		//cout << "fluideN seul" << endl;
		alpha_i1 = 0. ;
		dfdz_i1 = 0. ;
		int i =0;
		interpol_herm(mu1_N, press_N, n_n_N,
		  y,  i, f_i1 , dfdy_i1 ) ;
		if (f_i1 < 0.) { 
		  cout << " INTERPOLATION FLUID N i+1 --> negative pressure " << endl ;
		  abort();
		  // f_i1 = 0.;
		}
 		if (dfdy_i1 < 0.) { 
		  cout << " INTERPOLATION FLUID N i+1--> negative density " << endl ;
		  abort();
		  // dfdy_i1 = 0.;
		}
  }
  else if (Placei1 == 2 ) {	
		//cout << "fluideP seul" << endl;
		alpha_i1 = 0.;
		dfdy_i1 = 0. ;
		int i =0;
		interpol_herm( mu2_P, press_P, n_p_P,
		   z,  i, f_i1,  dfdz_i1) ;
		if (f_i1 < 0.) { 
		  cout << " INTERPOLATION FLUID P i+1--> negative pressure " << endl ;
		  abort();
		  // f_i1 = 0.;
		}
 		if (dfdz_i1 < 0.) { 
		  cout << " INTERPOLATION FLUID P i+1 --> negative density " << endl ;
		  abort();
		  // dfdz_i1= 0.;
		}
}
else if (Placei1 == 0 ) {
  /*
	// les deux fluides sont presents
	// on regarde alors si on a des termes nuls montrant qu'on est a la surface
		
	// ----------pour savoir si on est proche ou non d'une surface


	if ( (dfdytab( i1, j_near, k_near) > 0. ) && (dfdztab( i1, j_near, k_near) > 0. ) &&
		(dfdytab( i1, j1, k_near) > 0. ) && (dfdztab( i1, j1, k_near) > 0. ) && 
		(dfdytab( i1, j_near, k1) > 0. ) && (dfdztab( i1, j_near, k1) > 0. ) && 
 		(dfdytab( i1, j1, k1) > 0. ) && (dfdztab( i1, j1, k1) > 0. ) ) {

		//cout << "DEBUG : LOIN DES COURBES LIMITES" << endl;
		interpol_mixed_3d(xtab, ytab, ztab, ftab, 
		      dfdytab, dfdztab, d2fdydztab,
		      xtab(i1, j_near, k_near), y, z, f_i1, dfdy_i1 , dfdz_i1) ; 
	
		if (f_i1 < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS i+1--> negative pressure " << endl ;
//		  abort();
		}
 		if (dfdy_i1 < 0.) { 
//		   cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS i+1--> negative density for fluid N " << endl ;
//		  abort();
		}
		if (dfdz_i1 < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS i+1--> negative density for fluid P " << endl ;
//		  abort();
		}      
		      
		double der1 = 0., der2 = 0.;
		interpol_mixed_3d_mod(xtab, ytab, ztab, dlpsddelta_car, 
		     d2lpsdlent1ddelta_car, d2lpsdlent2ddelta_car, 
		      xtab(i1, j_near, k_near), y, z, alpha_i1,  der1 , der2 ) ;
		alpha_i1 = - alpha_i1 ;
		}
//
 	else if ( (dfdytab( i1, j_near, k_near) == 0. ) && (dfdztab( i1, j_near, k_near) == 0. ) &&
		(dfdytab( i1, j1, k_near) > 0. ) && (dfdztab( i1, j1, k_near) == 0. ) && 
		(dfdytab( i1, j_near, k1) == 0. ) && (dfdztab( i1, j_near, k1) > 0. ) && 
 		(dfdytab( i1, j1, k1) > 0. ) && (dfdztab( i1, j1, k1) > 0. ) ) {
//	
//
		//cout << "Croisement des courbes limites - Tranche i + 1 " << endl;
//		//********** point pris en haut a gauche
//
//
		double aP_sup, bP_sup;
		aP_sup = ( mu2_p0(i_nearP_i1,j_nearP_i1+1) - mu2_p0(i_nearP_i1,j_nearP_i1) ) / 
			 ( mu1_p0(i_nearP_i1,j_nearP_i1+1) - mu1_p0(i_nearP_i1,j_nearP_i1) ) ;
		bP_sup = mu2_p0(i_nearP_i1,j_nearP_i1)  - aP_sup * mu1_p0(i_nearP_i1,j_nearP_i1) ;

     		double mu2_nul_sup_Left = ztab(i1, j1, k1) ;
		double mu1_nul_sup_Left = (mu2_nul_sup_Left - bP_sup) / aP_sup ;	
	

		//********** point pris en bas a droite
	
		double aN_sup, bN_sup;
		aN_sup = (mu2_n0(i_nearN_i1,j_nearN_i1 +1) -  mu2_n0(i_nearN_i1,j_nearN_i1) ) / 
			 (mu1_n0(i_nearN_i1,j_nearN_i1 +1) -  mu1_n0(i_nearN_i1,j_nearN_i1) ) ;
		bN_sup = mu2_n0(i_nearN_i1,j_nearN_i1 )  - aN_sup * mu1_n0(i_nearN_i1,j_nearN_i1 ) ;

		double mu1_nul_sup_Right = ytab(i1, j1, k1) ;
      		double mu2_nul_sup_Right= aN_sup * mu1_nul_sup_Right + bN_sup;		

		double mu1_11, mu1_21, mu2_11, mu2_12 ; 
		double p_11,p_21,p_12,p_22 ; 
		double n1_11,n1_21,n1_12, n1_22, n2_11 ,n2_21, n2_12, n2_22 ; 
		double cross_11 , cross_21, cross_12, cross_22 ;
		double dp2sddelta_car_11,dp2sddelta_car_21,dp2sddelta_car_12, dp2sddelta_car_22 ; 
		double d2psdent1ddelta_car_11 ,d2psdent1ddelta_car_21, d2psdent1ddelta_car_12, d2psdent1ddelta_car_22 ;
		double	d2psdent2ddelta_car_11,d2psdent2ddelta_car_21, d2psdent2ddelta_car_12, d2psdent2ddelta_car_22 ; 
		
		mu1_11 = mu1_nul_sup_Left ;
		mu1_21 = mu1_nul_sup_Right ;
		mu2_11 = mu2_nul_sup_Right ;
		mu2_12 = mu2_nul_sup_Left ; 
		p_21 =  ftab(i1,j1, k_near) ;
		p_12 =  ftab(i1,j_near, k1) ;
		p_22 =  ftab(i1,j1, k1) ; 
		n1_21 =  dfdytab(i1,j1, k_near) ;
		n1_12 =  dfdytab(i1,j_near, k1) ;
		n1_22 =  dfdytab(i1,j1, k1) ; 
		n2_21 =  dfdztab(i1,j1, k_near) ;
		n2_12 =   dfdztab(i1,j_near, k1) ;
		n2_22 =  dfdztab(i1,j1, k1) ; 
		cross_21 =  d2fdydztab(i1,j1, k_near) ; 
		cross_12 =  d2fdydztab(i1,j_near, k1) ; 
		cross_22 =  d2fdydztab(i1,j1, k1) ; 
		dp2sddelta_car_21 = dlpsddelta_car(i1,j1, k_near) ;
		dp2sddelta_car_12 = dlpsddelta_car(i1,j_near, k1) ;
		dp2sddelta_car_22 = dlpsddelta_car(i1,j1, k1) ; 
		d2psdent1ddelta_car_21 = d2lpsdlent1ddelta_car(i1,j1, k_near) ;
		d2psdent1ddelta_car_12 = d2lpsdlent1ddelta_car(i1,j_near, k1) ;
		d2psdent1ddelta_car_22 = d2lpsdlent1ddelta_car(i1,j1, k1) ; 	
		d2psdent2ddelta_car_21 = d2lpsdlent2ddelta_car(i1,j1, k_near) ;
		d2psdent2ddelta_car_12 = d2lpsdlent2ddelta_car(i1,j_near, k1) ;
		d2psdent2ddelta_car_22 = d2lpsdlent2ddelta_car(i1,j1, k1) ; 
		p_11 = 0.;
		n1_11 = 0. ;
		n2_11 =  0.;
		cross_11 = 0.;
		dp2sddelta_car_11 = 0. ;
		d2psdent1ddelta_car_11 = 0. ;	
		d2psdent2ddelta_car_11 = 0. ;
		
	
		// changement
		interpol_herm_2d_new_avec( y, z, 
				mu1_11, mu1_21,  mu2_11, mu2_12, 
				p_11,p_21,  p_12,  p_22,
				n1_11, n1_21,  n1_12,  n1_22, 
				n2_11, n2_21,  n2_12,  n2_22, 
				cross_11, cross_21, cross_12, cross_22,
				f_i1, dfdy_i1, dfdz_i1) ;
		
		double der1= 0., der2=0.;
		interpol_herm_2d_new_sans( y, z, 
				mu1_11, mu1_21,  mu2_11, mu2_12, 
				dp2sddelta_car_11, dp2sddelta_car_21, dp2sddelta_car_12, dp2sddelta_car_22,
				d2psdent1ddelta_car_11, d2psdent1ddelta_car_21, d2psdent1ddelta_car_12,	d2psdent1ddelta_car_22,	
				d2psdent2ddelta_car_11, d2psdent2ddelta_car_21, d2psdent2ddelta_car_12, d2psdent2ddelta_car_22,
				alpha_i1, der1, der2) ;
		alpha_i1 = - alpha_i1 ;
		
		if (f_i1 < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS CAS 2 --> negative pressure " << endl ;
//		  abort();
		}
 		if (dfdy_i1 < 0.) { 
//		   cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS 2 --> negative density for fluid N " << endl ;
//		  abort();
		}
		if (dfdz_i1 < 0.) { 

//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS 2 --> negative density for fluid P " << endl ;
//		  abort();
		}
		
		
	
	}
//
	else if ( (dfdytab( i_near+1, j_near, k_near) == 0. ) && (dfdztab( i_near+1, j_near, k_near) > 0. ) &&
		(dfdytab( i_near+1, j_near, k1) == 0. ) && (dfdztab( i_near+1, j_near, k1) > 0. )  &&
		(dfdytab( i_near+1, j1, k_near) > 0. ) && (dfdztab( i_near+1, j1, k_near) > 0. ) && 
 		(dfdytab( i_near+1, j1, k1) > 0. ) && (dfdztab( i_near+1, j1, k1) > 0. ) ) {

		//cout << "cas Fluide de Neutrons seul / Courbe limite - Tranche i +1 " << endl;
	
		double aP_sup, bP_sup;
		aP_sup = (mu2_p0(i_nearP_i+1,j_nearP_i+1) - mu2_p0(i_nearP_i+1,j_nearP_i) ) / 
			( mu1_p0(i_nearP_i+1,j_nearP_i+1) - mu1_p0(i_nearP_i+1,j_nearP_i) ) ;
		bP_sup = mu2_p0(i_nearP_i+1,j_nearP_i)  - aP_sup * mu1_p0(i_nearP_i+1,j_nearP_i) ;

		
		double mu2_nul_sup = ztab(i1, j1, k1) ;				
		double mu1_nul_sup = (mu2_nul_sup - bP_sup)/aP_sup;
		double mu1_11, mu1_21, mu2_11, mu2_12 ; 
		double p_11,p_21,p_12,p_22 ; 
		double n1_11,n1_21,n1_12, n1_22, n2_11 ,n2_21, n2_12, n2_22 ; 
		double cross_11 , cross_21, cross_12, cross_22 ;
		double dp2sddelta_car_11,dp2sddelta_car_21,dp2sddelta_car_12, dp2sddelta_car_22 ; 
		double d2psdent1ddelta_car_11 ,d2psdent1ddelta_car_21, d2psdent1ddelta_car_12, d2psdent1ddelta_car_22 ;
		double	d2psdent2ddelta_car_11,d2psdent2ddelta_car_21, d2psdent2ddelta_car_12, d2psdent2ddelta_car_22 ; 
		
		
		mu1_11 = mu1_nul_sup ;
		mu1_21 = ytab(i1,j1, k_near) ;
		mu2_11 = ztab(i1,j1, k_near) ;
		mu2_12 = mu2_nul_sup ; 
		p_21 =  ftab(i1,j1, k_near) ;
		p_12 =  ftab(i1,j_near, k1) ;
		p_22 =  ftab(i1,j1, k1) ; 
		n1_21 =  dfdytab(i1,j1, k_near) ;
		n1_12 =  dfdytab(i1,j_near, k1) ;
		n1_22 =  dfdytab(i1,j1, k1) ; 
		n2_21 =  dfdztab(i1,j1, k_near) ;
		n2_12 =  dfdztab(i1,j_near, k1) ;
		n2_22 =  dfdztab(i1,j1, k1) ; 
		cross_21 =  d2fdydztab(i1,j1, k_near) ;
		cross_12 =  d2fdydztab(i1,j_near, k1) ;
		cross_22 =  d2fdydztab(i1,j1, k1) ; 
		p_11 =  ftab(i1,j_near, k_near) ;	
		n1_11 = 0.;// dfdytab(i1,j_near, k_near) ;
		n2_11 = dfdztab(i1,j_near, k_near) ;
		cross_11 = 0.; //d2fdydztab(i1,j_near, k_near) ;
		dp2sddelta_car_21 = dlpsddelta_car(i1,j1, k_near) ;
		dp2sddelta_car_12 = dlpsddelta_car(i1,j_near, k1) ;
		dp2sddelta_car_22 = dlpsddelta_car(i1,j1, k1) ; 
	        d2psdent1ddelta_car_21 = d2lpsdlent1ddelta_car(i1,j1, k_near) ;
		d2psdent1ddelta_car_12 = d2lpsdlent1ddelta_car(i1,j_near, k1) ;
		d2psdent1ddelta_car_22 = d2lpsdlent1ddelta_car(i1,j1, k1) ; 	
		d2psdent2ddelta_car_21 = d2lpsdlent2ddelta_car(i1,j1, k_near) ;
		d2psdent2ddelta_car_12 = d2lpsdlent2ddelta_car(i1,j_near, k1) ;
		d2psdent2ddelta_car_22 = d2lpsdlent2ddelta_car(i1,j1, k1) ;
		dp2sddelta_car_11 = 0.; //dlpsddelta_car(i1,j_near, k_near) ;
		d2psdent1ddelta_car_11 = 0.;//d2lpsdlent1ddelta_car(i1,j_near, k_near) ;	
		d2psdent2ddelta_car_11 = 0.;//d2lpsdlent2ddelta_car(i1,j_near, k_near) ;	
		 
		interpol_herm_2d_new_avec( y, z, 
				mu1_11, mu1_21,  mu2_11,mu2_12, 
				p_11,p_21,  p_12,  p_22,
				n1_11, n1_21,  n1_12,  n1_22, 
				n2_11, n2_21,  n2_12,  n2_22, 
				cross_11, cross_21, cross_12, cross_22,
				f_i1, dfdy_i1, dfdz_i1) ;
		
		double der1= 0., der2=0.;
		interpol_herm_2d_new_sans( y, z, 
				mu1_11, mu1_21,  mu2_11, mu2_12, 
				dp2sddelta_car_11, dp2sddelta_car_21, dp2sddelta_car_12, dp2sddelta_car_22,
				d2psdent1ddelta_car_11, d2psdent1ddelta_car_21, d2psdent1ddelta_car_12,	d2psdent1ddelta_car_22,	
				d2psdent2ddelta_car_11, d2psdent2ddelta_car_21, d2psdent2ddelta_car_12, d2psdent2ddelta_car_22,
				alpha_i1, der1, der2) ;
		alpha_i1 = - alpha_i1 ;

		if (f_i1 < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS CAS 3 --> negative pressure " << endl ;
//		  abort();
		}
 		if (dfdy_i1 < 0.) { 
//		   cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS 3 --> negative density for fluid N " << endl ;
//		  abort();
		}
		if (dfdz_i1 < 0.) { 

//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS 3 --> negative density for fluid P " << endl ;
//		  abort();
		}
		
}

	else { */
	  

	  
		//cout << "TRAITEMENT DE LA SURFACE PAS ENCORE REALISE - Tranche i + 1 " << endl;
		interpol_mixed_3d(xtab, ytab, ztab, ftab, 
		      dfdytab, dfdztab, d2fdydztab,
		      xtab(i1, j_near, k_near), y, z, f_i1, dfdy_i1 , dfdz_i1 ) ; 
	
		double der1b = 0., der2b = 0.;
		interpol_mixed_3d_mod(xtab, ytab, ztab, dlpsddelta_car, 
		     d2lpsdlent1ddelta_car, d2lpsdlent2ddelta_car, 
		      xtab(i1, j_near, k_near), y, z, alpha_i1,  der1b , der2b ) ;
		alpha_i1 = - alpha_i1 ;
			//cout << "surface non 6" << endl;
		// Zone en plein milieu de la table a deux fluides
/*		
		if (f_i1 < 0.) { 
//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS CAS NON TRAITE --> negative pressure " << endl ;
//		  abort();
		}
 		if (dfdy_i1 < 0.) { 
//		   cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS NON TRAITE --> negative density for fluid N " << endl ;
//		  abort();
		}
		if (dfdz_i1 < 0.) { 

//		  cout << " DEBUG : FAR FROM LIMIT CURVES " << endl;
		  cout << " INTERPOLATION 2 FLUIDS  CAS NON TRAITE --> negative density for fluid P " << endl ;
//		  abort();
		}
*/

 /*
	}*/
//
}


 //--------------------------------------------
 // ----- Linear Interpolation in Delta2 ------
 // -------------------------------------------
 
 double x1  = xtab(i_near, 0, 0) ;
 double x2  = xtab(i1, 0, 0) ;
 double x12 = x1-x2 ;
   
 //for f
 double y1 = f_i_near;
 double y2 = f_i1;
 double a  = (y1-y2)/x12 ;
 double b  = (x1*y2-y1*x2)/x12 ;

 f  = x*a+b ; 
  /*cout << " x1 " << x1 << endl ;
  cout << " x2 " << x2 << endl ;
  cout << " x " << x << endl ;
  cout << " y1 " << y1 << endl ;
  cout << " y2 " << y2 << endl ;
  cout << " y " <<f << endl ;*/

  //for df/dy
  double y1_y = dfdy_i_near;
  double y2_y = dfdy_i1;
  double a_y  = (y1_y-y2_y)/x12 ;
  double b_y  = (x1*y2_y-y1_y*x2)/x12 ;
	
  dfdy  = x*a_y+b_y ; 

  /*cout << " x1 " << x1 << endl ;
  cout << " x2 " << x2 << endl ;
  cout << " x " << x << endl ;
  cout << " y1_y " << y1_y << endl ;
  cout << " y2_y " << y2_y  << endl ;
  cout << " dfdy " << dfdy << endl ;*/
  
  //for df/dz
  double y1_z = dfdz_i_near;
  double y2_z = dfdz_i1;
  double a_z  = (y1_z-y2_z)/x12 ;
  double b_z  = (x1*y2_z-y1_z*x2)/x12 ;
	
  dfdz  = x*a_z+b_z ; 
  
  /*cout << " x1 " << x1 << endl ;
  cout << " x2 " << x2 << endl ;
  cout << " x " << x << endl ;
  cout << " y1_z " << y1_z  << endl ;
  cout << " y2_z " << y2_z << endl ;
  cout << " dfdz " <<dfdz  << endl ;*/
  
  // for alpha 
  double y1_alpha = alpha_i_near;
  double y2_alpha = alpha_i1;
  double a_alpha  = (y1_alpha-y2_alpha)/x12 ;
  double b_alpha  = (x1*y2_alpha-y1_alpha*x2)/x12 ;
	
  alpha  = x*a_alpha+b_alpha ; 

  /*cout << " x1 " << x1 << endl ;
  cout << " x2 " << x2 << endl ;
  cout << " x " << x << endl ;
  cout << " y1_alpha " <<y1_alpha<< endl ;
  cout << " y2_alpha " << y2_alpha  << endl ;
  cout << " alpha " << alpha<< endl ;*/
  
  // A vérifier mais on devrait pouvoir enlever ca :
/*if (( y <= m_1)&& ( z <= m_2)) {
alpha = 0. ;
f = 0.;
dfdy = 0.;
dfdz = 0.;


} */

       return;

}


  //--------------------------------------------------------------------
  // Quintic Hermite interpolation using data from the second derivative
  //--------------------------------------------------------------------
void interpol_herm_2nd_der(const Tbl& xtab, const Tbl& ytab, const Tbl& dytab,
			   const Tbl& d2ytab, double x, int& i, double& y, 
			   double& dy) {

	assert(ytab.dim == xtab.dim) ;
	assert(dytab.dim == xtab.dim) ;	
	assert(d2ytab.dim == xtab.dim) ;	
	
	huntm(xtab, x, i) ;
	
	int i1 = i + 1 ;
	
	double dx = xtab(i1) - xtab(i) ;

	double u = (x - xtab(i)) / dx ;
	double u2 = u*u ;
	double u3 = u2*u ;
	double u4 = u2*u2 ;
	double u5 = u3*u2 ;

	double v = 1. - u ;
	double v2 = v*v ;
	double v3 = v2*v ;
	double v4 = v2*v2 ;
	double v5 = v3*v2 ;

	y =   ytab(i) * ( -6.*u5 + 15.*u4 - 10.*u3 + 1. )
	  + ytab(i1) * ( -6.*v5 + 15.*v4 - 10.*v3 + 1. )
	  + dytab(i) * dx * ( -3.*u5 + 8.*u4 -6.*u3 + u )
	  - dytab(i1) * dx * ( -3.*v5 + 8.*v4 -6.*v3 + v ) 
	  + d2ytab(i) * dx*dx * ( -0.5*u5 + 1.5*u4 - 1.5*u3 + 0.5*u2 )
	  + d2ytab(i1) * dx*dx * ( -0.5*v5 + 1.5*v4 - 1.5*v3 + 0.5*v2 ) ; 
     	
 	dy = 30.*( ytab(i) / dx * ( -u4 + 2.*u3 - u2 ) 
		   - ytab(i1) / dx * ( -v4 + 2.*v3 - v2 ) )  
	  + dytab(i) * ( -15.*u4 + 32.*u3 - 18.*u2 + 1. ) 
	  + dytab(i1) * ( -15.*v4 + 32.*v3 - 18.*v2 + 1. ) 
	  + d2ytab(i) * dx * ( -2.5*u4 + 6.*u3 -4.5*u2 + u ) 
	  - d2ytab(i1) * dx * ( -2.5*v4 + 6.*v3 -4.5*v2 + v ) ;

}

} // End of namespace Lorene

