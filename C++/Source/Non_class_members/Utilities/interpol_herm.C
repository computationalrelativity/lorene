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


 

/*
 * $Id$
 * $Log$
 * Revision 1.19  2024/01/26 17:44:25  g_servignat
 * Updated the Pseudopolytrope_1D class to be consistent with the paper (i.e. with a GPP in the middle)
 *
 * Revision 1.18  2022/03/22 13:19:53  g_servignat
 * Modified row and column extraction of 2D Tbl's
 *
 * Revision 1.17  2022/03/01 10:03:38  g_servignat
 * Added subtbl extraction for interpol_linear_2D purposes
 *
 * Revision 1.16  2022/02/25 10:45:21  g_servignat
 * Added 2D linear interpolation
 *
 * Revision 1.15  2021/05/14 15:39:23  g_servignat
 * Added sound speed computation from enthalpy to Eos class and tabulated+polytropic derived classes
 *
 * Revision 1.14  2016/12/05 16:18:11  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.13  2015/06/15 15:08:22  j_novak
 * New file interpol_bifluid for interpolation of 2-fluid EoSs
 *
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

  //---------------------
  // Linear interpolation 
  //---------------------
  void interpol_linear_2d(const Tbl& x1tab, const Tbl& x2tab, const Tbl& ytab, 
		       double x1, double x2, int& i, int& j, double& y) {
    
    assert(ytab.dim.ndim == 2) ;
    assert(x1tab.dim.ndim == 1) ;
    assert(x2tab.dim.ndim == 1) ;
    assert(ytab.dim.dim[1] == x2tab.dim.dim[0]) ;
    assert(ytab.dim.dim[0] == x1tab.dim.dim[0]) ;
    
    huntm(x1tab, x1, i) ;
    huntm(x2tab, x2, j) ;
    
    int i1 = i + 1 ;
    int j1 = j + 1 ;
    
    // double dx  = xtab(i1) - xtab(i) ;
    double y11  = ytab(j,i) ;
    double y21  = ytab(j1,i) ;
    double y12  = ytab(j,i1) ;
    
    double x11  = x1tab(i) ;
    double x12  = x1tab(i1) ;
    double x21  = x2tab(j) ;
    double x22  = x2tab(j1) ;
    double x1diff = x12-x11 ;
    double x2diff = x22-x21 ;
    
    double a  = (y21-y11)/x2diff ;
    double b  = (y12-y11)/x1diff ;
	
    y  = y11 + (x1-x11)*b + (x2-x21)*a ; 
    
  }
  
  //------------------------
  // Quadratic interpolation
  //------------------------
  void interpol_quad(const Tbl& xtab, const Tbl& ytab, 
		       double x, int& i, double& y) {
    
    assert(ytab.dim == xtab.dim) ;
    //assert(dytab.dim == xtab.dim) ;	
    
    huntm(xtab, x, i) ;
    
    int i1 = i - 1 ;
    int i2 = i + 1 ;
    
    double y0=ytab(i1) ;
    double y1=ytab(i)  ;
    double y2=ytab(i2) ;
    
    double x0=xtab(i1) ;
    double x1=xtab(i)  ;
    double x2=xtab(i2) ;
    double x01=x0-x1 ;
    double x02=x0-x2 ;
    double x12=x1-x2 ;
    
    y = y0*(x-x1)*(x-x2)/(x01*x02) - y1*(x-x0)*(x-x2)/(x01*x12) + y2*(x-x0)*(x-x1)/(x02*x12) ;
    
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

    dy =   6. * ( ytab(i) / dx * ( u2 - u )
		  - ytab(i1) / dx * ( u2 - u ) )
      + dytab(i) * ( 3.*u2 - 4.*u + 1. )
      + dytab(i1) * ( 3.*u2 - 2.*u ) ;
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
			const Tbl& dfdxtab, const Tbl& dfdytab, const Tbl& 
			d2fdxdytab, double x, double y, double& f, double& 
			dfdx, double& dfdy) {

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
			     const Tbl& dfdxtab, const Tbl& dfdytab, double x, 
			     double y, double& f, double& dfdx, double& dfdy) {

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

  //--------------------------------------------------------------------
  // Quintic Hermite interpolation using data from the second derivative
  //--------------------------------------------------------------------
  void interpol_herm_2nd_der(const Tbl& xtab, const Tbl& ytab, const Tbl& dytab,
			     const Tbl& d2ytab, double x, int& i, double& y, 
			     double& dy, double& d2y) {
    
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
    
    d2y = 60.*(ytab(i)/(dx*dx) * (-2.*u3 + 3*u2 - u)
              +ytab(i1)/(dx*dx) *(-2.*v3 + 3*v2 - v))
      + 12.*dytab(i)/dx  *(-5.*u3 + 8.*u2 - 3.*u)
      - 12.*dytab(i1)/dx *(-5.*v3 + 8.*v2 - 3.*v)
      +    d2ytab(i)  *(-10.*u3 + 18.*u2 - 9.*u + 1.)
      +    d2ytab(i1) *(-10.*v3 + 18.*v2 - 9.*v + 1.) ;
  }
  
  
  ///------------------------------------------------------------
  /// Extraction of a column of a 2D Tbl 
  ///------------------------------------------------------------
  Tbl extract_column(const Tbl& xytab, const Tbl& ytab, double yy){
    assert(xytab.get_ndim() == 2) ;
    assert(ytab.get_ndim() == 1) ;
    int i_low = ytab.get_taille()/2;
    huntm(ytab,yy,i_low) ;
    int n_h = xytab.get_dim(0) ;
    Tbl resu(n_h) ; resu.set_etat_qcq() ;
    for (int i=0 ; i<n_h ; i++){
      resu.set(i) = xytab(i_low,i) ;
    }
    return resu;
  }
  
  ///------------------------------------------------------------
  /// Extraction of a row of a 2D Tbl 
  ///------------------------------------------------------------
  Tbl extract_row(const Tbl& xytab, const Tbl& xtab, double xx){
    assert(xytab.get_ndim() == 2) ;
    assert(xtab.get_ndim() == 1) ;
    int i_low = xtab.get_taille()/2;
    huntm(xtab,xx,i_low) ;
    int n_h = xytab.get_dim(1) ;
    Tbl resu(n_h) ; resu.set_etat_qcq() ;
    for (int i=0 ; i<n_h ; i++){
      resu.set(i) = xytab(i,i_low) ;
    }
    return resu;
  }
  
  
} // End of namespace Lorene

