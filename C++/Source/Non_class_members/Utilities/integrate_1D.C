/*
 * Integration of f(x) in the interval [xx(0), xx(n-1)], with non-equally spaced 
 * n-size xx grid.
 *
 * The function f is approximated by piecewise parabolae, The integral of f
 * is set to 0 at xx(0).
 */

/*
 *   Copyright (c) 2015 Jerome Novak
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


char integrate_1D_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2015/01/09 15:28:52  j_novak
 * New integration function for general non-equally-spaced grids.
 *
 *
 */

// Headers Lorene
#include "tbl.h"

namespace Lorene {

Tbl integ1D(const Tbl& xx, const Tbl& ff) {

  Tbl resu(ff) ;
  if (ff.get_etat() != ETATZERO) {

    assert (xx.get_etat() == ETATQCQ) ;
    assert (ff.get_etat() == ETATQCQ) ;
    int nx = xx.get_taille() ;
    assert(nx > 2) ;
    assert (ff.get_taille() == nx) ;
    
    resu.set(0) = 0. ;
    double x0 = xx(0) ;
    double x1(0), x2(0), x3(0);
    double a1(0), a2(0), a3(0);
    double b1(0), b2(0), b3(0);
    double c1(0), c2(0), c3(0) ;

    for (int i=1; i<nx-1; i++) {
      x1 = xx(i-1) ;
      x2 = xx(i) ;
      x3 = xx(i+1) ;
      a1 = ff(i-1) / ( (x1 - x2)*(x1 - x3) ) ;
      a2 = ff(i) / ( (x2 - x1)*(x2 - x3) ) ;
      a3 = ff(i+1) / ( (x3 - x1)*(x3 - x2) ) ;
      b1 = a1 + a2 + a3 ;
      b2 = -(x2 + x3)*a1 - (x1 + x3)*a2 - (x1 + x2)*a3 ;
      b3 = x2*x3*a1 + x1*x3*a2 + x1*x2*a3 ;
      if (i==1) {
	c1 = b1 ;
	c2 = b2 ;
	c3 = b3 ;
      }
      else {
	c1 = 0.5*(b1 + c1) ;
	c2 = 0.5*(b2 + c2) ;
	c3 = 0.5*(b3 + c3) ;
      }
      resu.set(i) = resu(i-1) + c1*(x2*x2*x2 - x0*x0*x0)/3.
	+ 0.5*c2*(x2*x2 - x0*x0) + c3*(x2 - x0) ;
      c1 = b1 ;
      c2 = b2 ;
      c3 = b3 ;
      x0 = x2 ;
    }
    
    x2 = xx(nx-1) ;
    resu.set(nx-1) = resu(nx-2) + c1*(x2*x2*x2 - x0*x0*x0)/3.
	+ 0.5*c2*(x2*x2 - x0*x0) + c3*(x2 - x0) ;
      
  }
  return resu ;
}

} // End of namespace Lorene
