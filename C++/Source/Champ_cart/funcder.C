/*
 * Methods of class FuncSpec to compute partial derivatives.
 *
 *  (see file funcspec.h for documentation)
 *
 */

/*
 *   Copyright (c) 2019 Jerome Novak
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


#include "funcspec.h"

namespace Lorene {

  void FuncSpec::del_deriv() const {

    if (p_dfdx != 0x0) {
      delete p_dfdx ;
      p_dfdx = 0x0 ;
    }

    if (p_dfdy != 0x0) {
      delete p_dfdy ;
      p_dfdy = 0x0 ;
    }

    if (p_dfdz != 0x0) {
      delete p_dfdz ;
      p_dfdz = 0x0 ;
    }

  }

  FuncSpec FuncSpec::get_partial_x() const {

    if (p_dfdx == 0x0) {
    
      if (nx >= 5) {
	if (!coefs_up_to_date) compute_coefs() ;
	p_dfdx = new FuncSpec(nx, ny, nz) ;
	p_dfdx->set_grids(xmin, xmax, ymin, ymax, zmin, zmax) ;
	p_dfdx->coefs = 0 ;
      
	int ntrans = ny*nz ;
	int nskip = ny*nz ;
	//      int next = 1 ;
	double* xi = coefs.tableau ;
	double* xo = (p_dfdx->coefs).tableau ;

	for (int n=0; n<ntrans; n++) {
	  //Starting index
	  int ind = n ;
	  double* xci = xi + ind ;
	  double* xco = xo + ind ;
	  double som ;
	  xco[(nx-1)*nskip] = 0 ;
	  som = 2*(nx-1) * xci[(nx-1)*nskip] ;
	  xco[(nx-2)*nskip] = som ;
	  for (int i= nx -4; i>=0; i -=2) {
	    som += 2*(i+1) * xci[(i+1)*nskip] ;
	    xco[i*nskip] = som ;
	  }
	  som = 2*(nx-2) * xci[(nx-2)*nskip] ;
	  xco[(nx-3)*nskip] = som ;
	  for (int i=nx-5; i >=0; i -= 2) {
	    som += 2*(i+1) * xci[(i+1)*nskip] ;
	    xco[i*nskip] = som ;
	  }
	  xco[0] *= 0.5 ;
	}
	p_dfdx->coefs = (2./(xmax - xmin)) * p_dfdx->coefs ;
      }
      else {
	p_dfdx = new FuncSpec(nx, ny, nz) ;
	p_dfdx->set_grids(xmin, xmax, ymin, ymax, zmin, zmax) ;
	p_dfdx->coefs = 0 ;
      }
      p_dfdx->coefs_up_to_date = true ;
      p_dfdx->values_up_to_date = false ;
    }
    return *p_dfdx ;
  }


  FuncSpec FuncSpec::get_partial_y() const {

    if (p_dfdy == 0x0) {
    
      if (ny >= 5) {
	if (!coefs_up_to_date) compute_coefs() ;
	p_dfdy = new FuncSpec(nx, ny, nz) ;
	p_dfdy->set_grids(xmin, xmax, ymin, ymax, zmin, zmax) ;
	p_dfdy->coefs = 0 ;
      
	int ntrans = nx*nz ;
	int nskip = nz ;
	int next = ny*nz ;
	int ntot = nx*ny*nz ;
	double* xi = coefs.tableau ;
	double* xo = p_dfdy->coefs.tableau ;

	for (int n=0; n<ntrans; n++) {
	  //Starting index
	  int ind = (n*next) % ntot + (n*next/ntot) ;
	  double* xci = xi + ind ;
	  double* xco = xo + ind ;
	  double som ;
	  xco[(ny-1)*nskip] = 0 ;
	  som = 2*(ny-1) * xci[(ny-1)*nskip] ;
	  xco[(ny-2)*nskip] = som ;
	  for (int i= ny -4; i>=0; i -=2) {
	    som += 2*(i+1) * xci[(i+1)*nskip] ;
	    xco[i*nskip] = som ;
	  }
	  som = 2*(ny-2) * xci[(ny-2)*nskip] ;
	  xco[(ny-3)*nskip] = som ;
	  for (int i=ny-5; i >=0; i -= 2) {
	    som += 2*(i+1) * xci[(i+1)*nskip] ;
	    xco[i*nskip] = som ;
	  }
	  xco[0] *= 0.5 ;
	}
	p_dfdy->coefs = (2./(ymax - ymin)) * p_dfdy->coefs ;
      }
      else {
	p_dfdy = new FuncSpec(nx, ny, nz) ;
	p_dfdy->set_grids(xmin, xmax, ymin, ymax, zmin, zmax) ;
	p_dfdy->coefs = 0 ;
      }
      p_dfdy->coefs_up_to_date = true ;
      p_dfdy->values_up_to_date = false ;
    }
    return *p_dfdy ;
  }

  FuncSpec FuncSpec::get_partial_z() const {

    if (p_dfdz == 0x0) {
    
      if (nz >= 5) {
	if (!coefs_up_to_date) compute_coefs() ;
	p_dfdz = new FuncSpec(nx, ny, nz) ;
	p_dfdz->set_grids(xmin, xmax, ymin, ymax, zmin, zmax) ;
	p_dfdz->coefs = 0 ;

	int ntrans = nx*ny ;
	int next = nz ;
	int ntot = nx*ny*nz ;
	double* xi = coefs.tableau ;
	double* xo = p_dfdz->coefs.tableau ;

	for (int n=0; n<ntrans; n++) {
	  //Starting index
	  int ind = (n*next) % ntot + (n*next/ntot) ;
	  double* xci = xi + ind ;
	  double* xco = xo + ind ;
	  double som ;
	  xco[nz-1] = 0 ;
	  som = 2*(nz-1) * xci[nz-1] ;
	  xco[nz-2] = som ;
	  for (int i= nz -4; i>=0; i -=2) {
	    som += 2*(i+1) * xci[i+1] ;
	    xco[i] = som ;
	  }
	  som = 2*(nz-2) * xci[nz-2] ;
	  xco[nz-3] = som ;
	  for (int i=nz-5; i >=0; i -= 2) {
	    som += 2*(i+1) * xci[i+1] ;
	    xco[i] = som ;
	  }
	  xco[0] *= 0.5 ;
	}
	p_dfdz->coefs = (2./(zmax - zmin)) * p_dfdz->coefs ;
      }
      else {
	p_dfdz = new FuncSpec(nx, ny, nz) ;
	p_dfdz->set_grids(xmin, xmax, ymin, ymax, zmin, zmax) ;
	p_dfdz->coefs = 0 ;
      }
      p_dfdz->coefs_up_to_date = true ;
      p_dfdz->values_up_to_date = false ;
    }
    return *p_dfdz ;
  }


} // end namespace Lorene
