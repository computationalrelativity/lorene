/*
 * Methods of class FuncSpec for interpolation.
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


#include<cmath>
#include<fstream>
#include<iomanip>
#include "funcspec.h"

namespace Lorene {

  bool check_grid(const TabSpec& x_from, const TabSpec& y_from, const TabSpec& z_from, const TabSpec& x_to, const TabSpec& y_to, const TabSpec& z_to) {
  
    bool resu = true ;

    int nx_from = x_from.get_sizex() ;
    int nx_to = x_to.get_sizex() ;
    double xf_min = x_from(0) ;
    double xf_max = x_from(nx_from - 1) ;
    double xt_min = x_to(0) ;
    double xt_max = x_to(nx_to - 1) ;

    if (nx_to > 1) {
      resu = (resu && (xf_min <= xt_min) ) ;
      resu = (resu && (xf_max >= xt_max) ) ;
    }
  
    int ny_from = y_from.get_sizex() ;
    int ny_to = y_to.get_sizex() ;
    double yf_min = y_from(0) ;
    double yf_max = y_from(ny_from - 1) ;
    double yt_min = y_to(0) ;
    double yt_max = y_to(ny_to - 1) ;

    if (ny_to > 1) {
      resu = (resu && (yf_min <= yt_min) ) ;
      resu = (resu && (yf_max >= yt_max) ) ;
    }
    
    int nz_from = z_from.get_sizex() ;
    int nz_to = z_to.get_sizex() ;
    double zf_min = z_from(0) ;
    double zf_max = z_from(nz_from - 1) ;
    double zt_min = z_to(0) ;
    double zt_max = z_to(nz_to - 1) ;

    if (nz_to > 1) {
      resu = (resu && (zf_min <= zt_min) ) ;
      resu = (resu && (zf_max >= zt_max) ) ;
    }

    return resu ;
  }

  void interpol_1d(const double* t_from, double* x_from, int ninterp_f, int next_f,
		   int, double* t_to, double* x_to, int ninterp_t, int next_t,
		   int, int type_inter){

    switch(type_inter) {
    case 1: {
      int im1 = 0 ;
      int i0 = 1 ;
      for (int i=0; i<ninterp_t; i++) {
	while(x_from[i0] < x_to[i]) i0++ ;
	if (i0 >= ninterp_f)
	  throw(std::range_error("interpol_1d: destination interval not contained in source interval."));
	im1 = i0 - 1 ;
	t_to[i*next_t] = (t_from[i0*next_f]*(x_from[im1] - x_to[i]) +
			  t_from[im1*next_f]*(x_to[i] - x_from[i0])) /
	  (x_from[im1] - x_from[i0]) ;
      }
      break;
    }
    
    case 2: {
      int i1, i2, i3 ;
      double xr, x1, x2, x3, y1, y2, y3 ;
      i2 = 0 ;
      i3 = 1 ;
      for (int i=0; i<ninterp_t; i++) {
	xr = x_to[i] ;
	while(x_from[i3] < xr) i3++ ;
	if (i3 >= ninterp_f)
	  throw(std::range_error("interpol_1d: destination interval not contained in source interval."));
	if (i3 == 1) {
	  i1 = 0 ;
	  i2 = 1 ;
	  i3 = 2 ;
	}
	else {
	  i2 = i3 - 1 ;
	  i1 = i2 - 1 ;
	}
	x1 = x_from[i1] ;
	x2 = x_from[i2] ;
	x3 = x_from[i3] ;
	y1 = t_from[i1*next_f] ;
	y2 = t_from[i2*next_f] ;
	y3 = t_from[i3*next_f] ;
	double c = y1 ;
	double b = (y2 - y1) / (x2 - x1) ;
	double a = ( (y3 - y2)/(x3 - x2) - b ) / (x3 - x1) ;
	t_to[i*next_t] = c + (xr - x1)*( b + a*(xr - x2)) ;
      }
      break;
    }
    
    default: {
      throw(std::invalid_argument("Unknown type of interpolation for FuncSpec"));
    }
    }
  }

  void FuncSpec::interpolate_from_Tab(const TabSpec& val_from, const TabSpec& x_from,
				      const TabSpec& y_from, const TabSpec& z_from) {

    bool boundaries_ok = check_grid(x_from, y_from, z_from, xx, yy, zz) ;
    if (!boundaries_ok) {
      throw(std::range_error("interpolate_from_tab() : source grid is not larger than destination one")) ;
    }
  
    int nx_from = x_from.get_sizex() ;
    int ny_from = y_from.get_sizex() ;
    int nz_from = z_from.get_sizex() ;

    int one = 1 ;
    int type_inter = 2 ; // ## Move this choice to another place ?

    // First interpolation in x-direction
    TabSpec interm1(nx, ny_from, nz_from) ;
    //#pragma omp parallel for collapse(2)
    for (int i=0; i<ny_from; i++) {
      for (int j=0; j<nz_from; j++) {
	double* xci = &(val_from.tableau[i*nz_from + j]) ;
	double* xco = &(interm1.tableau[i*nz_from + j]) ;
	interpol_1d(xci, x_from.tableau, nx_from, ny_from*nz_from,
		    nx_from*ny_from*nz_from, xco, xx.tableau, nx, ny_from*nz_from,
		    nx*ny_from*nz_from, type_inter) ;
      }
    }

    // Second interpolation in y-direction
    TabSpec interm2(nx, ny, nz_from) ;
    //#pragma omp parallel for collapse(2)
    for (int i=0; i<nx; i++) {
      for (int j=0; j<nz_from; j++) {
	double* xci = &(interm1.tableau[i*nz_from*ny_from + j]) ;
	double* xco = &(interm2.tableau[i*nz_from*ny + j]) ;
	interpol_1d(xci, y_from.tableau, ny_from, nz_from, nx*ny_from*nz_from, xco,
		    yy.tableau, ny, nz_from, nx*ny*nz_from, type_inter) ;
      }
    }

    // Last interpolation in z-direction
    //#pragma omp parallel for collapse(2)
    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
	double* xci = &(interm2.tableau[(i*ny + j)*nz_from]) ;
	double* xco = &(values.tableau[(i*ny + j)*nz]) ;
	interpol_1d(xci, z_from.tableau, nz_from, one, nx*ny*nz_from, xco,
		    zz.tableau, nz, one, nx*ny*nz, type_inter) ;
      }
    }
    values_up_to_date = true ;
    coefs_up_to_date = false ;

  }


} // namespace Lorene
