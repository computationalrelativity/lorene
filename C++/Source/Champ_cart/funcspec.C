/*
 * Methods of class FuncSpec.
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
  
  FuncSpec::FuncSpec(int nbx, int nby, int nbz):
    nx(nbx), ny(nby), nz(nbz), xmin(-1.), xmax(1.), ymin(-1.), ymax(1.),
    zmin(-1.), zmax(1.), xx(nbx), yy(nby), zz(nbz), values(nbx, nby, nbz),
    coefs(nbx, nby, nbz), coefs_up_to_date(false), values_up_to_date(false),
    p_dfdx(0x0), p_dfdy(0x0), p_dfdz(0x0)
  {
    if ( (nx<=1) || (ny<1) || (nz<1) ) {
      throw(std::invalid_argument("Too small size in FuncSpec constructor")) ;
    }
  }
  
  FuncSpec::FuncSpec(const TabSpec& ts):
    nx(ts.get_sizex()), ny(ts.get_sizey()), nz(ts.get_sizez()), xmin(-1.),
    xmax(1.), ymin(-1.), ymax(1.), zmin(-1.), zmax(1.), xx(nx), yy(ny), zz(nz),
    values(ts), coefs(nx, ny, nz), coefs_up_to_date(false), values_up_to_date(true),
    p_dfdx(0x0), p_dfdy(0x0), p_dfdz(0x0)
  {
    if ( (nx<=1) || (ny<1) || (nz<1) ) {
      throw(std::invalid_argument("Too small size in FuncSpec constructor")) ;
    }
  }

  FuncSpec::FuncSpec(const string& coord_name, const string& field)
    :xx(1), yy(1), zz(1), values(1), coefs(1), coefs_up_to_date(false),
     values_up_to_date(false), p_dfdx(0x0), p_dfdy(0x0), p_dfdz(0x0)
  {
    string tmp_str = coord_name + "x.dat" ;
    try {
      TabSpec tmp_tab(tmp_str) ;
      nx = tmp_tab.get_sizex() ;
      if (nx > 0) {
	xx.resize(nx, 1, 1) ;
	xx = tmp_tab ;
      }
      else
	throw std::ios::failure( "Bad value for nx!" ) ;
    }
    catch (const ifstream::failure& e) {
      cerr << "Cannot read x-coordinate in file " << tmp_str ;
      cerr << " (FuncSpec constructor)" << endl ;
      throw ;
    }
    xmin = xx(0) ;
    xmax = xx(nx-1) ;
  
    tmp_str = coord_name + "y.dat" ;
    try {
      TabSpec tmp_tab(tmp_str) ;
      ny = tmp_tab.get_sizex() ;
      if (ny > 0) {
	yy.resize(ny, 1, 1) ;
	yy = tmp_tab ;
      }
      else
	throw std::ios::failure( "Bad value for ny!" ) ;
    }
    catch (const ifstream::failure& e) {
      cerr << "Cannot read y-coordinate in file " << tmp_str ;
      cerr << " (FuncSpec constructor)" << endl ;
      throw ;
    }
    ymin = yy(0) ;
    ymax = yy(ny-1) ;
  
    tmp_str = coord_name + "z.dat" ;
    try {
      TabSpec tmp_tab(tmp_str) ;
      nz = tmp_tab.get_sizex() ;
      if (nz > 0) {
	zz.resize(nz, 1, 1) ;
	zz = tmp_tab ;
      }
      else
	throw std::ios::failure( "Bad value for nz!" ) ;
    }
    catch (const ifstream::failure& e) {
      cerr << "Cannot read z-coordinate in file " << tmp_str ;
      cerr << " (FuncSpec constructor)" << endl ;
      throw ;
    }
    zmin = zz(0) ;
    zmax = zz(nz-1) ;

    values.resize(nx, ny, nz) ;
    try {
      TabSpec tmp_tab(field) ;
      values = tmp_tab ;
    }
    catch (const ifstream::failure& e) {
      cerr << "Cannot read field values in file " << tmp_str ;
      cerr << " (FuncSpec constructor)" << endl ;
      throw ;
    }
    values_up_to_date = true ;
  
    coefs.resize(nx, ny, nz) ;
    coefs_up_to_date = false ;

  }

  FuncSpec::FuncSpec(const FuncSpec& fs):
    nx(fs.nx), ny(fs.ny), nz(fs.nz), xmin(fs.xmin), xmax(fs.xmax), ymin(fs.ymin),
    ymax(fs.ymax), zmin(fs.zmin), zmax(fs.zmax), xx(fs.xx), yy(fs.yy), zz(fs.zz),
    values(fs.values), coefs(fs.coefs), coefs_up_to_date(fs.coefs_up_to_date),
    values_up_to_date(fs.values_up_to_date), p_dfdx(0x0), p_dfdy(0x0), p_dfdz(0x0)
  {}

  void FuncSpec::set_grids(double x1, double x2, double y1, double y2, double z1,
			   double z2) {
    xmin = x1 ;
    xmax = x2 ;
    ymin = y1 ;
    ymax = y2 ;
    zmin = z1 ;
    zmax = z2 ;
    for (int i=0; i<nx; i++) {
      xx.set(i) = 0.5*( (xmin - xmax)*std::cos(i*M_PI/double(nx-1)) + xmin + xmax) ;
    }
    if (ny > 1)
      for (int i=0; i<ny; i++) {
	yy.set(i) = 0.5*( (ymin - ymax)*std::cos(i*M_PI/double(ny-1)) + ymin + ymax) ;
      }
    else 
      yy.set(0) = 0. ;
    if (nz > 1)
      for (int i=0; i<nz; i++) {
	zz.set(i) = 0.5*( (zmin - zmax)*std::cos(i*M_PI/double(nz-1)) + zmin + zmax) ;
      }
    else 
      zz.set(0) = 0. ;
    coefs_up_to_date = false ;
    del_deriv() ;
  }

  FuncSpec::~FuncSpec()
  {
    del_deriv() ;
  }

  void FuncSpec::operator=(const FuncSpec& fs)
  {
    nx = fs.nx ;
    ny = fs.ny ;
    nz = fs.nz ;
    xmin = fs.xmin ; xmax = fs.xmax ;
    ymin = fs.ymin ; ymax = fs.ymax ;
    zmin = fs.zmin ; zmax = fs.zmax ;

    xx.resize(nx, 1, 1) ;
    yy.resize(ny, 1, 1) ;
    zz.resize(nz, 1, 1) ;

    xx = fs.xx ;
    yy = fs.yy ;
    zz = fs.zz ;
  
    values.resize(nx, ny, nz) ;
    coefs.resize(nx, ny, nz) ;

    values = fs.values ;
    coefs = fs.coefs ;

    coefs_up_to_date = fs.coefs_up_to_date ; 
    values_up_to_date = fs.values_up_to_date ;
    del_deriv() ;
  }

  void FuncSpec::operator=(const TabSpec& ts)
  {
    bool new_grid = false ;
    if (nx != ts.get_sizex()) {
      nx = ts.get_sizex() ;
      xx.resize(nx) ;
      new_grid = true ;
    }
    if (ny != ts.get_sizey()) {
      ny = ts.get_sizey() ;
      yy.resize(ny) ;
      new_grid = true ;
    }
    if (nz != ts.get_sizez()) {
      nz = ts.get_sizez() ;
      zz.resize(nz) ;
      new_grid = true ;
    }

    if (new_grid) set_grids(xmin, xmax, ymin, ymax, zmin, zmax) ;

    values.resize(nx, ny, nz) ;
    values = ts ;
    values_up_to_date = true ;

    coefs.resize(nx, ny, nz) ;
    coefs_up_to_date = false ;
    del_deriv() ;
  }

  void FuncSpec::operator=(double x)
  {
    values = x ;
    coefs_up_to_date = false ;
    values_up_to_date = true ;
    del_deriv() ;
  }

  void FuncSpec::set_coefs(const TabSpec& ts)
  {
    bool new_grid = false ;
    if (nx != ts.get_sizex()) {
      nx = ts.get_sizex() ;
      xx.resize(nx) ;
      new_grid = true ;
    }
    if (ny != ts.get_sizey()) {
      ny = ts.get_sizey() ;
      yy.resize(ny) ;
      new_grid = true ;
    }
    if (nz != ts.get_sizez()) {
      nz = ts.get_sizez() ;
      zz.resize(nz) ;
      new_grid = true ;
    }

    if (new_grid) set_grids(xmin, xmax, ymin, ymax, zmin, zmax) ;

    coefs.resize(nx, ny, nz) ;
    coefs = ts ;
    coefs_up_to_date = true ;

    values.resize(nx, ny, nz) ;
    values_up_to_date = false ;
    del_deriv() ;
  }

  double FuncSpec::compute_in_xyz(double x, double y, double z) const {

#ifndef NDEBUG
    if (!coefs_up_to_date)
      throw(std::runtime_error("Coefficients not up to date in FuncSpec::compute_in_xyz")) ;

    if ( (x<xmin) || (x>xmax) )
      throw(std::range_error("x variable out of range in FuncSpec::compute_in_xyz")) ;
    if ( (y<ymin) || (y>ymax) )
      throw(std::range_error("y variable out of range in FuncSpec::compute_in_xyz")) ;
    if ( (ny == 1) && (y != 0.) )
      throw(std::range_error("y must be zero in FuncSpec::compute_in_xyz if only one gridpoint is set")) ;
    if ( (z<zmin) || (z>zmax) )
      throw(std::range_error("z variable out of range in FuncSpec::compute_in_xyz")) ;
    if ( (nz == 1) && (z != 0.) )
      throw(std::range_error("z must be zero in FuncSpec::compute_in_xyz if only one gridpoint is set")) ;
#endif

    double xi = (2*x - xmax - xmin) / (xmax - xmin) ;
    double yi = (2*y - ymax - ymin) / (ymax - ymin) ;
    double zi = (2*z - zmax - zmin) / (zmax - zmin) ;

    double resu = 0. ;
    const double* cf = coefs.tableau ;

    // ====================================
    // Special cases: 1D or 2D with nz = 1
    // ====================================
    if (nz == 1) {
      if (ny == 1) { // 1D case, only x-variable
	double tnm1x = 1. ;
	resu += (*cf)*tnm1x ;
	cf++ ;
	double tnx = xi ;
	resu += (*cf)*tnx ;
	cf++ ;
	for (int i=2; i<nx; i++) {
	  double newx = 2*xi*tnx - tnm1x ;
	  tnm1x = tnx ;
	  tnx = newx ;
	  resu += (*cf)*tnx ;
	  cf++ ;
	}
	return resu ;
      }// end of 1D case

      // 2D case, x, y variables
      double tnm1x = 1. ;
      double tnm1y = 1. ;
      resu += (*cf)*tnm1y*tnm1x ;
      cf++ ;
      double tny = yi ;
      resu += (*cf)*tny*tnm1x ;
      cf++ ;
      for (int j=2; j<ny; j++) {
	double newy = 2*yi*tny - tnm1y ;
	tnm1y = tny ;
	tny = newy ;
	resu += (*cf)*tny*tnm1x ;
	cf++ ;
      }
      double tnx = xi ;
      tnm1y = 1. ;
      resu += (*cf)*tnm1y*tnx ;
      cf++ ;
      tny = yi ;
      resu += (*cf)*tny*tnx ;
      cf++ ;
      for (int j=2; j<ny; j++) {
	double newy = 2*yi*tny - tnm1y ;
	tnm1y = tny ;
	tny = newy ;
	resu += (*cf)*tny*tnx ;
	cf++ ;
      }
      for (int i=2; i<nx; i++) {
	double newx = 2*xi*tnx - tnm1x ;
	tnm1x = tnx ;
	tnx = newx ;
	tnm1y = 1. ;
	resu += (*cf)*tnm1y*tnx ;
	cf++ ;
	tny = yi ;
	resu += (*cf)*tny*tnx ;
	cf++ ;
	for (int j=2; j<ny; j++) {
	  double newy = 2*yi*tny - tnm1y ;
	  tnm1y = tny ;
	  tny = newy ;
	  resu += (*cf)*tny*tnx ;
	  cf++ ;
	}
      }
      return resu ;
    } // end of 2D case, with nz = 1

    // 2D case, x, z variables
    if (ny == 1) {
      double tnm1x = 1. ;
      double tnm1z = 1. ;
      resu += (*cf)*tnm1z*tnm1x ;
      cf++ ;
      double tnz = zi ;
      resu += (*cf)*tnz*tnm1x ;
      cf++ ;
      for (int k=2; k<nz; k++) {
	double newz = 2*zi*tnz - tnm1z ;
	tnm1z = tnz ;
	tnz = newz ;
	resu += (*cf)*tnz*tnm1x ;
	cf++ ;
      }
      double tnx = xi ;
      tnm1z = 1. ;
      resu += (*cf)*tnm1z*tnx ;
      cf++ ;
      tnz = zi ;
      resu += (*cf)*tnz*tnx ;
      cf++ ;
      for (int k=2; k<nz; k++) {
	double newz = 2*zi*tnz - tnm1z ;
	tnm1z = tnz ;
	tnz = newz ;
	resu += (*cf)*tnz*tnx ;
	cf++ ;
      }
      for (int i=2; i<nx; i++) {
	double newx = 2*xi*tnx - tnm1x ;
	tnm1x = tnx ;
	tnx = newx ;
	tnm1z = 1. ;
	resu += (*cf)*tnm1z*tnx ;
	cf++ ;
	tnz = zi ;
	resu += (*cf)*tnz*tnx ;
	cf++ ;
	for (int k=2; k<nz; k++) {
	  double newz = 2*zi*tnz - tnm1z ;
	  tnm1z = tnz ;
	  tnz = newz ;
	  resu += (*cf)*tnz*tnx ;
	  cf++ ;
	}
      }    
      return resu ;
    }

    // ====================================
    // General case: 3D => 3 loops
    // ====================================
    // Initialization in x, y and z
    //--------------------------------
    double tnm1x = 1. ;
    double tnm1y = 1. ;
    double tnm1z = 1. ;
    resu += (*cf)*tnm1z*tnm1y*tnm1x ;
    cf++ ;
    double tnz = zi ; // second step in z
    resu += (*cf)*tnz*tnm1y*tnm1x ;
    cf++ ;
    for (int k=2; k<nz; k++) { // loop in z (first step in x and y)
      double newz = 2*zi*tnz - tnm1z ;
      tnm1z = tnz ;
      tnz = newz ;
      resu += (*cf)*tnz*tnm1y*tnm1x ;
      cf++ ;
    }
    double tny = yi ; //second step in y
    tnm1z = 1. ; // first step in z
    resu += (*cf)*tnm1z*tny*tnm1x ;
    cf++ ;
    tnz = zi ; // second step in z 
    resu += (*cf)*tnz*tny*tnm1x ;
    cf++ ;
    for (int k=2; k<nz; k++) { // loop in z (first step in x, 2nd step in y)
      double newz = 2*zi*tnz - tnm1z ;
      tnm1z = tnz ;
      tnz = newz ;
      resu += (*cf)*tnz*tny*tnm1x ;
      cf++ ;
    }
    for (int j=2; j<ny; j++) { // loop in y (first step in x)
      double newy = 2*yi*tny - tnm1y ;
      tnm1y = tny ;
      tny = newy ;
      tnm1z = 1. ;
      resu += (*cf)*tnm1z*tny*tnm1x ;
      cf++ ;
      tnz = zi ;
      resu += (*cf)*tnz*tny*tnm1x ;
      cf++ ;
      for (int k=2; k<nz; k++) {
	double newz = 2*zi*tnz - tnm1z ;
	tnm1z = tnz ;
	tnz = newz ;
	resu += (*cf)*tnz*tny*tnm1x ;
	cf++ ;
      }
    }

    // Initializing: second step in x
    //--------------------------------
    double tnx = xi ;
    tnm1y = 1. ; // first step in y
    tnm1z = 1. ; // first step in z
    resu += (*cf)*tnm1z*tnm1y*tnx ;
    cf++ ;
    tnz = zi ; // second step in z
    resu += (*cf)*tnz*tnm1y*tnx ;
    cf++ ;
    for (int k=2; k<nz; k++) { // loop in z (second step in x, first step in y)
      double newz = 2*zi*tnz - tnm1z ;
      tnm1z = tnz ;
      tnz = newz ;
      resu += (*cf)*tnz*tnm1y*tnx ;
      cf++ ;
    }
    tny = yi ; // second step in y
    tnm1z = 1. ; // first step in z 
    resu += (*cf)*tnm1z*tny*tnx ;
    cf++ ;
    tnz = zi ; //second step in z
    resu += (*cf)*tnz*tny*tnx ;
    cf++ ;
    for (int k=2; k<nz; k++) { // loop in z (second step in x, 2nd in y)
      double newz = 2*zi*tnz - tnm1z ;
      tnm1z = tnz ;
      tnz = newz ;
      resu += (*cf)*tnz*tny*tnx ;
      cf++ ;
    }
    for (int j=2; j<ny; j++) { // loop in y (second step in x)
      double newy = 2*yi*tny - tnm1y ;
      tnm1y = tny ;
      tny = newy ;
      tnm1z = 1. ;
      resu += (*cf)*tnm1z*tny*tnx ;
      cf++ ;
      tnz = zi ;
      resu += (*cf)*tnz*tny*tnx ;
      cf++ ;
      for (int k=2; k<nz; k++) {
	double newz = 2*zi*tnz - tnm1z ;
	tnm1z = tnz ;
	tnz = newz ;
	resu += (*cf)*tnz*tny*tnx ;
	cf++ ;
      }
    }

    // Loop in x
    //-----------
    for (int i=2; i<nx; i++) {
      double newx = 2*xi*tnx -tnm1x ;
      tnm1x = tnx ;
      tnx = newx ;
      tnm1y = 1. ;
      tnm1z = 1. ;
      resu += (*cf)*tnm1z*tnm1y*tnx ;
      cf++ ;
      tnz = zi ;
      resu += (*cf)*tnz*tnm1y*tnx ;
      cf++ ;
      for (int k=2; k<nz; k++) {
	double newz = 2*zi*tnz - tnm1z ;
	tnm1z = tnz ;
	tnz = newz ;
	resu += (*cf)*tnz*tnm1y*tnx ;
	cf++ ;
      }
      tny = yi ;
      tnm1z = 1. ;
      resu += (*cf)*tnm1z*tny*tnx ;
      cf++ ;
      tnz = zi ;
      resu += (*cf)*tnz*tny*tnx ;
      cf++ ;
      for (int k=2; k<nz; k++) {
	double newz = 2*zi*tnz - tnm1z ;
	tnm1z = tnz ;
	tnz = newz ;
	resu += (*cf)*tnz*tny*tnx ;
	cf++ ;
      }
      for (int j=2; j<ny; j++) {
	double newy = 2*yi*tny - tnm1y ;
	tnm1y = tny ;
	tny = newy ;
	tnm1z = 1. ;
	resu += (*cf)*tnm1z*tny*tnx ;
	cf++ ;
	tnz = zi ;
	resu += (*cf)*tnz*tny*tnx ;
	cf++ ;
	for (int k=2; k<nz; k++) {
	  double newz = 2*zi*tnz - tnm1z ;
	  tnm1z = tnz ;
	  tnz = newz ;
	  resu += (*cf)*tnz*tny*tnx ;
	  cf++ ;
	}
      }
    }
    return resu ;

  }

  TabSpec FuncSpec::grid_x() const {

    TabSpec resu(nx, ny, nz) ;
    for (int i=0; i<nx; i++)
      for (int j=0; j<ny; j++)
	for (int k=0; k<nz; k++)
	  resu.set(i,j,k) = xx(i) ;
    return resu ;
  }

  TabSpec FuncSpec::grid_y() const {

    TabSpec resu(nx, ny, nz) ;
    for (int i=0; i<nx; i++)
      for (int j=0; j<ny; j++)
	for (int k=0; k<nz; k++)
	  resu.set(i,j,k) = yy(j) ;
    return resu ;
  }

  TabSpec FuncSpec::grid_z() const {

    TabSpec resu(nx, ny, nz) ;
    for (int i=0; i<nx; i++)
      for (int j=0; j<ny; j++)
	for (int k=0; k<nz; k++)
	  resu.set(i,j,k) = zz(k) ;
    return resu ;
  }


  ostream& operator<<(ostream& o, const FuncSpec& f) {

    o << "3D spectral function on " << f.nx << "x" << f.ny << "x" << f.nz
      << " points" << endl ;
    o << "x-interval : [" << f.xmin << ", " << f.xmax << "]" << endl ; 
    o << "y-interval : [" << f.ymin << ", " << f.ymax << "]" << endl ; 
    o << "z-interval : [" << f.zmin << ", " << f.zmax << "]" << endl ;

    if (!f.values_up_to_date) {
      o << "Values not known." << endl ;
    }
    else {
      o << "Values at grid points: " << endl ;
      o << f.values ;
    }
    if (!f.coefs_up_to_date)
      o << "Coefficients not known." << endl ;
    else {
      o << "Coefficients:" << endl ;
      o << f.coefs ;
    }
    return o ;
  }


  void FuncSpec::write_grids(const string& str) {

    string tmp_str = str + "x.dat" ;
    xx.write_file(tmp_str) ;

    tmp_str = str + "y.dat" ;
    yy.write_file(tmp_str) ;

    tmp_str = str + "z.dat" ;
    zz.write_file(tmp_str) ;
  }

  
  void FuncSpec::write_values(const string& str) {
    
    if (!values_up_to_date) {
      if (!coefs_up_to_date) {
	throw(std::runtime_error("Both values & coefficients not up to date in FuncSpec::write_values")) ;
      }
      compute_values() ;
      values_up_to_date = true ;
    }
    
    values.write_file(str) ;
    
  }

} // end namespace Lorene 


  
