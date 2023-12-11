 /*
 *  Definition of Lorene class FuncSpec
 *
 */

/*
 *   Copyright (c) 2019 Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

#ifndef __FUNCSPEC__
#define __FUNCSPEC__
#include "tabspec.h"

namespace Lorene {

  /** 
   * Class for representing functions of 3 variables, supposed to be Cartesian
   * coordinates (x, y, z). \ingroup (map_cart)
   */
  class FuncSpec {
  
    //Data
    //----
  protected:
    int nx, ny, nz ;
    double xmin, xmax, ymin, ymax, zmin, zmax ;
    TabSpec xx, yy, zz ;
    mutable TabSpec values, coefs ;
    mutable bool coefs_up_to_date ;
    mutable bool values_up_to_date ;

    mutable FuncSpec* p_dfdx ;
    mutable FuncSpec* p_dfdy ;
    mutable FuncSpec* p_dfdz ;
  
    //Constructors - Destructor
    //-------------------------
  public:
    explicit FuncSpec(int=2, int=2, int=2) ;
    explicit FuncSpec(const TabSpec&) ;
    /** Constructor from files.'coord_name' is a prefix for files containing 
     *  coordinate data (coord_namex.dat, coord_namey.dat, coord_namez.dat);
     *  'field', is the complete name of the fiel containing the field data.*/
    explicit FuncSpec(const string& coord_name, const string& field) ; 
    FuncSpec(const FuncSpec&) ; ///< Copy constructor

    void set_grids(double, double, double, double, double, double) ;
  
    virtual ~FuncSpec() ; ///< Destructor

  protected:
    void del_deriv() const ;

    // Assignments
    //------------
  public:
    void operator=(const FuncSpec&) ;
    void operator=(const TabSpec&) ;
    void operator=(double) ;

    void set_coefs(const TabSpec&) ;

    // Computing
    //----------
    void compute_coefs() const ;
    double compute_in_xyz(double, double, double) const ;
    void compute_values() const ;
    FuncSpec get_partial_x() const ;
    FuncSpec get_partial_y() const ;
    FuncSpec get_partial_z() const ;

    FuncSpec primitive_x() const ;

    // Saving into files
    //-------------------
    void write_grids(const string&) ;
    void write_values(const string&) ;

    // Access to data
    //---------------
    TabSpec grid_x() const ;
    TabSpec grid_y() const ;
    TabSpec grid_z() const ;

    // Interpolate from a 3D TabSpec
    //------------------------------
    void interpolate_from_Tab(const TabSpec& values, const TabSpec& x_coord,
			      const TabSpec& y_coord, const TabSpec& z_coord) ;

    friend ostream& operator<<(ostream&, const FuncSpec&) ;
    friend FuncSpec operator+(const FuncSpec&, const FuncSpec&) ;
    friend FuncSpec operator-(const FuncSpec&) ;
    friend FuncSpec operator+(const FuncSpec&, double) ;
    friend FuncSpec operator+(double, const FuncSpec&) ;
    friend FuncSpec operator-(const FuncSpec&, const FuncSpec&) ;
    friend FuncSpec operator-(const FuncSpec&, double) ;
    friend FuncSpec operator-(double, const FuncSpec&) ;
    friend FuncSpec operator*(const FuncSpec&, const FuncSpec&) ;
    friend FuncSpec operator*(const FuncSpec&, double) ;
    friend FuncSpec operator*(double, const FuncSpec&) ;
    friend FuncSpec operator/(const FuncSpec&, const FuncSpec&) ;
    friend FuncSpec operator/(const FuncSpec&, double) ;
    friend FuncSpec operator/(double, const FuncSpec&) ;

    friend FuncSpec sin(const FuncSpec&) ;
    friend FuncSpec cos(const FuncSpec&) ;
    friend FuncSpec tan(const FuncSpec&) ;
    friend FuncSpec exp(const FuncSpec&) ;
    friend FuncSpec log(const FuncSpec&) ;
    friend FuncSpec sqrt(const FuncSpec&) ;
    friend FuncSpec pow(const FuncSpec&, double) ;
    friend FuncSpec abs(const FuncSpec&) ;

    friend double max(const FuncSpec&) ;

  };

  ostream& operator<<(ostream&, const FuncSpec&) ;

  FuncSpec operator-(const FuncSpec&) ;
  FuncSpec operator+(const FuncSpec&, const FuncSpec&) ;
  FuncSpec operator+(const FuncSpec&, double) ;
  FuncSpec operator+(double, const FuncSpec&) ;
  FuncSpec operator-(const FuncSpec&, const FuncSpec&) ;
  FuncSpec operator-(const FuncSpec&, double) ;
  FuncSpec operator-(double, const FuncSpec&) ;
  FuncSpec operator*(const FuncSpec&, const FuncSpec&) ;
  FuncSpec operator*(const FuncSpec&, double) ;
  FuncSpec operator*(double, const FuncSpec&) ;
  FuncSpec operator/(const FuncSpec&, const FuncSpec&) ;
  FuncSpec operator/(const FuncSpec&, double) ;
  FuncSpec operator/(double, const FuncSpec&) ;

  FuncSpec sin(const FuncSpec&) ;
  FuncSpec cos(const FuncSpec&) ;
  FuncSpec tan(const FuncSpec&) ;
  FuncSpec exp(const FuncSpec&) ;
  FuncSpec log(const FuncSpec&) ;
  FuncSpec sqrt(const FuncSpec&) ;
  FuncSpec pow(const FuncSpec&, double) ;
  FuncSpec abs(const FuncSpec&) ;
  double max(const FuncSpec&) ;
}
#endif
