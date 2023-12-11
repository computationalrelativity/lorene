/*
 * Mathematics for the class FuncSpec.
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
#include <cmath>

namespace Lorene {

  // Operators+
  //-------------
  FuncSpec operator+(const FuncSpec& t1, const FuncSpec& t2) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      if (t2.values_up_to_date) {
	resu = t1.values + t2.values ;
	if (t1.coefs_up_to_date && t2.coefs_up_to_date) {
	  resu.coefs = t1.coefs + t2.coefs ;
	  resu.coefs_up_to_date = true ;
	}
      }
      else {
	if (!t2.coefs_up_to_date)
	  throw invalid_argument("Ill-formed second argument in operator+(FuncSpec)") ;
	if (t1.coefs_up_to_date) {
	  resu.set_coefs(t1.coefs + t2.coefs) ;
	}
	else {
	  resu.set_coefs(t2.coefs) ;
	  resu.compute_values() ;
	  resu = t1.values + resu.values ;
	}
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator+(FuncSpec)") ;
      if (t2.coefs_up_to_date) {
	resu.set_coefs(t1.coefs + t2.coefs) ;
      }
      else {
	if (!t2.values_up_to_date)
	  throw invalid_argument("Ill-formed second argument in operator+(FuncSpec)") ;
	resu.set_coefs(t1.coefs) ;
	resu.compute_values() ;
	resu = resu.values + t2.values ;
      }
    }
    return resu ;

  }

  FuncSpec operator+(const FuncSpec& t1, double x) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      resu.values = t1.values + x ;
      resu.values_up_to_date = true ;
      if (t1.coefs_up_to_date) {
	resu.coefs = t1.coefs ;
	resu.coefs.set(0) += x ;
	resu.coefs_up_to_date = true ;
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator+(FuncSpec)") ;
      resu.set_coefs(t1.coefs) ;
      resu.coefs.set(0) += x ;
    }
    return resu ;
  }

  FuncSpec operator+(double x, const FuncSpec& t1) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      resu.values = t1.values + x ;
      resu.values_up_to_date = true ;
      if (t1.coefs_up_to_date) {
	resu.coefs = t1.coefs ;
	resu.coefs.set(0) += x ;
	resu.coefs_up_to_date = true ;
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator+(FuncSpec)") ;
      resu.set_coefs(t1.coefs) ;
      resu.coefs.set(0) += x ;
    }
    return resu ;
  }

  // Operators-
  //-------------

  FuncSpec operator-(const FuncSpec& t1) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      resu.values = -t1.values ;
      resu.values_up_to_date = true ;
      if (t1.coefs_up_to_date) {
	resu.coefs = -t1.coefs ;
	resu.coefs_up_to_date = true ;
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator-(FuncSpec)") ;
      resu.set_coefs(-t1.coefs) ;
    }
    return resu ;
  }

  FuncSpec operator-(const FuncSpec& t1, const FuncSpec& t2) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      if (t2.values_up_to_date) {
	resu = t1.values - t2.values ;
	if (t1.coefs_up_to_date && t2.coefs_up_to_date) {
	  resu.coefs = t1.coefs - t2.coefs ;
	  resu.coefs_up_to_date = true ;
	}
      }
      else {
	if (!t2.coefs_up_to_date)
	  throw invalid_argument("Ill-formed second argument in operator-(FuncSpec)") ;
	if (t1.coefs_up_to_date) {
	  resu.set_coefs(t1.coefs - t2.coefs) ;
	}
	else {
	  resu.set_coefs(t2.coefs) ;
	  resu.compute_values() ;
	  resu = t1.values - resu.values ;
	}
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator-(FuncSpec)") ;
      if (t2.coefs_up_to_date) {
	resu.set_coefs(t1.coefs - t2.coefs) ;
      }
      else {
	if (!t2.values_up_to_date)
	  throw invalid_argument("Ill-formed second argument in operator-(FuncSpec)") ;
	resu.set_coefs(t1.coefs) ;
	resu.compute_values() ;
	resu = resu.values - t2.values ;
      }
    }
    return resu ;

  }

  FuncSpec operator-(const FuncSpec& t1, double x) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      resu.values = t1.values - x ;
      resu.values_up_to_date = true ;
      if (t1.coefs_up_to_date) {
	resu.coefs = t1.coefs ;
	resu.coefs.set(0) -= x ;
	resu.coefs_up_to_date = true ;
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator-(FuncSpec)") ;
      resu.set_coefs(t1.coefs) ;
      resu.coefs.set(0) -= x ;
    }
    return resu ;
  }

  FuncSpec operator-(double x, const FuncSpec& t1) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      resu.values = x - t1.values ;
      resu.values_up_to_date = true ;
      if (t1.coefs_up_to_date) {
	resu.coefs = -t1.coefs ;
	resu.coefs.set(0) += x ;
	resu.coefs_up_to_date = true ;
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator-(FuncSpec)") ;
      resu.set_coefs(-t1.coefs) ;
      resu.coefs.set(0) += x ;
    }
    return resu ;
  }

  // Operators*
  //------------

  FuncSpec operator*(const FuncSpec& t1, const FuncSpec& t2) {

    FuncSpec resu(t1) ;
    if (!resu.values_up_to_date) resu.compute_values() ;
    if (t2.values_up_to_date) {
      resu.values = resu.values*t2.values ;
      resu.values_up_to_date = true ;
      resu.coefs_up_to_date = false ;
    }
    else {
      FuncSpec tmp(t2) ;
      tmp.compute_values() ;
      resu.values = resu.values*tmp.values ;
      resu.values_up_to_date = true ;
      resu.coefs_up_to_date = false ;
    }
    return resu ;
  }

  FuncSpec operator*(const FuncSpec& t1, double x) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      resu.values = t1.values * x ;
      resu.values_up_to_date = true ;
      if (t1.coefs_up_to_date) {
	resu.coefs = t1.coefs * x ;
	resu.coefs_up_to_date = true ;
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator*(FuncSpec)") ;
      resu.set_coefs(x*t1.coefs) ;
    }
    return resu ;
  }

  FuncSpec operator*(double x, const FuncSpec& t1) {

    FuncSpec resu(t1.nx, t1.ny, t1.nz) ;
  
    if (t1.values_up_to_date) {
      resu.values = t1.values * x ;
      resu.values_up_to_date = true ;
      if (t1.coefs_up_to_date) {
	resu.coefs = t1.coefs * x ;
	resu.coefs_up_to_date = true ;
      }
    }
    else {
      if (!t1.coefs_up_to_date)
	throw invalid_argument("Ill-formed first argument in operator*(FuncSpec)") ;
      resu.set_coefs(x*t1.coefs) ;
    }
    return resu ;
  }

  // Function abs
  //--------------
  FuncSpec abs(const FuncSpec& t) {

    FuncSpec resu(t.nx, t.ny, t.nz) ;

    if (t.values_up_to_date) {
      resu.values = abs(t.values) ;
      resu.values_up_to_date = true ;
      resu.coefs_up_to_date = false ;
    }
    else {
      if (!t.coefs_up_to_date)
	throw invalid_argument("Ill-formed argument in abs(FuncSpec)") ;
      FuncSpec tmp = t ;
      tmp.compute_values() ;
      resu.values = abs(tmp.values) ;
      resu.values_up_to_date = true ;
      resu.coefs_up_to_date = false ;
    }
    return resu ;  
  }

  // Function max
  //---------------

  double max(const FuncSpec& fs) {

    if (!fs.values_up_to_date)
      throw(std::runtime_error("Values not up to date in max(FuncSpec)")) ;  
    return max(fs.values) ;

  }


} // end namespace Lorene
