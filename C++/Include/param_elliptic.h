/*
 *  Definition of Lorene class Param_elliptic
 *
 */

/*
 *   Copyright (c) 2003 Philippe Grandclement
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

#ifndef __PARAM_ELLIPTIC_H_ 
#define __PARAM_ELLIPTIC_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2004/05/10 15:28:21  j_novak
 * First version of functions for the solution of the r-component of the
 * vector Poisson equation.
 *
 * Revision 1.6  2004/03/23 14:54:45  j_novak
 * More documentation
 *
 * Revision 1.5  2004/03/17 15:58:47  p_grandclement
 * Slight modification of sol_elliptic_no_zec
 *
 * Revision 1.4  2004/03/05 09:18:48  p_grandclement
 * Addition of operator sec_order_r2
 *
 * Revision 1.3  2004/02/11 09:47:44  p_grandclement
 * Addition of a new elliptic solver, matching with the homogeneous solution
 * at the outer shell and not solving in the external domain (more details
 * coming soon ; check your local Lorene dealer...)
 *
 * Revision 1.2  2004/01/28 16:46:22  p_grandclement
 * Addition of the sol_elliptic_fixe_der_zero stuff
 *
 * Revision 1.1  2003/12/11 14:57:00  p_grandclement
 * I had forgotten the .h (sorry folks...)
 *
 * Revision 1.2  2001/12/11 06:44:41  e_gourgoulhon
 *
 * $Header$
 *
 */

#include "map.h"
#include "ope_elementary.h"
#include "change_var.h"
#include "scalar.h"


/**
 * This class contains the parameters needed to call the general
 * elliptic solver.
 * 
 * For every domain and every spherical harmonics, it contains the
 * appropriate operator of type \c Ope_elementary and the appropriate
 * variable given by a \c Change_var . \ingroup (ellip)
 *
 * This class is only defined on an affine mapping \c Map_af .
 * 
 **/

class Param_elliptic {

 protected:
  const Map_af* mp ; ///< The affine mapping.
  Ope_elementary** operateurs ; ///< Array on the elementary operators.
  Change_var** variables ; ///< Array on the variable changes.

 public:
  /**
   * Standard constructor from a \c Scalar 
   * @param so [parameter] type
   * of the source of the elliptic equation. The actual values are not
   * used but \c *this will be constructed using the same number of
   * points, domains and symetry than \c so .
   * 
   * This constructor initializes everything to solve a Poisson
   * equation with non variable changes from domains to another.
   **/
  Param_elliptic (const Scalar&) ;
  ~Param_elliptic() ; ///< Destructor.
  
  /// Returns the affine mapping.
  const Map_af& get_mp() const {return *mp ;} ;

 public:  
  /**
   * Set the operator to \f$\left(\Delta - m^2\right)\f$ in one domain
   * (not in the nucleus).
   *
   * @param zone [input] : the domain.
   * @param mas [input] : the masse \f$m\f$.
   **/
  void set_helmholtz_minus (int zone, double mas) ;
   /**
    * Set the operator to \f$\left(\Delta + m^2\right)\f$ in one
    * domain (only in the shells).
    *
    * @param zone [input] : the domain.
    * @param mas [input] : the masse \f$m\f$.
    **/
  void set_helmholtz_plus (int zone, double mas) ; 

  /**
    * Set the operator to \f$a r^2 \partial^2/\partial r^2 + 
    * b r \partial /\partial r + c\f$ in one domain (only in the shells).
    *
    * @param zone [input] : the domain.
    * @param a [input] : the parameter \f$a\f$.
    * @param b [input] : the parameter \f$b\f$.
    * @param c [input] : the parameter \f$c\f$.
    **/
  void set_sec_order_r2 (int zone, double a, double b, double c) ;
  
   /**
    * Sets the operator to \f$\Delta + \frac{2}{r} \frac{\partial}{\partial r} 
    * + \frac{2}{r^2} \f$ in all domains, for \f$ l \not= 0 \f$; and to 
    * \f$\frac{\partial^2}{\partial r^2} + \frac{2}{r} 
    * \frac{\partial}{\partial r} - \frac{2}{r^2} \f$ 
    * in all domains otherwise.
    *
    * @param zone [input] : the domain.
    **/
  void set_poisson_vect_r(int zone) ; 

  /**
   * Increases the quantum number \e l in the domain \c zone .
   **/
  void inc_l_quant (int zone) ;
  
  /**
   * Changes the variable to the new type \c tipe  in the domain \c zone .
   **/
  void set_variable (int tipe, int zone) ;

  friend Mtbl_cf elliptic_solver  (const Param_elliptic&, const Mtbl_cf&) ;
  friend Mtbl_cf elliptic_solver_no_zec  
    (const Param_elliptic&, const Mtbl_cf&, double) ;
  friend Mtbl_cf elliptic_solver_sin_zec  
    (const Param_elliptic&, const Mtbl_cf&, double, int, double&, double&) ;
  friend Mtbl_cf elliptic_solver_fixe_der_zero  
    (double, const Param_elliptic&, const Mtbl_cf&) ;
} ;

#endif
