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
 * Revision 1.18  2007/01/16 15:05:59  n_vasset
 * New constructor (taking a Scalar in mono-domain angular grid for
 * boundary) for function sol_elliptic_boundary
 *
 * Revision 1.17  2005/11/30 11:09:03  p_grandclement
 * Changes for the Bin_ns_bh project
 *
 * Revision 1.16  2005/08/26 14:02:38  p_grandclement
 * Modification of the elliptic solver that matches with an oscillatory exterior solution
 * small correction in Poisson tau also...
 *
 * Revision 1.15  2005/06/09 07:56:25  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 * Revision 1.14  2005/02/15 15:43:16  j_novak
 * First version of the block inversion for the vector Poisson equation (method 6).
 *
 * Revision 1.13  2004/12/23 16:30:14  j_novak
 * New files and class for the solution of the rr component of the tensor Poisson
 * equation.
 *
 * Revision 1.12  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.11  2004/06/22 08:49:57  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.10  2004/06/14 15:23:07  j_novak
 * Modif. comments.
 *
 * Revision 1.9  2004/06/14 15:07:10  j_novak
 * New methods for the construction of the elliptic operator appearing in
 * the vector Poisson equation (acting on eta).
 *
 * Revision 1.8  2004/05/14 08:51:00  p_grandclement
 * *** empty log message ***
 *
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
#include "scalar.h"

#define MAP_AFF 0 
#define MAP_LOG 1 

/**
 * This class contains the parameters needed to call the general
 * elliptic solver.
 * 
 * For every domain and every spherical harmonics, it contains the
 * appropriate operator of type \c Ope_elementary and the appropriate
 * variable given by a \c Change_var . \ingroup (ellip)
 *
 * This class is only defined on an affine mapping \c Map_af or 
 * a logarithmic one \c Map_log
 * 
 **/
class Param_elliptic {

 protected:
  int type_map ; ///< Type of mapping either MAP_AFF or MAP_LOG
  const Map_af* mp_af ; ///< The mapping, if affine.
  const Map_log* mp_log ; ///< The mapping if log type.
 
  Ope_elementary** operateurs ; ///< Array on the elementary operators.
  
  Scalar var_F ; ///< Additive variable change function.
  Scalar var_G ; ///< Multiplicative variable change that must be sphericaly symetric !

  mutable Itbl done_F ; ///< Stores what has been computed for \c F
  mutable Itbl done_G ; ///< Stores what has been computed for \c G
  mutable Tbl val_F_plus ; ///< Values of F at the outer boundaries of the various domains.
  mutable Tbl val_F_minus ; ///< Values of F at the inner boundaries of the various domains.
  mutable Tbl val_dF_plus ; ///< Values of the derivative of F at the outer boundaries of the various domains.
  mutable Tbl val_dF_minus ; ///< Values of the derivative of F at the inner boundaries of the various domains.
  mutable Tbl val_G_plus ; ///< Values of G at the outer boundaries of the various domains.
  mutable Tbl val_G_minus ; ///< Values of G at the inner boundaries of the various domains.
  mutable Tbl val_dG_plus ; ///< Values of the derivative of G at the outer boundaries of the various domains.
  mutable Tbl val_dG_minus ; ///< Values of the derivative of G at the inner boundaries of the various domains.

 private:
  void compute_val_F(int, int, int) const ; ///< Computes the various values of \c F
  void compute_val_G(int) const ; ///< Computes the various values of \c G

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
  
  /// Returns the mapping.
  const Map_radial& get_mp() const ;
  double get_alpha (int) const ;
  double get_beta (int) const ;
  int get_type (int) const ;

 public:  
  /**
   * Set the operator to \f$\left(\Delta - m^2\right)\f$ in one domain
   * (not in the nucleus).
   *
   * @param zone [input] : the domain.
   * @param mas [input] : the masse \f$m\f$.
   * @param so [input] : the source (used only to get the right basis).
   **/
  void set_helmholtz_minus (int zone, double mas, Scalar& so) ;
  /**
   * Set the operator to \f$\left(\Delta - 2 \partial_r / r\right)\f$ 
   * everywhere but in the compactified domain.
   *
   * @param so [input] : the source (used only to get the right basis).
   **/
  void set_poisson_pseudo_1d (Scalar& so) ;
  /**
   * Set the operator to \f$\left(\Delta - 2 \partial_r / r - m^2\right)\f$ in one domain
   *
   * @param zone [input] : the domain.
   * @param mas [input] : the masse \f$m\f$.
   * @param so [input] : the source (used only to get the right basis).
   **/
  void set_helmholtz_minus_pseudo_1d (int zone, double mas, Scalar& so) ;

   /**
    * Set the operator to \f$\left(\Delta + m^2\right)\f$ in one
    * domain (only in the shells).
    *
    * @param zone [input] : the domain.
    * @param mas [input] : the masse \f$m\f$.
    **/
  void set_helmholtz_plus (int zone, double mas) ; 
 
  /**
   * Set everything to do a 2d-Poisson, with or without l=0 (not put by default...)
   **/
  void set_poisson_2d (const Scalar &, bool indic = false) ;  
  
  /**
   * Set the 2D Helmholtz operator (with minus sign)
   *
   *  @param zone [input] : the domain.
   * @param mas [input] : the masse parameter.
   **/
  void set_helmholtz_minus_2d (int zone, double mas, const Scalar&) ; 

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
    * Set the operator to \f$a \partial^2/\partial r^2 + 
    * b \partial /\partial r + c\f$ in one domain (only in the shells).
    *
    * @param zone [input] : the domain.
    * @param a [input] : the parameter \f$a\f$.
    * @param b [input] : the parameter \f$b\f$.
    * @param c [input] : the parameter \f$c\f$.
    **/
  void set_sec_order (int zone, double a, double b, double c) ;
  
   /**
    * Sets the operator to \f$\Delta + \frac{2}{r} \frac{\partial}{\partial r} 
    * + \frac{2 - l(l+1)}{r^2} \f$ in all domains, for \f$ l \not= 0 \f$; 
    * and to \f$\frac{\partial^2}{\partial r^2} + \frac{2}{r} 
    * \frac{\partial}{\partial r} - \frac{2}{r^2} \f$ in all domains otherwise.
    *
    * @param zone [input] : the domain.
    * @param only_l_zero [input] : the operator in built only for l=0 
    **/
  void set_poisson_vect_r(int zone, bool only_l_zero = false) ; 

   /**
    * Sets the operator to be a regular elliptic operator to solve for the
    * \f$\eta \f$ component of the vector Poisson equation. The operator is
    * \f$\frac{\partial^2}{\partial r^2} + 
    * \frac{2}{r} \frac{\partial}{\partial r} - \frac{l(l-1)}{r^2} \f$ 
    * (Poisson with the decrease of \e l by one unit)
    * in all domains but the CED, for \f$ l \not= 0 \f$; it is not defined 
    * for \e l = 0. In the CED, the operator is also the Laplace one, but 
    * with \e l increased by one unit: \f$\frac{\partial^2}{\partial r^2} + 
    * \frac{2}{r} \frac{\partial}{\partial r} - \frac{(l+1)(l+2)}{r^2} \f$. 
    * This is intended to solve the equation for \f$ \eta \f$ arising in 
    * the decomposition of the vector Poisson equation.
    *
    * @param zone [input] : the domain.
    **/
  void set_poisson_vect_eta(int zone) ; 

   /**
    * Sets the operator to \f$\Delta + \frac{4}{r} \frac{\partial}{\partial r} 
    * + \frac{6 - l(l+1)}{r^2} \f$ in all domains, for \f$l \geq 2\f$ only. 
    *
    * @param zone [input] : the domain.
    **/
  void set_poisson_tens_rr(int zone) ; 

  /**
   * Increases the quantum number \e l in the domain \c zone .
   **/
  void inc_l_quant (int zone) ;
  
  /**
   * Changes the variable function F
   **/
  void set_variable_F (const Scalar&) ;

  /**
   * Changes the variable function G
   **/
  void set_variable_G (const Scalar&) ;

  /**
   * Returns the value of F, for a given angular point, at the outer boundary of 
   * the domain \c zone ;
   **/
  double F_plus (int zone, int k, int j) const ;

   /**
   * Returns the value of F, for a given angular point, at the inner boundary of 
   * the domain \c zone ;
   **/
  double F_minus (int zone, int k, int j) const ;

  /**
   * Returns the value of the radial derivative of F, 
   * for a given angular point, at the outer boundary of 
   * the domain \c zone ;
   **/
  double dF_plus (int zone, int k, int j) const ;
 
  /**
   * Returns the value of the radial derivative of F, 
   * for a given angular point, at the inner boundary of 
   * the domain \c zone ;
   **/
  double dF_minus (int zone, int k, int j) const ;

  /**
   * Returns the value of G, for a given angular point, at the outer boundary of 
   * the domain \c zone ;
   **/
  double G_plus (int zone) const ;

   /**
   * Returns the value of G, for a given angular point, at the inner boundary of 
   * the domain \c zone ;
   **/
  double G_minus (int zone) const ;

  /**
   * Returns the value of the radial derivative of G, 
   * for a given angular point, at the outer boundary of 
   * the domain \c zone ;
   **/
  double dG_plus (int zone) const ;
 
  /**
   * Returns the value of the radial derivative of G, 
   * for a given angular point, at the inner boundary of 
   * the domain \c zone ;
   **/
  double dG_minus (int zone) const ;

  // A lot of friend functions... Possiblement pas toutes utiles...
  friend Mtbl_cf elliptic_solver  (const Param_elliptic&, const Mtbl_cf&) ;
  friend Mtbl_cf elliptic_solver_boundary (const Param_elliptic&, const Mtbl_cf&, const Mtbl_cf& bound, double fact_dir, double fact_neu   ) ;
  friend Mtbl_cf elliptic_solver_no_zec  
    (const Param_elliptic&, const Mtbl_cf&, double) ;
  friend Mtbl_cf elliptic_solver_only_zec  
    (const Param_elliptic&, const Mtbl_cf&, double) ;
  friend Mtbl_cf elliptic_solver_sin_zec  
    (const Param_elliptic&, const Mtbl_cf&, double, double&, double) ;
  friend Mtbl_cf elliptic_solver_fixe_der_zero  
    (double, const Param_elliptic&, const Mtbl_cf&) ;

  friend void Map_af::sol_elliptic(Param_elliptic&, const Scalar&, Scalar&) const ;
  friend void Map_af::sol_elliptic_boundary(Param_elliptic&, const Scalar&, Scalar&, const Mtbl_cf& , 
					    double , double ) const ;
  friend void Map_af::sol_elliptic_boundary(Param_elliptic&, const Scalar&, Scalar&, const Scalar& , 
					    double , double ) const ;


  friend void Map_af::sol_elliptic_no_zec(Param_elliptic&, const Scalar&, Scalar&, double) const ;
  friend void Map_af::sol_elliptic_only_zec(Param_elliptic&, const Scalar&, Scalar&, double) const ;
  friend void Map_af::sol_elliptic_sin_zec(Param_elliptic&, const Scalar&, Scalar&, double, double&, double) const ;
  friend void Map_af::sol_elliptic_fixe_der_zero(double, Param_elliptic&, const Scalar&, Scalar&) const ;

  friend void Map_af::sol_elliptic_2d(Param_elliptic&, const Scalar&, Scalar&) const ;
  friend void Map_af::sol_elliptic_pseudo_1d(Param_elliptic&, const Scalar&, Scalar&) const ;

  friend void Map_log::sol_elliptic(Param_elliptic&, const Scalar&, Scalar&) const ;
  friend void Map_log::sol_elliptic_boundary(Param_elliptic&, const Scalar&, Scalar&, const Mtbl_cf&, 
					     double, double) const ;
 
 friend void Map_log::sol_elliptic_boundary(Param_elliptic&, const Scalar&, Scalar&, const Scalar&, 
					     double, double) const ;


  friend void Map_log::sol_elliptic_no_zec(Param_elliptic&, const Scalar&, Scalar&, double) const ;
} ;

#endif
