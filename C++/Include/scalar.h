/*
 *  Definition of Lorene class Scalar
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
 *   
 *   Copyright (c) 1999-2000 Jean-Alain Marck (for previous class Cmp)
 *   Copyright (c) 1999-2002 Eric Gourgoulhon (for previous class Cmp)
 *   Copyright (c) 1999-2001 Philippe Grandclement (for previous class Cmp)
 *   Copyright (c) 2000-2002 Jerome Novak (for previous class Cmp)
 *   Copyright (c) 2000-2001 Keisuke Taniguchi (for previous class Cmp)
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


#ifndef __SCALAR_H_ 
#define __SCALAR_H_ 


/*
 * $Id$
 * $Log$
 * Revision 1.74  2005/11/17 15:29:46  e_gourgoulhon
 * Added arithmetics with Mtbl.
 *
 * Revision 1.73  2005/10/25 08:56:34  p_grandclement
 * addition of std_spectral_base in the case of odd functions near the origin
 *
 * Revision 1.72  2005/09/07 13:10:47  j_novak
 * Added a filter setting to zero all mulitpoles in a given range.
 *
 * Revision 1.71  2005/08/26 14:02:38  p_grandclement
 * Modification of the elliptic solver that matches with an oscillatory exterior solution
 * small correction in Poisson tau also...
 *
 * Revision 1.70  2005/08/25 12:14:07  p_grandclement
 * Addition of a new method to solve the scalar Poisson equation, based on a multi-domain Tau-method
 *
 * Revision 1.69  2005/06/09 07:56:25  f_limousin
 * Implement a new function sol_elliptic_boundary() and
 * Vector::poisson_boundary(...) which solve the vectorial poisson
 * equation (method 6) with an inner boundary condition.
 *
 * Revision 1.68  2005/06/08 12:35:18  j_novak
 * New method for solving divergence-like ODEs.
 *
 * Revision 1.67  2005/05/20 14:42:30  j_novak
 * Added the method Scalar::get_spectral_base().
 *
 * Revision 1.66  2005/04/04 21:28:57  e_gourgoulhon
 * Added argument lambda to method poisson_angu
 * to deal with the generalized angular Poisson equation:
 *    Lap_ang u + lambda u = source.
 *
 * Revision 1.65  2004/12/14 09:09:39  f_limousin
 * Modif. comments.
 *
 * Revision 1.64  2004/11/23 12:41:53  f_limousin
 * Intoduce function poisson_dir_neu(...) to solve a scalar poisson
 * equation with a mixed boundary condition (Dirichlet + Neumann).
 *
 * Revision 1.63  2004/10/11 15:09:00  j_novak
 * The radial manipulation functions take Scalar as arguments, instead of Cmp.
 * Added a conversion operator from Scalar to Cmp.
 * The Cmp radial manipulation function make conversion to Scalar, call to the
 * Map_radial version with a Scalar argument and back.
 *
 * Revision 1.62  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.61  2004/07/27 08:24:26  j_novak
 * Modif. comments
 *
 * Revision 1.60  2004/07/26 16:02:21  j_novak
 * Added a flag to specify whether the primitive should be zero either at r=0
 * or at r going to infinity.
 *
 * Revision 1.59  2004/07/06 13:36:27  j_novak
 * Added methods for desaliased product (operator |) only in r direction.
 *
 * Revision 1.58  2004/06/22 08:49:57  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.57  2004/06/14 15:24:23  e_gourgoulhon
 * Added method primr (radial primitive).
 *
 * Revision 1.56  2004/06/11 14:29:56  j_novak
 * Scalar::multipole_spectrum() is now a const method.
 *
 * Revision 1.55  2004/05/24 14:07:31  e_gourgoulhon
 * Method set_domain now includes a call to del_deriv() for safety.
 *
 * Revision 1.54  2004/05/07 11:26:10  f_limousin
 * New method filtre_r(int*)
 *
 * Revision 1.53  2004/03/22 13:12:43  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.52  2004/03/17 15:58:47  p_grandclement
 * Slight modification of sol_elliptic_no_zec
 *
 * Revision 1.51  2004/03/11 12:07:30  e_gourgoulhon
 * Added method visu_section_anim.
 *
 * Revision 1.50  2004/03/08 15:45:38  j_novak
 * Modif. comment
 *
 * Revision 1.49  2004/03/05 15:09:40  e_gourgoulhon
 * Added method smooth_decay.
 *
 * Revision 1.48  2004/03/01 09:57:02  j_novak
 * the wave equation is solved with Scalars. It now accepts a grid with a
 * compactified external domain, which the solver ignores and where it copies
 * the values of the field from one time-step to the next.
 *
 * Revision 1.47  2004/02/27 09:43:58  f_limousin
 * New methods filtre_phi(int) and filtre_theta(int).
 *
 * Revision 1.46  2004/02/26 22:46:26  e_gourgoulhon
 * Added methods derive_cov, derive_con and derive_lie.
 *
 * Revision 1.45  2004/02/21 17:03:49  e_gourgoulhon
 * -- Method "point" renamed "val_grid_point".
 * -- Method "set_point" renamed "set_grid_point".
 *
 * Revision 1.44  2004/02/19 22:07:35  e_gourgoulhon
 * Added argument "comment" in method spectral_display.
 *
 * Revision 1.43  2004/02/11 09:47:44  p_grandclement
 * Addition of a new elliptic solver, matching with the homogeneous solution
 * at the outer shell and not solving in the external domain (more details
 * coming soon ; check your local Lorene dealer...)
 *
 * Revision 1.42  2004/01/28 16:46:22  p_grandclement
 * Addition of the sol_elliptic_fixe_der_zero stuff
 *
 * Revision 1.41  2004/01/28 13:25:40  j_novak
 * The ced_mult_r arguments have been suppressed from the Scalar::*dsd* methods.
 * In the div/mult _r_dzpuis, there is no more default value.
 *
 * Revision 1.40  2004/01/28 10:39:17  j_novak
 * Comments modified.
 *
 * Revision 1.39  2004/01/27 15:10:01  j_novak
 * New methods Scalar::div_r_dzpuis(int) and Scalar_mult_r_dzpuis(int)
 * which replace div_r_inc*. Tried to clean the dzpuis handling.
 * WARNING: no testing at this point!!
 *
 * Revision 1.38  2004/01/23 13:25:44  e_gourgoulhon
 * Added methods set_inner_boundary and set_outer_boundary.
 * Methods set_val_inf and set_val_hor, which are particular cases of
 * the above, have been suppressed.
 *
 * Revision 1.37  2004/01/22 16:10:09  e_gourgoulhon
 * Added (provisory) method div_r_inc().
 *
 * Revision 1.36  2003/12/16 06:32:20  e_gourgoulhon
 * Added method visu_box.
 *
 * Revision 1.35  2003/12/14 21:46:35  e_gourgoulhon
 * Added argument start_dx in visu_section.
 *
 * Revision 1.34  2003/12/11 16:19:38  e_gourgoulhon
 * Added method visu_section for visualization with OpenDX.
 *
 * Revision 1.33  2003/12/11 14:48:47  p_grandclement
 * Addition of ALL (and that is a lot !) the files needed for the general elliptic solver ... UNDER DEVELOPEMENT...
 *
 * Revision 1.32  2003/11/13 13:43:53  p_grandclement
 * Addition of things needed for Bhole::update_metric (const Etoile_bin&, double, double)
 *
 * Revision 1.31  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.30  2003/11/04 22:55:50  e_gourgoulhon
 * Added new methods mult_cost(), mult_sint() and div_sint().
 *
 * Revision 1.29  2003/10/29 13:09:11  e_gourgoulhon
 * -- Added integer argument to derivative functions dsdr, etc...
 *    so that one can choose the dzpuis of the result (default=2).
 * -- Change of method name: laplacien --> laplacian.
 *
 * Revision 1.28  2003/10/29 11:00:42  e_gourgoulhon
 * Virtual functions dec_dzpuis and inc_dzpuis have now an integer argument to
 *  specify by which amount dzpuis is to be increased.
 * Accordingly virtual methods dec2_dzpuis and inc2_dzpuis have been suppressed.
 *
 * Revision 1.27  2003/10/20 14:26:02  j_novak
 * New assignement operators.
 *
 * Revision 1.26  2003/10/19 19:46:33  e_gourgoulhon
 * -- Method spectral_display now virtual (from Tensor), list of argument
 *    changed.
 *
 * Revision 1.25  2003/10/17 13:46:14  j_novak
 * The argument is now between 1 and 3 (instead of 0->2)
 *
 * Revision 1.24  2003/10/16 15:23:41  e_gourgoulhon
 * Name of method div_r_ced() changed to div_r_inc2().
 * Name of method div_rsint_ced() changed to div_rsint_inc2().
 *
 * Revision 1.23  2003/10/15 21:10:11  e_gourgoulhon
 * Added method poisson_angu().
 *
 * Revision 1.22  2003/10/15 16:03:35  j_novak
 * Added the angular Laplace operator for Scalar.
 *
 * Revision 1.21  2003/10/15 10:29:05  e_gourgoulhon
 * Added new members p_dsdt and p_stdsdp.
 * Added new methods dsdt(), stdsdp() and div_tant().
 *
 * Revision 1.20  2003/10/13 13:52:39  j_novak
 * Better managment of derived quantities.
 *
 * Revision 1.19  2003/10/10 15:57:27  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.18  2003/10/08 14:24:08  j_novak
 * replaced mult_r_zec with mult_r_ced
 *
 * Revision 1.17  2003/10/06 16:16:02  j_novak
 * New constructor from a Tensor.
 *
 * Revision 1.16  2003/10/06 13:58:45  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.15  2003/10/05 21:06:31  e_gourgoulhon
 * - Added new methods div_r_ced() and div_rsint_ced().
 * - Added new virtual method std_spectral_base()
 * - Removed method std_spectral_base_scal()
 *
 * Revision 1.14  2003/10/01 13:02:58  e_gourgoulhon
 * Suppressed the constructor from Map* .
 *
 * Revision 1.13  2003/09/29 12:52:56  j_novak
 * Methods for changing the triad are implemented.
 *
 * Revision 1.12  2003/09/25 09:33:36  j_novak
 * Added methods for integral calculation and various manipulations
 *
 * Revision 1.11  2003/09/25 09:11:21  e_gourgoulhon
 * Added functions for radial operations (divr, etc...)
 *
 * Revision 1.10  2003/09/25 08:55:23  e_gourgoulhon
 * Added members raccord*.
 *
 * Revision 1.9  2003/09/25 08:50:11  j_novak
 * Added the members import
 *
 * Revision 1.8  2003/09/25 08:13:51  j_novak
 * Added method for calculating derivatives
 *
 * Revision 1.7  2003/09/25 07:59:26  e_gourgoulhon
 * Added prototypes for PDE resolutions.
 *
 * Revision 1.6  2003/09/25 07:17:58  j_novak
 * Method asymptot implemented.
 *
 * Revision 1.5  2003/09/24 20:53:38  e_gourgoulhon
 * Added  -- constructor by conversion from a Cmp
 *        -- assignment from Cmp
 *
 * Revision 1.4  2003/09/24 15:10:54  j_novak
 * Suppression of the etat flag in class Tensor (still present in Scalar)
 *
 * Revision 1.3  2003/09/24 12:01:44  j_novak
 * Added friend functions for math.
 *
 * Revision 1.2  2003/09/24 10:22:01  e_gourgoulhon
 * still in progress...
 *
 * Revision 1.1  2003/09/22 12:50:47  e_gourgoulhon
 * First version: not ready yet!
 *
 *
 * $Header$
 *
 */

// Headers Lorene 
#include "valeur.h"
#include "tensor.h"

class Param ; 
class Cmp ;
class Param_elliptic ;

/**
 * Tensor field of valence 0 (or component of a tensorial field).
 * \ingroup (tensor)
 * 
 * 
 */

class Scalar : public Tensor {
  
  // Data : 
  // -----
 protected:
  
  /** The logical state \c ETATNONDEF  (undefined), \c ETATZERO  (null),
   *  \c ETATUN  (one), or \c ETATQCQ (ordinary).
   */
  int etat ; 
  
  /**
   * Power of \e r  by which the quantity represented by \c this  
   * must be divided in the  compactified external domain (CED) in order 
   * to get the correct physical values
   */
  int dzpuis ;	
  
  Valeur va ;		///< The numerical value of the \c Scalar    
  
  // Derived data : 
  // ------------
 protected:
  /// Pointer on \f$\partial/\partial r\f$ of \c *this  (0x0 if not up to date)
  mutable Scalar* p_dsdr ;	

  /** Pointer on \f$1/r \partial/\partial \theta\f$ of \c *this  
   *  (0x0 if not up to date)
   */
  mutable Scalar* p_srdsdt ;	

  /** Pointer on \f$1/(r\sin\theta) \partial/\partial \phi\f$ of \c *this 
   *  (0x0 if not up to date)
   */
  mutable Scalar* p_srstdsdp ;

  /// Pointer on \f$\partial/\partial \theta\f$ of \c *this  (0x0 if not up to date)
  mutable Scalar* p_dsdt ;	

  /** Pointer on \f$1/\sin\theta \partial/\partial \phi\f$ of \c *this 
   *  (0x0 if not up to date)
   */
  mutable Scalar* p_stdsdp ;	
  
  /** Pointer on \f$\partial/\partial x\f$ of \c *this ,
   *  where \f$x=r\sin\theta \cos\phi\f$ (0x0 if not up to date)
   */
  mutable Scalar* p_dsdx ;	
  
  /** Pointer on \f$\partial/\partial y\f$ of \c *this ,
   *  where \f$y=r\sin\theta \sin\phi\f$(0x0 if not up to date)
   */
  mutable Scalar* p_dsdy ;	

  /** Pointer on \f$\partial/\partial z\f$ of \c *this ,
   *  where \f$z=r\cos\theta\f$ (0x0 if not up to date)
   */
  mutable Scalar* p_dsdz ;	
  
  /** Pointer on the Laplacian of \c *this  (0x0 if not up to date)
   */
  mutable Scalar* p_lap ;	
  
  /** Pointer on the Laplacian of \c *this  (0x0 if not up to date)
   */
  mutable Scalar* p_lapang ;	
   
  /// Pointer on \f$\partial/\partial radial \f$ of \c *this  
  mutable Scalar* p_dsdradial ;	

  /// Pointer on \f$\partial/\partial \rho \f$ of \c *this  
  mutable Scalar* p_dsdrho ;	

  /** Power of \e r  by which the last computed Laplacian has been 
   *  multiplied in the compactified external domain.  
   */
  mutable int ind_lap ; 

  /** Pointer on the space integral of \c *this  (values in each 
   *  domain) (0x0 if not up to date)
   */
  mutable Tbl* p_integ ; 
  
  // Constructors - Destructor
  // -------------------------
  
 public:
  
  explicit Scalar(const Map& mpi) ;	///< Constructor from mapping

  /// Constructor from a Tensor (must be of valence 0)
  Scalar(const Tensor& a) ;     

  Scalar(const Scalar& a) ;		///< Copy constructor
  
  Scalar(const Cmp& a) ;	///< Constructor by conversion of a Cmp
  
  /// Constructor from a file (see \c sauve(FILE*) )
  Scalar(const Map&, const Mg3d&, FILE* ) ;    		
  
  virtual ~Scalar() ;			///< Destructor
  
  
  // Memory management
  // -----------------
 protected:
  void del_t() ;		    ///< Logical destructor
  virtual void del_deriv() const;	    ///< Logical destructor of the derivatives
  void set_der_0x0() const;	    ///< Sets the pointers for derivatives to 0x0
  
 public:
  
  /**
   * Sets the logical state to \c ETATNONDEF  (undefined). 
   * Calls the logical destructor of the \c Valeur \c va   * deallocates the memory occupied by all the derivatives. 
   */
  virtual void set_etat_nondef() ;   
  
  /**
   * Sets the logical state to \c ETATZERO  (zero). 
   * Calls the logical destructor of the \c Valeur \c va  and
   * deallocates the memory occupied by all the derivatives. 
   */
  virtual void set_etat_zero() ;	    
  
  /**
   * Sets the logical state to \c ETATQCQ  (ordinary state).
   * If the state is already \c ETATQCQ , this function does nothing.
   * Otherwise, it calls the logical destructor of the \c Valeur \c va  and
   * deallocates the memory occupied by all the derivatives.
   */
  virtual void set_etat_qcq() ;	    
  
  /**
   * Sets the logical state to \c ETATUN  (one). 
   * Fills the \c Valeur \c va  with ones and
   * deallocates the memory occupied by all the derivatives. 
   */
  void set_etat_one() ;	    
  
  /**
   * Sets the logical state to \c ETATQCQ  (ordinary state)
   *  and performs the memory allocation of all the 
   *  elements, down to the \c double  arrays of the \c Tbl s. 
   *  This function performs in fact recursive calls to \c set_etat_qcq() 
   *  on each element of the chain \c Scalar ->
   *  \c Valeur  -> \c Mtbl  -> \c Tbl . 
   */
  virtual void allocate_all() ; 
  
  /**
   * Sets the \c Scalar to zero in a hard way. 
   * 1/ Sets the logical state to \c ETATQCQ , i.e. to an ordinary state.
   * 2/ Fills the \c Valeur \c va  with zeros. 
   * NB: this function must be used for debugging purposes only.
   * For other operations, the functions \c set_etat_zero() 
   * or \c annule(int,int) must be perferred. 
   */
  void annule_hard() ;
  
  // Extraction of information
  // -------------------------
    public:
  /** Returns the logical state \c ETATNONDEF  (undefined), 
   * \c ETATZERO (null) or \c ETATQCQ (ordinary).
   */
  int get_etat() const {return etat;} ; 
  
  /// Returns \c dzpuis 
  int get_dzpuis() const {return dzpuis;} ; 
  
  /** Returns \c true  if the last domain is compactified and
   *  \c *this  is not zero in this domain
   */
  bool dz_nonzero() const ; 
	
  /** Returns \c false  if the last domain is compactified 
   *  and \c *this  is not zero in this domain and \c dzpuis 
   *  is not equal to \c dzi , otherwise return true. 
   */
  bool check_dzpuis(int dzi) const ; 
  
  // Assignment
  // -----------
 public: 
  /// Assignment to another \c Scalar defined on the same mapping
  void operator=(const Scalar& a) ;	
  
  /// Assignment to a \c Tensor  (of valence 0)
  virtual void operator=(const Tensor& a) ; 

  void operator=(const Cmp& a) ; 	///< Assignment to a \c Cmp 
  void operator=(const Valeur& a) ; ///< Assignment to a \c Valeur 
  void operator=(const Mtbl& a) ;	 ///< Assignment to a \c Mtbl 
  void operator=(double ) ;	 ///< Assignment to a \c double 
  void operator=(int ) ;		 ///< Assignment to an \c int 
  
  // Conversion oprators
  //---------------------
  operator Cmp() const ; ///< Conversion to a \c Cmp

  // Access to individual elements
  // -----------------------------
    public:
  
  /// Returns \c va  (read only version)
  const Valeur& get_spectral_va() const {return va;} ; 
  
  /// Returns \c va  (read/write version)
  Valeur& set_spectral_va() {return va;} ; 
  
  /** Read/write of the value in a given domain.
   *  This method should be used only to set the value in a given
   *  domain (it performs a call to \c del_deriv); for reading the
   *  value in a domain without changing it, the method \c domain(int )
   *  is preferable.
   *
   * @param l [input] domain index
   * @return writable \c Tbl containing the value of the field in domain \c l .
   */ 
  Tbl& set_domain(int l) {
    assert(etat == ETATQCQ) ;
    del_deriv() ; 
    return va.set(l) ;
  };
  
  /** Read-only of the value in a given domain.
   * @param l [input] domain index
   * @return \c Tbl containing the value of the field in domain \c l .
   */ 
  const Tbl& domain(int l) const {
    assert( (etat == ETATQCQ) || (etat == ETATUN) ) ;
    return va(l) ;
  };
  
  
  /** Returns the value of the field at a specified grid point.
   * @param l [input] domain index
   * @param k [input] \f$\phi\f$ index
   * @param j [input] \f$\theta\f$ index
   * @param i [input] \e r  (\f$\xi\f$) index
   */ 
  double val_grid_point(int l, int k, int j, int i) const {
    assert(etat != ETATNONDEF) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else {
      if (etat == ETATUN) {
      double one = 1. ;
      return one ;
      }
      else{ 	    
	return va(l, k, j, i) ;
      }
    }
  };
  
  /** Computes the value of the field at an
   *   arbitrary point \f$(r, \theta, \phi)\f$, by means of the spectral 
   *   expansion.
   *   NB: if \f$(r, \theta, \phi)\f$ is a point of the spectral grid, 
   *     the method \c val_grid_point  is to be preferred, 
   *     being much more efficient. 
   *	 @param r [input] value of the coordinate \e r 
   *	 @param theta [input] value of the coordinate \f$\theta\f$
   *	 @param phi [input] value of the coordinate \f$\phi\f$
   *	 @return value at the point \f$(r, \theta, \phi)\f$ 
   *		 of the field represented by \c *this . 
   */
  double val_point(double r, double theta, double phi) const ; 
  
  
  /** Setting the value of the field at a given grid point.
   * CAUTION: to gain in efficiency (especially when this method is
   *  invoqued inside a loop), the method \c del_deriv()  (to delete
   *     the derived members) is not called by \c set_grid_point . 
   *     It must thus be invoqued by the user, after all the calls
   *     to  \c set_grid_point  have been performed.   
   *     
   * @param l [input] domain index
   * @param k [input] \f$\phi\f$ index
   * @param j [input] \f$\theta\f$ index
   * @param i [input] \e r  (\f$\xi\f$) index
   * @return writable value of the field at the specified grid point
   */ 
  double& set_grid_point(int l, int k, int j, int i) {
    assert(etat == ETATQCQ) ;
    return va.set(l, k, j, i) ;
  };
  
	
  /**
   * Sets the \c Scalar to zero in several domains.
   *	@param l_min [input] The \c Scalar will be set (logically) to zero
   *			     in the domains whose indices are in the range
   *			     \c [l_min,l_max] .
   *	@param l_max [input] see the comments for \c l_min .
   * 
   * Note that \c annule(0,va.mg->get_nzone()-1) is equivalent to
   *	 \c set_etat_zero() .
   */
  virtual void annule(int l_min, int l_max) ; 
  
  /** Sets the value of the \c Scalar at the inner boundary of a given 
   * domain. 
   * @param l [input] domain index
   * @param x [input] (constant) value at the inner boundary of domain no. \c l 
   */
  void set_inner_boundary(int l, double x) ;
    
  /** Sets the value of the \c Scalar at the outer boundary of a given 
   * domain. 
   * @param l [input] domain index
   * @param x [input] (constant) value at the outer boundary of domain no. \c l 
   */
  void set_outer_boundary(int l, double x) ;

  /**
   * Gives the spectrum in terms of multipolar modes \e l .
   *  @return a \c Tbl  of size (nzone, lmax), where lmax is the
   *  maximal multipolar momentum over all domains. The \e l -th
   *  element contains the L1 norm of the \e l -th multipole 
   *  (\e i.e. a sum over all \e m of the norms (coefficient space)
   *  of the component of a given \f$Y_l^m\f$.
   */
  Tbl multipole_spectrum () const ;
  
  // Differential operators and others
  // ---------------------------------
 public:
  /** Returns \f$\partial / \partial r\f$ of \c *this .
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdr() const ; 
 
  /** Returns \f$1/r \partial / \partial \theta\f$ of \c *this .
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   */
  const Scalar& srdsdt() const ; 
  
  /** Returns \f$1/(r\sin\theta) \partial / \partial \phi\f$ of \c *this .
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   */
  const Scalar& srstdsdp() const ; 
  
  /** Returns \f$\partial / \partial \theta\f$ of \c *this .
   */
  const Scalar& dsdt() const ; 
  
  /** Returns \f$\partial / \partial r\f$ of \c *this if the mapping is
   * affine (class \c Map_af) and \f$\partial / \partial \ln r\f$ 
   *  of \c *this if the mapping
   * is logarithmic (class \c Map_log).
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdradial() const ; 
  
  /** Returns \f$\partial / \partial \rho \f$ of \c *this .
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdrho() const ;
 
  /** Returns \f$1/\sin\theta \partial / \partial \phi\f$ of \c *this .
   */
  const Scalar& stdsdp() const ; 
  
  /** Returns \f$\partial/\partial x\f$ of \c *this ,
   *  where \f$x=r\sin\theta \cos\phi\f$.
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdx() const ;	
  
  /** Returns \f$\partial/\partial y\f$ of \c *this ,
   *  where \f$y=r\sin\theta \sin\phi\f$.
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdy() const ;	
  
  /** Returns \f$\partial/\partial z\f$ of \c *this ,
   *  where \f$z=r\cos\theta\f$.
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdz() const ;	
  
  /** Returns \f$\partial/\partial x_i\f$ of \c *this ,
   *  where \f$x_i = (x, y, z)\f$.
   *  If \c dzpuis  is zero, then the returned \c Scalar has 
   *  \c dzpuis  = 2. It is increased by 1 otherwise.
   *  @param i [input] i=1 for \e x ,  i=2 for \e y , i=3 for \e z .
   */
  const Scalar& deriv(int i) const ;	
  
    /** Returns the gradient (1-form = covariant vector) of \c *this  
     *  @param gam metric components only used to get the triad with 
     *    respect to which the components of the result are defined        
     */
    const Vector& derive_cov(const Metric& gam) const ; 


    /** Returns the "contravariant" derivative of \c *this  with respect 
     * to some metric \f$\gamma\f$, by raising the index of the
     * gradient (cf. method \c derive_cov() ) with 
     * \f$\gamma\f$.
     */
    const Vector& derive_con(const Metric& gam) const ; 

    /// Computes the derivative of \c this  along a vector field \c v 
    Scalar derive_lie(const Vector& v) const ; 


  /** Returns the Laplacian of \c *this 
   *   @param ced_mult_r [input] Determines the quantity computed in
   *			 the  compactified external domain (CED) 
   *		(\e u  in the field represented by \c *this ) :  
   *		    \li ced_mult_r = 0 : \f$\Delta u\f$	
   *		    \li ced_mult_r = 2 : \f$r^2 \,  \Delta u\f$	
   *		    \li ced_mult_r = 4 (default) : \f$r^4 \, \Delta u\f$	
   */
  const Scalar& laplacian(int ced_mult_r = 4) const ; 
  
  /** Returns the angular Laplacian \f$\Delta_{\theta\varphi}\f$ of \c *this ,
   *  where \f$\Delta_{\theta\varphi} f = \frac{\partial^2 f}
   *  {\partial \theta^2} + \frac{1}{\tan \theta} \frac{\partial f}
   *  {\partial \theta} +\frac{1}{\sin^2 \theta}\frac{\partial^2 f}
   *  {\partial \varphi^2}\f$
   * 
   */
  const Scalar& lapang() const ; 
  
  /// Division by \e r  everywhere; \c dzpuis  is not changed.
  void div_r() ;    
 
  /** Division by \e r  everywhere but with the output flag \c dzpuis  
   *  set to \c ced_mult_r .
   *  @param  ced_mult_r [input] value of \c dzpuis  of the result.
   */
  void div_r_dzpuis(int ced_mult_r) ; 
  
  /** Division by \e r  in the compactified external domain (CED), the 
   * \c dzpuis  flag is not changed.
   */
  void div_r_ced() ;

  /// Multiplication by \e r  everywhere; \c dzpuis  is not changed.
  void mult_r() ;  
  
  /** Multiplication by \e r  everywhere but with the output flag \c dzpuis  
   *  set to \c ced_mult_r .
   *  @param  ced_mult_r [input] value of \c dzpuis  of the result. 
   */
  void mult_r_dzpuis(int ced_mult_r) ; 
  
  /** Multiplication by \e r  in the compactified external domain (CED), the 
   * \c dzpuis  flag is not changed.
   */
  void mult_r_ced() ;
  
  /// Multiplication by \f$r\sin\theta\f$ everywhere; \c dzpuis  is not changed.
  void mult_rsint() ;   
  
  /** Multiplication by \f$r\sin\theta\f$ but with the output flag \c dzpuis 
   *  set to \c ced_mult_r .
   *  @param  ced_mult_r [input] value of \c dzpuis  of the result. 
   */
  void mult_rsint_dzpuis(int ced_mult_r) ; 

  /// Division by \f$r\sin\theta\f$ everywhere; \c dzpuis  is not changed.
  void div_rsint() ;    
  
  /** Division by \f$r\sin\theta\f$ but with the output flag \c dzpuis 
   *  set to \c ced_mult_r .
   *  @param  ced_mult_r [input] value of \c dzpuis  of the result. 
   */
  void div_rsint_dzpuis(int ced_mult_r) ; 
  
  void mult_cost() ;   ///< Multiplication by \f$\cos\theta\f$

  void mult_sint() ;   ///< Multiplication by \f$\sin\theta\f$

  void div_sint() ;    ///< Division by \f$\sin\theta\f$

  void div_tant() ;    ///< Division by \f$\tan\theta\f$
  
  /** Computes the radial primitive which vanishes for \f$r\to \infty\f$.
   *  i.e. the function 
   *      \f$ F(r,\theta,\varphi) = \int_r^\infty f(r',\theta,\varphi) \, dr' \f$
   *  where \e f is the function represented by \c *this 
   * (and must have a \c dzpuis = 2). 
   *      
   *      @param null_infty if true (default), the primitive is null
   *      at infinity (or on the grid boundary). \e F is null at the
   *      center otherwise
   *      @return function \e F
   */ 
  Scalar primr(bool null_infty = true) const  ;  	    

  /** Computes the integral over all space of \c *this .
   *  The computed quantity is (\e u  being the field represented by
   *   \c *this )
   *    \f$\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi\f$.
   *  Note that in the compactified external domain (CED), \c dzpuis  
   *  must be 4 for the computation to take place. 
   */
  double integrale() const ; 
	
  /** Computes the integral in each domain of \c *this .
   *  The computed quantity is (\e u  being the field represented by
   *   \c *this )
   *    \f$\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi\f$
   *  in each domain. The result is returned a \c Tbl  on the 
   *  various domains. 
   *  Note that in the compactified external domain (CED), \c dzpuis  
   *  must be 4 for the computation to take place. 
   */
  const Tbl& integrale_domains() const ; 
	
  /** Decreases by \c dec  units the value of \c dzpuis  and 
   *  changes accordingly the values of the \c Scalar in the 
   *  compactified external domain (CED).
   */
  virtual void dec_dzpuis(int dec = 1) ; 

  /** Increases by \c inc  units the value of \c dzpuis  and 
   *  changes accordingly the values of the \c Scalar in the 
   *  compactified external domain (CED).
   */
  virtual void inc_dzpuis(int inc = 1) ; 
	
  /** Sets a new vectorial basis (triad) of decomposition and modifies
   *  the components accordingly. 
   */
  virtual void change_triad(const Base_vect& new_triad) ; 
    
  /**
   * Sets the \c n  lasts coefficients in \e r  to 0 in the 
   *  external domain.
   */
  void filtre (int n) ;
    
  /**
   * Sets the \c n  lasts coefficients in \e r  to 0 in all
   *  domains.
   */
  void filtre_r (int* nn) ;

  /**
   * Sets all the multipolar components between \c l_min and \c l_max
   * to zero. This is done for [ \c l_min , \c l_max ] and all relevant 
   * \c m in the spherical harmonics expansion basis. If \c ylm_output 
   * is set to \c true , then the spectral expansion basis of \c this is 
   * left to be that of spherical harmonics.
   */
  void annule_l (int l_min, int l_max, bool ylm_output= false ) ;

  /**
   * Sets the \c n  lasts coefficients in \f$\Phi\f$ to 0 in the 
   * domain \c zone .
   */
  void filtre_phi (int n, int zone) ;
 
  /**
   * Sets the \c n  lasts coefficients in \f$\theta\f$ to 0 in the 
   * domain \c nz1 to \c nz2 when expressed in spherical harmonics.
   */
  void filtre_tp(int nn, int nz1, int nz2) ;
  
   
  /**
   * Substracts all the components behaving like \f$r^{-n}\f$ in the external 
   * domain, with \e n  strictly lower than \c puis , so that \c *this  
   * decreases at least like \f$r^{\tt puis} \f$ at infinity.
   */
  void fixe_decroissance (int puis) ;

    /** Performs a \f$C^k\f$ matching of the last non-compactified shell with
     * a decaying function \f$\sum_{j=0}^k {\alpha_j \over r^{\ell+n+j}}\f$ where
     * \f$\ell\f$ is the spherical harmonic index and \e n  is some 
     * specifiable parameter. 
     */
    void smooth_decay(int k, int n) ; 

  /**
   * Performs the \f$C^n\f$ matching of the nucleus with respect to the 
   * first shell.
   */
  void raccord(int n) ;
	
  /**
   * Performs the \f$C^1\f$ matching of the external domain with respect to
   * the last shell using function like \f$\frac{1}{r^i}\f$ with 
   * \f${\tt puis}  \leq i \leq {\tt puis+nbre}\f$ for each spherical harmonics 
   * with \f$l \leq {\tt lmax}\f$.
   */
  void raccord_c1_zec(int puis, int nbre, int lmax) ;

  /**
   * Matching of the external domain with the outermost shell
   */
  void raccord_externe(int puis, int nbre, int lmax) ;

  // Outputs
  // -------
 public:
  virtual void sauve(FILE *) const ;	    ///< Save in a file
    
	/** Displays the spectral coefficients and the associated
	 *  basis functions. This function shows only the values greater than a 
	 *  given threshold.
         *   @param comment comment to be printed at top of the display
         *      (default: 0x0 = nothing printed)
	 *   @param threshold [input] Value above which a coefficient is printed
	 *    (default: 1.e-7)
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param ostr [input] Output stream used for the printing (default: cout)
	 */
	virtual void spectral_display(const char* comment = 0x0, 
                            double threshold = 1.e-7, int precision = 4, 
			    ostream& ostr = cout) const ;

  /// Display
  friend ostream& operator<<(ostream& , const Scalar & ) ;	
  
  /** 3D visualization via a plane section.
   * Prepares files for visualization by OpenDX of the values of the field in
   * a plane x=const, y=const or z=const
   *
   * @param section_type [input] defines the type of section : 
   *    \li 'x' for a plane x = a with a = const (parameter \c aa ) 
   *    \li 'y' for a plane y = a with a = const (parameter \c aa )
   *    \li 'z' for a plane z = a with a = const (parameter \c aa )
   * @param aa [input] constant a defining the section plane
   * @param umin [input] defines with \c umax  the range of the plane coordinate u 
   * @param umax [input] defines with \c umin  the range of the plane coordinate u 
   * @param vmin [input] defines with \c vmax  the range of the plane coordinate v 
   * @param vmax [input] defines with \c vmin  the range of the plane coordinate v 
   * @param title [input] title for the graph (for OpenDX legend)
   * @param filename [input] name for the file which will be the input for 
   *    OpenDX; the default 0x0 is transformed into "scalar_section"
   * @param start_dx [input] determines whether OpenDX must be launched (as a
   *     subprocess) to view the field; if set to \c false , only input files
   *     for future usage of OpenDX are created 
   * @param nu [input] number of points in the u direction (uniform sampling)   
   * @param nv [input] number of points in the v direction (uniform sampling)   
   *
   */
    void visu_section(const char section_type, double aa, double umin, double umax, double vmin,
        double vmax, const char* title = 0x0, const char* filename = 0x0,
        bool start_dx = true, int nu = 200, int nv = 200) const ;   

  /** 3D visualization via a plane section.
   * Prepares files for visualization by OpenDX of the values of the field in
   * any given plane.
   *
   * @param plane [input] : 2D \c Tbl  defining the section plane: \c plane 
   *    must of dimension 3x3 with the following content: 
   *    \li \c plane(0,i) : absolute Cartesian coordinates (xa0,ya0,za0) of some
   *    point in the plane considered as the origin for the plane coordinates
   *    (u,v): \c plane(0,0) = xa0 , \c plane(0,1) = ya0 , 
   *    \li \c plane(0,2) = za0   
   *    \li \c plane(1,i) : components w.r.t. absolute Cartesian coordinates 
   *        of the u-coordinate unit vector in the section plane 
   *    \li \c plane(2,i) : components w.r.t. absolute Cartesian coordinates 
   *        of the v-coordinate unit vector in the section plane
   * @param umin [input] defines with \c umax  the range of the plane coordinate u 
   * @param umax [input] defines with \c umin  the range of the plane coordinate u 
   * @param vmin [input] defines with \c vmax  the range of the plane coordinate v 
   * @param vmax [input] defines with \c vmin  the range of the plane coordinate v 
   * @param title [input] title for the graph (for OpenDX legend)
   * @param filename [input] name for the file which will be the input for 
   *    OpenDX; the default 0x0 is transformed into "scalar_section"
   * @param start_dx [input] determines whether OpenDX must be launched (as a
   *     subprocess) to view the field; if set to \c false , only input files
   *     for future usage of OpenDX are created 
   * @param nu [input] number of points in the u direction (uniform sampling)   
   * @param nv [input] number of points in the v direction (uniform sampling)   
   *
   */
    void visu_section(const Tbl& plane, double umin, double umax, double vmin,
        double vmax, const char* title = 0x0, const char* filename = 0x0,
        bool start_dx = true, int nu = 200, int nv = 200) const ;   

  /** 3D visualization via time evolving plane section (animation).
   * Prepares files for visualization by OpenDX of the values of the field in
   * a plane x=const, y=const or z=const at successive time steps
   *
   * @param section_type [input] defines the type of section : 
   *    \li 'x' for a plane x = a with a = const (parameter \c aa ) 
   *    \li 'y' for a plane y = a with a = const (parameter \c aa )
   *    \li 'z' for a plane z = a with a = const (parameter \c aa )
   * @param aa [input] constant a defining the section plane
   * @param umin [input] defines with \c umax  the range of the plane coordinate u 
   * @param umax [input] defines with \c umin  the range of the plane coordinate u 
   * @param vmin [input] defines with \c vmax  the range of the plane coordinate v 
   * @param vmax [input] defines with \c vmin  the range of the plane coordinate v 
   * @param jtime [input] time step label
   * @param ttime [input] time t corresponding to \c jtime 
   * @param jgraph [input] number of time steps between two graphs: the graph
   *    will be generated only if \c jtime  is a multiple of \c jgraph 
   * @param title [input] title for the graph (for OpenDX legend)
   * @param filename_root [input] beginning of the names for the files which will 
   *    be the input for OpenDX (the end of names will be automatically generated
   *    from the time steps); the default 0x0 is transformed into "anim"
   * @param start_dx [input] determines whether OpenDX must be launched (as a
   *     subprocess) to view the field; if set to \c false , only input files
   *     for future usage of OpenDX are created 
   * @param nu [input] number of points in the u direction (uniform sampling)   
   * @param nv [input] number of points in the v direction (uniform sampling)   
   *
   */
   void visu_section_anim(const char section_type, double aa, double umin, 
        double umax, double vmin, double vmax, int jtime, double ttime, 
        int jgraph = 1, const char* title = 0x0, const char* filename_root = 0x0, 
        bool start_dx = false, int nu = 200, int nv = 200) const ;   

  /** 3D visualization (volume rendering) via OpenDX.
   * Prepares files for visualization by OpenDX of the values of the field in
   * some rectangular box.
   *
   * @param xmin [input] defines with \c xmax  the x range of the visualization box 
   * @param xmax [input] defines with \c xmin  the x range of the visualization box 
   * @param ymin [input] defines with \c ymax  the y range of the visualization box 
   * @param ymax [input] defines with \c ymin  the y range of the visualization box 
   * @param zmin [input] defines with \c zmax  the z range of the visualization box 
   * @param zmax [input] defines with \c zmin  the z range of the visualization box 
   * @param title [input] title for the graph (for OpenDX legend)
   * @param filename [input] name for the file which will be the input for 
   *    OpenDX; the default 0x0 is transformed into "scalar_box"
   * @param start_dx [input] determines whether OpenDX must be launched (as a
   *     subprocess) to view the field; if set to \c false , only input files
   *     for future usage of OpenDX are created 
   * @param nx [input] number of points in the x direction (uniform sampling)   
   * @param ny [input] number of points in the y direction (uniform sampling)   
   * @param nz [input] number of points in the z direction (uniform sampling)   
   *
   */
    void visu_box(double xmin, double xmax, double ymin, double ymax,
        double zmin, double zmax, const char* title0 = 0x0, 
        const char* filename0 = 0x0, bool start_dx = true, int nx = 40, int ny = 40, 
        int nz = 40) const ;      
        


  // Member arithmetics
  // ------------------
 public:
  void operator+=(const Scalar &) ;		    ///< += Scalar
  void operator-=(const Scalar &) ;		    ///< -= Scalar
  void operator*=(const Scalar &) ;		    ///< *= Scalar

  // Manipulation of spectral bases
  // ------------------------------    
  /** Sets the spectral bases of the \c Valeur \c va  to the standard ones 
   *  for a scalar field
   */
  virtual void std_spectral_base() ;	 
  
 /** Sets the spectral bases of the \c Valeur \c va  to the standard odd ones 
   *  for a scalar field
   */
  virtual void std_spectral_base_odd() ;	 
  
  /** Sets the spectral bases of the \c Valeur \c va  
   */
  void set_spectral_base(const Base_val& ) ;	 

  /// Returns the spectral bases of the \c Valeur \c va  
  const Base_val& get_spectral_base( ) const {return va.base ;} ;	 

  /** Modifies the \c dzpuis  flag.
   *  NB: this method does not change the field values stored in
   *  the compactified external domain (use methods \c dec_dzpuis() ,
   *  etc... for this purpose).  
   */
  void set_dzpuis(int ) ; 

  /** Asymptotic expansion at r = infinity. 
   * 
   *  Determines the coefficients \f$a_k(\theta, \phi)\f$ of the expansion
   *  \f[
   *	\sum_{k=0}^n {a_k(\theta, \phi) \over r^k}
   *  \f]
   *  of \c *this  when \f$r \rightarrow \infty\f$. 
   *
   *	@param n order of the expansion
   *	@param flag : output
   *	@return Array of \c n +1 \c Valeur s on \c mg->angu  
   *		describing the coefficients \f$a_k(\theta, \phi)\f$. 
   *		This array is allocated by the routine. 
   * 
   */
  Valeur** asymptot(int n, const int flag = 0) const ; 
	

  // PDE resolution 
  // --------------
 public:
  /** Solves the scalar Poisson equation with \c *this  as a source.
   *   The source \f$\sigma\f$ of the equation \f$\Delta u = \sigma\f$ is 
   *   represented by the \c Scalar \c *this . 
   *   Note that \c dzpuis  must be equal to 2 or 4, i.e. that the
   *   quantity stored in \c *this  is in fact \f$r^2 \sigma\f$ or
   *   \f$r^4 \sigma\f$ in the compactified external domain. 
   *   The solution \e u  with the boundary condition \e u =0 at spatial
   *   infinity is the returned \c Scalar. 
   */
  Scalar poisson() const ;

  /** Solves the scalar Poisson equation with \c *this  as a source
   *   (version with parameters to control the resolution).
   *   The source \f$\sigma\f$ of the equation \f$\Delta u = \sigma\f$ is 
   *   represented by the \c Scalar \c *this . 
   *   Note that \c dzpuis  must be equal to 2 or 4, i.e. that the
   *   quantity stored in \c *this  is in fact \f$r^2 \sigma\f$ or
   *   \f$r^4 \sigma\f$ in the compactified external domain. 
   *   @param par [input/output] possible parameters
   *   @param uu [input/output] solution \e u  with the boundary condition 
   *   \e u =0 at spatial infinity. 
   */
  void poisson(Param& par, Scalar& uu) const ;
	  
  /** Solves the scalar Poisson equation with \c *this  as a source using a real Tau method
   *   The source \f$\sigma\f$ of the equation \f$\Delta u = \sigma\f$ is 
   *   represented by the \c Scalar \c *this . 
   *   Note that \c dzpuis  must be equal to 2, 3 or 4, i.e. that the
   *   quantity stored in \c *this  is in fact \f$r^2 \sigma\f$, \f$r^3 \sigma\f$ or
   *   \f$r^4 \sigma\f$ in the compactified external domain. 
   *   The solution \e u  with the boundary condition \e u =0 at spatial
   *   infinity is the returned \c Scalar. 
   */
  Scalar poisson_tau() const ;
   
   /** Solves the scalar Poisson equation with \c *this  as a source using a real Tau method
   * (version with parameters to control the resolution)
   *   The source \f$\sigma\f$ of the equation \f$\Delta u = \sigma\f$ is 
   *   represented by the \c Scalar \c *this . 
   *   Note that \c dzpuis  must be equal to 2, 3 or 4, i.e. that the
   *   quantity stored in \c *this  is in fact \f$r^2 \sigma\f$, \f$r^3 \sigma\f$ or
   *   \f$r^4 \sigma\f$ in the compactified external domain. 
   *   The solution \e u  with the boundary condition \e u =0 at spatial
   *   infinity is the returned \c Scalar. 
   */
  void poisson_tau(Param& par, Scalar& uu) const ;

  /**
   * Is identicall to \c Scalar::poisson() . The regularity condition at the 
   * origin is replace by a boundary condition of the Dirichlet type.
   * 
   * @param limite [input] : angular function. The boundary condition is 
   * given by \c limite[num] .
   * @param num [input] : index of the boudary at which the condition is to 
   * be fullfilled.
   * 
   * More precisely we impose the solution is equal to \c limite[num]  at the
   * boundary between the domains \c num  and \c num+1  (the latter one being 
   * a shell).
   * 
   */
  Scalar poisson_dirichlet (const Valeur& limite, int num) const ;
	
  /**
   * Idem as \c Scalar::poisson_dirichlet , the boundary condition being on 
   * the radial derivative of the solution.
   */
  Scalar poisson_neumann   (const Valeur&, int) const ;


  /**
   * Is identicall to \c Scalar::poisson() . The regularity condition at the 
   * origin is replace by a mixed boundary condition (Dirichlet + Neumann).
   * 
   * @param limite [input] : angular function. The boundary condition is 
   * given by \c limite[num] .
   * @param num [input] : index of the boudary at which the condition is to 
   * be fullfilled.
   * @param fact_dir [input] : double in front of \f$\Psi\f$ (if \f$\Psi\f$
   * is the variable solved).
   * @param fact_neu [input] : double in front of the radial derivative
   * of \f$\Psi\f$.
   *
   * More precisely we impose \f$ fact\_dir.\Psi + fact\_neu.\frac{\partial 
   * \Phi}{\partial r}\f$ is equal to the source at the
   * boundary between the domains \c num  and \c num+1  (the latter one being 
   * a shell).
   */
  Scalar poisson_dir_neu  (const Valeur& limite , int num, 
			   double fact_dir, double fact_neu) const ;

  /**
   * Idem as \c Scalar::poisson_dirichlet , the boundary condition being on 
   * both the function and its radial derivative. The boundary condition 
   * at infinity is relaxed.
   */

  Scalar poisson_frontiere_double   (const Valeur&, const Valeur&, int) const ;

  /** Solves the scalar Poisson equation with \c *this  as a source
   *   (version with parameters to control the resolution).
   *   The source \f$\sigma\f$ of the equation \f$\Delta u = \sigma\f$ is 
   *   represented by the \c Scalar \c *this . 
   *   The regularized source
   *   \f$\sigma_{\rm regu} = \sigma - \sigma_{\rm div}\f$
   *   is constructed and solved.
   *   Note that \c dzpuis  must be equal to 2 or 4, i.e. that the
   *   quantity stored in \c *this  is in fact \f$r^2 \sigma\f$ or
   *   \f$r^4 \sigma\f$ in the compactified external domain.
   *   @param k_div [input] regularization degree of the procedure
   *   @param nzet [input] number of domains covering the star
   *   @param unsgam1 [input] parameter \f$1/(\gamma-1)\f$ where \f$\gamma\f$
   *          denotes the adiabatic index
   *   @param par [input/output] possible parameters
   @param uu [input/output] solution
   *   @param uu_regu [output] solution of the regular part of
   *          the source.
   *   @param uu_div [output] solution of the diverging part of
   *          the source.
   *   @param duu_div [output] derivative of the diverging potential.
   *   @param source_regu [output] regularized source
   *   @param source_div [output] diverging part of the source
   */
  void poisson_regular(int k_div, int nzet, double unsgam1, Param& par,
		       Scalar& uu, Scalar& uu_regu, Scalar& uu_div,
		       Tensor& duu_div,
		       Scalar& source_regu, Scalar& source_div) const ;

  /** Checks if a Poisson equation with \c *this  as a source
   *  has been correctly solved.
   * 
   *  @param uu [input] Solution \e u  of the Poisson equation
   *		      \f$\Delta u = \sigma\f$,  \f$\sigma\f$ being 
   *		      represented by the \c Scalar \c *this .
   * 
   *  @param ostr [input/output] Output stream used for displaying
   *		\c err .
   *
   *  @param detail [input] \li if \c true  displays \c err(0,*) , 
   *		    \c err(1,*) and \c err(2,*) 
   *		\li if \c false (default),  displays only 
   *		the relative error \c err(0,*) . 
   *  
   *  @return 2-D \c Tbl  \c err  decribing the errors in each 
   *	    domain: 
   *	\li \c err(0,l) : Relative error in domain no. \c l , 
   *	    defined as the maximum value of 
   *	    \f$|\Delta u - \sigma|\f$ in that domain divided by \e M , 
   *	    where \e M  is the maximum value of \f$|\sigma|\f$ 
   *	    over all domains if \c dzpuis = 0  or \f$\sigma\f$ is
   *	    zero in the compactified external domain (CED). If 
   *	    \c dzpuis != 0  and \f$\sigma\f$ does not vanish in the 
   *	    CED, the value of \e M  used in the
   *	    non-compactified domains is the maximum value over
   *	    these domains, whereas the value of \e M  used in the
   *	    compactified external domain is the maximum value
   *	    on that particular domain. 
   *	\li \c err(1,l) :    Maximum value of the absolute error
   *			\f$|\Delta u - \sigma|\f$ in domain no. \c l  
   *	\li \c err(2,l) :    Maximum value of \f$|\sigma|\f$ in domain 
   *			    no. \c l  
   */
  Tbl test_poisson(const Scalar& uu, ostream& ostr, 
		   bool detail = false) const ;  

	/** Solves the (generalized) angular Poisson equation with \c *this  
         * as source. 
	 * The generalized angular Poisson equation is 
         * \f$\Delta_{\theta\varphi} u + \lambda u = \sigma\f$,
	 * where \f$\Delta_{\theta\varphi} u := \frac{\partial^2 u}
	 *  {\partial \theta^2} + \frac{1}{\tan \theta} \frac{\partial u}
	 *  {\partial \theta} +\frac{1}{\sin^2 \theta}\frac{\partial^2 u}
	 *  {\partial \varphi^2}\f$.
	 * 
         *    @param lambda [input] coefficient \f$\lambda\f$ in the above equation
         *      (default value = 0)
	 *   @return solution \e u . 
	 */
	Scalar poisson_angu(double lambda =0) const ;

  /** Performs one time-step integration (from \f$t=J \to J+1\f$) of the 
   *   scalar d'Alembert equation with \c *this  being the value of 
   *   the function \e f  at time \e J .
   *
   *   Works only with an affine mapping (class \c Map_af ) and,
   *   if the last domain is a compactified one, it simply copies
   *   the value of the field in this last domain at the time-step \e J 
   *   to the last domain of the returned solution.
   *   @param par [input/output] possible parameters to control the
   *   resolution of the d'Alembert equation: 
   *   \li \c par.get_double(0)  : [input] the time step \e dt ,
   *   \li \c par.get_int(0)  : [input] the type of boundary conditions
   *   set at the outer boundary (0 : reflexion, 1 : Sommerfeld 
   *   outgoing wave, valid only for \e l=0  components, 2 : Bayliss 
   *   \& Turkel outgoing wave, valid for \e l=0, 1, 2  components)
   *   \li \c par.get_int_mod(0)  : [input/output] set to 0 at first
   *   call, is used as a working flag after (must not be modified after
   *   first call)
   *   \li \c par.get_tensor_mod(0)  : [input] (optional) if the wave 
   *   equation is on a curved space-time, this is the potential in front
   *   of the Laplace operator. It has to be a Scalar and updated at 
   *    every time-step (for a potential depending on time).
   *   Note: there are many other working objects attached to this
   *   \c Param , so one should not modify it.
   *   There should be exactly one \c Param  for each wave equation to be 
   *   solved. 
   *   @param fJm1 [input] solution \f$f^{J-1}\f$ at time \e J-1 
   *   @param source [input] source \f$\sigma\f$ of the d'Alembert equation 
   *	    \f$\diamond u = \sigma\f$.
   *   @return solution \f$f^{J+1}\f$ at time \e J+1
   *   with boundary conditions defined by \c par.get_int(0) .
   */
  Scalar avance_dalembert(Param& par, const Scalar& fJm1, const Scalar& source) 
    const ;

 /**
   * Resolution of a general elliptic equation, putting zero at infinity.
   * @param params [input] the operators and variables to be used.
   **/
  Scalar sol_elliptic(Param_elliptic& params) const ;
 
/**
   * Resolution of a general elliptic equation, putting zero at infinity
   * and with inner boundary conditions.
   * @param params [input] the operators and variables to be used.
   * @param bound [input] : the boundary condition
   * @param fact_dir : 1 Dirchlet condition, 0 Neumann condition
   * @param fact_neu : 0 Dirchlet condition, 1 Neumann condition

   **/
  Scalar sol_elliptic_boundary(Param_elliptic& params, const Mtbl_cf& bound,
			       double fact_dir, double fact_neu) const ;
 

  /** Solves the scalar 2-dimensional elliptic equation 
   *   with \c *this  as a source.
   *   Note that \c dzpuis  must be equal to 2, 3 or 4, i.e. 
   *   The solution \e u  with the boundary condition \e u =0 at spatial
   *   infinity is the returned \c Scalar. 
   */
  Scalar sol_elliptic_2d(Param_elliptic&) const ;
 
  /** Solves a pseudo-1d elliptic equation 
   *   with \c *this  as a source.
   *  
   */
  Scalar sol_elliptic_pseudo_1d(Param_elliptic&) const ;

   /**
   * Resolution of a general elliptic equation, putting a given value 
   * at the outermost 
   * shell and not solving in the compactified domain.
   * @param params [input] the operators and variables to be used.
   * @param val [input] value at the last shell.
   **/
  Scalar sol_elliptic_no_zec(Param_elliptic& params, double val = 0) const ;
  
  /**
   * Resolution of a general elliptic equation solving in the 
   * compactified domain and putting a given value at the inner boundary.
   * @param params [input] the operators and variables to be used.
   * @param val [input] value at the inner boundary of the external domain.
   **/
  Scalar sol_elliptic_only_zec(Param_elliptic& params, double val) const ;

  /**
   * General elliptic solver.
   * The equation is not solved in the compactified domain and the 
   * matching is done with a solution of the type 
   * \f$\frac{\sin \left( f r + \pi/2 \right)}{r}\f$
   * @param params [input] the operators and variables to be used.
   * @param freq [input] : the frequency \f$f\f$.
   * @param coef [output] : coefficient of the oscillatory solution 
   * in the external domain.
   * @param phase [output] : phase that minimizes \c coef .
   **/
  Scalar sol_elliptic_sin_zec(Param_elliptic& params, double freq,
			      double& coef) const ;

   
  /**
   * Resolution of a general elliptic equation fixing the dericative at 
   * the origin and relaxing one continuity  condition.
   * 
   * @param val [input] value of the derivative.
   * @param params [input] the operators and variables to be used.
   **/
  Scalar sol_elliptic_fixe_der_zero(double val, 
				    Param_elliptic& params) const ;
  

  /**
   * Resolution of a divergence-like equation.
   * The equation solved reads: \f$ \frac{\partial \phi}{\partial r} + \frac{n}{r} 
   * \phi = \sigma\f$, with \f$\phi\f$ the unknown and \f$\sigma\f$ the source
   * represented by \c this.
   *
   * @param n [input] the coefficient in front of the 1/r term.
   *
   * @return the solution to the equation.
   */
  Scalar sol_divergence(int n) const ;
	
  // Import from other mapping 
  // -------------------------

  /** Assignment to another \c Scalar defined on a different mapping.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import(const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping.
   *  Case where the \c Scalar is symmetric with respect to the plane y=0.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_symy(const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping.
   *  Case where the \c Scalar is antisymmetric with respect to the 
   *  plane y=0.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_asymy(const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import(int nzet, const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping.
   *  Case where the \c Scalar is symmetric with respect to the plane y=0.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_symy(int nzet, const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping.
   *  Case where the \c Scalar is antisymmetric with respect to the 
   *  plane y=0.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_asymy(int nzet, const Scalar& ci) ;	 

 protected:
  /** Assignment to another \c Scalar defined on a different mapping,
   *  when the two mappings do not have a particular relative orientation.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_gal(int nzet, const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping,
   *  when the two mappings have aligned Cartesian axis. 
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_align(int nzet, const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping,
   *  when the two mappings have anti-aligned Cartesian axis (i.e.
   *  \f$x_1 = - x_2\f$,  \f$y_1 = - y_2\f$,  \f$z_1 = z_2\f$). 
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_anti(int nzet, const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping,
   *  when the two mappings have aligned Cartesian axis. 
   *  Case where the \c Scalar is symmetric with respect to the plane y=0.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_align_symy(int nzet, const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping,
   *  when the two mappings have anti-aligned Cartesian axis (i.e.
   *  \f$x_1 = - x_2\f$,  \f$y_1 = - y_2\f$,  \f$z_1 = z_2\f$). 
   *  Case where the \c Scalar is symmetric with respect to the plane y=0.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_anti_symy(int nzet, const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping,
   *  when the two mappings have aligned Cartesian axis. 
   *  Case where the \c Scalar is antisymmetric with respect to the 
   *  plane y=0.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_align_asymy(int nzet, const Scalar& ci) ;	 

  /** Assignment to another \c Scalar defined on a different mapping,
   *  when the two mappings have anti-aligned Cartesian axis (i.e.
   *  \f$x_1 = - x_2\f$,  \f$y_1 = - y_2\f$,  \f$z_1 = z_2\f$). 
   *  Case where the \c Scalar is antisymmetric with respect to the 
   *  plane y=0.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original \c Scalar. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. \c this->mp ) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and \c nzet-1 . In the other
   *			    domains, \c *this  is set to zero. 
   *	@param ci [input] \c Scalar to be imported.
   */
  void import_anti_asymy(int nzet, const Scalar& ci) ;	 
	

  friend Scalar operator-(const Scalar& ) ;			
  friend Scalar operator+(const Scalar&, const Scalar &) ;	
  friend Scalar operator+(const Scalar&, const Mtbl&) ;	
  friend Scalar operator+(const Scalar&, double ) ;		
  friend Scalar operator-(const Scalar &, const Scalar &) ;
  friend Scalar operator-(const Scalar&, const Mtbl&) ;	
  friend Scalar operator-(const Scalar&, double ) ;		
  friend Scalar operator*(const Scalar &, const Scalar &) ;
  friend Scalar operator%(const Scalar &, const Scalar &) ;
  friend Scalar operator|(const Scalar &, const Scalar &) ;
  friend Scalar operator*(const Mtbl&, const Scalar &) ;		
  friend Scalar operator*(double, const Scalar &) ;		
  friend Scalar operator/(const Scalar &, const Scalar &) ;
  friend Scalar operator/(const Scalar &, const Mtbl &) ;
  friend Scalar operator/(const Mtbl &, const Scalar &) ;
  friend Scalar operator/(const Scalar&, double ) ;	       
  friend Scalar operator/(double, const Scalar &) ;

  friend Scalar sin(const Scalar& ) ;
  friend Scalar cos(const Scalar& ) ;
  friend Scalar tan(const Scalar& ) ;
  friend Scalar asin(const Scalar& ) ;
  friend Scalar acos(const Scalar& ) ;
  friend Scalar atan(const Scalar& ) ;
  friend Scalar exp(const Scalar& ) ;	
  friend Scalar log(const Scalar& ) ;	
  friend Scalar log10(const Scalar& ) ;	
  friend Scalar sqrt(const Scalar& ) ;	
  friend Scalar racine_cubique (const Scalar& ) ;
  friend Scalar pow(const Scalar& , int ) ;	
  friend Scalar pow(const Scalar& , double ) ; 
  friend Scalar abs(const Scalar& ) ;	

  friend Tbl max(const Scalar& ) ;   
  friend Tbl min(const Scalar& ) ;   
  friend Tbl norme(const Scalar& ) ;   
  friend Tbl diffrel(const Scalar& a, const Scalar& b) ; 
  friend Tbl diffrelmax(const Scalar& a, const Scalar& b) ; 

};

ostream& operator<<(ostream& , const Scalar & ) ;	

// Prototypage de l'arithmetique
/**
 * \defgroup sal_mat Scalar mathematics
 * \ingroup (tensor)
 * @{
 */

Scalar operator+(const Scalar& ) ;			///< + Scalar
Scalar operator-(const Scalar& ) ;			///< \c - Scalar
Scalar operator+(const Scalar&, const Scalar &) ;	///< Scalar + Scalar
Scalar operator+(const Scalar&, const Mtbl&) ;	///< Scalar + Mbtl
Scalar operator+(const Mtbl&, const Scalar&) ;	///< Mtbl + Scalar
Scalar operator+(const Scalar&, double ) ;		///< Scalar + double
Scalar operator+(double, const Scalar& ) ;		///< double + Scalar 
Scalar operator+(const Scalar&, int ) ;		///< Scalar + int
Scalar operator+(int, const Scalar& ) ;		///< int + Scalar 
Scalar operator-(const Scalar &, const Scalar &) ;	///< Scalar - Scalar
Scalar operator-(const Scalar&, const Mtbl&) ;	///< Scalar - Mbtl
Scalar operator-(const Mtbl&, const Scalar&) ;	///< Mtbl - Scalar
Scalar operator-(const Scalar&, double ) ;		///< Scalar - double
Scalar operator-(double, const Scalar& ) ;		///< double - Scalar 
Scalar operator-(const Scalar&, int ) ;		///< Scalar - int
Scalar operator-(int, const Scalar& ) ;		///< int - Scalar 
Scalar operator*(const Scalar &, const Scalar &) ;	///< Scalar * Scalar

/// Scalar * Scalar with desaliasing
Scalar operator%(const Scalar &, const Scalar &) ;	

/// Scalar * Scalar with desaliasing only in \e r
Scalar operator|(const Scalar &, const Scalar &) ;

Scalar operator*(const Mtbl&, const Scalar&) ; ///< Mtbl * Scalar
Scalar operator*(const Scalar&, const Mtbl&) ; ///< Scalar * Mtbl
	
Scalar operator*(const Scalar&, double ) ;		///< Scalar * double
Scalar operator*(double, const Scalar &) ;		///< double * Scalar
Scalar operator*(const Scalar&, int ) ;		///< Scalar * int
Scalar operator*(int, const Scalar& ) ;		///< int * Scalar 
Scalar operator/(const Scalar &, const Scalar &) ;	///< Scalar / Scalar
Scalar operator/(const Scalar&, double ) ;		///< Scalar / double
Scalar operator/(double, const Scalar &) ;		///< double / Scalar
Scalar operator/(const Scalar&, int ) ;		///< Scalar / int
Scalar operator/(int, const Scalar &) ;		///< int / Scalar
Scalar operator/(const Scalar &, const Mtbl&) ; ///< Scalar / Mtbl  
Scalar operator/(const Mtbl&, const Scalar &) ; ///< Mtbl / Scalar
	

Scalar sin(const Scalar& ) ;		///< Sine
Scalar cos(const Scalar& ) ;		///< Cosine
Scalar tan(const Scalar& ) ;		///< Tangent
Scalar asin(const Scalar& ) ;		///< Arcsine
Scalar acos(const Scalar& ) ;		///< Arccosine
Scalar atan(const Scalar& ) ;		///< Arctangent
Scalar exp(const Scalar& ) ;		///< Exponential
Scalar log(const Scalar& ) ;		///< Neperian logarithm
Scalar log10(const Scalar& ) ;	///< Basis 10 logarithm
Scalar sqrt(const Scalar& ) ;		///< Square root
Scalar racine_cubique (const Scalar& ) ;		///< Cube root
Scalar pow(const Scalar& , int ) ;	///< Power \f${\tt Scalar}^{\tt int}\f$
Scalar pow(const Scalar& , double ) ; ///< Power \f${\tt Scalar}^{\tt double}\f$
Scalar abs(const Scalar& ) ;		///< Absolute value

/**
 * Maximum values of a \c Scalar in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Scalar& ) ;   

/**
 * Minimum values of a \c Scalar in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Scalar& ) ;   

/**
 * Sums of the absolute values of all the values of the \c Scalar 
 * in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Scalar& ) ;   

/**
 * Relative difference between two \c Scalar (norme version).
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are \c norme[a(l)-b(l)]/norme[b(l)]  if \c b(l)!=0  and
 *	   \c norme[a(l)-b(l)]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c a  and \c b  
 *	   in domain no. \c l . 
 */
Tbl diffrel(const Scalar& a, const Scalar& b) ; 

/**
 * Relative difference between two \c Scalar (max version).
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are \c max[abs(a(l)-b(l))]/max[abs(b(l))]  if \c b(l)!=0  and
 *	   \c max[abs(a(l)-b(l))]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c a  and \c b  
 *	   in domain no. \c l . 
 */
Tbl diffrelmax(const Scalar& a, const Scalar& b) ; 

/** @} */
#endif
