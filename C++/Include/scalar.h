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
 * 
 * @version #$Id$#
 * 
 */

class Scalar : public Tensor {
  
  // Data : 
  // -----
 protected:
  
  /** The logical state {\tt ETATNONDEF} (undefined), {\tt ETATZERO} (null),
   *  {\tt ETATUN} (one), or {\tt ETATQCQ}(ordinary).
   */
  int etat ; 
  
  /**
   * Power of {\it r} by which the quantity represented by {\tt this} 
   * must be divided in the  compactified external domain (CED) in order 
   * to get the correct physical values
   */
  int dzpuis ;	
  
  Valeur va ;		/// The numerical value of the {\tt Scalar}    
  
  // Derived data : 
  // ------------
 protected:
  /// Pointer on $\partial/\partial r$ of {\tt *this} (0x0 if not up to date)
  mutable Scalar* p_dsdr ;	

  /** Pointer on $1/r \partial/\partial \theta$ of {\tt *this} 
   *  (0x0 if not up to date)
   */
  mutable Scalar* p_srdsdt ;	

  /** Pointer on $1/(r\sin\theta) \partial/\partial \phi$ of {\tt *this}
   *  (0x0 if not up to date)
   */
  mutable Scalar* p_srstdsdp ;

  /// Pointer on $\partial/\partial \theta$ of {\tt *this} (0x0 if not up to date)
  mutable Scalar* p_dsdt ;	

  /** Pointer on $1/\sin\theta \partial/\partial \phi$ of {\tt *this}
   *  (0x0 if not up to date)
   */
  mutable Scalar* p_stdsdp ;	
  
  /** Pointer on $\partial/\partial x$ of {\tt *this},
   *  where $x=r\sin\theta \cos\phi$ (0x0 if not up to date)
   */
  mutable Scalar* p_dsdx ;	
  
  /** Pointer on $\partial/\partial y$ of {\tt *this},
   *  where $y=r\sin\theta \sin\phi$(0x0 if not up to date)
   */
  mutable Scalar* p_dsdy ;	

  /** Pointer on $\partial/\partial z$ of {\tt *this},
   *  where $z=r\cos\theta$ (0x0 if not up to date)
   */
  mutable Scalar* p_dsdz ;	
  
  /** Pointer on the Laplacian of {\tt *this} (0x0 if not up to date)
   */
  mutable Scalar* p_lap ;	
  
  /** Pointer on the Laplacian of {\tt *this} (0x0 if not up to date)
   */
  mutable Scalar* p_lapang ;	
  
  /** Power of {\it r} by which the last computed Laplacian has been 
   *  multiplied in the compactified external domain.  
   */
  mutable int ind_lap ; 

  /** Pointer on the space integral of {\tt *this} (values in each 
   *  domain) (0x0 if not up to date)
   */
  mutable Tbl* p_integ ; 
  
  // Constructors - Destructor
  // -------------------------
  
 public:
  
  explicit Scalar(const Map& mpi) ;	/// Constructor from mapping

  /// Constructor from a Tensor (must be of valence 0)
  Scalar(const Tensor& a) ;     

  Scalar(const Scalar& a) ;		/// Copy constructor
  
  explicit Scalar(const Cmp& a) ;	/// Constructor by conversion of a Cmp
  
  /// Constructor from a file (see {\tt sauve(FILE* )})
  Scalar(const Map&, const Mg3d&, FILE* ) ;    		
  
  virtual ~Scalar() ;			/// Destructor
  
  
  // Memory management
  // -----------------
 protected:
  void del_t() ;		    /// Logical destructor
  virtual void del_deriv() const;	    /// Logical destructor of the derivatives
  void set_der_0x0() const;	    /// Sets the pointers for derivatives to 0x0
  
 public:
  
  /**
   * Sets the logical state to {\tt ETATNONDEF} (undefined). 
   * Calls the logical destructor of the {\tt Valeur va} and
   * deallocates the memory occupied by all the derivatives. 
   */
  virtual void set_etat_nondef() ;   
  
  /**
   * Sets the logical state to {\tt ETATZERO} (zero). 
   * Calls the logical destructor of the {\tt Valeur va} and
   * deallocates the memory occupied by all the derivatives. 
   */
  virtual void set_etat_zero() ;	    
  
  /**
   * Sets the logical state to {\tt ETATQCQ} (ordinary state).
   * If the state is already {\tt ETATQCQ}, this function does nothing.
   * Otherwise, it calls the logical destructor of the {\tt Valeur va} and
   * deallocates the memory occupied by all the derivatives.
   */
  virtual void set_etat_qcq() ;	    
  
  /**
   * Sets the logical state to {\tt ETATUN} (one). 
   * Fills the {\tt Valeur va} with ones and
   * deallocates the memory occupied by all the derivatives. 
   */
  void set_etat_one() ;	    
  
  /**
   * Sets the logical state to {\tt ETATQCQ} (ordinary state)
   *  and performs the memory allocation of all the 
   *  elements, down to the {\tt double} arrays of the {\tt Tbl}s. 
   *  This function performs in fact recursive calls to {\tt set\_etat\_qcq()}
   *  on each element of the chain {\tt Scalar} ->
   *  {\tt Valeur} -> {\tt Mtbl} -> {\tt Tbl}. 
   */
  virtual void allocate_all() ; 
  
  /**
   * Sets the {\tt Scalar} to zero in a hard way. 
   * 1/ Sets the logical state to {\tt ETATQCQ}, i.e. to an ordinary state.
   * 2/ Fills the {\tt Valeur va} with zeros. 
   * NB: this function must be used for debugging purposes only.
   * For other operations, the functions {\tt set\_etat\_zero()}
   * or {\tt annule(int, int)} must be perferred. 
   */
  void annule_hard() ;
  
  // Extraction of information
  // -------------------------
    public:
  /** Returns the logical state {\tt ETATNONDEF} (undefined), 
   * {\tt ETATZERO}(null) or {\tt ETATQCQ}(ordinary).
   */
  int get_etat() const {return etat;} ; 
  
  int get_dzpuis() const {return dzpuis;} ; /// Returns {\tt dzpuis}
  
  /** Returns {\tt true} if the last domain is compactified and
   *  {\tt *this} is not zero in this domain
   */
  bool dz_nonzero() const ; 
	
  /** Returns {\tt false} if the last domain is compactified 
   *  and {\tt *this} is not zero in this domain and {\tt dzpuis}
   *  is not equal to {\tt dzi}, otherwise return true. 
   */
  bool check_dzpuis(int dzi) const ; 
  
  // Assignment
  // -----------
 public: 
  /// Assignment to another {\tt Scalar} defined on the same mapping
  void operator=(const Scalar& a) ;	
  
  /// Assignment to a {\tt Tensor} (of valence 0)
  virtual void operator=(const Tensor& a) ; 

  void operator=(const Cmp& a) ; 	/// Assignment to a {\tt Cmp}
  void operator=(const Valeur& a) ; /// Assignment to a {\tt Valeur}
  void operator=(const Mtbl& a) ;	 /// Assignment to a {\tt Mtbl}
  void operator=(double ) ;	 /// Assignment to a {\tt double}
  void operator=(int ) ;		 /// Assignment to an {\tt int}
  
  // Access to individual elements
  // -----------------------------
    public:
  
  /// Returns {\tt va} (read only version)
  const Valeur& get_spectral_va() const {return va;} ; 
  
  /// Returns {\tt va} (read/write version)
  Valeur& set_spectral_va() {return va;} ; 
  
  /** Read/write of the value in a given domain.
   * CAUTION: to gain in efficiency, the method {\tt del\_deriv()} (to delete
   *     the derived members) is not called by this function. It must
   *     thus be invoqued by the user.  
   *
   * @param l [input] domain index
   * @return Tbl containing the value of the field in domain {\tt l}.
   */ 
  Tbl& set_domain(int l) {
    assert(etat == ETATQCQ) ;
    return va.set(l) ;
  };
  
  /** Read-only of the value in a given domain.
   * @param l [input] domain index
   * @return Tbl containing the value of the field in domain {\tt l}.
   */ 
  const Tbl& domain(int l) const {
    assert( (etat == ETATQCQ) || (etat == ETATUN) ) ;
    return va(l) ;
  };
  
  
  /** Returns the value of the field at a specified grid point.
   * @param l [input] domain index
   * @param k [input] $\phi$ index
   * @param j [input] $\theta$ index
   * @param i [input] {\it r} ($\xi$) index
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
   *   arbitrary point $(r, \theta, \phi)$, by means of the spectral 
   *   expansion.
   *   NB: if $(r, \theta, \phi)$ is a point of the spectral grid, 
   *     the method {\tt val\_grid\_point} is to be preferred, 
   *     being much more efficient. 
   *	 @param r [input] value of the coordinate {\it r}
   *	 @param theta [input] value of the coordinate $\theta$
   *	 @param phi [input] value of the coordinate $\phi$
   *	 @return value at the point $(r, \theta, \phi)$ 
   *		 of the field represented by {\tt *this}. 
   */
  double val_point(double r, double theta, double phi) const ; 
  
  
  /** Setting the value of the field at a given grid point.
   * CAUTION: to gain in efficiency (especially when this method is
   *  invoqued inside a loop), the method {\tt del\_deriv()} (to delete
   *     the derived members) is not called by {\tt set\_grid\_point}. 
   *     It must thus be invoqued by the user, after all the calls
   *     to  {\tt set\_grid\_point} have been performed.   
   *     
   * @param l [input] domain index
   * @param k [input] $\phi$ index
   * @param j [input] $\theta$ index
   * @param i [input] {\it r} ($\xi$) index
   * @return writable value of the field at the specified grid point
   */ 
  double& set_grid_point(int l, int k, int j, int i) {
    assert(etat == ETATQCQ) ;
    return va.set(l, k, j, i) ;
  };
  
	
  /**
   * Sets the {\tt Scalar} to zero in several domains.
   *	@param l_min [input] The {\tt Scalar} will be set (logically) to zero
   *			     in the domains whose indices are in the range
   *			     {\tt [l\_min, l\_max]}.
   *	@param l_max [input] see the comments for {\tt l\_min}.
   * 
   * Note that {\tt annule(0, va.mg->get\_nzone()-1)} is equivalent to
   *	 {\tt set\_etat\_zero()}.
   */
  virtual void annule(int l_min, int l_max) ; 
  
  /** Sets the value of the {\tt Scalar} at the inner boundary of a given 
   * domain. 
   * @param l [input] domain index
   * @param x [input] (constant) value at the inner boundary of domain no. {\tt l}
   */
  void set_inner_boundary(int l, double x) ;
    
  /** Sets the value of the {\tt Scalar} at the outer boundary of a given 
   * domain. 
   * @param l [input] domain index
   * @param x [input] (constant) value at the outer boundary of domain no. {\tt l}
   */
  void set_outer_boundary(int l, double x) ;

  /**
   * Gives the spectrum in terms of multipolar modes {\it l}.
   *  @return a {\tt Tbl} of size (nzone, lmax), where lmax is the
   *  maximal multipolar momentum over all domains. The {\it l}-th
   *  element contains the L1 norm of the {\it l}-th multipole 
   *  ({\it i.e.} a sum over all {\it m} of the norms (coefficient space)
   *  of the component of a given $Y_l^m$.
   */
  Tbl multipole_spectrum () ;
  
  // Differential operators and others
  // ---------------------------------
 public:
  /** Returns $\partial / \partial r$ of {\tt *this}.
   *  If {\tt dzpuis} is zero, then the returned {\tt Scalar} has 
   *  {\tt dzpuis} = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdr() const ; 
  
  /** Returns $1/r \partial / \partial \theta$ of {\tt *this}.
   *  If {\tt dzpuis} is zero, then the returned {\tt Scalar} has 
   *  {\tt dzpuis} = 2. It is increased by 1 otherwise.
   */
  const Scalar& srdsdt() const ; 
  
  /** Returns $1/(r\sin\theta) \partial / \partial \phi$ of {\tt *this}.
   *  If {\tt dzpuis} is zero, then the returned {\tt Scalar} has 
   *  {\tt dzpuis} = 2. It is increased by 1 otherwise.
   */
  const Scalar& srstdsdp() const ; 
  
  /** Returns $\partial / \partial \theta$ of {\tt *this}.
   */
  const Scalar& dsdt() const ; 
  
  /** Returns $1/\sin\theta \partial / \partial \phi$ of {\tt *this}.
   */
  const Scalar& stdsdp() const ; 
  
  /** Returns $\partial/\partial x$ of {\tt *this},
   *  where $x=r\sin\theta \cos\phi$.
   *  If {\tt dzpuis} is zero, then the returned {\tt Scalar} has 
   *  {\tt dzpuis} = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdx() const ;	
  
  /** Returns $\partial/\partial y$ of {\tt *this},
   *  where $y=r\sin\theta \sin\phi$.
   *  If {\tt dzpuis} is zero, then the returned {\tt Scalar} has 
   *  {\tt dzpuis} = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdy() const ;	
  
  /** Returns $\partial/\partial z$ of {\tt *this},
   *  where $z=r\cos\theta$.
   *  If {\tt dzpuis} is zero, then the returned {\tt Scalar} has 
   *  {\tt dzpuis} = 2. It is increased by 1 otherwise.
   */
  const Scalar& dsdz() const ;	
  
  /** Returns $\partial/\partial x_i$ of {\tt *this},
   *  where $x_i = (x, y, z)$.
   *  If {\tt dzpuis} is zero, then the returned {\tt Scalar} has 
   *  {\tt dzpuis} = 2. It is increased by 1 otherwise.
   *  @param i [input] i=1 for {\it x},  i=2 for {\it y}, i=3 for {\it z}.
   */
  const Scalar& deriv(int i) const ;	
  
    /** Returns the gradient (1-form = covariant vector) of {\tt *this} 
     *  @param gam metric components only used to get the triad with 
     *    respect to which the components of the result are defined        
     */
    const Vector& derive_cov(const Metric& gam) const ; 


    /** Returns the "contravariant" derivative of {\tt *this} with respect 
     * to some metric $\gamma$, by raising the index of the
     * gradient (cf. method {\tt derive\_cov()}) with 
     * $\gamma$.
     */
    const Vector& derive_con(const Metric& gam) const ; 

    /// Computes the derivative of {\tt this} along a vector field {\tt v}
    Scalar derive_lie(const Vector& v) const ; 


  /** Returns the Laplacian of {\tt *this}
   *   @param ced_mult_r [input] Determines the quantity computed in
   *			 the  compactified external domain (CED) 
   *		({\it u} in the field represented by {\tt *this}) :  \\
   *		    ced\_mult\_r = 0 : $\Delta u$	\\
   *		    ced\_mult\_r = 2 : $r^2 \,  \Delta u$	\\
   *		    ced\_mult\_r = 4 (default) : $r^4 \, \Delta u$	
   */
  const Scalar& laplacian(int ced_mult_r = 4) const ; 
  
  /** Returns the angular Laplacian $\Delta_{\theta\varphi}$ of {\tt *this},
   *  where $\Delta_{\theta\varphi} f = \frac{\partial^2 f}
   *  {\partial \theta^2} + \frac{1}{\tan \theta} \frac{\partial f}
   *  {\partial \theta} +\frac{1}{\sin^2 \theta}\frac{\partial^2 f}
   *  {\partial \varphi^2}$
   * 
   */
  const Scalar& lapang() const ; 
  
  /// Division by {\it r} everywhere; {\tt dzpuis} is not changed.
  void div_r() ;    
 
  /** Division by {\it r} everywhere but with the output flag {\tt dzpuis} 
   *  set to {\tt ced\_mult\_r}.
   *  @param  ced_mult_r [input] value of {\tt dzpuis} of the result.
   */
  void div_r_dzpuis(int ced_mult_r) ; 
  
  /** Division by {\it r} in the compactified external domain (CED), the 
   * {\tt dzpuis} flag is not changed.
   */
  void div_r_ced() ;

  /// Multiplication by {\it r} everywhere; {\tt dzpuis} is not changed.
  void mult_r() ;  
  
  /** Multiplication by {\it r} everywhere but with the output flag {\tt dzpuis} 
   *  set to {\tt ced\_mult\_r}.
   *  @param  ced_mult_r [input] value of {\tt dzpuis} of the result. 
   */
  void mult_r_dzpuis(int ced_mult_r) ; 
  
  /** Multiplication by {\it r} in the compactified external domain (CED), the 
   * {\tt dzpuis} flag is not changed.
   */
  void mult_r_ced() ;
  
  /// Multiplication by $r\sin\theta$ everywhere; {\tt dzpuis} is not changed.
  void mult_rsint() ;   
  
  /** Multiplication by $r\sin\theta$ but with the output flag {\tt dzpuis}
   *  set to {\tt ced\_mult\_r}.
   *  @param  ced_mult_r [input] value of {\tt dzpuis} of the result. 
   */
  void mult_rsint_dzpuis(int ced_mult_r) ; 

  /// Division by $r\sin\theta$ everywhere; {\tt dzpuis} is not changed.
  void div_rsint() ;    
  
  /** Division by $r\sin\theta$ but with the output flag {\tt dzpuis}
   *  set to {\tt ced\_mult\_r}.
   *  @param  ced_mult_r [input] value of {\tt dzpuis} of the result. 
   */
  void div_rsint_dzpuis(int ced_mult_r) ; 
  
  void mult_cost() ;   /// Multiplication by $\cos\theta$

  void mult_sint() ;   /// Multiplication by $\sin\theta$

  void div_sint() ;    /// Division by $\sin\theta$

  void div_tant() ;    /// Division by $\tan\theta$
  
  /** Computes the integral over all space of {\tt *this}.
   *  The computed quantity is ({\it u} being the field represented by
   *   {\tt *this})
   *    $\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi$.
   *  Note that in the compactified external domain (CED), {\tt dzpuis} 
   *  must be 4 for the computation to take place. 
   */
  double integrale() const ; 
	
  /** Computes the integral in each domain of {\tt *this}.
   *  The computed quantity is ({\it u} being the field represented by
   *   {\tt *this})
   *    $\int u \, r^2 \sin\theta \,  dr\, d\theta \, d\phi$
   *  in each domain. The result is returned a {\tt Tbl} on the 
   *  various domains. 
   *  Note that in the compactified external domain (CED), {\tt dzpuis} 
   *  must be 4 for the computation to take place. 
   */
  const Tbl& integrale_domains() const ; 
	
  /** Decreases by {\tt dec} units the value of {\tt dzpuis} and 
   *  changes accordingly the values of the {\tt Scalar} in the 
   *  compactified external domain (CED).
   */
  virtual void dec_dzpuis(int dec = 1) ; 

  /** Increases by {\tt inc} units the value of {\tt dzpuis} and 
   *  changes accordingly the values of the {\tt Scalar} in the 
   *  compactified external domain (CED).
   */
  virtual void inc_dzpuis(int inc = 1) ; 
	
  /** Sets a new vectorial basis (triad) of decomposition and modifies
   *  the components accordingly. 
   */
  virtual void change_triad(const Base_vect& new_triad) ; 
    
  /**
   * Sets the {\tt n} lasts coefficients in {\it r} to 0 in the 
   *  external domain.
   */
  void filtre (int n) ;
    
  /**
   * Sets the {\tt n} lasts coefficients in $\Phi$ to 0 in the 
   * domain {\tt zone}.
   */
  void filtre_phi (int n, int zone) ;
    
  /**
   * Sets the {\tt n} lasts coefficients in $\Phi$ to 0 in all domains
   */
  void filtre_phi (int n) ;

   /**
   * Sets the {\tt n} lasts coefficients in $\Theta$ to 0 in all domains
   */
  void filtre_theta (int n) ;
   
  /**
   * Substracts all the components behaving like $r^{-n}$ in the external 
   * domain, with {\it n} strictly lower than {\tt puis}, so that {\tt *this} 
   * decreases at least like $r^{\tt puis}$ at infinity.
   */
  void fixe_decroissance (int puis) ;

    /** Performs a $C^k$ matching of the last non-compactified shell with
     * a decaying function $\sum_{j=0}^k {\alpha_j \over r^{\ell+n+j}$ where
     * $\ell$ is the spherical harmonic index and {\it n} is some 
     * specifiable parameter. 
     */
    void smooth_decay(int k, int n) ; 

  /**
   * Performs the $C^n$ matching of the nucleus with respect to the 
   * first shell.
   */
  void raccord(int n) ;
	
  /**
   * Performs the $C^1$ matching of the external domain with respect to
   * the last shell using function like $\frac{1}{r^i}$ with 
   * ${\tt puis} \leq i \leq {\tt puis+nbre}$ for each spherical harmonics 
   * with $l \leq {\tt lmax}$.
   */
  void raccord_c1_zec(int puis, int nbre, int lmax) ;

  /**
   * Matching of the external domain with the outermost shell
   */
  void raccord_externe(int puis, int nbre, int lmax) ;

  // Outputs
  // -------
 public:
  virtual void sauve(FILE *) const ;	    /// Save in a file
    
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
   * @param section_type [input] defines the type of section : \\
   *    'x' for a plane x = a with a = const (parameter {\tt aa}) \\
   *    'y' for a plane y = a with a = const (parameter {\tt aa})\\
   *    'z' for a plane z = a with a = const (parameter {\tt aa})
   * @param aa [input] constant a defining the section plane
   * @param umin [input] defines with {\tt umax} the range of the plane coordinate u 
   * @param umax [input] defines with {\tt umin} the range of the plane coordinate u 
   * @param vmin [input] defines with {\tt vmax} the range of the plane coordinate v 
   * @param vmax [input] defines with {\tt vmin} the range of the plane coordinate v 
   * @param title [input] title for the graph (for OpenDX legend)
   * @param filename [input] name for the file which will be the input for 
   *    OpenDX; the default 0x0 is transformed into "scalar_section"
   * @param start_dx [input] determines whether OpenDX must be launched (as a
   *     subprocess) to view the field; if set to {\tt false}, only input files
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
   * @param plane [input] : 2D {\tt Tbl} defining the section plane: {\tt plane}
   *    must of dimension 3x3 with the following content: \\
   *    {\tt plane(0,i)}: absolute Cartesian coordinates (xa0,ya0,za0) of some
   *    point in the plane considered as the origin for the plane coordinates
   *    (u,v): {\tt plane(0,0) = xa0}, {\tt plane(0,1) = ya0}, 
   *    {\tt plane(0,2) = za0}  \\
   *    {\tt plane(1,i)}: components w.r.t. absolute Cartesian coordinates 
   *        of the u-coordinate unit vector in the section plane \\
   *    {\tt plane(2,i)}: components w.r.t. absolute Cartesian coordinates 
   *        of the v-coordinate unit vector in the section plane
   * @param umin [input] defines with {\tt umax} the range of the plane coordinate u 
   * @param umax [input] defines with {\tt umin} the range of the plane coordinate u 
   * @param vmin [input] defines with {\tt vmax} the range of the plane coordinate v 
   * @param vmax [input] defines with {\tt vmin} the range of the plane coordinate v 
   * @param title [input] title for the graph (for OpenDX legend)
   * @param filename [input] name for the file which will be the input for 
   *    OpenDX; the default 0x0 is transformed into "scalar_section"
   * @param start_dx [input] determines whether OpenDX must be launched (as a
   *     subprocess) to view the field; if set to {\tt false}, only input files
   *     for future usage of OpenDX are created 
   * @param nu [input] number of points in the u direction (uniform sampling)   
   * @param nv [input] number of points in the v direction (uniform sampling)   
   *
   */
    void visu_section(const Tbl& plane, double umin, double umax, double vmin,
        double vmax, const char* title = 0x0, const char* filename = 0x0,
        bool start_dx = true, int nu = 200, int nv = 200) const ;   

  /** 3D visualization (volume rendering) via OpenDX.
   * Prepares files for visualization by OpenDX of the values of the field in
   * some rectangular box.
   *
   * @param xmin [input] defines with {\tt xmax} the x range of the visualization box 
   * @param xmax [input] defines with {\tt xmin} the x range of the visualization box 
   * @param ymin [input] defines with {\tt ymax} the y range of the visualization box 
   * @param ymax [input] defines with {\tt ymin} the y range of the visualization box 
   * @param zmin [input] defines with {\tt zmax} the z range of the visualization box 
   * @param zmax [input] defines with {\tt zmin} the z range of the visualization box 
   * @param title [input] title for the graph (for OpenDX legend)
   * @param filename [input] name for the file which will be the input for 
   *    OpenDX; the default 0x0 is transformed into "scalar_box"
   * @param start_dx [input] determines whether OpenDX must be launched (as a
   *     subprocess) to view the field; if set to {\tt false}, only input files
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
  void operator+=(const Scalar &) ;		    /// += Scalar
  void operator-=(const Scalar &) ;		    /// -= Scalar
  void operator*=(const Scalar &) ;		    /// *= Scalar

  // Manipulation of spectral bases
  // ------------------------------    
  /** Sets the spectral bases of the {\tt Valeur va} to the standard ones 
   *  for a scalar field
   */
  virtual void std_spectral_base() ;	 

  /** Sets the spectral bases of the {\tt Valeur va} 
   */
  void set_spectral_base(const Base_val& ) ;	 

  /** Modifies the {\tt dzpuis} flag.
   *  NB: this method does not change the field values stored in
   *  the compactified external domain (use methods {\tt dec\_dzpuis()},
   *  etc... for this purpose).  
   */
  void set_dzpuis(int ) ; 

  /** Asymptotic expansion at r = infinity. 
   * 
   *  Determines the coefficients $a_k(\theta, \phi)$ of the expansion
   *  \begin{equation}
   *	\sum_{k=0}^n {a_k(\theta, \phi) \over r^k}
   *  \end{equation} 
   *  of {\tt *this} when $r \rightarrow \infty$. 
   *
   *	@param n order of the expansion
   *	@param flag : output
   *	@return Array of {\tt n}+1 {\tt Valeur}s on {\tt mg->angu} 
   *		describing the coefficients $a_k(\theta, \phi)$. 
   *		This array is allocated by the routine. 
   * 
   */
  Valeur** asymptot(int n, const int flag = 0) const ; 
	

  // PDE resolution 
  // --------------
 public:
  /** Solves the scalar Poisson equation with {\tt *this} as a source.
   *   The source $\sigma$ of the equation $\Delta u = \sigma$ is 
   *   represented by the {\tt Scalar} {\tt *this}. 
   *   Note that {\tt dzpuis} must be equal to 2 or 4, i.e. that the
   *   quantity stored in {\tt *this} is in fact $r^2 \sigma$ or
   *   $r^4 \sigma$ in the compactified external domain. 
   *   The solution {\it u} with the boundary condition {\it u}=0 at spatial
   *   infinity is the returned {\tt Scalar}. 
   */
  Scalar poisson() const ;

  /** Solves the scalar Poisson equation with {\tt *this} as a source
   *   (version with parameters to control the resolution).
   *   The source $\sigma$ of the equation $\Delta u = \sigma$ is 
   *   represented by the {\tt Scalar} {\tt *this}. 
   *   Note that {\tt dzpuis} must be equal to 2 or 4, i.e. that the
   *   quantity stored in {\tt *this} is in fact $r^2 \sigma$ or
   *   $r^4 \sigma$ in the compactified external domain. 
   *   @param par [input/output] possible parameters
   *   @param uu [input/output] solution {\it u} with the boundary condition 
   *   {\it u}=0 at spatial infinity. 
   */
  void poisson(Param& par, Scalar& uu) const ;
	
  /**
   * Is identicall to {\tt Scalar::poisson()}. The regularity condition at the 
   * origin is replace by a boundary condition of the Dirichlet type.
   * 
   * @param limite [input] : angular function. The boundary condition is 
   * given by {\tt limite[num]}.
   * @param num [input] : index of the boudary at which the condition is to 
   * be fullfilled.
   * 
   * More precisely we impose the solution is equal to {\tt limite[num]} at the
   * boundary between the domains {\tt num} and {\tt num+1} (the latter one being 
   * a shell).
   * 
   */
  Scalar poisson_dirichlet (const Valeur& limite, int num) const ;
	
  /**
   * Idem as {\tt Scalar::poisson\_dirichlet}, the boundary condition being on 
   * the radial derivative of the solution.
   */
  Scalar poisson_neumann   (const Valeur&, int) const ;

  /**
   * Idem as {\tt Scalar::poisson\_dirichlet}, the boundary condition being on 
   * both the function and its radial derivative. The boundary condition 
   * at infinity is relaxed.
   */

  Scalar poisson_frontiere_double   (const Valeur&, const Valeur&, int) const ;

  /** Solves the scalar Poisson equation with {\tt *this} as a source
   *   (version with parameters to control the resolution).
   *   The source $\sigma$ of the equation $\Delta u = \sigma$ is 
   *   represented by the {\tt Scalar} {\tt *this}. 
   *   The regularized source
   *   $\sigma_{\rm regu} = \sigma - \sigma_{\rm div}$
   *   is constructed and solved.
   *   Note that {\tt dzpuis} must be equal to 2 or 4, i.e. that the
   *   quantity stored in {\tt *this} is in fact $r^2 \sigma$ or
   *   $r^4 \sigma$ in the compactified external domain.
   *   @param k_div [input] regularization degree of the procedure
   *   @param nzet [input] number of domains covering the star
   *   @param unsgam1 [input] parameter $1/(\gamma-1)$ where $\gamma$
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

  /** Checks if a Poisson equation with {\tt *this} as a source
   *  has been correctly solved.
   * 
   *  @param uu [input] Solution {\it u} of the Poisson equation
   *		      $\Delta u = \sigma$,  $\sigma$ being 
   *		      represented by the {\tt Scalar} {\tt *this}.
   * 
   *  @param ostr [input/output] Output stream used for displaying
   *		{\tt err}.
   *
   *  @param detail [input] if {\tt true} displays {\tt err(0, *)}, 
   *		    {\tt err(1, *)} and {\tt err(2, *)} \\
   *		if {\tt false} (default),  displays only 
   *		the relative error {\tt err(0, *)}. 
   *  
   *  @return 2-D {\tt Tbl} {\tt err} decribing the errors in each 
   *	    domain: \\
   *	{\tt err(0, l) : } Relative error in domain no. {\tt l}, 
   *	    defined as the maximum value of 
   *	    $|\Delta u - \sigma|$ in that domain divided by {\it M}, 
   *	    where {\it M} is the maximum value of $|\sigma|$ 
   *	    over all domains if {\tt dzpuis = 0} or $\sigma$ is
   *	    zero in the compactified external domain (CED). If 
   *	    {\tt dzpuis != 0} and $\sigma$ does not vanish in the 
   *	    CED, the value of {\it M} used in the
   *	    non-compactified domains is the maximum value over
   *	    these domains, whereas the value of {\it M} used in the
   *	    compactified external domain is the maximum value
   *	    on that particular domain. \\
   *	{\tt err(1, l) : }  Maximum value of the absolute error
   *			$|\Delta u - \sigma|$ in domain no. {\tt l} \\
   *	{\tt err(2, l) : }  Maximum value of $|\sigma|$ in domain 
   *			    no. {\tt l} 
   */
  Tbl test_poisson(const Scalar& uu, ostream& ostr, 
		   bool detail = false) const ;  

	/** Solves the angular Poisson equation with {\tt *this} as source. 
	 * The angular Poisson equation is $\Delta_{\theta\varphi} u = \sigma$,
	 * where $\Delta_{\theta\varphi} u := \frac{\partial^2 u}
	 *  {\partial \theta^2} + \frac{1}{\tan \theta} \frac{\partial u}
	 *  {\partial \theta} +\frac{1}{\sin^2 \theta}\frac{\partial^2 u}
	 *  {\partial \varphi^2}$.
	 * 
	 *   @return solution {\it u}. 
	 */
	Scalar poisson_angu() const ;


  /** Performs one time-step integration (from $t=J \to J+1$) of the 
   *   scalar d'Alembert equation with {\tt *this} being the value of 
   *   the function {\it f} at time {\it J}.
   *
   *   Works only with an affine mapping (class {\tt Map\_af}) and,
   *   if the last domain is a compactified one, it simply copies
   *   the value of the field in this last domain at the time-step {\it J}
   *   to the last domain of the returned solution.
   *   @param par [input/output] possible parameters to control the
   *   resolution of the d'Alembert equation: \\
   *   {\tt par.get\_double(0)} : [input] the time step {\it dt},\\
   *   {\tt par.get\_int(0)} : [input] the type of boundary conditions
   *   set at the outer boundary (0 : reflexion, 1 : Sommerfeld 
   *   outgoing wave, valid only for {\it l=0} components, 2 : Bayliss 
   *   \& Turkel outgoing wave, valid for {\it l=0, 1, 2} components)\\
   *   {\tt par.get\_int\_mod(0)} : [input/output] set to 0 at first
   *   call, is used as a working flag after (must not be modified after
   *   first call)\\
   *   {\tt par.get\_tensor\_mod(0)} : [input] (optional) if the wave 
   *   equation is on a curved space-time, this is the potential in front
   *   of the Laplace operator. It has to be a Scalar and updated at 
   *    every time-step (for a potential depending on time).\\
   *   Note: there are many other working objects attached to this
   *   {\tt Param}, so one should not modify it.\\
   *   There should be exactly one {\tt Param} for each wave equation to be 
   *   solved. 
   *   @param fJm1 [input] solution $f^{J-1}$ at time {\it J-1}
   *   @param source [input] source $\sigma$ of the d'Alembert equation 
   *	    $\diamond u = \sigma$.
   *   @return solution $f^{J+1}$ at time {\it J+1}
   *   with boundary conditions defined by {\tt par.get\_int(0)}.
   */
  Scalar avance_dalembert(Param& par, const Scalar& fJm1, const Scalar& source) 
    const ;

 /**
   * Resolution of a general elliptic equation, putting zero at infinity.
   * @param params [input] the operators and variables to be used.
   **/
  Scalar sol_elliptic(const Param_elliptic& params) const ;

   /**
   * Resolution of a general elliptic equation, putting zero at the outermost 
   * shell and not solving in the compactified domain.
   * @param params [input] the operators and variables to be used.
   **/
  Scalar sol_elliptic_no_zec(const Param_elliptic& params) const ;

  /**
   * General elliptic solver.
   * The equation is not solved in the compactified domain and the 
   * matching is done with a solution of the type 
   * $\frac{\sin \left( f r + \phi \right)}{r}$ , minimizing the coefficient 
   * with respect to $\phi$.
   * @param params [input] the operators and variables to be used.
   * @param freq [input] : the frequency $f$.
   * @param nbr_phase [input] : number of points for the maximization 
   * over $\phi$.
   * @param error [output] : coefficient of the oscillatory solution 
   * in the external domain.
   * @param phase [output] : phase that minimizes {\tt coef}.
   **/
  Scalar sol_elliptic_sin_zec(const Param_elliptic& params, double freq,
			      int nbr_phase, double& coef, 
			      double& phase) const ;

   
  /**
   * Resolution of a general elliptic equation fixing the dericative at 
   * the origin and relaxing one continuity  condition.
   * 
   * @param val [input] value of the derivative.
   * @param params [input] the operators and variables to be used.
   **/
  Scalar sol_elliptic_fixe_der_zero(double val, 
				    const Param_elliptic& params) const ;
  
	
  // Import from other mapping 
  // -------------------------

  /** Assignment to another {\tt Scalar} defined on a different mapping.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import(const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping.
   *  Case where the {\tt Scalar} is symmetric with respect to the plane y=0.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_symy(const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping.
   *  Case where the {\tt Scalar} is antisymmetric with respect to the 
   *  plane y=0.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_asymy(const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import(int nzet, const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping.
   *  Case where the {\tt Scalar} is symmetric with respect to the plane y=0.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_symy(int nzet, const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping.
   *  Case where the {\tt Scalar} is antisymmetric with respect to the 
   *  plane y=0.
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_asymy(int nzet, const Scalar& ci) ;	 

 protected:
  /** Assignment to another {\tt Scalar} defined on a different mapping,
   *  when the two mappings do not have a particular relative orientation.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_gal(int nzet, const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping,
   *  when the two mappings have aligned Cartesian axis. 
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_align(int nzet, const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping,
   *  when the two mappings have anti-aligned Cartesian axis (i.e.
   *  $x_1 = - x_2$,  $y_1 = - y_2$,  $z_1 = z_2$). 
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_anti(int nzet, const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping,
   *  when the two mappings have aligned Cartesian axis. 
   *  Case where the {\tt Scalar} is symmetric with respect to the plane y=0.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_align_symy(int nzet, const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping,
   *  when the two mappings have anti-aligned Cartesian axis (i.e.
   *  $x_1 = - x_2$,  $y_1 = - y_2$,  $z_1 = z_2$). 
   *  Case where the {\tt Scalar} is symmetric with respect to the plane y=0.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_anti_symy(int nzet, const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping,
   *  when the two mappings have aligned Cartesian axis. 
   *  Case where the {\tt Scalar} is antisymmetric with respect to the 
   *  plane y=0.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_align_asymy(int nzet, const Scalar& ci) ;	 

  /** Assignment to another {\tt Scalar} defined on a different mapping,
   *  when the two mappings have anti-aligned Cartesian axis (i.e.
   *  $x_1 = - x_2$,  $y_1 = - y_2$,  $z_1 = z_2$). 
   *  Case where the {\tt Scalar} is antisymmetric with respect to the 
   *  plane y=0.
   *
   *  This assignment is performed point to point by means of the
   *  spectral expansion of the original {\tt Scalar}. 
   *	@param nzet [input] Number of domains of the destination
   *			    mapping (i.e. {\tt this->mp}) where the 
   *			    importation is performed: the assignment
   *			    is done for the domains whose indices are
   *			    between 0 and {\tt nzet-1}. In the other
   *			    domains, {\tt *this} is set to zero. 
   *	@param ci [input] {\tt Scalar} to be imported.
   */
  void import_anti_asymy(int nzet, const Scalar& ci) ;	 
	

  friend Scalar operator-(const Scalar& ) ;			
  friend Scalar operator+(const Scalar&, const Scalar &) ;	
  friend Scalar operator+(const Scalar&, double ) ;		
  friend Scalar operator-(const Scalar &, const Scalar &) ;
  friend Scalar operator-(const Scalar&, double ) ;		
  friend Scalar operator*(const Scalar &, const Scalar &) ;
  friend Scalar operator%(const Scalar &, const Scalar &) ;
  friend Scalar operator*(double, const Scalar &) ;		
  friend Scalar operator/(const Scalar &, const Scalar &) ;
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
 * @name Scalar mathematics
 */
//@{
Scalar operator+(const Scalar& ) ;			/// + Scalar
Scalar operator-(const Scalar& ) ;			/// - Scalar
Scalar operator+(const Scalar&, const Scalar &) ;	/// Scalar + Scalar
Scalar operator+(const Scalar&, double ) ;		/// Scalar + double
Scalar operator+(double, const Scalar& ) ;		/// double + Scalar 
Scalar operator+(const Scalar&, int ) ;		/// Scalar + int
Scalar operator+(int, const Scalar& ) ;		/// int + Scalar 
Scalar operator-(const Scalar &, const Scalar &) ;	/// Scalar - Scalar
Scalar operator-(const Scalar&, double ) ;		/// Scalar - double
Scalar operator-(double, const Scalar& ) ;		/// double - Scalar 
Scalar operator-(const Scalar&, int ) ;		/// Scalar - int
Scalar operator-(int, const Scalar& ) ;		/// int - Scalar 
Scalar operator*(const Scalar &, const Scalar &) ;	/// Scalar * Scalar
Scalar operator%(const Scalar &, const Scalar &) ;	/// Scalar * Scalar with desaliasing
Scalar operator*(const Scalar&, double ) ;		/// Scalar * double
Scalar operator*(double, const Scalar &) ;		/// double * Scalar
Scalar operator*(const Scalar&, int ) ;		/// Scalar * int
Scalar operator*(int, const Scalar& ) ;		/// int * Scalar 
Scalar operator/(const Scalar &, const Scalar &) ;	/// Scalar / Scalar
Scalar operator/(const Scalar&, double ) ;		/// Scalar / double
Scalar operator/(double, const Scalar &) ;		/// double / Scalar
Scalar operator/(const Scalar&, int ) ;		/// Scalar / int
Scalar operator/(int, const Scalar &) ;		/// int / Scalar

Scalar sin(const Scalar& ) ;		/// Sine
Scalar cos(const Scalar& ) ;		/// Cosine
Scalar tan(const Scalar& ) ;		/// Tangent
Scalar asin(const Scalar& ) ;		/// Arcsine
Scalar acos(const Scalar& ) ;		/// Arccosine
Scalar atan(const Scalar& ) ;		/// Arctangent
Scalar exp(const Scalar& ) ;		/// Exponential
Scalar log(const Scalar& ) ;		/// Neperian logarithm
Scalar log10(const Scalar& ) ;	/// Basis 10 logarithm
Scalar sqrt(const Scalar& ) ;		/// Square root
Scalar racine_cubique (const Scalar& ) ;		/// Cube root
Scalar pow(const Scalar& , int ) ;	/// Power ${\tt Scalar}^{\tt int}$
Scalar pow(const Scalar& , double ) ; /// Power ${\tt Scalar}^{\tt double}$
Scalar abs(const Scalar& ) ;		/// Absolute value

/**
 * Maximum values of a {\tt Scalar} in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Scalar& ) ;   

/**
 * Minimum values of a {\tt Scalar} in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Scalar& ) ;   

/**
 * Sums of the absolute values of all the values of the {\tt Scalar} 
 * in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Scalar& ) ;   

/**
 * Relative difference between two {\tt Scalar} (norme version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt norme[a(l)-b(l)]/norme[b(l)]} if {\tt b(l)!=0} and
 *	   {\tt norme[a(l)-b(l)]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrel(const Scalar& a, const Scalar& b) ; 

/**
 * Relative difference between two {\tt Scalar} (max version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt max[abs(a(l)-b(l))]/max[abs(b(l))]} if {\tt b(l)!=0} and
 *	   {\tt max[abs(a(l)-b(l))]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrelmax(const Scalar& a, const Scalar& b) ; 

//@}
#endif
