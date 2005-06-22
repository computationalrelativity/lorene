/*
 *  Definition of Lorene class Tbl_val
 *
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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


#ifndef	__TBL_VAL_H_
#define	__TBL_VAL_H_

/*
 * $Id$
 * $Log$
 * Revision 1.7  2005/06/22 09:09:38  lm_lin
 *
 * Grid wedding: convert from the old C++ object "Cmp" to "Scalar".
 *
 * Revision 1.6  2004/11/26 17:02:18  j_novak
 * Added a function giving a smooth transition to the atmosphere.
 *
 * Revision 1.5  2004/03/22 13:12:43  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.4  2002/11/13 11:22:57  j_novak
 * Version "provisoire" de l'interpolation (sommation depuis la grille
 * spectrale) aux interfaces de la grille de Valence.
 *
 * Revision 1.3  2002/11/12 10:03:53  j_novak
 * The method "Tbl_val::get_gval" has been changed to "get_grid".
 *
 * Revision 1.2  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1  2001/11/22 13:38:09  j_novak
 * added Include files for Valencia objects: tbl_val.h and grille_val.h
 *
 *
 * $Header$
 *
 */

// Fichiers includes
#include <assert.h>
#include <stdlib.h>

#include "grille_val.h"
#include "tensor.h"

class Grille_val ; 

/**
 * Finite-difference array intended to store field values.\ingroup (mdm)
 *
 * Class defined on a cartesian (\c Gval_cart ) or spherical 
 * (\c Gval_spher ) grid, in order to represent 
 * Godunov-type arrays in 1,2 or 3D. 
 *
 */
class Tbl_val {

  // Data : 
  // -----
 private:
  /// logical state (\c ETATNONDEF , \c ETATQCQ  or \c ETATZERO ).
  int etat ;
  /**
   * The \c Dim_tbl  giving the dimensions and number of points (without 
   * the hidden cells).
   */
  const Dim_tbl* dim ;
/// The \c Grille_val (cartesian or spherical) on which the array is defined
  const Grille_val* gval ;  
  
 public:
  /// The array of \c double at the nodes
  double* t ;	    
  /// The array at z (or r) interfaces
  double* tzri ;    
  /// The array at x (or \f$\theta\f$) interfaces
  double* txti ;    
  /// The array at y (or \f$\phi\f$) interfaces
  double* typi ;    
  
  // Constructors - Destructor
  // -------------------------
  
 public:
  /// Constructor from a 3D grid
  explicit Tbl_val(const Grille_val* ) ; 
  /// Constructor from a file (see \c sauve(FILE*) )
  explicit Tbl_val(const Grille_val*, FILE* ) ;	
  /// Copy constructor
  Tbl_val(const Tbl_val& ) ;		
  
  /// Destructor
  ~Tbl_val() ;			
  
  // Assignement
  // -----------
  /// Assignment to another \c Tbl_val
  void operator=(const Tbl_val& ) ;	
  /// Assignment to a \c double
  void operator=(double ) ; 
  /// Assignment to a \c int
  void operator=(int ) ;	 

  // Memory management
  // -----------------
 private:
  /** Logical destructor: dellocates the memory occupied by the array
   *  \c t  and sets the logical state to ETATNONDEF. 
   */
  void del_t() ;		
  
 public:
  
  /**
   * Sets the logical state to \c ETATNONDEF  (undefined). 
   * Deallocates the memory occupied by the \c double array \c t .
   */
  void set_etat_nondef() ;	
  
  /**
   * Sets the logical state to \c ETATZERO  (zero). 
   * Deallocates the memory occupied by the \c double array \c t .
   */
  void set_etat_zero() ;	    	
  
  /**
   * Sets the logical state to \c ETATQCQ  (ordinary state).
   * If the state (member \c etat ) is already \c ETATQCQ , this 
   * function does nothing. Otherwise, it performs the memory allocation
   * for the \c double array \c t .  
   */
  void set_etat_qcq() ;	    	
    
  /**
   * Sets the \c Tbl_val to zero in a hard way. 
   * 1/ Sets the logical state to \c ETATQCQ , i.e. to an ordinary state.
   * 2/ Allocates the memory of the \c double array \c t , and fills it
   * with zeros. NB: this function must be used for debugging purposes only.
   * For other operations, the function \c set_etat_zero()  must
   * be perferred. 
   */
  void annule_hard() ;			
  
  // Access to individual elements
  // -----------------------------
 public:
  /// Read/write of a particular element (index \c i )  (1D case)
  double& set(int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 1 ) ;
    int fant = gval->get_fantome() ; 
    assert( i >= - fant) ; 
    assert( i < dim->dim[0] + fant) ;
    return t[i + fant] ;
  } ;
  
  /// Read/write of a particular element on the interface (index \c i )  (1D case)
  double& set_zri(int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 1 ) ; 
    int fant = gval->get_fantome() ; 
    assert( i >= -fant ) ; 
    assert( i < dim->dim[0] + fant + 1) ;
    return tzri[i+fant] ;
  } ;
  
  /// Read-only of a particular element (index \c i ) (1D case)
  double operator()(int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 1 ) ; 
    int fant = gval->get_fantome() ; 
    assert( i >= -fant ) ; 
    assert( i < dim->dim[0] + fant ) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return t[i+fant] ;
  };
  
  /// Read-only of a particular element on the interface (index \c i ) (1D case)
  double get_zri(int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 1 ) ; 
    int fant = gval->get_fantome() ; 
    assert( i >= -fant ) ; 
    assert( i < dim->dim[0] + fant +1) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return tzri[i+fant] ;
  };
  
  /// Read/write of a particular element (index \c (j,i) ) (2D case)
  double& set(int j, int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 2 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0]+fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1]+fant) ) ;
    return t[(dim->dim[0] +2*fant)* (j+fant) + i + fant] ;
  };

  /**
   * Read/write of a particular element on the x (or \f$\theta\f$) 
   * interface (index \c (j,i) ) (2D case)
   */
  double& set_xti(int j, int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 2 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0]+fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1]+fant+1) ) ;
    return txti[(dim->dim[0] +2*fant)*(j+fant) + i + fant] ;
  };
  
  /**
   * Read/write of a particular element on the z (or r) 
   * interface (index \c (j,i) ) (2D case)
   */
  double& set_zri(int j, int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 2 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant+1) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    return tzri[(dim->dim[0] +2*fant+1)*(j+fant) + i + fant] ;
  };
  
  /// Read-only of a particular element (index \c (j,i) ) (2D case)
  double operator()(int j, int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 2 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return t[(dim->dim[0] + 2*fant) *(j+fant) + i + fant] ;
  };
  
  /**
   * Read-only of a particular element on the x (or \f$\theta\f$) interface 
   * (index \c (j,i) ) (2D case)
   */
  double get_xti(int j, int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 2 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant + 1) ) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return txti[(dim->dim[0] + 2*fant) *(j+fant) + i + fant] ;
  };
  
  /**
   * Read-only of a particular element on the z (or r) interface 
   * (index \c (j,i) ) (2D case)
   */
  double get_zri(int j, int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 2 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant + 1) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return tzri[(dim->dim[0] + 2*fant + 1) *(j+fant) + i + fant] ;
  };
  
  /// Read/write of a particular element (index \c (k,j,i) ) (3D case)
  double& set(int k, int j, int i) {	
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 3 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    assert( (k>=-fant) && (k<dim->dim[2] + fant) ) ;
    return t[(dim->dim[1]+2*fant)*(dim->dim[0]+2*fant)*(k+fant) + 
	    (dim->dim[0]+2*fant)*(j+fant) + i +fant] ;
  };
  
  /**
   * Read/write of a particular element on the y (or \f$\phi\f$) 
   * interface (index \c (k,j,i) ) (3D case)
   */
  double& set_ypi(int k, int j, int i) {	
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 3 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    assert( (k>=-fant) && (k<dim->dim[2] + fant + 1) ) ;
    return typi[(dim->dim[1]+2*fant)*(dim->dim[0]+2*fant)*(k+fant) + 
	      (dim->dim[0]+2*fant)*(j+fant) + i +fant] ;
  };

  /**
   * Read/write of a particular element on the x (or \f$\theta\f$) 
   * interface (index \c (k,j,i) ) (3D case)
   */
  double& set_xti(int k, int j, int i) {	
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 3 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant + 1) ) ;
    assert( (k>=-fant) && (k<dim->dim[2] + fant) ) ;
    return txti[(dim->dim[1]+2*fant+1)*(dim->dim[0]+2*fant)*(k+fant) + 
	      (dim->dim[0]+2*fant)*(j+fant) + i +fant] ;
  };

  /**
   * Read/write of a particular element on the z (or r) 
   * interface (index \c (k,j,i) ) (3D case)
   */
  double& set_zri(int k, int j, int i) {	
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 3 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant + 1) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    assert( (k>=-fant) && (k<dim->dim[2] + fant) ) ;
    return tzri[(dim->dim[1]+2*fant)*(dim->dim[0]+2*fant+1)*(k+fant) + 
	      (dim->dim[0]+2*fant+1)*(j+fant) + i +fant] ;
  };
  
  /// Read-only of a particular element (index \c (k,j,i) ) (3D case)
  double operator()(int k, int j, int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 3 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    assert( (k>=-fant) && (k<dim->dim[2] + fant) ) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return t[(dim->dim[1]+2*fant)*(dim->dim[0]+2*fant)*(k+fant) 
		 + (dim->dim[0]+2*fant)*(j+fant) + i +fant] ;
  };

  /**
   * Read-only of a particular element on the y (or \f$\phi\f$) interface 
   * (index \c (k,j,i) ) (3D case)
   */
  double get_ypi(int k, int j, int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 3 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    assert( (k>=-fant) && (k<dim->dim[2] + fant + 1) ) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return typi[(dim->dim[1]+2*fant)*(dim->dim[0]+2*fant)*(k+fant) 
		   + (dim->dim[0]+2*fant)*(j+fant) + i +fant] ;
  };
  
  /**
   * Read-only of a particular element on the x (or \f$\theta\f$) interface 
   * (index \c (k,j,i) ) (3D case)
   */
  double get_xti(int k, int j, int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 3 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant + 1) ) ;
    assert( (k>=-fant) && (k<dim->dim[2] + fant) ) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return txti[(dim->dim[1]+2*fant+1)*(dim->dim[0]+2*fant)*(k+fant) 
		   + (dim->dim[0]+2*fant)*(j+fant) + i +fant] ;
  };
  
  /**
   * Read-only of a particular element on the z (or r) interface 
   * (index \c (k,j,i) ) (3D case)
   */
  double get_zri(int k, int j, int i) const {
    assert(etat != ETATNONDEF) ;
    assert( dim->ndim == 3 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant + 1) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    assert( (k>=-fant) && (k<dim->dim[2] + fant) ) ;
    if (etat == ETATZERO) {
      double zero = 0. ;
      return zero ; 
    }
    else return tzri[(dim->dim[1]+2*fant)*(dim->dim[0]+2*fant+1)*(k+fant) 
		   + (dim->dim[0]+2*fant+1)*(j+fant) + i +fant] ;
  };

  // Extraction of information
  // -------------------------
/// Gives the logical state
  int get_etat() const { return etat ; };	    
  
  /// Gives the size of the node array (including the hidden cells)
  int get_taille() const { 
    int resu = 1 ;
    for (int i=0; i<dim->ndim; i++) 
      resu *= dim->dim[i] + 2*(gval->get_fantome()) ;
    return resu ; }; 
  
  /// Gives the size of the interface arrays (including the hidden cells)
  int get_taille_i(int i) const {
    assert (i<dim->ndim) ; 
    int resu = 1 ;
    for (int j=0; j<dim->ndim; j++) 
      if (j!=i) {
	resu *= dim->dim[j] + 2*gval->get_fantome() ;
      }
      else {
	resu *= dim->dim[j] + 2*gval->get_fantome() + 1 ;
      }
    return resu ; }; 
  
  /// Gives the number of dimensions (ie \c dim->ndim )
  int get_ndim() const { return dim->ndim ; };	
  
  /// Gives the \c i th dimension (ie \c dim->dim[i] , without hidden cells)
  int get_dim(int i) const {	
    assert( (i>=0) && (i<dim->ndim) ) ;
    return dim->dim[i] ;
  };

  /// Returns a pointer on the grid on which the \c Tbl_val is defined
  const Grille_val* get_grille() const { return gval ; } ;
  
  // Outputs
  // -------
 public:
/// Save in a file
  void sauve(FILE* ) const ;	
  
  /** Prints only the values greater than a given threshold.
   *   @param ostr [input] Output stream used for the printing 
   *   @param precision [input] Number of printed digits (default: 4)
   *   @param threshold [input] Value above which an array element is printed
   *    (default: 1.e-7)
   */
  void affiche_seuil(ostream& ostr, int precision = 4, 
		     double threshold = 1.e-7) const ;
  /// Display   
  friend ostream& operator<<(ostream& , const Tbl_val& ) ;	
  
  // Member arithmetics
  // ------------------
 public:
  
/// Addition of a \c Tbl_val to \c this
  void operator+=(const Tbl_val &) ;	
/// Addition of a \c double to \c this
  void operator+=(double) ;	
/// Subtraction of a \c Tbl_val to \c this
  void operator-=(const Tbl_val &) ;	
/// Subtraction of a \c double to \c this
  void operator-=(double) ;	
/// Multiplication of \c this by a \c Tbl_val
  void operator*=(const Tbl_val &) ;	
/// Multiplication of \c this by a \c double
  void operator*=(double) ;	
/// Division of \c this by a \c Tbl_val
  void operator/=(const Tbl_val &) ;	
/// Division of \c this by a \c double
  void operator/=(double) ;	

  /**
   * Interpolation from a \c Tbl_val to a \c Cmp . The 
   * \c Cmp  is evaluated only in zones [lmin, lmax[. 
   * @param map [input] The \c Mapping  to which the \c Tbl_val is 
   *                    interpolated. The symetries of both grids must be
   *                    the same (see \c Mg3d  and \c Grille_val 
   *                    documentation), and the spectral grid (between lmin
   *                    and lmax-1) must be included in the Godunov one.
   *                    The number of points in \f$\theta\f$ and \f$\phi\f$ of the
   *                    spectral grid may be different/domain. Still, the
   *                    domain with the highest number of points in \f$\theta\f$
   *                    (resp.\f$\phi\f$) must contain the collocation points
   * @param lmax [input] index of the outer zone \b +1
   * @param lmin [input] index of the inner zone 
   * @param type_inter [input] type of interpolation: \\
   *    0 -> uses the \c INSMTS  routine of second derivative minimization\\
   *    1 -> linear interpolation\\
   *    2 -> parabolic interpolation\\
   *    3 -> spline interpolation (not implemented yet)\\
   * @return Cmp containing the value of the field at spectral collocation
   * points.
   */
  Scalar to_spectral(const Map& map, const int lmax, const int lmin=0, 
		      int type_inter = 2) const ;
  
  /**
   * Interpolation from a \c Cmp  to a \c Tbl_val (spectral
   * summation). The \c Cmp  is considered only in zones [lmin,lmax[.
   * @param meudon [input] The \c Cmp  from which the interpolation is done
   * @param lmax [input] index of the outer zone \b +1
   * @param lmin [input] index of the inner zone 
   */
  void from_spectral(const Scalar& meudon, int lmax, int lmin=0,
		     bool interfr = false, bool interft = false) ;


  void smooth_atmosphere(double atmosphere_thr) ;
} ;


/**
 * \defgroup tbl_val_m Tbl_val Mathematics
 * \ingroup (mdm)
 * @{
 */

/// + Tbl_val
Tbl_val operator+(const Tbl_val&) ;			
/// \c - Tbl_val
Tbl_val operator-(const Tbl_val&) ;			
/// Tbl_val + Tbl_val
Tbl_val operator+(const Tbl_val&, const Tbl_val&) ;	
/// Tbl_val + double
Tbl_val operator+(const Tbl_val&, double) ;		
/// double + Tbl_val
Tbl_val operator+(double, const Tbl_val&) ;		
/// Tbl_val + int
Tbl_val operator+(const Tbl_val&, int) ;		
/// int + Tbl_val
Tbl_val operator+(int, const Tbl_val&) ;		
/// Tbl_val - Tbl_val
Tbl_val operator-(const Tbl_val&, const Tbl_val&) ;	
/// Tbl_val - double
Tbl_val operator-(const Tbl_val&, double) ;		
/// double - Tbl_val
Tbl_val operator-(double, const Tbl_val&) ;		
/// Tbl_val - int
Tbl_val operator-(const Tbl_val&, int) ;		
/// int - Tbl_val
Tbl_val operator-(int, const Tbl_val&) ;		
/// Tbl_val * Tbl_val
Tbl_val operator*(const Tbl_val&, const Tbl_val&) ;	
/// Tbl_val * double
Tbl_val operator*(const Tbl_val&, double) ;		
/// double * Tbl_val
Tbl_val operator*(double, const Tbl_val&) ;		
/// Tbl_val * int
Tbl_val operator*(const Tbl_val&, int) ;		
/// int * Tbl_val
Tbl_val operator*(int, const Tbl_val&) ;		
/// Tbl_val / Tbl_val
Tbl_val operator/(const Tbl_val&, const Tbl_val&) ;	
/// Tbl_val / double
Tbl_val operator/(const Tbl_val&, double) ;		
/// double / Tbl_val
Tbl_val operator/(double, const Tbl_val&) ;		
/// Tbl_val / int
Tbl_val operator/(const Tbl_val&, int) ;		
/// int / Tbl_val
Tbl_val operator/(int, const Tbl_val&) ;		

/// Sine
Tbl_val sin(const Tbl_val& ) ;	    
/// Cosine
Tbl_val cos(const Tbl_val& ) ;	    
/// Tangent
Tbl_val tan(const Tbl_val& ) ;	    
/// Arcsine
Tbl_val asin(const Tbl_val& ) ;	    
/// Arccosine
Tbl_val acos(const Tbl_val& ) ;	    
/// Arctangent
Tbl_val atan(const Tbl_val& ) ;	    
/// Exponential
Tbl_val exp(const Tbl_val& ) ;	    
/// Neperian logarithm
Tbl_val log(const Tbl_val& ) ;	    
/// Basis 10 logarithm
Tbl_val log10(const Tbl_val& ) ;    
/// Square root
Tbl_val sqrt(const Tbl_val& ) ;	    
/// cube root
Tbl_val racine_cubique (const Tbl_val&) ; 
/// Power \f${\tt Tbl_val}^{\tt int}\f$
Tbl_val pow(const Tbl_val& , int ) ;  
/// Power \f${\tt Tbl_val}^{\tt double}\f$
Tbl_val pow(const Tbl_val& , double ) ; 
/// Absolute value
Tbl_val abs(const Tbl_val& ) ;	    
/// Maximum value of the \c Tbl_val elements
double max(const Tbl_val& ) ;   
/// Minimum value of the \c Tbl_val elements
double min(const Tbl_val& ) ;   

/// Sum of the absolute values of all the \c Tbl_val elements
double norme(const Tbl_val& ) ;   

/**
 * Relative difference between two \c Tbl_val (norme version).
 * Returns \c norme(a-b)/norme(b)  unless \c b=0 , in which
 * case it returns \c norme(a-b) .
 */
double diffrel(const Tbl_val& a, const Tbl_val& b) ; 

/**
 * Relative difference between two \c Tbl_val (max version).
 * Returns \c max(abs(a-b))/max(abs(b))  unless \c b=0 , in which
 * case it returns \c max(abs(a-b)) .
 */
double diffrelmax(const Tbl_val& a, const Tbl_val& b) ; 

/** @} */

#endif

