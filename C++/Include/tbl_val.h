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
#include "tenseur.h"

class Grille_val ; 

/**
 * Tbl_val.
 *
 * Class defined on a cartesian ({\tt Gval_cart}) or spherical 
 * ({\tt Gval_spher}) grid, in order to represent 
 * Godunov-type arrays in 1,2 or 3D. 
 *
 * @version #$Id$#
 */
class Tbl_val {

  // Data : 
  // -----
 private:
  /// logical state ({\tt ETATNONDEF}, {\tt ETATQCQ} or {\tt ETATZERO}).
  int etat ;
  /**
   * The {\tt Dim_tbl} giving the dimensions and number of points (without 
   * the hidden cells).
   */
  const Dim_tbl* dim ;
/// The {\tt Grille_val} (cartesian or spherical) on which the array is defined
  const Grille_val* gval ;  
  
 public:
  /// The array of {\tt double} at the nodes
  double* t ;	    
  /// The array at z (or r) interfaces
  double* tzri ;    
  /// The array at x (or $\theta$) interfaces
  double* txti ;    
  /// The array at y (or $\phi$) interfaces
  double* typi ;    
  
  // Constructors - Destructor
  // -------------------------
  
 public:
  /// Constructor from a 3D grid
  explicit Tbl_val(const Grille_val* ) ; 
  /// Constructor from a file (see {\tt sauve(FILE* )})
  explicit Tbl_val(const Grille_val*, FILE* ) ;	
  /// Copy constructor
  Tbl_val(const Tbl_val& ) ;		
  
  /// Destructor
  ~Tbl_val() ;			
  
  // Assignement
  // -----------
  /// Assignment to another {\tt Tbl_val}
  void operator=(const Tbl_val& ) ;	
  /// Assignment to a {\tt double}
  void operator=(double ) ; 
  /// Assignment to a {\tt int}
  void operator=(int ) ;	 

  // Memory management
  // -----------------
 private:
  /** Logical destructor: dellocates the memory occupied by the array
   *  {\tt t} and sets the logical state to ETATNONDEF. 
   */
  void del_t() ;		
  
 public:
  
  /**
   * Sets the logical state to {\tt ETATNONDEF} (undefined). 
   * Deallocates the memory occupied by the {\tt double} array {\tt t}.
   */
  void set_etat_nondef() ;	
  
  /**
   * Sets the logical state to {\tt ETATZERO} (zero). 
   * Deallocates the memory occupied by the {\tt double} array {\tt t}.
   */
  void set_etat_zero() ;	    	
  
  /**
   * Sets the logical state to {\tt ETATQCQ} (ordinary state).
   * If the state (member {\tt etat}) is already {\tt ETATQCQ}, this 
   * function does nothing. Otherwise, it performs the memory allocation
   * for the {\tt double} array {\tt t}.  
   */
  void set_etat_qcq() ;	    	
    
  /**
   * Sets the {\tt Tbl_val} to zero in a hard way. 
   * 1/ Sets the logical state to {\tt ETATQCQ}, i.e. to an ordinary state.
   * 2/ Allocates the memory of the {\tt double} array {\tt t}, and fills it
   * with zeros. NB: this function must be used for debugging purposes only.
   * For other operations, the function {\tt set\_etat\_zero()} must
   * be perferred. 
   */
  void annule_hard() ;			
  
  // Access to individual elements
  // -----------------------------
 public:
  /// Read/write of a particular element (index {\tt i})  (1D case)
  double& set(int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 1 ) ;
    int fant = gval->get_fantome() ; 
    assert( i >= - fant) ; 
    assert( i < dim->dim[0] + fant) ;
    return t[i + fant] ;
  } ;
  
  /// Read/write of a particular element on the interface (index {\tt i})  (1D case)
  double& set_zri(int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 1 ) ; 
    int fant = gval->get_fantome() ; 
    assert( i >= -fant ) ; 
    assert( i < dim->dim[0] + fant + 1) ;
    return tzri[i+fant] ;
  } ;
  
  /// Read-only of a particular element (index {\tt i}) (1D case)
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
  
  /// Read-only of a particular element on the interface (index {\tt i}) (1D case)
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
  
  /// Read/write of a particular element (index {\tt (j,i)}) (2D case)
  double& set(int j, int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 2 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0]+fant) ) ;
    assert( (j>=-fant) && (j<dim->dim[1]+fant) ) ;
    return t[(dim->dim[0] +2*fant)* (j+fant) + i + fant] ;
  };

  /**
   * Read/write of a particular element on the x (or $\theta$) 
   * interface (index {\tt (j,i)}) (2D case)
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
   * interface (index {\tt (j,i)}) (2D case)
   */
  double& set_zri(int j, int i) {
    assert (etat == ETATQCQ) ;
    assert( dim->ndim == 2 ) ;
    int fant = gval->get_fantome() ; 
    assert( (i>=-fant) && (i<dim->dim[0] + fant+1) ) ;
    assert( (j>=-fant) && (j<dim->dim[1] + fant) ) ;
    return tzri[(dim->dim[0] +2*fant+1)*(j+fant) + i + fant] ;
  };
  
  /// Read-only of a particular element (index {\tt (j,i)}) (2D case)
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
   * Read-only of a particular element on the x (or $\theta$) interface 
   * (index {\tt (j,i)}) (2D case)
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
   * (index {\tt (j,i)}) (2D case)
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
  
  /// Read/write of a particular element (index {\tt (k,j,i)}) (3D case)
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
   * Read/write of a particular element on the y (or $\phi$) 
   * interface (index {\tt (k,j,i)}) (3D case)
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
   * Read/write of a particular element on the x (or $\theta$) 
   * interface (index {\tt (k,j,i)}) (3D case)
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
   * interface (index {\tt (k,j,i)}) (3D case)
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
  
  /// Read-only of a particular element (index {\tt (k,j,i)}) (3D case)
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
   * Read-only of a particular element on the y (or $\phi$) interface 
   * (index {\tt (k,j,i)}) (3D case)
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
   * Read-only of a particular element on the x (or $\theta$) interface 
   * (index {\tt (k,j,i)}) (3D case)
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
   * (index {\tt (k,j,i)}) (3D case)
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
  
  /// Gives the number of dimensions (ie {\tt dim->ndim})
  int get_ndim() const { return dim->ndim ; };	
  
  /// Gives the {\tt i}th dimension (ie {tt dim->dim[i]}, without hidden cells)
  int get_dim(int i) const {	
    assert( (i>=0) && (i<dim->ndim) ) ;
    return dim->dim[i] ;
  };

  /// Returns a pointer on the grid on which the {\tt Tbl\_val} is defined
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
  
/// Addition of a {\tt Tbl_val} to {\tt this}
  void operator+=(const Tbl_val &) ;	
/// Addition of a {\tt double} to {\tt this}
  void operator+=(double) ;	
/// Subtraction of a {\tt Tbl_val} to {\tt this}
  void operator-=(const Tbl_val &) ;	
/// Subtraction of a {\tt double} to {\tt this}
  void operator-=(double) ;	
/// Multiplication of {\tt this} by a {\tt Tbl_val}
  void operator*=(const Tbl_val &) ;	
/// Multiplication of {\tt this} by a {\tt double}
  void operator*=(double) ;	
/// Division of {\tt this} by a {\tt Tbl_val}
  void operator/=(const Tbl_val &) ;	
/// Division of {\tt this} by a {\tt double}
  void operator/=(double) ;	

  /**
   * Interpolation from a {\tt Tbl\_val} to a {\tt Cmp}. The 
   * {\tt Cmp} is evaluated only in zones [lmin, lmax[. 
   * @param map [input] The {\tt Mapping} to which the {\tt Tbl_val} is 
   *                    interpolated. The symetries of both grids must be
   *                    the same (see {\tt Mg3d} and {\tt Grille_val} 
   *                    documentation), and the spectral grid (between lmin
   *                    and lmax-1) must be included in the Godunov one.
   *                    The number of points in $\theta$ and $\phi$ of the
   *                    spectral grid may be different/domain. Still, the
   *                    domain with the highest number of points in $\theta$
   *                    (resp.$\phi$) must contain the collocation points
   * @param lmax [input] index of the outer zone {\bf +1}
   * @param lmin [input] index of the inner zone 
   * @param type_inter [input] type of interpolation: \\
   *    0 -> uses the {\tt INSMTS} routine of second derivative minimization\\
   *    1 -> linear interpolation\\
   *    2 -> parabolic interpolation\\
   *    3 -> spline interpolation (not implemented yet)\\
   * @return Cmp containing the value of the field at spectral collocation
   * points.
   */
  Cmp to_spectral(const Map& map, const int lmax, const int lmin=0, 
		      int type_inter = 2) const ;
  
  /**
   * Interpolation from a {\tt Cmp} to a {\tt Tbl\_val} (spectral
   * summation). The {\tt Cmp} is considered only in zones [lmin,lmax[.
   * @param meudon [input] The {\tt Cmp} from which the interpolation is done
   * @param lmax [input] index of the outer zone {\bf +1}
   * @param lmin [input] index of the inner zone 
   */
  void from_spectral(const Cmp& meudon, int lmax, int lmin=0,
		     bool interfr = false, bool interft = false) ;
} ;


/**
 * @name Tbl_val Mathematics
 */
//@{
/// + Tbl_val
Tbl_val operator+(const Tbl_val&) ;			
/// - Tbl_val
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
/// Power ${\tt Tbl_val}^{\tt int}$
Tbl_val pow(const Tbl_val& , int ) ;  
/// Power ${\tt Tbl_val}^{\tt double}$
Tbl_val pow(const Tbl_val& , double ) ; 
/// Absolute value
Tbl_val abs(const Tbl_val& ) ;	    
/// Maximum value of the {\tt Tbl_val} elements
double max(const Tbl_val& ) ;   
/// Minimum value of the {\tt Tbl_val} elements
double min(const Tbl_val& ) ;   

/// Sum of the absolute values of all the {\tt Tbl_val} elements
double norme(const Tbl_val& ) ;   

/**
 * Relative difference between two {\tt Tbl_val} (norme version).
 * Returns {\tt norme(a-b)/norme(b)} unless {\tt b=0}, in which
 * case it returns {\tt norme(a-b)}.
 */
double diffrel(const Tbl_val& a, const Tbl_val& b) ; 

/**
 * Relative difference between two {\tt Tbl_val} (max version).
 * Returns {\tt max(abs(a-b))/max(abs(b))} unless {\tt b=0}, in which
 * case it returns {\tt max(abs(a-b))}.
 */
double diffrelmax(const Tbl_val& a, const Tbl_val& b) ; 

//@}

#endif

