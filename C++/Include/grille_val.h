/*
 *  Definition of Lorene class Grille_val, Gval_cart and Gval_spher
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


#ifndef __GRILLE_VAL_H_
#define __GRILLE_VAL_H_

/*
 * $Id$
 * $Log$
 * Revision 1.11  2014/10/13 08:52:35  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/10/06 15:09:39  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2008/02/18 13:53:37  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.8  2007/11/02 16:49:12  j_novak
 * Suppression of intermediate array for spectral summation.
 *
 * Revision 1.7  2005/06/22 09:09:38  lm_lin
 *
 * Grid wedding: convert from the old C++ object "Cmp" to "Scalar".
 *
 * Revision 1.6  2004/05/07 12:32:12  j_novak
 * New summation from spectral to FD grid. Much faster!
 *
 * Revision 1.5  2004/03/22 13:12:41  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.4  2002/11/13 11:22:57  j_novak
 * Version "provisoire" de l'interpolation (sommation depuis la grille
 * spectrale) aux interfaces de la grille de Valence.
 *
 * Revision 1.3  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.1  2001/11/22 13:38:09  j_novak
 * added Include files for Valencia objects: tbl_val.h and grille_val.h
 *
 *
 * $Header$
 *
 */


#include <cassert>
#include <cmath>
#include "tensor.h"

namespace Lorene {
/**
 * Base class for Godunov-type grids. \ingroup (mdm)
 * 
 * Takes into account hidden cells and interfaces between cells, 
 * in 1,2 or 3D. 
 * The derived classes stand for spherical and cartesian grids.
 * this base class settles quantities for 1D grid, which is minimal
 * and independent from the geometry (cartesian or spherical).
 * Therefore, the coordinates \e z and \e r are identified.
 *
 * This is an abstract class (cannot be instanciated).
 *
 **/
class Grille_val {
  /// Arrays defined on Godunov-type grids.
  friend class Tbl_val ;
  
  // Data : 
  // -----
 protected:
  /// The dimensions of the grid.
  Dim_tbl dim ;
  /// The number of hidden cells (same on each side)
  int nfantome ;
  /**
   * Type of symmetry in \f$\theta\f$:\li \c SYM -> \f$\theta \to -\theta\f$
   *                               \li \c NONSYM no symmetry
   */                  
  int type_t ;
  /**
   * Type of symmetry in \f$\phi\f$: \li \c SYM -> \f$(x,y) \to (-x,-y)\f$
   *                             \li \c NONSYM no symmetry
   */                  
  int type_p ;

  /// Lower boundary for \e z (or \e r ) direction  
  double *zrmin;

  /// Higher boundary for \e z (or \e r ) direction  
  double *zrmax ;
  
 public:
  /// Arrays containing the values of coordinate \e z (or \e r) on the nodes  
  Tbl *zr;
  /// Arrays containing the values of coordinate \e z (or \e r) on the interfaces
  Tbl *zri ;

  // Constructors - Destructor
  // -------------------------
 protected: 
  /// Auxilliary function used to allocate memory and construct 1D grid
  Tbl* fait_grille1D(const double rmin, const double rmax, const int n) ;

  /// Standard 1D constructor (the size is to be given without hidden cells)
  Grille_val(const double, const double, const int n1,
	     const int fantome = 2) ; 
  
  /// Standard 2D constructor (the sizes are to be given without hidden cells)
  Grille_val(const double, const double, const int n2, const int n1, 
	     const int itype_t, const int fantome = 2) ; 

  /// Standard 3D constructor (the sizes are to be given without hidden cells)
  Grille_val(const double, const double, const int n3, const int n2, const
	     int n1, const int itype_t, const int itype_p, const int
	     fantome = 2); 
  

  /// Copy constructor
  Grille_val(const Grille_val& ) ;		
  
  /// Constructor from a file (see \c sauve(FILE*) )
  Grille_val(FILE* ) ;    		
  
  /// Destructor
  virtual ~Grille_val() ;			
  
  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another Grille_val
  void operator=(const Grille_val&) ;	
  
  // Accessors
  // ---------
 public:
  /// Returns the number of hidden cells
  int get_fantome() const {
    return nfantome ;
  } ;
  
  /// Returns the type of symmetry in \f$\theta\f$
  int get_type_t() const {
    return type_t ;
  } ;
  
  /// Returns the type of symmetry in \f$\phi\f$
  int get_type_p() const {
    return type_p ;
  } ;
  
  /// Returns the number of dimensions
  int get_ndim() const {
    return dim.ndim ;
  } ;
  
  /// Returns the size (without hidden cells)
  int get_dim(const int i) const {
    assert ( (i>=0) && (i<dim.ndim) ) ;
    return dim.dim[i] ;
  } ;
  
  /// Returns the \c Dim_tbl  associated with the grid.
  const Dim_tbl* get_dim_tbl() const {
    return &dim ;
  } ;

  /// Read-only of a particular value of the coordinate \e z (or \e r ) at the nodes
  double get_zr(const int i) const {
    assert (i>= -nfantome) ;
    assert (i<dim.dim[0]+nfantome) ;
    
    return zr->t[i+nfantome] ;
  } ;
  
  /// Read-only of a particular value of the coordinate \e z (or \e r ) at the interfaces
  double get_zri(const int i) const {
    assert (i>= -nfantome) ;
    assert (i<dim.dim[0]+nfantome+1) ;
    
    return zri->t[i+nfantome] ;
  } ;
  
  
  // Outputs
  // -------
 public:
  /// Save in a file
  virtual void sauve(FILE *) const ;	    
  
  /// Display
  friend ostream& operator<<(ostream& , const Grille_val& ) ;	

 protected:
  /// Operator >> (virtual function called by the operator <<). 
  virtual ostream& operator>>(ostream& ) const ;    
  
  // Interpolation
  // -------------
 public:
  /**
   * Checks if the spectral grid and mapping are compatible with the
   * \c Grille_val  caracteristics for the interpolation to be done.
   * It checks wether the spectral grid is included in the Godunov one,
   * if the numbers of dimensions are the same (1,2 or 3D), and if the
   * spectral collocation points in \f$\theta\f$ and \f$\phi\f$ are well
   * defined across all the domains (see the documentation of \c
   * Tbl_val ).
   **/
  virtual bool compatible(const Map* mp, const int lmax, const int lmin=0) 
    const = 0 ;

  /**
   * Performs 1D interpolation.
   * @param rdep [input] the coordinates \e r of the source points
   * @param rarr [input] the coordinates \e r of the destination points
   * @param fdep [input] values of the function at the source points
   * @param flag [input] = 1 used for \c INSMTS  -- ought to disappear
   * @param type_inter [input] type of interpolation (see \c Tbl_val )
   * @return \c Tbl 1D of the same size as rarr, containing the values
   * of the function at destination points
   **/
  Tbl interpol1(const Tbl& rdep, const Tbl& rarr, const Tbl& fdep, 
	       int flag, const int type_inter) const ;

  /** 
   * Performs 2D interpolation.
   * @param fdep [input] values of the function at the source points
   * @param rarr [input] the coordinates \e r of the destination points
   * @param tetarr [input] the coordinates \f$\theta\f$ of the destination points
   * @param type_inter [input] type of interpolation (see \c Tbl_val )
   * @return Tbl 2D size1:that of rarr, size2: that of tetarr,
   * containing the values of the function at destination points.
   **/
  virtual Tbl interpol2(const Tbl& fdep, const Tbl& rarr, 
		       const Tbl& tetarr, const int type_inter) const = 0 ;
  /** 
   * Performs 3D interpolation.
   * @param fdep [input] values of the function at the source points
   * @param rarr [input] the coordinates \e r of the destination points
   * @param tetarr [input] the coordinates \f$\theta\f$ of the destination points
   * @param phiarr [input] the coordinates \f$\phi\f$ of the destination points
   * @param type_inter [input] type of interpolation (see \c Tbl_val )
   * @return Tbl 3D size1:that of rarr, size2: that of tetarr,
   * size3: that of phiarr, containing the values of the function at 
   * destination points.
   **/
  virtual Tbl interpol3(const Tbl& fdep, const Tbl& rarr, 
		       const Tbl& tetarr, const Tbl& phiarr, 
		       const int type_inter) const = 0 ;

  /**
   * Checks if \c Grille_val  is contained inside the spectral grid/mapping
   * within the domains [lmin, lmax[, if the numbers of dimensions are 
   * the same (1,2 or 3D), and if the symmetries are compatible.
   **/
  virtual bool contenue_dans(const Map& mp, const int lmax, const int lmin=0) 
    const = 0 ;

 protected:
  /**
   * Makes the sommation of the spectral basis functions to know
   * the values of the function described by the \c Scalar meudon 
   * at the points of the 1D Godunov grid \c this .
   * The result is an array \c t of size \c taille of all the values of 
   * the function at \c this grid points. 
   */
  void somme_spectrale1(const Scalar& meudon, double* t, int taille) const ;

  /// Same as before but for the 2D case
  virtual void somme_spectrale2(const Scalar& meudon, double* t, int taille) const = 0 ;

  /// Same as before but for the 3D case
  virtual void somme_spectrale3(const Scalar& meudon, double* t, int taille) const = 0 ;
  
};
ostream& operator<<(ostream& , const Grille_val& ) ;	


		    //------------------------------------//
		    //		class Gval_cart		  //
		    //------------------------------------//

/**
 * Class for cartesian Godunov-type grids.\ingroup (mdm)
 *         
 * Can be used for 1D (only z-coordinate), 2D (x and z) or 3D (y,x and z)
 * grids. The coordinates of the nodes are stored in \c Tbl 's zr
 * (derived from \c Grille_val ),x and y. 
 * The coordinates of the interfaces are stored in \c Tbl 's zri 
 * (derived from \c Grille_val ),xi and yi.
 * The standard constructors only allow for equally-spaced nodes.
 *
 **/
class Gval_cart : public Grille_val {
  /// Arrays defined on Godunov-type grids.
  friend class Tbl_val ;

  // Data : 
  // -----
 protected:
  /// Lower boundary for \e x dimension
  double *xmin ;
  /// Higher boundary for \e x dimension
  double *xmax ;
  /// Lower boundary for \e y dimension
  double *ymin ;
  /// Higher boundary for \e y dimension
  double *ymax ;
  
 public:
  /// Arrays containing the values of coordinate \e x on the nodes
  Tbl *x ;
  /// Arrays containing the values of coordinate \e x on the interfaces
  Tbl *xi ;
  /// Arrays containing the values of coordinate \e y on the nodes
  Tbl *y ; 
  /// Arrays containing the values of coordinate \e y on the interfaces
  Tbl *yi ;

  // Constructors - Destructor
  // -------------------------

  /**
   * Standard 1D constructor.
   * @param izmin [input] lower \e z boundary
   * @param izmax [input] higher \e z boundary
   * @param n1 [input] the number of cells (without the hidden ones)
   * @param fantome [input] the number of hidden cells on each side
   */
  Gval_cart(const double izmin, const double izmax, const int n1,
	    const int fantome = 2) ; 
  
  /**
   * Standard 2D constructor.
   * @param ixmin [input] lower \e x boundary
   * @param ixmax [input] higher \e x boundary
   * @param izmin [input] lower \e z boundary
   * @param izmax [input] higher \e z boundary
   * @param nx [input] the number of cells in \e x direction 
   * (without the hidden ones)
   * @param nz [input] the number of cells in \e z direction 
   * (without the hidden ones)
   * @param type_t [input] the type of symmetry in \f$\theta\f$ (SYM, NONSYM,
   * see base class documentation)
   * @param fantome [input] the number of hidden cells on each side
   */
  Gval_cart(const double ixmin, const double ixmax, const double izmin, 
	    const double izmax, const int nx, const int nz, const int type_t, 
	    const int fantome = 2) ; 

  /**
   * Standard 3D constructor.
   * @param iymin [input] lower \e y boundary
   * @param iymax [input] higher \e y boundary
   * @param ixmin [input] lower \e x boundary
   * @param ixmax [input] higher \e x boundary
   * @param izmin [input] lower \e z boundary
   * @param izmax [input] higher \e z boundary
   * @param ny [input] the number of cells in \e y direction 
   * (without the hidden ones)
   * @param nx [input] the number of cells in \e x direction 
   * (without the hidden ones)
   * @param nz [input] the number of cells in \e z direction 
   * (without the hidden ones)
   * @param type_t [input] the type of symmetry in \f$\theta\f$ (SYM, NONSYM,
   * see base class documentation)
   * @param type_p [input] the type of symmetry in \f$\phi\f$ (SYM, NONSYM,
   * see base class documentation)
   * @param fantome [input] the number of hidden cells on each side
   */
  Gval_cart(const double iymin, const double iymax, const double ixmin, 
	    const double ixmax, const double izmin, const double izmax, 
	    const int ny, const int nx, const int nz, const int itype_t, 
	    const int itype_p, const int fantome = 2); 
  
  /// Copy constructor
  Gval_cart(const Gval_cart& ) ;		
  
  /// Constructor from a file (see \c sauve(FILE*) )
  Gval_cart(FILE* ) ;    		
  
  /// Destructor
  virtual ~Gval_cart() ;			
  
  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another Gval_cart
  void operator=(const Gval_cart&) ;	
  
  // Accessors
  // ---------
 public:
  /// Read-only of a particular value of the coordinate \e x at the nodes
  double get_x(const int i) const {
    assert (i>= -nfantome) ;
    assert (i<dim.dim[1]+nfantome) ;
    assert (dim.ndim >= 2) ;
    
    return (*x)(i+nfantome) ;
  } ;
  
  /// Read-only of a particular value of the coordinate \e y at the nodes
  double get_y(const int i) const {
    assert (i>= -nfantome) ;
    assert (i<dim.dim[2]+nfantome) ;
    assert (dim.ndim == 3) ;
    
    return (*y)(i+nfantome) ;
  } ;
  
  /// Read-only of a particular value of the coordinate \e x at the interfaces
  double get_xi(const int i) const {
    assert (i>= -nfantome) ;
    assert (i<dim.dim[1]+nfantome+1) ;
    assert (dim.ndim >= 2) ;
    
    return (*xi)(i+nfantome) ;
  } ;
  
  /// Read-only of a particular value of the coordinate \e y at the interfaces
  double get_yi(const int i) const {
    assert (i>= -nfantome) ;
    assert (i<dim.dim[2]+nfantome+1) ;
    assert (dim.ndim == 3) ;
    
    return (*yi)(i+nfantome) ;
  } ;	
  
  /// Returns the lower boundary for x
  double get_xmin() const {
    return *xmin ;
  } ;
  
  /// Returns the higher boundary for x
  double get_xmax() const {
    return *xmax ;
  } ;
  
  /// Returns the lower boundary for y
  double get_ymin() const {
    return *ymin ;
  } ;

  /// Returns the higher boundary for x
  double get_ymax() const {
    return *ymax ;
  } ;
  
  // Outputs
  // -------
 public:
  /// Save in a file
  virtual void sauve(FILE *) const ;	    
  
 protected:
  /// Operator >> (virtual function called by the operator <<). 
  virtual ostream& operator>>(ostream& ) const ;    
  
  // Interpolation
  // -------------
 public:
  /**
   * Checks if the spectral grid and mapping are compatible with the
   * \c Grille_val  caracteristics for the interpolation to be done.
   *
   * It checks wether the spectral grid is included in the Godunov one,
   * if the numbers of dimensions are the same (1,2 or 3D), and if the
   * spectral collocation points in \f$\theta\f$ and \f$\phi\f$ are well
   * defined across all the domains (see the documentation of \c Tbl_val .
   **/
  virtual bool compatible(const Map* mp, const int lmax, const int lmin=0) 
    const ;

   /** 
   * Performs 2D interpolation.
   * @param fdep [input] values of the function at the source points,
   * defined as the nodes of the Godunov grid 
   * @param rarr [input] the coordinates \e r of the destination points
   * @param tetarr [input] the coordinates \f$\theta\f$ of the destination points
   * @param type_inter [input] type of interpolation (see \c Tbl_val )
   * @return Tbl 2D size1:that of rarr, size2: that of tetarr,
   * containing the values of the function at destination points.
   **/
  virtual Tbl interpol2(const Tbl& fdep, const Tbl& rarr, 
		       const Tbl& tetarr, const int type_inter) const ;

  /**
   * Same as before, but the coordinates of source points are passed
   * explicitly (xdep, zdep).
   */
  Tbl interpol2c(const Tbl& xdep, const Tbl& zdep, const Tbl& fdep, 
		const Tbl& rarr, const Tbl& tetarr, 
		const int type_inter) const ;
  /** 
   * Performs 3D interpolation.
   * @param fdep [input] values of the function at the source points
   * @param rarr [input] the coordinates \e r of the destination points
   * @param tetarr [input] the coordinates \f$\theta\f$ of the destination points
   * @param phiarr [input] the coordinates \f$\phi\f$ of the destination points
   * @param type_inter [input] type of interpolation (see \c Tbl_val )
   * @return Tbl 3D size1:that of rarr, size2: that of tetarr,
   * size3: that of phiarr, containing the values of the function at 
   * destination points.
   **/
  virtual Tbl interpol3(const Tbl& fdep, const Tbl& rarr, 
		       const Tbl& tetarr, const Tbl& phiarr, 
		       const int type_inter) const ;
  /**
   * Checks if \c Gval_cart  is contained inside the spectral grid/mapping
   * within the domains [lmin, lmax[, if the numbers of dimensions are 
   * the same (1,2 or 3D), and if the symmetries are compatible.
   **/
  virtual bool contenue_dans(const Map& mp, const int lmax, const int lmin=0) 
    const ;

 protected:
  /**
   * Makes the sommation of the spectral basis functions to know
   * the values of the function described by the \c Scalar meudon 
   * at the points of the 2D Godunov grid \c this .
   * The result is an array \c t of size \c taille of all the values 
   * of the function at \c this grid points.
   */
  virtual void somme_spectrale2(const Scalar& meudon, double* t, int taille) const  ;

  /// Same as before but for the 3D case
  virtual void somme_spectrale3(const Scalar& meudon, double* t, int taille) const  ;
};

		    //------------------------------------//
		    //		class Gval_spher		  //
		    //------------------------------------//

/**
 * Class for spherical Godunov-type grids.\ingroup (mdm)
 *         
 * Can be used for 1D (only r-coordinate), 2D (\f$\theta\f$ and \e r ) 
 * or 3D (\f$\phi\f$, \f$\theta\f$ and \e r )
 * grids. The coordinates of the nodes are stored in \c Tbl 's zr
 * (derived from \c Grille_val ),tet and phi. 
 * The coordinates of the interfaces are stored in \c Tbl 's zri 
 * (derived from \c Grille_val ),teti and phii.
 * The standard constructors only allow for equally-spaced nodes.
 *
 **/
class Gval_spher : public Grille_val {
  /// Arrays defined on Godunov-type grids.
  friend class Tbl_val ;
  
  // Data : 
  // -----
 public:
  /// Arrays containing the values of coordinate \f$\theta\f$ on the nodes
  Tbl *tet ;
  /// Arrays containing the values of coordinate \f$\theta\f$ on the interfaces
  Tbl *teti ;
  /// Arrays containing the values of coordinate \f$\phi\f$ on the nodes
  Tbl *phi ;
  /// Arrays containing the values of coordinate \f$\phi\f$ on the interfaces
  Tbl *phii ;

  // Constructors - Destructor
  // -------------------------

  /**
   * Standard 1D constructor.
   * @param irmin [input] lower \e r boundary
   * @param irmax [input] higher \e r boundary
   * @param nr [input] the number of cells (without the hidden ones)
   * @param fantome [input] the number of hidden cells on each side
   */
  Gval_spher(const double irmin, const double irmax, const int nr,
	     const int fantome = 2) ; 
  
  /**
   * Standard 2D constructor.
   * @param irmin [input] lower \e r boundary
   * @param irmax [input] higher \e r boundary
   * @param nt [input] the number of cells in \f$\theta\f$ (without hidden ones)
   * @param nr [input] the number of cells in r(without the hidden ones)
   * @param type_t [input] the type of symmetry in \f$\theta\f$ (SYM, NONSYM,
   * see base class documentation)
   * @param fantome [input] the number of hidden cells on each side
   */
  Gval_spher(const double irmin, const double irmax, const int nt, const 
	     int nr, const int type_t, const int fantome = 2) ; 

  /**
   * Standard 3D constructor.
   * @param irmin [input] lower \e r boundary
   * @param irmax [input] higher \e r boundary
   * @param np [input] the number of cells in \f$\phi\f$ (without the hidden ones)
   * @param nt [input] the number of cells in \f$\theta\f$ (without hidden ones)
   * @param nr [input] the number of cells in r(without the hidden ones)
   * @param type_t [input] the type of symmetry in \f$\theta\f$ (SYM, NONSYM,
   * see base class documentation)
   * @param type_p [input] the type of symmetry in \f$\phi\f$ (SYM, NONSYM,
   * see base class documentation)
   * @param fantome [input] the number of hidden cells on each side
   */
  Gval_spher(const double irmin, const double irmax, const int np, const 
	     int nt, const int nr, const int itype_t, const int itype_p, 
	     const int fantome = 2); 
  
  /// Copy constructor
  Gval_spher(const Gval_spher& ) ;		
  
  /// Constructor from a file (see \c sauve(FILE*) )
  Gval_spher(FILE* ) ;    		
  
  /// Destructor
  virtual ~Gval_spher() ;			
  
  // Mutators / assignment
  // ---------------------
 public:

  /// Assignment to another Gval_spher
  void operator=(const Gval_spher&) ;	
  
  // Accessors
  // ---------
 public:

  /// Read-only of a particular value of the coordinate \f$\theta\f$ at the nodes
  double get_tet(const int i) const {
    assert (i>= -nfantome) ;
    assert (dim.ndim >= 2) ;
    assert (i<dim.dim[1]+nfantome) ;
    
    return tet->t[i+nfantome] ;
  } ;
  
  /// Read-only of a particular value of the coordinate \f$\phi\f$ at the nodes
  double get_phi(const int i) const {
    assert (i>= -nfantome) ;
    assert (dim.ndim == 3) ;
    assert (i<dim.dim[2]+nfantome) ;
    
    return phi->t[i+nfantome] ;
  } ;
  
  /// Read-only of a particular value of coordinate \f$\theta\f$ at the interfaces
  double get_teti(const int i) const {
    assert (i>= -nfantome) ;
    assert (dim.ndim >= 2) ;
    assert (i<dim.dim[1]+nfantome+1) ;
    
    return teti->t[i+nfantome] ;
  } ;
  
  /// Read-only of a particular value of coordinate \f$\phi\f$ at the interfaces
  double get_phii(const int i) const {
    assert (i>= -nfantome) ;
    assert (dim.ndim == 3) ;
    assert (i<dim.dim[2]+nfantome+1) ;
    
    return phii->t[i+nfantome] ;
  } ;	
  
  // Outputs
  // -------
 public:
  /// Save in a file
  virtual void sauve(FILE *) const ;	    
  
 protected:
  /// Operator >> (virtual function called by the operator <<). 
  virtual ostream& operator>>(ostream& ) const ;    
  
  // Interpolation
  // -------------
 public:
  /**
   * Checks if the spectral grid and mapping are compatible with the
   * \c Grille_val  caracteristics for the interpolation to be done.
   *
   * It checks wether the spectral grid is included in the Godunov one,
   * if the numbers of dimensions are the same (1,2 or 3D), and if the
   * spectral collocation points in \f$\theta\f$ and \f$\phi\f$ are well
   * defined across all the domains (see the documentation of \c Tbl_val .
   **/
  virtual bool compatible(const Map* mp, const int lmax, const int lmin=0) 
    const ;

   /** 
   * Performs 2D interpolation.
   * @param fdep [input] values of the function at the source points,
   * defined as the nodes of the Godunov grid 
   * @param rarr [input] the coordinates \e r of the destination points
   * @param tetarr [input] the coordinates \f$\theta\f$ of the destination points
   * @param type_inter [input] type of interpolation (see \c Tbl_val )
   * @return Tbl 2D size1:that of rarr, size2: that of tetarr,
   * containing the values of the function at destination points.
   **/
  virtual Tbl interpol2(const Tbl& fdep, const Tbl& rarr, const Tbl& tetarr, 
		       const int type_inter) const ;
   /** 
   * Performs 3D interpolation.
   * @param fdep [input] values of the function at the source points
   * defined as the nodes of the Godunov grid \c this 
   * @param rarr [input] the coordinates \e r of the destination points
   * @param tetarr [input] the coordinates \f$\theta\f$ of the destination points
   * @param phiarr [input] the coordinates \f$\phi\f$ of the destination points
   * @param type_inter [input] type of interpolation (see \c Tbl_val )
   * @return Tbl 3D size1:that of rarr, size2: that of tetarr,
   * size3: that of phiarr, containing the values of the function at 
   * destination points.
   **/
 virtual Tbl interpol3(const Tbl& fdep, const Tbl& rarr, const Tbl& tetarr, 
		       const Tbl& phiarr, const int type_inter) const ;

  /**
   * Checks if \c Gval_spher  is contained inside the spectral grid/mapping
   * within the domains [lmin, lmax[, if the numbers of dimensions are 
   * the same (1,2 or 3D), and if the symmetries are compatible.
   **/
 virtual bool contenue_dans(const Map& mp, const int lmax, const int lmin=0) 
    const ;

 protected:
  /**
   * Makes the sommation of the spectral basis functions to know
   * the values of the function described by the \c Scalar meudon 
   * at the points of the 2D Godunov grid \c this .
   * The result is an array \c t of size \c taille of all the values 
   * of the function at \c this grid points.
   */
  virtual void somme_spectrale2(const Scalar& meudon, double* t, int taille) const  ;

  /// Same as before but at radial grid interfaces
  double* somme_spectrale2ri(const Scalar& meudon) const ;

  /// Same as before but at angular grid interfaces
  double* somme_spectrale2ti(const Scalar& meudon) const ;

  /// Same as before but for the 3D case
  virtual void somme_spectrale3(const Scalar& meudon, double* t, int taille) const  ;

  void initialize_spectral_r(const Map& mp, const Base_val& base, int*& idom, 
			     double*& chebnri) const ;
  void initialize_spectral_theta(const Map& mp, const Base_val& base, 
				 double*& tetlj) const ;
  void initialize_spectral_phi(const Map& mp, const Base_val& base,
			       double*& expmk) const ;

 };

}
#endif


