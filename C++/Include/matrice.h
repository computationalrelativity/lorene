/*
 *  Definition of Lorene class Matrice
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


#ifndef __MATRICE_H_
#define __MATRICE_H_

/*
 * $Id$
 * $Log$
 * Revision 1.9  2004/12/29 12:27:35  j_novak
 * permute is now a Itbl* which array is sent directly to the LAPACK routines.
 * It is now possible to solve a general system (i.e. even if the Matrice
 * is not in a banded form).
 *
 * Revision 1.8  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.7  2004/03/22 13:12:42  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.6  2002/09/24 10:49:41  e_gourgoulhon
 *
 * Modif commentaires.
 *
 * Revision 1.5  2002/09/24 08:34:12  e_gourgoulhon
 *
 * Added member function transpose()
 * and matrix multiplication
 *
 * Revision 1.4  2002/09/13 09:17:33  j_novak
 * Modif. commentaires
 *
 * Revision 1.3  2002/06/17 14:05:17  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.2  2002/01/03 13:18:40  j_novak
 * Optimization: the members set(i,j) and operator(i,j) of class Matrice are
 * now defined inline. Matrice is a friend class of Tbl.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.12  1999/11/30  17:45:28  phil
 * changement de prototypage
 *
 * Revision 2.11  1999/10/25  08:02:11  eric
 * Changements commentaires.
 *
 * Revision 2.10  1999/10/12  16:04:06  phil
 * Doc
 *
 * Revision 2.9  1999/10/12  15:59:09  phil
 * Documentation
 *
 * Revision 2.8  1999/10/12  09:42:07  phil
 * retour
 *
 * Revision 2.7  1999/10/12  09:39:38  phil
 * passage en const
 *
 * Revision 2.6  1999/10/11  09:35:40  phil
 * changement prototypage arithmetique
 *
 * Revision 2.5  1999/10/05  17:02:32  phil
 * ajout de determinant et val_propre
 *
 * Revision 2.4  1999/09/20  11:25:25  phil
 * passage en Doc++
 *
 * Revision 2.3  1999/09/20  11:18:53  phil
 * *** empty log message ***
 *
 * Revision 2.2  1999/04/13  13:56:11  phil
 * suppression de proto.h
 *
 * Revision 2.1  1999/04/07  14:53:38  phil
 * Changement de prototypage
 *
 * Revision 2.0  1999/04/07  14:03:59  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//fichiers includes
#include <stdio.h>

#include "type_parite.h"
#include "tbl.h"
#include "itbl.h"

/**
 * Matrix handling.
 * The matrix can be stored in the usual way in \c std,  in a band-way by
 * \c band and on a LU-decomposition by the two arrays \c lu and 
 * \c permute. All the storage conventions are those af \b LAPACK which is
 * used to make the LU-decomposition,  the inversion and to compute the
 * eigenvalues of the matrix. All those representations are redondant, that means
 * that doing the LU-decomposition, for example,  does NOT destroy 
 * previously calculated type of storage.
 * 
 * \ingroup (util)
 */

class Matrice {
    //Elements
    private:
    
	int etat ; ///< logical state \c (ETATZERO, ETATQCQ or ETATNONDEF)
	
	Tbl* std ; ///< Pointer on the array of the standard representation.
	
	
	mutable int ku ;    ///< Number of upper-diagonals in the band representation.
	mutable int kl ;    ///< Number of lower-diagonals in the band representation.
	
	/**
	 * Pointer on the array of the band representation of a square matrix.
	 * To be precise, \f$ A(i, j)\f$ is stored in \c band
	 * \f$[ku+1+i-j, j]\f$ for \f$\mathrm {max}(1, j-ku) \leq i \leq
	 * \mathrm{min} (n, j+kl)\f$, \e n being the size of the matrix.
	 */
	mutable Tbl* band ;
	
	
	mutable Tbl* lu ;   ///< Pointer on the first array of the LU-representation.
	mutable Itbl* permute ;	///< Pointer on the second array of the LU-representation.
	
    // Constructeurs destructeurs
    public:
	/**
	 * Standard constructor.
	 * All the representations are set to \c ETATNONDEF.
	 * @param size1 [input] number of lines.
	 * @param size2 [input] number of columns.
	 */
	Matrice (int size1, int size2 ) ;
	
	Matrice (const Matrice& ) ; ///< Constructor by copy.
	
	/**
	 * Constructor from a \c Tbl.
	 * @param tab [input]  2-dimension or 1-dimension array
	 *
        * If \c tab is a 1-dimension \c Tbl, a single-column matrix is created,
	 *  otherwise \c *std is simply constructed by a \c Tbl copy of \c tab.
	 *
	 */
	Matrice (const Tbl& tab) ;

	~Matrice() ; ///< Destructor
    
        //Gestion memoire
	/**
	 * Logical destructor : dellocates the memory of the various used 
	 * representations.
	 * 
	 */
    private:
	void del_t() ;

    // manipulation des etats
    public:
	/// Returns the logical state.
	int get_etat() const { return etat ; }; 

	/**
	 * Sets the logical state to \c ETATQCQ (ordinary state).
	 * The state of \c *std is now \c ETATQCQ and the one of all the
	 * other representations is \c ETATNONDEF.
	 */
	void set_etat_qcq()  ;

	/**
	 * Sets the logical state to \c ETATZERO (zero).
	 * The state of \c *std is now \c ETATZERO and the one of all the
	 * other representations is \c ETATNONDEF.
	 */
	void set_etat_zero() ;

	/**
	 * Sets the logical state to \c ETATNONDEF (undefined state).
	 * The state of of all the representations is now \c ETATNONDEF.
	 */
	 void set_etat_nondef() ;

    public:
	/**
	 * Returns the dimension of the matrix.
	 * @param i [input] if \e i =0 returns the number of lines and if \e i =2 
	 * returns the number of columns.
	 */
	int get_dim(int i) const ;
	
	/// Returns the array of matrix elements
	Tbl get_array() const {return *std; } ;	
	
    // affectation
    public:	
	/**
	 * Sets all the element of \c *std to \e x.
	 * The other representations are set to \c ETATNONDEF.
	 */
	void operator=(double x) ;

	void operator=(const Matrice& ) ; ///< Assignement to another \c Matrice.
 
    //Impression
	friend ostream& operator<<(ostream& , const Matrice& ) ; ///< Display
    	
    // extraction d'un element :
    public:
	/**
	 * Read/write of a particuliar element.
	 * This is done in \c *std and all the other representations are no
	 * longer valid.
	 * @param j [input] line coordinate.
	 * @param i [input] column coordinate.
	 * 
	 */
	double& set(int j, int i) {
	  assert (etat == ETATQCQ) ;
	  assert ((i>=0) && (i<std->dim.dim[0])) ;
	  assert( (j>=0) && (j<std->dim.dim[1]) ) ;
	  if (band != 0x0) {
	    delete band ;
	    band = 0x0 ;
	  }
	  
	  if (lu != 0x0) {
	    delete lu ;
	    delete permute ;
	    lu = 0x0 ;
	    permute = 0x0 ;
	  }
	  return std->t[std->dim.dim[0] * j + i] ;
	} ;
	
	/**
	 * Read-only of a particuliar element.
	 * @param j [input] line coordinate.
	 * @param i [input] column coordinate. 
	 */	
	double operator()(int j , int i) const {
	  assert(etat != ETATNONDEF) ;
	  assert( (i>=0) && (i<std->dim.dim[0]) ) ;
	  assert( (j>=0) && (j<std->dim.dim[1]) ) ;
	  if (etat == ETATZERO) {
		double zero = 0. ;
		return zero ; 
	  }
	  else return std->t[std->dim.dim[0] * j + i] ;
	};
	
    // Passage matrice a bande
	/**
	 * Calculate the band storage of \c *std.
	 * Please note that this function does NOT check if \c *std 
	 * represents a real band-matrix.
	 * @param up [input] number of upper-diagonals.
	 * @param low [input] number of lower-diagonals.
	 */
	void set_band (int up, int low) const ;
	
    // Decomposition LU
	/**
	 * Calculate the LU-representation,  assuming the band-storage has been
	 * done. The calculus is done using \b LAPACK.
	 */
	void set_lu () const ;
        
    // Inversion de la matrice
	/**
	 * Solves the linear system represented by the matrix.
	 * The calculus assumes the the LU-decomposition has been done and is
	 * conducted using \b LAPACK.
	 * @param sec_membre [input] the right-hand side of the system.
	 */
	Tbl inverse (const Tbl& sec_membre) const ;

    // Les valeurs propres :
	/**
	 * Returns the eigenvalues of the matrix, calculated using \b LAPACK.
	 * @return contains the real and the imaginary parts of the 
	 * eigenvalues. The real parts are in \c Tbl \e (0, *) and 
	 * the imaginary parts in \c Tbl \e (1, *).
	 */
	Tbl val_propre() const ;	
	/**
	 * Returns the eigenvectors of the matrix, calculated using \b LAPACK.
	 *
	 */
	Matrice vect_propre() const ;
	
	/**
	 * Computes the determinant of the matrix, using \b LAPACK and the
	 * standard decomposition.
	 */
	double determinant() const ;

	/**
	 * Computes the transpose matrix
	 */
	Matrice transpose() const ;
	
    // Operateurs amis
	friend Matrice operator+ (const Matrice&, const Matrice& ) ;
	friend Matrice operator- (const Matrice&, const Matrice& ) ;
	friend Matrice operator* (const Matrice&, double ) ;
	friend Matrice operator* (double, const Matrice& ) ;
	friend Matrice operator* (const Matrice&, const Matrice& ) ;
	friend Matrice operator/ (const Matrice&,  double ) ;
} ;
ostream& operator<<(ostream& , const Matrice& ) ; 

/**
 * \defgroup mat_mat Matrice arithmetics
 * \ingroup (util)
 * @{
 */
Matrice operator+ (const Matrice&, const Matrice& ) ; ///< \c Matrice + \c Matrice
Matrice operator- (const Matrice&, const Matrice& ) ; ///< \c Matrice - \c Matrice
Matrice operator* (const Matrice&, double ) ;///< \c Matrice * \c double
Matrice operator* (double, const Matrice& ) ;///< \c double * \c Matrice
Matrice operator* (const Matrice&, const Matrice& ) ; ///< Matrix product
Matrice operator/ (const Matrice&,  double ) ; ///< \c Matrice / \c double

/**@} */

#endif	
