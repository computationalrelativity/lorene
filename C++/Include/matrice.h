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
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
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

/**
 * Matrix handling.
 * The matrix can be stored in the usual way in {\tt std},  in a band-way by
 * {\tt band} and on a LU-decomposition by the two arrays {\tt lu} and 
 * {\tt permute}. All the storage conventions are those af {\bf LAPACK} which is
 * used to make the LU-decomposition,  the inversion and to compute the
 * eigenvalues of the matrix. All those representations are redondant, that means
 * that doing the LU-decomposition, for example,  does NOT destroy 
 * previously calculated type of storage.
 * 
 * @version #$Id$#
 * 
 */

class Matrice {
    //Elements
    private:
    
	int etat ; /// logical state {\tt (ETATZERO, ETATQCQ or ETATNONDEF)}
	
	Tbl* std ; /// Pointer on the array of the standard representation.
	
	
	mutable int ku ;    /// Number of upper-diagonals in the band representation.
	mutable int kl ;    /// Number of lower-diagonals in the band representation.
	
	/**
	 * Pointer on the array of the band representation of a square matrix.
	 * To be precise, {$ A(i, j)$} is stored in {\tt band}
	 * $[ku+1+i-j, j]$ for $\mathrm {max}(1, j-ku) \leq i \leq
	 * \mathrm{min} (n, j+kl)$, $n$ being the size of the matrix.
	 */
	mutable Tbl* band ;
	
	
	mutable Tbl* lu ;   /// Pointer on the first array of the LU-representation.
	mutable Tbl* permute ;	/// Pointer on the second array of the LU-representation.
	
    // Constructeurs destructeurs
    public:
	/**
	 * Standard constructor.
	 * All the representations are set to {\tt ETATNONDEF}.
	 * @param size1 [input] number of lines.
	 * @param size2 [input] number of columns.
	 */
	Matrice (int size1, int size2 ) ;
	
	Matrice (const Matrice& ) ; /// Constructor by copy.
	
	/**
	 * Constructor from a {\tt Tbl}.
	 * @param tab [input] a 2-dimension array.
	 * 
	 * {\tt *std} is set to {\tt tab} and all the other representations to 
	 * {\tt ETATNONDEF}
	 * 
	 */
	Matrice (const Tbl& tab) ;

	~Matrice() ; /// Destructor
    
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
	int get_etat() const { return etat ; }; /// Returns the logical state.
	/**
	 * Sets the logical state to {\tt ETATQCQ} (ordinary state).
	 * The state of {\tt *std} is now {\tt ETATQCQ} and the one of all the
	 * other representations is {\tt ETATNONDEF}.
	 */
	
	void set_etat_qcq()  ;
	/**
	 * Sets the logical state to {\tt ETATZERO} (zero).
	 * The state of {\tt *std} is now {\tt ETATZERO} and the one of all the
	 * other representations is {\tt ETATNONDEF}.
	 */
	void set_etat_zero() ;

	/**
	 * Sets the logical state to {\tt ETATNONDEF} (undefined state).
	 * The state of of all the representations is now {\tt ETATNONDEF}.
	 */
	 void set_etat_nondef() ;


    // void dimensions :
    public:
	/**
	 * Returns the dimension of the matrix.
	 * @param i [input] if $i=0$ returns the number of lines and if $i=2$ 
	 * returns the number of columns.
	 */
	int get_dim(int i) const ;
	
    // affectation
    public:	
	/**
	 * Sets all the element of {\tt *std} to $x$.
	 * The other representations are set to {\tt ETATNONDEF}.
	 */
	void operator=(double x) ;

	void operator=(const Matrice& ) ; /// Assignement to another {\tt Matrice}.
 
    //Impression
	friend ostream& operator<<(ostream& , const Matrice& ) ; /// Display
    	
    // extraction d'un element :
    public:
	/**
	 * Read/write of a particuliar element.
	 * This is done in {\tt *std} and all the other representations are no
	 * longer valid.
	 * @param i [input] line coordinate.
	 * @param j [input] column coordinate.
	 * 
	 */
	double& set(int i, int j) ;
	
	/**
	 * Read-only of a particuliar element.
	 * @param i [input] line coordinate.
	 * @param j [input] column coordinate. 
	 */	
	double operator()(int i , int j) const;
	
    // Passage matrice a bande
	/**
	 * Calculate the band storage of {\tt *std}.
	 * Please note that this function does NOT check if {\tt *std} 
	 * represents a real band-matrix.
	 * @param up [input] number of upper-diagonals.
	 * @param low [input] number of lower-diagonals.
	 */
	void set_band (int up, int low) const ;
	
    // Decomposition LU
	/**
	 * Calculate the LU-representation,  assuming the band-storage has been
	 * done. The calculus is done using {\bf LAPACK}.
	 */
	void set_lu () const ;
        
    // Inversion de la matrice
	/**
	 * Solves the linear system represented by the matrix.
	 * The calculus assumes the the LU-decomposition has been done and is
	 * conducted using {\bf LAPACK}.
	 * @param sec_membre [input] the right-hand side of the system.
	 */
	Tbl inverse (const Tbl& sec_membre) const ;

    // Les valeurs propres :
	/**
	 * Returns the eigenvalues of the matrix, calculated using {\bf LAPACK}.
	 * @return contains the real and the imaginary parts of the 
	 * eigenvalues. The real parts are in {\tt Tbl$(0, *)$} and 
	 * the imaginary parts in {\tt Tbl$(1, *)$}.
	 */
	Tbl val_propre() const ;
	
	/**
	 * Computes the determinant of the matrix, using {\bf LAPACK} and the
	 * standard decomposition.
	 */
	double determinant() const ;
    
    // Operateurs d'arithmetique
	friend Matrice operator+ (const Matrice&, const Matrice& ) ; /// {\tt Matrice $+$ Matrice}
	friend Matrice operator- (const Matrice&, const Matrice& ) ; /// {\tt Matrice $-$ Matrice}
	friend Matrice operator* (const Matrice&, double ) ;/// {\tt Matrice $*$ double}
	friend Matrice operator* (double, const Matrice& ) ;/// {\tt double $*$ Matrice}
	friend Matrice operator/ (const Matrice&,  double ) ; /// {\tt Matrice $/$ double}
} ;


#endif	
