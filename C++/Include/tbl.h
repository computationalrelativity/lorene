/*
 *  Definition of Lorene class Tbl
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


#ifndef	__TBL_H_
#define	__TBL_H_

/*
 * Classe de tableaux de dimension determinee (actuellement 1, 2 et 3)
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.7  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.6  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.5  2002/09/24 10:49:42  e_gourgoulhon
 *
 * Modif commentaires.
 *
 * Revision 1.4  2002/09/24 08:33:35  e_gourgoulhon
 *
 * Added constructor from Matrice
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
 * Revision 2.11  1999/12/02  17:52:59  phil
 * *** empty log message ***
 *
 * Revision 2.10  1999/11/23  13:31:27  eric
 * Le constructeur Tbl::Tbl(const Dim_tbl ) devient Tbl::Tbl(const Dim_tbl& ).
 * Le constructeur Tbl::Tbl(const Grille3d* ) devient
 *   Tbl(const Grille3d& ).
 *
 * Revision 2.9  1999/11/15  16:36:08  eric
 * Le membre dim est desormais un Dim_tbl et non plus un pointeur sur un
 * Dim_tbl.
 *
 * Revision 2.8  1999/10/29  15:04:07  eric
 * Suppression des fonctions membres min() et max():
 * elles deviennent des fonctions externes.
 * Ajout de fonctions mathematiques (abs, norme, etc...).
 *
 * Revision 2.7  1999/10/18  15:05:48  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 *
 * Revision 2.6  1999/10/01  10:18:18  eric
 * Amelioration des commentaires.
 *
 * Revision 2.5  1999/09/30  12:50:59  eric
 *  Constructeur a 1 parametre rendu explicit.
 * Amelioration des commentaires.
 *
 * Revision 2.4  1999/09/24  14:23:01  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.3  1999/03/02  18:54:36  eric
 * Ajout de la fonction affiche_seuil.
 *
 * Revision 2.2  1999/02/23  10:40:55  hyc
 * *** empty log message ***
 *
 *
 * Revision 2.0  1998/12/02  11:11:16  hyc
 * Version 2
 *
 * $Header$
 *
 */


// Fichiers includes
#include <assert.h>
#include <stdlib.h>

#include "type_parite.h"
#include "dim_tbl.h"

class Grille3d ;
class Matrice ;

/**
 * Basic array class.
 * 
 * This class is essentially an {\tt double} array class. The elements of the 
 * array are stored continuously using the C convention.
 * A {\tt Tbl} is initialy created with a {\sl logical} state 
 * {\tt ETATNONDEF} (i.e. undefined),  except by the copy constructor and
 * the constructor from a file. 
 * The general logical state of an initialized 
 * {\tt Tbl} is {\tt ETATQCQ}; it is the only state for which the memory
 * allocation is performed for the {\tt double} array {\tt t}. 
 * The value zero is treated as a special logical state ({\tt ETATZERO}), 
 * without any memory allocation. 
 * Arithmetic operations are provided with the usual meaning (see 
 * below).
 * 
 * @version #$Id$#
 */
class Tbl {

  friend class Matrice ;

    // Data : 
    // -----
    private:
	/// logical state ({\tt ETATNONDEF}, {\tt ETATQCQ} or {\tt ETATZERO}).
	int etat ;  

    public:
    	Dim_tbl dim ;	    /// Number of dimensions, size,...
	double*	t ;	    /// The array of {\tt double}

    // Constructors - Destructor
    // -------------------------
	
    public:
	/**
	 * 1D constructor 
	 * @param size0  [input] Number of elements of the array {\tt t}.
	 *		  Will be assigned to {\tt dim.dim[0]}.  
	 */
	explicit Tbl(int size0) ;			

	/**
	 * 2D constructor
	 * @param size1  [input] Defines the range [0, size1-1] of the 
	 *		  outermost index in the storage of the array {\tt t}. 
	 *		  Will be assigned to {\tt dim.dim[1]}.
	 * @param size0  [input] Defines the range [0, size0-1] of the 
	 *		  innermost index in the storage of the array {\tt t}. 
	 *		  Will be assigned to {\tt dim.dim[0]}.
	 */ 
	Tbl(int size1, int size0) ;		
	
	/**
	 * 3D constructor
	 * @param size2  [input] Defines the range [0, size2-1]  of the 
	 *		  outermost index in the storage of the array {\tt t}. 
	 *		  Will be assigned to {\tt dim.dim[2]}.
	 * @param size1  [input] Defines the range [0, size1-1] of the 
	 *		  intermediate index in the storage of the array {\tt t}. 
	 *		  Will be assigned to {\tt dim.dim[1]}.
	 * @param size0  [input] Defines the range [0, size0-1] of the 
	 *		  innermost index in the storage of the array {\tt t}. 
	 *		  Will be assigned to {\tt dim.dim[0]}.
	 */ 
	Tbl(int size2, int size1, int size0) ;		
	
	explicit Tbl(const Grille3d& grid) ; /// Constructor from a 3D grid
	explicit Tbl(const Dim_tbl& dim) ;	 /// Constructor from a {\tt Dim\_tbl}
	/// Constructor from a file (see {\tt sauve(FILE* )})
	explicit Tbl(FILE* ) ;	
	Tbl(const Tbl& a) ;		/// Copy constructor
	
	/** Constructor from a matrix.
	 *  If the matrix has only one row or one column, the {\tt Tbl}
	 *  is 1D, otherwise it is 2D.
	 */
	explicit Tbl(const Matrice& mat) ;
	

	~Tbl() ;			/// Destructor

    // Assignement
    // -----------
	void operator=(const Tbl& ) ;	/// Assignment to another {\tt Tbl}
	void operator=(double ) ; /// Assignment to a {\tt double}
	void operator=(int ) ;	 /// Assignment to a {\tt int}

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
     * Sets the {\tt Tbl} to zero in a hard way. 
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
	    assert( dim.ndim == 1 ) ; 
	    assert( i >= 0 ) ; 
	    assert( i < dim.dim[0] ) ;
	    return t[i] ;
	} ;
	
	/// Read-only of a particular element (index {\tt i}) (1D case)
	double operator()(int i) const {
	    assert(etat != ETATNONDEF) ;
	    assert( dim.ndim == 1 ) ; 
	    assert( i >= 0 ) ; 
	    assert( i < dim.dim[0] ) ;
	    if (etat == ETATZERO) {
		double zero = 0. ;
		return zero ; 
	    }
	    else return t[i] ;
	};

	/// Read/write of a particular element (index {\tt (j,i)}) (2D case)
	double& set(int j, int i) {
	    assert (etat == ETATQCQ) ;
	    assert( dim.ndim == 2 ) ;
	    assert( (i>=0) && (i<dim.dim[0]) ) ;
	    assert( (j>=0) && (j<dim.dim[1]) ) ;
	    return t[dim.dim[0] * j + i] ;
	};

	/// Read-only of a particular element (index {\tt (j,i)}) (2D case)
	double operator()(int j, int i) const {
	    assert(etat != ETATNONDEF) ;
	    assert( dim.ndim == 2 ) ;
	    assert( (i>=0) && (i<dim.dim[0]) ) ;
	    assert( (j>=0) && (j<dim.dim[1]) ) ;
	    if (etat == ETATZERO) {
		double zero = 0. ;
		return zero ; 
	    }
	    else return t[dim.dim[0] * j + i] ;
	};

	/// Read/write of a particular element (index {\tt (k,j,i)}) (3D case)
	double& set(int k, int j, int i) {	
	    assert (etat == ETATQCQ) ;
	    assert( dim.ndim == 3 ) ;
	    assert( (i>=0) && (i<dim.dim[0]) ) ;
	    assert( (j>=0) && (j<dim.dim[1]) ) ;
	    assert( (k>=0) && (k<dim.dim[2]) ) ;
	    return t[dim.dim[1]*dim.dim[0]*k + dim.dim[0]*j + i] ;
	};

	/// Read-only of a particular element (index {\tt (k,j,i)}) (3D case)
	double operator()(int k, int j, int i) const {
	    assert(etat != ETATNONDEF) ;
	    assert( dim.ndim == 3 ) ;
	    assert( (i>=0) && (i<dim.dim[0]) ) ;
	    assert( (j>=0) && (j<dim.dim[1]) ) ;
	    assert( (k>=0) && (k<dim.dim[2]) ) ;
	    if (etat == ETATZERO) {
		double zero = 0. ;
		return zero ; 
	    }
	    else return t[dim.dim[1]*dim.dim[0]*k + dim.dim[0]*j + i] ;
	};

    // Extraction of information
    // -------------------------
	int get_etat() const { return etat ; };	    /// Gives the logical state

	/// Gives the total size (ie {\tt dim.taille})
	int get_taille() const { return dim.taille ; };

	/// Gives the number of dimensions (ie {\tt dim.ndim})
	int get_ndim() const { return dim.ndim ; };
	
	/// Gives the {\tt i}th dimension (ie {\tt dim.dim[i]})
	int get_dim(int i) const {	
	    assert( (i>=0) && (i<dim.ndim) ) ;
	    return dim.dim[i] ;
	};
	
    // Outputs
    // -------
    public:
	void sauve(FILE* ) const ;	/// Save in a file

	/** Prints only the values greater than a given threshold.
	 *   @param ostr [input] Output stream used for the printing 
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param threshold [input] Value above which an array element is printed
	 *    (default: 1.e-7)
	 */
	void affiche_seuil(ostream& ostr, int precision = 4, 
			   double threshold = 1.e-7) const ;
	/// Display   
	friend ostream& operator<<(ostream& , const Tbl& ) ;	

    // Member arithmetics
    // ------------------
    public:
	
	void operator+=(const Tbl &) ;	/// Addition of a {\tt Tbl} to {\tt this}
	void operator+=(double) ;	/// Addition of a {\tt double} to {\tt this}
	void operator-=(const Tbl &) ;	/// Subtraction of a {\tt Tbl} to {\tt this}
	void operator-=(double) ;	/// Subtraction of a {\tt double} to {\tt this}
	void operator*=(const Tbl &) ;	/// Multiplication of {\tt this} by a {\tt Tbl}
	void operator*=(double) ;	/// Multiplication of {\tt this} by a {\tt double}
	void operator/=(const Tbl &) ;	/// Division of {\tt this} by a {\tt Tbl}
	void operator/=(double) ;	/// Division of {\tt this} by a {\tt double}

} ;
ostream& operator<<(ostream& , const Tbl& ) ;	


/**
 * @name Tbl Mathematics
 */
    //@{
Tbl operator+(const Tbl&) ;			/// + Tbl
Tbl operator-(const Tbl&) ;			/// - Tbl
Tbl operator+(const Tbl&, const Tbl&) ;		/// Tbl + Tbl
Tbl operator+(const Tbl&, double) ;		/// Tbl + double
Tbl operator+(double, const Tbl&) ;		/// double + Tbl
Tbl operator+(const Tbl&, int) ;		/// Tbl + int
Tbl operator+(int, const Tbl&) ;		/// int + Tbl
Tbl operator-(const Tbl&, const Tbl&) ;		/// Tbl - Tbl
Tbl operator-(const Tbl&, double) ;		/// Tbl - double
Tbl operator-(double, const Tbl&) ;		/// double - Tbl
Tbl operator-(const Tbl&, int) ;		/// Tbl - int
Tbl operator-(int, const Tbl&) ;		/// int - Tbl
Tbl operator*(const Tbl&, const Tbl&) ;		/// Tbl * Tbl
Tbl operator*(const Tbl&, double) ;		/// Tbl * double
Tbl operator*(double, const Tbl&) ;		/// double * Tbl
Tbl operator*(const Tbl&, int) ;		/// Tbl * int
Tbl operator*(int, const Tbl&) ;		/// int * Tbl
Tbl operator/(const Tbl&, const Tbl&) ;		/// Tbl / Tbl
Tbl operator/(const Tbl&, double) ;		/// Tbl / double
Tbl operator/(double, const Tbl&) ;		/// double / Tbl
Tbl operator/(const Tbl&, int) ;		/// Tbl / int
Tbl operator/(int, const Tbl&) ;		/// int / Tbl

Tbl sin(const Tbl& ) ;	    /// Sine
Tbl cos(const Tbl& ) ;	    /// Cosine
Tbl tan(const Tbl& ) ;	    /// Tangent
Tbl asin(const Tbl& ) ;	    /// Arcsine
Tbl acos(const Tbl& ) ;	    /// Arccosine
Tbl atan(const Tbl& ) ;	    /// Arctangent
Tbl exp(const Tbl& ) ;	    /// Exponential
Tbl log(const Tbl& ) ;	    /// Neperian logarithm
Tbl log10(const Tbl& ) ;    /// Basis 10 logarithm
Tbl sqrt(const Tbl& ) ;	    /// Square root
Tbl racine_cubique (const Tbl&) ; /// cube root
Tbl pow(const Tbl& , int ) ;  /// Power ${\tt Tbl}^{\tt int}$
Tbl pow(const Tbl& , double ) ; /// Power ${\tt Tbl}^{\tt double}$
Tbl abs(const Tbl& ) ;	    /// Absolute value
double max(const Tbl& ) ;   /// Maximum value of the {\tt Tbl} elements
double min(const Tbl& ) ;   /// Minimum value of the {\tt Tbl} elements

/// Sum of the absolute values of all the {\tt Tbl} elements
double norme(const Tbl& ) ;   

/**
 * Relative difference between two {\tt Tbl} (norme version).
 * Returns {\tt norme(a-b)/norme(b)} unless {\tt b=0}, in which
 * case it returns {\tt norme(a-b)}.
 */
double diffrel(const Tbl& a, const Tbl& b) ; 

/**
 * Relative difference between two {\tt Tbl} (max version).
 * Returns {\tt max(abs(a-b))/max(abs(b))} unless {\tt b=0}, in which
 * case it returns {\tt max(abs(a-b))}.
 */
double diffrelmax(const Tbl& a, const Tbl& b) ; 

    //@}


#endif
