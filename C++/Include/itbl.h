/*
 *  Definition of Lorene class Itbl
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


#ifndef __ITBL_H_
#define __ITBL_H_

/*
 * Classe de tableaux d'entiers de dimension determinee (actuellement 1, 2 et 3)
 *
 */

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/10/11 16:43:05  e_gourgoulhon
 *
 * IMPORTANT CHANGE: the standard constructors sets now the logical state
 *   to ETATQCQ, and no longer to ETATNONDEF.
 *
 * Revision 1.3  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  1999/11/23  13:16:12  eric
 * Le constructeur Itbl::Itbl(const Dim_tbl ) devient desormais
 *   Itbl::Itbl(const Dim_tbl& ).
 * Modif commentaires (taille 0 autorisee).
 *
 * Revision 2.1  1999/11/19  12:02:17  phil
 * /
 *
 * Revision 2.0  1999/11/17  16:04:27  phil
 * *** empty log message ***
 *
 * Revision 1.1  1999/11/17  16:03:06  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */
 
// Fichiers includes
#include <assert.h>
#include <stdlib.h>

#include "type_parite.h"
#include "dim_tbl.h"


/**
 * Basic integer array class.
 * 
 * This class is essentially an {\tt int} array class. The elements of the 
 * array are stored continuously using the C convention.
 * 
 * The general logical state of an initialized 
 * {\tt Itbl} is {\tt ETATQCQ}; it is the only state for which the memory
 * allocation is performed for the {\tt int} array {\tt t}. 
 * The value zero is treated as a special logical state ({\tt ETATZERO}), 
 * without any memory allocation. 
 *
 * Contrary to a {\tt Tbl}, an {\tt Itbl} is initialy created in a
 * logical state {\tt ETATQCQ} (i.e. ordinary state, with memory allocated
 * for the array {\tt t}).
 *
 * Arithmetic operations are provided with the usual meaning (see 
 * below).
 * 
 * A 1D Itbl can be of dimension 0 (size 0). Its logical state is then
 * {\tt ETATZERO} (by convention). 
 * 
 * @version #$Id$#.
 */
class Itbl {

    // Data : 
    // -----
    private:
	/// logical state ({\tt ETATNONDEF}, {\tt ETATQCQ} or {\tt ETATZERO}).
	int etat ;  
    
    public: 
    Dim_tbl dim ;	    /// Number of dimensions, size,...
	int*	t ;	    /// The array of {\tt int}'s

    // Constructors - Destructor
    // -------------------------
	
    public:
	/**
	 * 1D constructor.
	 * This constructor sets the {\tt Itbl} in the logical state {\tt ETATQCQ},
	 * so that it is ready for initialization via the method {\tt set(int )}.
	 * 
	 * @param size0  [input] Number of elements of the array {\tt t}.
	 *		  Will be assigned to {\tt dim.dim[0]}.
	 *		  The size 0 is allowed (the corresponding 
	 *		  logical state will be then ETATZERO).
	 */
	explicit Itbl(int size0) ;			

	/**
	 * 2D constructor.
	 * This constructor sets the {\tt Itbl} in the logical state {\tt ETATQCQ},
	 * so that it is ready for initialization via the method 
	 * {\tt set(int, int )}.
	 *
	 * @param size1  [input] Defines the range [0, size1-1] of the 
	 *		  outermost index in the storage of the array {\tt t}. 
	 *		  Will be assigned to {\tt dim.dim[1]}.
	 * @param size0  [input] Defines the range [0, size0-1] of the 
	 *		  innermost index in the storage of the array {\tt t}. 
	 *		  Will be assigned to {\tt dim.dim[0]}.
	 */ 
	Itbl(int size1, int size0) ;		
	
	/**
	 * 3D constructor
	 * This constructor sets the {\tt Itbl} in the logical state {\tt ETATQCQ},
	 * so that it is ready for initialization via the method 
	 * {\tt set(int, int, int )}.
	 *
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
	Itbl(int size2, int size1, int size0) ;		
	

	/** Constructor from a {\tt Dim\_tbl}.
	 * This constructor sets the {\tt Itbl} in the logical state {\tt ETATQCQ},
	 * so that it is ready for initialization via the method 
	 * {\tt set}.
	 */
	explicit Itbl(const Dim_tbl& ) ; 
	
	/// Constructor from a file (see {\tt sauve(FILE* )})
	explicit Itbl(FILE* ) ;	
	
	Itbl(const Itbl& ) ;		/// Copy constructor

	~Itbl() ;			/// Destructor

    // Assignement
    // -----------
	void operator=(const Itbl& ) ;	/// Assignment to another {\tt Itbl}
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
     * Deallocates the memory occupied by the {\tt int} array {\tt t}.
     */
	void set_etat_nondef() ;	
	
    /**
     * Sets the logical state to {\tt ETATZERO} (zero). 
     * Deallocates the memory occupied by the {\tt int} array {\tt t}.
     */
	void set_etat_zero() ;	    	

    /**
     * Sets the logical state to {\tt ETATQCQ} (ordinary state).
     * If the state (member {\tt etat}) is already {\tt ETATQCQ}, this 
     * function does nothing. Otherwise, it performs the memory allocation
     * for the {\tt int} array {\tt t}.  
     */
	void set_etat_qcq() ;	    	
    
    /**
     * Sets the {\tt Itbl} to zero in a hard way. 
     * 1/ Sets the logical state to {\tt ETATQCQ}, i.e. to an ordinary state.
     * 2/ Allocates the memory of the {\tt int} array {\tt t}, and fills it
     * with zeros. NB: this function must be used for debugging purposes only.
     * For other operations, the function {\tt set\_etat\_zero()} must
     * be perferred. 
     */
	void annule_hard() ;			
	
    // Access to individual elements
    // -----------------------------
    public:
	/// Read/write of a particular element (index {\tt i})  (1D case)
	int& set(int i) {
	    assert (etat == ETATQCQ) ;
	    assert( dim.ndim == 1 ) ; 
	    assert( i >= 0 ) ; 
	    assert( i < dim.dim[0] ) ;
	    return t[i] ;
	} ;
	
	/// Read-only of a particular element (index {\tt i}) (1D case)
	int operator()(int i) const {
	    assert(etat != ETATNONDEF) ;
	    assert( dim.ndim == 1 ) ; 
	    assert( i >= 0 ) ; 
	    assert( i < dim.dim[0] ) ;
	    if (etat == ETATZERO) {
		int zero = 0 ;
		return zero ; 
	    }
	    else return t[i] ;
	};

	/// Read/write of a particular element (index {\tt (j,i)}) (2D case)
	int& set(int j, int i) {
	    assert (etat == ETATQCQ) ;
	    assert( dim.ndim == 2 ) ;
	    assert( (i>=0) && (i<dim.dim[0]) ) ;
	    assert( (j>=0) && (j<dim.dim[1]) ) ;
	    return t[dim.dim[0] * j + i] ;
	};

	/// Read-only of a particular element (index {\tt (j,i)}) (2D case)
	int operator()(int j, int i) const {
	    assert(etat != ETATNONDEF) ;
	    assert( dim.ndim == 2 ) ;
	    assert( (i>=0) && (i<dim.dim[0]) ) ;
	    assert( (j>=0) && (j<dim.dim[1]) ) ;
	    if (etat == ETATZERO) {
		int zero = 0 ;
		return zero ; 
	    }
	    else return t[dim.dim[0] * j + i] ;
	};

	/// Read/write of a particular element (index {\tt (k,j,i)}) (3D case)
	int& set(int k, int j, int i) {	
	    assert (etat == ETATQCQ) ;
	    assert( dim.ndim == 3 ) ;
	    assert( (i>=0) && (i<dim.dim[0]) ) ;
	    assert( (j>=0) && (j<dim.dim[1]) ) ;
	    assert( (k>=0) && (k<dim.dim[2]) ) ;
	    return t[dim.dim[1]*dim.dim[0]*k + dim.dim[0]*j + i] ;
	};

	/// Read-only of a particular element (index {\tt (k,j,i)}) (3D case)
	int operator()(int k, int j, int i) const {
	    assert(etat != ETATNONDEF) ;
	    assert( dim.ndim == 3 ) ;
	    assert( (i>=0) && (i<dim.dim[0]) ) ;
	    assert( (j>=0) && (j<dim.dim[1]) ) ;
	    assert( (k>=0) && (k<dim.dim[2]) ) ;
	    if (etat == ETATZERO) {
		int zero = 0 ;
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
	
	/// Gives the {\tt i}th dimension (ie {tt dim.dim[i]})
	int get_dim(int i) const {	
	    assert( (i>=0) && (i<dim.ndim) ) ;
	    return dim.dim[i] ;
	};
	
    // Outputs
    // -------
    public:
	void sauve(FILE* ) const ;	/// Save in a file

	/// Display   
	friend ostream& operator<<(ostream& , const Itbl& ) ;	

    // Member arithmetics
    // ------------------
    public:
	
	void operator+=(const Itbl &) ;	/// Addition of a {\tt Itbl} to {\tt this}
	void operator+=(int) ;	/// Addition of a {\tt int} to {\tt this}
	void operator-=(const Itbl &) ;	/// Subtraction of a {\tt Itbl} to {\tt this}
	void operator-=(int) ;	/// Subtraction of a {\tt int} to {\tt this}
	void operator*=(const Itbl &) ;	/// Multiplication of {\tt this} by a {\tt Itbl}
	void operator*=(int) ;	/// Multiplication of {\tt this} by a {\tt int}
} ;
ostream& operator<<(ostream& , const Itbl& ) ;	


/**
 * @name Itbl Mathematics
 */
    //@{
Itbl operator+(const Itbl&) ;			/// + Ibl
Itbl operator-(const Itbl&) ;			/// - Itbl
Itbl operator+(const Itbl&, const Itbl&) ;	/// Itbl + Itbl
Itbl operator+(const Itbl&, int) ;		/// Itbl + int
Itbl operator+(int, const Itbl&) ;		/// int + Itbl
Itbl operator+(const Itbl&, int) ;		/// Itbl + int
Itbl operator-(const Itbl&, const Itbl&) ;	/// Itbl - Itbl
Itbl operator-(const Itbl&, int) ;		/// Itbl - int
Itbl operator-(int, const Itbl&) ;		/// int - Itbl
Itbl operator*(const Itbl&, const Itbl&) ;	/// Itbl * Itbl
Itbl operator*(const Itbl&, int) ;		/// Itbl * int
Itbl operator*(int, const Itbl&) ;		/// int * Itbl

Itbl abs(const Itbl& ) ;	    /// Absolute value
int max(const Itbl& ) ;   /// Maximum value of the {\tt Itbl} elements
int min(const Itbl& ) ;   /// Minimum value of the {\tt Itbl} elements

/// Sum of the absolute values of all the {\tt Itbl} elements
int norme(const Itbl& ) ;   

/**
 * Relative difference between two {\tt Itbl} (norme version).
 * Returns {\tt norme(a-b)/norme(b)} unless {\tt b=0}, in which
 * case it returns {\tt norme(a-b)}.
 */
double diffrel(const Itbl& a, const Itbl& b) ; 

/**
 * Relative difference between two {\tt Itbl} (max version).
 * Returns {\tt max(abs(a-b))/max(abs(b))} unless {\tt b=0}, in which
 * case it returns {\tt max(abs(a-b))}.
 */
double diffrelmax(const Itbl& a, const Itbl& b) ; 

    //@}
#endif
