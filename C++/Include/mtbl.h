/*
 *  Definition of Lorene class Mtbl
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


#ifndef __MTBL_H_
#define __MTBL_H_

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/06/17 14:05:17  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.9  2000/08/16  10:29:45  eric
 * Suppression du membre dzpuis.
 *
 * Revision 2.8  2000/08/04  11:40:58  eric
 * Ajout de l'operateur (int l) et de la fonction set(int l) pour l'acces
 * individuel aux Tbl.
 *
 * Revision 2.7  1999/12/02  17:55:07  phil
 * *** empty log message ***
 *
 * Revision 2.6  1999/10/29  15:05:43  eric
 * Suppression des fonctions membres min() et max():
 * elles deviennent des fonctions externes.
 * Ajout de fonctions mathematiques (abs, norme, etc...).
 *
 * Revision 2.5  1999/10/18  15:06:25  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 * Introduction de la fonction membre annule(int, int).
 *
 * Revision 2.4  1999/10/01  10:35:58  eric
 * Amelioration des commentaires.
 *
 * Revision 2.3  1999/10/01  10:08:25  eric
 * Depoussierage
 * Documentation.
 *
 * Revision 2.2  1999/03/02  18:54:50  eric
 * Ajout de la fonction affiche_seuil.
 *
 * Revision 2.1  1999/02/22  15:23:49  hyc
 * *** empty log message ***
 *
 *
 * Revision 2.0  1999/01/15  09:10:39  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */


// Headers Lorene 
#include "tbl.h"
#include "grilles.h"

class Coord ;

/**
 * Multi-domain array.
 * 
 * This class is essentially an array of {\tt Tbl}. It is intended to be 
 * used in conjunction with the class {\tt Mtbl\_cf}.
 * A {\tt Mtbl} is initialy created with a {\sl logical} state {\tt NONDEF}.
 * Arithmetic operations are provided with the usual meaning (see 
 * below).
 * 
 * @version #$Id$#
 *
 */
class Mtbl {

    // Data : 
    // -----
    private:
	const Mg3d* mg ;  /// Pointer on the multi-grid {\tt Mgd3} on which {\tt this} is defined
	int nzone ;	/// Number of domains (zones)
	/// Logical state ({\tt ETATNONDEF}, {\tt ETATQCQ} or {\tt ETATZERO}).
	int etat ;	

    public:
	Tbl** t;	/// Array (size {\tt nzone}) of pointers on the {\tt Tbl}'s

    // Constructors - Destructor
    // -------------------------
	
    public:
	explicit Mtbl(const Mg3d& ) ;	/// Constructor
	explicit Mtbl(const Mg3d* ) ;	/// Constructor
	/// Constructor from a file (see {\tt sauve(FILE* )})
	Mtbl(const Mg3d&, FILE* ) ;		    
	Mtbl(const Coord& ) ;		/// Constructor from a Coord
	Mtbl(const Mtbl& ) ;		    /// Copy constructor
	~Mtbl() ;			    /// Destructor

    // Assignement
    // -----------
    public:
	void operator=(const Mtbl& ) ;	    /// Assignement to another {\tt Mtbl}
	void operator=(double ) ;	    /// Assignement to a {\tt double}
	void operator=(int ) ;		    /// Assignement to a {\tt int}

    // Memory management
    // -----------------
    private:
	/** Logical destructor: dellocates the memory occupied by the {\tt Tbl}
	 * array {\tt t}. 
	 */
	void del_t() ;	
			
    public:

    /**
     * Sets the logical state to {\tt ETATNONDEF} (undefined). 
     * Deallocates the memory occupied by the {\tt Tbl} array {\tt t}.
     */
	void set_etat_nondef() ;
	
    /**
     * Sets the logical state to {\tt ETATZERO} (zero). 
     * Deallocates the memory occupied by the {\tt Tbl} array {\tt t}.
     */
	void set_etat_zero() ;	    	

    /**
     * Sets the logical state to {\tt ETATQCQ} (ordinary state).
     * If the state (member {\tt etat}) is already {\tt ETATQCQ}, this 
     * function does nothing. Otherwise, it performs the memory allocation
     * for the {\tt Tbl} array {\tt t}.  
     */
	void set_etat_qcq() ;	    	

    /**
     * Sets the {\tt Mtbl} to zero in a hard way. 
     * 1/ Sets the logical state to {\tt ETATQCQ}, i.e. to an ordinary state.
     * 2/ Allocates the memory of the {\tt Tbl} array {\tt t}, and fills it
     * with zeros. NB: this function must be used for debugging purposes only.
     * For other operations, the functions {\tt set\_etat\_zero()}
     * or {\tt annule(int, int)} must be perferred. 
     */
	void annule_hard() ;	
	
    /**
     * Sets the {\tt Mtbl} to zero in some domains.
     *	@param l_min [input] The {\tt Mtbl} will be set (logically) to zero
     *			     in the domains whose indices are in the range
     *			     {\tt [l\_min, l\_max]}.
     *	@param l_max [input] see the comments for {\tt l\_min}.
     * 
     * Note that {\tt annule(0, nzone-1)} is equivalent to
     *	 {\tt set\_etat\_zero()}.
     */
	void annule(int l_min, int l_max) ; 


    // Access to individual elements
    // -----------------------------
    public:
	/** 
	 * Read/write of the {\tt Tbl} in a given domain.
	 * @param l [input] domain index
	 */ 
	Tbl& set(int l) {
	    assert(l < nzone) ;
	    assert(etat == ETATQCQ) ;
	    return *(t[l]) ;
	};
	
	
	/** 
	 * Read-only of the {\tt Tbl} in a given domain.
	 * @param l [input] domain index
	 */ 
	const Tbl& operator()(int l) const {
	    assert(l < nzone) ;
	    assert(etat == ETATQCQ) ;
	    return *(t[l]) ;
	};

	/** 
	 * Read/write of a particular element.
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] $r$ ($\xi$) index
	 */ 
	double& set(int l, int k, int j, int i) {
	    assert(l < nzone) ;
	    assert(etat == ETATQCQ) ;
	    return (t[l])->set(k, j, i) ;
	};
	
	
	/** 
	 * Read-only of a particular element.
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] $r$ ($\xi$) index
	 */ 
	double operator()(int l, int k, int j, int i) const {
	    assert(etat != ETATNONDEF) ;
	    assert(l < nzone) ;
	    if (etat == ETATZERO) {
		double zero = 0. ;
		return zero ; 
	    }
	    else return (*t[l])(k, j, i) ;
	};


    // Extraction of information
    // -------------------------
    public:
	/// Gives the {\tt Mg3d} on which the {\tt Mtbl} is defined
	const Mg3d* get_mg() const { return mg ; }; 

	int get_etat() const { return etat ; };   /// Gives the logical state
	
	/// Gives the number of zones (domains)
	int get_nzone() const { return nzone ; } ; 
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file
	
	/** Prints only the values greater than a given threshold.
	 *   @param ostr [input] Output stream used for the printing 
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param threshold [input] Value above which an array element is printed
	 *    (default: 1.e-7)
	 */
	void affiche_seuil(ostream& ostr, int precision = 4, 
			   double threshold = 1.e-7) const ;
	/// Display
	friend ostream& operator<<(ostream& , const Mtbl & ) ;   

    // Member arithmetics
    // ------------------
    public:
	void operator+=(const Mtbl & ) ;	/// += Mtbl
	void operator+=(double ) ;		/// += double
	void operator-=(const Mtbl & ) ;	/// -= Mtbl
	void operator-=(double ) ;		/// -= double
	void operator*=(const Mtbl & ) ;	/// *= Mtbl
	void operator*=(double ) ;		/// *= double
	void operator/=(const Mtbl & ) ;	/// /= Mtbl
	void operator/=(double ) ;		/// /= double

} ;
ostream& operator<<(ostream& , const Mtbl & ) ;   

/**
 * @name Mtbl Mathematics
 */
    //@{
Mtbl operator+(const Mtbl& ) ;		    /// + Mtbl
Mtbl operator-(const Mtbl& ) ;		    /// - Mtbl
Mtbl operator+(const Mtbl&, const Mtbl& ) ; /// Mtbl + Mtbl
Mtbl operator+(const Mtbl&, double ) ;	    /// Mtbl + double
Mtbl operator+(double , const Mtbl& ) ;	    /// double + Mtbl
Mtbl operator+(const Mtbl&, int ) ;	    /// Mtbl + int
Mtbl operator+(int, const Mtbl& ) ;	    /// int + Mtbl
Mtbl operator-(const Mtbl&, const Mtbl& ) ; /// Mtbl - Mtbl
Mtbl operator-(const Mtbl&, double ) ;	    /// Mtbl - double
Mtbl operator-(double, const Mtbl& ) ;	    /// double - Mtbl
Mtbl operator-(const Mtbl&, int ) ;	    /// Mtbl - int
Mtbl operator-(int, const Mtbl& ) ;	    /// int - Mtbl
Mtbl operator*(const Mtbl&, const Mtbl& ) ; /// Mtbl * Mtbl
Mtbl operator*(const Mtbl&, double ) ;	    /// Mtbl * double
Mtbl operator*(double, const Mtbl& ) ;	    /// double * Mtbl
Mtbl operator*(const Mtbl&, int ) ;	    /// Mtbl * int
Mtbl operator*(int, const Mtbl& ) ;	    /// int * Mtbl
Mtbl operator/(const Mtbl&, const Mtbl& ) ; /// Mtbl / Mtbl
Mtbl operator/(const Mtbl&, double ) ;	    /// Mtbl / double
Mtbl operator/(double, const Mtbl& ) ;	    /// double / Mtbl
Mtbl operator/(const Mtbl&, int ) ;	    /// Mtbl / int
Mtbl operator/(int, const Mtbl& ) ;	    /// int / Mtbl

Mtbl sin(const Mtbl& ) ;	    /// Sine
Mtbl cos(const Mtbl& ) ;	    /// Cosine
Mtbl tan(const Mtbl& ) ;	    /// Tangent
Mtbl asin(const Mtbl& ) ;	    /// Arcsine
Mtbl acos(const Mtbl& ) ;	    /// Arccosine
Mtbl atan(const Mtbl& ) ;	    /// Arctangent
Mtbl exp(const Mtbl& ) ;	    /// Exponential
Mtbl log(const Mtbl& ) ;	    /// Neperian logarithm
Mtbl log10(const Mtbl& ) ;    /// Basis 10 logarithm
Mtbl sqrt(const Mtbl& ) ;	    /// Square root
Mtbl racine_cubique (const Mtbl&) ; /// Cube root
Mtbl pow(const Mtbl& , int ) ;  /// Power ${\tt Mtbl}^{\tt int}$
Mtbl pow(const Mtbl& , double ) ; /// Power ${\tt Mtbl}^{\tt double}$
Mtbl abs(const Mtbl& ) ;	    /// Absolute value

/**
 * Maximum values of the {\tt Mtbl} elements in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Mtbl& ) ;   

/**
 * Minimum values of the {\tt Mtbl} elements in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Mtbl& ) ;   

/**
 * Sums of the absolute values of all the {\tt Mtbl} elements in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Mtbl& ) ;   

/**
 * Relative difference between two {\tt Mtbl} (norme version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt norme[a(l)-b(l)]/norme[b(l)]} if {\tt b(l)!=0} and
 *	   {\tt norme[a(l)-b(l)]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrel(const Mtbl& a, const Mtbl& b) ; 

/**
 * Relative difference between two {\tt Mtbl} (max version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt max[abs(a(l)-b(l))]/max[abs(b(l))]} if {\tt b(l)!=0} and
 *	   {\tt max[abs(a(l)-b(l))]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrelmax(const Mtbl& a, const Mtbl& b) ; 

    //@}


#endif
