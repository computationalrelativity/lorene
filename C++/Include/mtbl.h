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
 * Revision 1.5  2004/03/22 13:12:42  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.4  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.3  2002/09/13 09:17:33  j_novak
 * Modif. commentaires
 *
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
 * Multi-domain array. \ingroup (spec)
 * 
 * This class is essentially an array of \c Tbl . It is intended to be 
 * used in conjunction with the class \c Mtbl_cf .
 * A \c Mtbl  is initialy created with a \a logical  state \c NONDEF .
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
/// Pointer on the multi-grid \c Mgd3  on which \c this  is defined
	const Mg3d* mg ;  
/// Number of domains (zones)
	int nzone ;	
	/// Logical state (\c ETATNONDEF , \c ETATQCQ  or \c ETATZERO ).
	int etat ;	

    public:
/// Array (size \c nzone ) of pointers on the \c Tbl 's
	Tbl** t;	

    // Constructors - Destructor
    // -------------------------
	
    public:
/// Constructor
	explicit Mtbl(const Mg3d& mgrid) ;	
/// Constructor
	explicit Mtbl(const Mg3d* p_mgrid) ;	
	/// Constructor from a file (see \c sauve(FILE*))
	Mtbl(const Mg3d&, FILE* ) ;		    
/// Constructor from a Coord
	Mtbl(const Coord& c) ;		
/// Copy constructor
	Mtbl(const Mtbl& a) ;		    
/// Destructor
	~Mtbl() ;			    

    // Assignement
    // -----------
    public:
/// Assignement to another \c Mtbl 
	void operator=(const Mtbl& ) ;	    
/// Assignement to a \c double 
	void operator=(double ) ;	    
/// Assignement to a \c int 
	void operator=(int ) ;		    

    // Memory management
    // -----------------
    private:
	/** Logical destructor: dellocates the memory occupied by the \c Tbl 
	 * array \c t . 
	 */
	void del_t() ;	
			
    public:

    /**
     * Sets the logical state to \c ETATNONDEF  (undefined). 
     * Deallocates the memory occupied by the \c Tbl  array \c t .
     */
	void set_etat_nondef() ;
	
    /**
     * Sets the logical state to \c ETATZERO  (zero). 
     * Deallocates the memory occupied by the \c Tbl  array \c t .
     */
	void set_etat_zero() ;	    	

    /**
     * Sets the logical state to \c ETATQCQ  (ordinary state).
     * If the state (member \c etat ) is already \c ETATQCQ , this 
     * function does nothing. Otherwise, it performs the memory allocation
     * for the \c Tbl  array \c t .  
     */
	void set_etat_qcq() ;	    	

    /**
     * Sets the \c Mtbl  to zero in a hard way. 
     * 1/ Sets the logical state to \c ETATQCQ , i.e. to an ordinary state.
     * 2/ Allocates the memory of the \c Tbl  array \c t , and fills it
     * with zeros. NB: this function must be used for debugging purposes only.
     * For other operations, the functions \c set_etat_zero() 
     * or \c annule(int, int)  must be perferred. 
     */
	void annule_hard() ;	
	
    /**
     * Sets the \c Mtbl  to zero in some domains.
     *	@param l_min [input] The \c Mtbl  will be set (logically) to zero
     *			     in the domains whose indices are in the range
     *			     \c [l_min,l_max] .
     *	@param l_max [input] see the comments for \c l_min .
     * 
     * Note that \c annule(0,nzone-1) is equivalent to
     *	 \c set_etat_zero() .
     */
	void annule(int l_min, int l_max) ; 


    // Access to individual elements
    // -----------------------------
    public:
	/** 
	 * Read/write of the \c Tbl  in a given domain.
	 * @param l [input] domain index
	 */ 
	Tbl& set(int l) {
	    assert(l < nzone) ;
	    assert(etat == ETATQCQ) ;
	    return *(t[l]) ;
	};
	
	
	/** 
	 * Read-only of the \c Tbl  in a given domain.
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
	 * @param k [input] \f$\phi\f$ index
	 * @param j [input] \f$\theta\f$ index
	 * @param i [input] \e r  (\f$\xi\f$) index
	 */ 
	double& set(int l, int k, int j, int i) {
	    assert(l < nzone) ;
	    assert(etat == ETATQCQ) ;
	    return (t[l])->set(k, j, i) ;
	};
	
	
	/** 
	 * Read-only of a particular element.
	 * @param l [input] domain index
	 * @param k [input] \f$\phi\f$ index
	 * @param j [input] \f$\theta\f$ index
	 * @param i [input] \e r  (\f$\xi\f$) index
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
	/// Gives the \c Mg3d  on which the \c Mtbl  is defined
	const Mg3d* get_mg() const { return mg ; }; 

/// Gives the logical state
	int get_etat() const { return etat ; };   
	
	/// Gives the number of zones (domains)
	int get_nzone() const { return nzone ; } ; 
	
    // Outputs
    // -------
    public:
/// Save in a file
	void sauve(FILE *) const ;	    
	
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
/// += Mtbl
	void operator+=(const Mtbl & ) ;	
/// += double
	void operator+=(double ) ;		
/// -= Mtbl
	void operator-=(const Mtbl & ) ;	
/// -= double
	void operator-=(double ) ;		
/// *= Mtbl
	void operator*=(const Mtbl & ) ;	
/// *= double
	void operator*=(double ) ;		
/// /= Mtbl
	void operator/=(const Mtbl & ) ;	
/// /= double
	void operator/=(double ) ;		

} ;
ostream& operator<<(ostream& , const Mtbl & ) ;   

/**
 * \defgroup mtbl_mat Mtbl Mathematics
 * \ingroup (spec)
 *
 * @{
 */
/// + Mtbl
Mtbl operator+(const Mtbl& ) ;		    
/// \c - Mtbl
Mtbl operator-(const Mtbl& ) ;		    
/// Mtbl + Mtbl
Mtbl operator+(const Mtbl&, const Mtbl& ) ; 
/// Mtbl + double
Mtbl operator+(const Mtbl&, double ) ;	    
/// double + Mtbl
Mtbl operator+(double , const Mtbl& ) ;	    
/// Mtbl + int
Mtbl operator+(const Mtbl&, int ) ;	    
/// int + Mtbl
Mtbl operator+(int, const Mtbl& ) ;	    
/// Mtbl - Mtbl
Mtbl operator-(const Mtbl&, const Mtbl& ) ; 
/// Mtbl - double
Mtbl operator-(const Mtbl&, double ) ;	    
/// double - Mtbl
Mtbl operator-(double, const Mtbl& ) ;	    
/// Mtbl - int
Mtbl operator-(const Mtbl&, int ) ;	    
/// int - Mtbl
Mtbl operator-(int, const Mtbl& ) ;	    
/// Mtbl * Mtbl
Mtbl operator*(const Mtbl&, const Mtbl& ) ; 
/// Mtbl * double
Mtbl operator*(const Mtbl&, double ) ;	    
/// double * Mtbl
Mtbl operator*(double, const Mtbl& ) ;	    
/// Mtbl * int
Mtbl operator*(const Mtbl&, int ) ;	    
/// int * Mtbl
Mtbl operator*(int, const Mtbl& ) ;	    
/// Mtbl / Mtbl
Mtbl operator/(const Mtbl&, const Mtbl& ) ; 
/// Mtbl / double
Mtbl operator/(const Mtbl&, double ) ;	    
/// double / Mtbl
Mtbl operator/(double, const Mtbl& ) ;	    
/// Mtbl / int
Mtbl operator/(const Mtbl&, int ) ;	    
/// int / Mtbl
Mtbl operator/(int, const Mtbl& ) ;	    

/// Sine
Mtbl sin(const Mtbl& ) ;	    
/// Cosine
Mtbl cos(const Mtbl& ) ;	    
/// Tangent
Mtbl tan(const Mtbl& ) ;	    
/// Arcsine
Mtbl asin(const Mtbl& ) ;	    
/// Arccosine
Mtbl acos(const Mtbl& ) ;	    
/// Arctangent
Mtbl atan(const Mtbl& ) ;	    
/// Exponential
Mtbl exp(const Mtbl& ) ;	    
/// Neperian logarithm
Mtbl log(const Mtbl& ) ;	    
/// Basis 10 logarithm
Mtbl log10(const Mtbl& ) ;    
/// Square root
Mtbl sqrt(const Mtbl& ) ;	    
/// Cube root
Mtbl racine_cubique (const Mtbl&) ; 
/// Power \f${\tt Mtbl}^{\tt int}\f$
Mtbl pow(const Mtbl& , int ) ;  
/// Power \f${\tt Mtbl}^{\tt double}\f$
Mtbl pow(const Mtbl& , double ) ; 
/// Absolute value
Mtbl abs(const Mtbl& ) ;	    

/**
 * Maximum values of the \c Mtbl  elements in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Mtbl& ) ;   

/**
 * Minimum values of the \c Mtbl  elements in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Mtbl& ) ;   

/**
 * Sums of the absolute values of all the \c Mtbl  elements in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Mtbl& ) ;   

/**
 * Relative difference between two \c Mtbl  (norme version).
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are \c norme[a(l)-b(l)]/norme[b(l)]  if \c b(l)!=0  and
 *	   \c norme[a(l)-b(l)]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c a  and \c b  
 *	   in domain no. \c l . 
 */
Tbl diffrel(const Mtbl& a, const Mtbl& b) ; 

/**
 * Relative difference between two \c Mtbl  (max version).
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are \c max[abs(a(l)-b(l))]/max[abs(b(l))]  if \c b(l)!=0  and
 *	   \c max[abs(a(l)-b(l))]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c a  and \c b  
 *	   in domain no. \c l . 
 */
Tbl diffrelmax(const Mtbl& a, const Mtbl& b) ; 

/**@} */


#endif
