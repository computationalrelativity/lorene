/*
 *  Definition of Lorene class Mtbl_cf
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 1999-2001 Jerome Novak
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


#ifndef __MTBL_CF_H_
#define __MTBL_CF_H_

/*
 * $Id$
 * $Log$
 * Revision 1.6  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.5  2003/10/19 19:44:41  e_gourgoulhon
 * Introduced new method display (to replace the old affiche_seuil).
 *
 * Revision 1.4  2003/10/15 21:09:22  e_gourgoulhon
 * Added method poisson_regu.
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
 * Revision 2.28  2000/09/11  13:52:21  eric
 * Ajout des methodes mult_cp() et mult_sp().
 *
 * Revision 2.27  2000/08/16  10:42:54  eric
 * Suppression du membre dzpuis.
 *
 * Revision 2.26  2000/08/04  11:41:52  eric
 * Ajout de l'operateur (int l) et de la fonction set(int l) pour l'acces
 * individuel aux Tbl.
 *
 * Revision 2.25  2000/03/06  10:26:24  eric
 * Ajout des fonctions val_point_symy, val_point_asymy.
 *
 * Revision 2.24  2000/02/25  13:53:44  eric
 * Suppression de la fonction nettoie().
 *
 * Revision 2.23  1999/12/29  13:11:34  eric
 * Ajout de la fonction val_point_jk.
 *
 * Revision 2.22  1999/12/07  14:51:52  eric
 * Changement ordre des arguments (phi,theta,xi) --> (xi,theta,phi)
 *  dans la routine val_point.
 *
 * Revision 2.21  1999/12/06  16:45:48  eric
 * Ajout de la fonction val_point.
 *
 * Revision 2.20  1999/11/23  14:30:12  novak
 * Ajout des membres mult_ct et mult_st
 *
 * Revision 2.19  1999/11/16  13:06:22  novak
 * Ajout de mult_x et scost
 *
 * Revision 2.18  1999/10/29  15:07:10  eric
 * Suppression des fonctions membres min() et max():
 * elles deviennent des fonctions externes.
 * Ajout de fonctions mathematiques (abs, norme, etc...).
 *
 * Revision 2.17  1999/10/21  13:42:00  eric
 * *** empty log message ***
 *
 * Revision 2.16  1999/10/21  12:49:13  eric
 * Ajout de la fonction membre nettoie().
 *
 * Revision 2.15  1999/10/18  15:06:50  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 * Introduction de la fonction membre annule(int, int).
 *
 * Revision 2.14  1999/10/18  13:39:13  eric
 * Suppression de l'argument base dans les routines de derivation.
 *
 * Revision 2.13  1999/10/13  15:49:22  eric
 * Ajout du membre base.
 * Modification des constructeurs (la base doit etre passee en argument).
 *
 * Revision 2.12  1999/10/01  14:49:19  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.11  1999/09/07  14:33:43  phil
 * ajout de la fonction ssint(int*)
 *
 * Revision 2.10  1999/04/26  17:02:40  phil
 * ajout de sx2(int*)
 *
 * Revision 2.9  1999/04/26  16:24:02  phil
 * ajout de mult2_xm1_zec(int*)
 *
 * Revision 2.8  1999/04/26  16:12:27  phil
 * ajout de mult_xm1_zec(int *)
 *
 * Revision 2.7  1999/04/26  15:48:43  phil
 * ajout de sxm1_zec(int*)
 *
 * Revision 2.6  1999/04/26  14:50:30  phil
 * ajout de sx(int*)
 *
 * Revision 2.5  1999/04/26  12:19:38  phil
 * ajout lapang
 *
 * Revision 2.4  1999/03/03  10:32:44  hyc
 * *** empty log message ***
 *
 * Revision 2.3  1999/03/02  18:55:00  eric
 * Ajout de la fonction affiche_seuil.
 *
 *
 * $Header$
 *
 */

// Headers Lorene 
#include "tbl.h"
#include "base_val.h"
#include "grilles.h"

class Mg3d ;

/**
 * Coefficients storage for the multi-domain spectral method.
 * 
 * This class is essentially an array (on the various physical domains)
 * of {\tt Tbl} specially designed for 
 * storage of the coefficients of the spectral expansions in each domain.  
 * It is intended to be 
 * used in conjunction with the class {\tt Mtbl} (see class {\tt Valeur}).
 * A difference between a {\tt Mtbl} and a {\tt Mtbl\_cf}, both defined one
 * the same grid {\tt Mg3d}, is that each {\tt Tbl} of the {\tt Mtbl\_cf}
 * has 2 more elements in the $\phi$-dimension (Dim\_tbl::dim[2]) than the
 * corresponding {\tt Tbl} of the {\tt Mtbl}. 
 * A {\tt Mbl\_cf} is initialy created with a {\it logical} state {\tt ETATZERO}.
 * Arithmetic operations are provided with the usual meaning (see 
 * below). 
 * 
 * @version #$Id$#
 *
 */
class Mtbl_cf {

    // Data : 
    // -----
    private:
	const Mg3d* mg ;  /// Pointer on the multi-grid {\tt Mgd3} on which {\tt this} is defined
	int nzone ;	/// Number of domains (zones)
	/// Logical state ({\tt ETATNONDEF}, {\tt ETATQCQ} or {\tt ETATZERO}).
	int etat ;	

    public:
	/// Bases of the spectral expansions
	Base_val base ;
	
	/** Array (size {\tt nzone}) of pointers on the {\tt Tbl}'s which 
	 * contain the spectral coefficients in each domain
	 */
	Tbl** t ;	

    // Constructors - Destructor
    // -------------------------
	
    public:
	Mtbl_cf(const Mg3d& mgrid, const Base_val& basis) ; /// Constructor 
	Mtbl_cf(const Mg3d* p_mgrid, const Base_val& basis) ; /// Constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	Mtbl_cf(const Mg3d&, FILE* ) ;		    

	Mtbl_cf(const Mtbl_cf& ) ;	    /// Copy constructor
	~Mtbl_cf() ;			    /// Destructor

    // Assignement
    // -----------
	void operator=(const Mtbl_cf& ) ;   /// Assignement to another {\tt Mtbl\_cf}
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
     * Sets the {\tt Mtbl\_cf} to zero in a hard way. 
     * 1/ Sets the logical state to {\tt ETATQCQ}, i.e. to an ordinary state.
     * 2/ Allocates the memory of the {\tt Tbl} array {\tt t}, and fills it
     * with zeros. NB: this function must be used for debugging purposes only.
     * For other operations, the functions {\tt set\_etat\_zero()}
     * or {\tt annule(int, int)} must be perferred. 
     */
	void annule_hard() ;	
	
    /**
     * Sets the {\tt Mtbl\_cf} to zero in some domains.
     *	@param l_min [input] The {\tt Mtbl\_cf} will be set (logically) to zero
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
	 * Read/write of the {\tt Tbl} containing the coefficients
	 * in a given domain.
	 * @param l [input] domain index
	 */ 
	Tbl& set(int l) {
	    assert(l < nzone) ;
	    assert(etat == ETATQCQ) ;
	    return *(t[l]) ;
	};
	
	
	/** 
	 * Read-only of the {\tt Tbl} containing the coefficients
	 * in a given domain.
	 * @param l [input] domain index
	 */ 
	const Tbl& operator()(int l) const {
	    assert(l < nzone) ;
	    assert(etat == ETATQCQ) ;
	    return *(t[l]) ;
	};


	/** Read/write of a particular element.
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] {\it r} ($\xi$) index
	 */ 
	double& set(int l, int k, int j, int i) {
	    assert(l < nzone) ;
	    assert(etat == ETATQCQ) ;
	    return (t[l])->set(k, j, i) ;
	};
	
	
	/** Read-only of a particular element.
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] {\it r} ($\xi$) index
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

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point, by means of the spectral expansion.
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate $\xi$
	*	 @param theta [input] value of the coordinate $\theta'$
	*	 @param phi [input] value of the coordinate $\phi'$
	*	 @return value at the point $(\xi, \theta', \phi')$ in
	*	    the domain no. {\it l} of the field whose spectral coefficients
	*	    are stored in {\tt *this}. 
	*/
	double val_point(int l, double x, double theta, double phi) const ; 

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point, by means of the spectral expansion.
	*   Case where the field is symmetric with respect to the y=0 plane. 
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate $\xi$
	*	 @param theta [input] value of the coordinate $\theta'$
	*	 @param phi [input] value of the coordinate $\phi'$
	*	 @return value at the point $(\xi, \theta', \phi')$ in
	*	    the domain no. {\it l} of the field whose spectral coefficients
	*	    are stored in {\tt *this}. 
	*/
	double val_point_symy(int l, double x, double theta, double phi) const ; 

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point, by means of the spectral expansion.
	*   Case where the field is antisymmetric with respect to the y=0 plane. 
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate $\xi$
	*	 @param theta [input] value of the coordinate $\theta'$
	*	 @param phi [input] value of the coordinate $\phi'$
	*	 @return value at the point $(\xi, \theta', \phi')$ in
	*	    the domain no. {\it l} of the field whose spectral coefficients
	*	    are stored in {\tt *this}. 
	*/
	double val_point_asymy(int l, double x, double theta, double phi) const ; 

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point in $\xi$, but collocation point in 
	*   $(\theta', \phi')$, by means of the spectral expansion.
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate $\xi$
	*	 @param j [input] index of the collocation point in $\theta'$
	*	 @param k [input] index of the collocation point in $\phi'$
	*	 @return value at the point 
	*		    $(\xi, {\theta'}_j, {\phi'}_k)$ in
	*	    the domain no. {\it l} of the field whose spectral coefficients
	*	    are stored in {\tt *this}. 
	*/
	double val_point_jk(int l, double x, int j, int k) const ; 

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point in $\xi$, but collocation point in 
	*   $(\theta', \phi')$, by means of the spectral expansion.
	*   Case where the field is symmetric with respect to the y=0 plane. 
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate $\xi$
	*	 @param j [input] index of the collocation point in $\theta'$
	*	 @param k [input] index of the collocation point in $\phi'$
	*	 @return value at the point 
	*		    $(\xi, {\theta'}_j, {\phi'}_k)$ in
	*	    the domain no. {\it l} of the field whose spectral coefficients
	*	    are stored in {\tt *this}. 
	*/
	double val_point_jk_symy(int l, double x, int j, int k) const ; 

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point in $\xi$, but collocation point in 
	*   $(\theta', \phi')$, by means of the spectral expansion.
	*   Case where the field is antisymmetric with respect to the y=0 plane. 
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate $\xi$
	*	 @param j [input] index of the collocation point in $\theta'$
	*	 @param k [input] index of the collocation point in $\phi'$
	*	 @return value at the point 
	*		    $(\xi, {\theta'}_j, {\phi'}_k)$ in
	*	    the domain no. {\it l} of the field whose spectral coefficients
	*	    are stored in {\tt *this}. 
	*/
	double val_point_jk_asymy(int l, double x, int j, int k) const ; 


    // Extraction of information
    // -------------------------
    public:
	/// Returns the {\tt Mg3d} on which the {\tt Mtbl\_cf} is defined
	const Mg3d* get_mg() const { return mg ; }; 

	int get_etat() const { return etat ; };   /// Returns the logical state
	
	/// Returns the number of zones (domains)
	int get_nzone() const { return nzone ; } ; 
	
		
    // Outputs
    // -------
    public:

	void sauve(FILE *) const ;	    /// Save in a file

	/** Prints the coefficients whose values are greater than a given threshold,
	 *  as well as the corresponding basis
	 *   @param threshold [input] Value above which a coefficient is printed
	 *    (default: 1.e-7)
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param ostr [input] Output stream used for the printing (default: cout)
	 */
	void display(double threshold = 1.e-7, int precision = 4, 
			   ostream& ostr = cout) const ;

	/** Prints only the values greater than a given threshold.
	 *   @param ostr [input] Output stream used for the printing 
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param threshold [input] Value above which an array element is printed
	 *    (default: 1.e-7)
	 */
	void affiche_seuil(ostream& ostr, int precision = 4, 
			   double threshold = 1.e-7) const ;
	/// Display
	friend ostream& operator<<(ostream& , const Mtbl_cf& ) ;   

    // Member arithmetics
    // ------------------
    public:
	void operator+=(const Mtbl_cf & ) ;	/// += Mtbl\_cf
	void operator-=(const Mtbl_cf & ) ;	/// -= Mtbl\_cf
	void operator*=(double ) ;		/// *= double
	void operator/=(double ) ;		/// /= double

    // Linear operators
    // ----------------
    public:
	/// ${\partial \over \partial \xi}$ 
	void dsdx() ;		    

	/// ${\partial^2\over \partial \xi^2}$ 
	void d2sdx2() ;		    

	/** ${1 \over \xi}$ ({\it r}-sampling = {\tt RARE}) \\
	 * Id ({\it r} sampling = {\tt FIN}) \\
	 * ${1 \over \xi-1}$ ({\it r}-sampling = {\tt UNSURR})
	 */
	void sx() ;		    

	/** ${1 \over \xi^2}$ ({\it r}-sampling = {\tt RARE}) \\
	 * Id ({\it r} sampling = {\tt FIN}) \\
	 * ${1 \over (\xi-1)^2}$ ({\it r}-sampling = {\tt UNSURR})
	 */
	void sx2() ;		    
	
	/** $\xi \, Id$ ({\it r}-sampling = {\tt RARE}) \\
	 * Id ({\it r} sampling = {\tt FIN}) \\
	 * $(\xi-1) \, Id $ ({\it r}-sampling = {\tt UNSURR})
	 */
	void mult_x() ;		    
	
	/** Id ({\it r} sampling = {\tt RARE, FIN}) \\
	 * ${1 \over (\xi-1)}$ ({\it r}-sampling = {\tt UNSURR})
	 */
	void sxm1_zec() ;		    

	/** Id ({\it r} sampling = {\tt RARE, FIN}) \\
	 * $(\xi-1) \, Id$ ({\it r}-sampling = {\tt UNSURR})
	 */
	void mult_xm1_zec() ;	    

	/** Id ({\it r} sampling = {\tt RARE, FIN}) \\
	 * $(\xi-1)^2 \, Id$ ({\it r}-sampling = {\tt UNSURR})
	 */
	void mult2_xm1_zec() ;	    

	/// ${\partial \over \partial \theta}$ 
	void dsdt() ;		   

	/// ${\partial^2 \over \partial \theta^2}$
	void d2sdt2() ;		    

	/// $Id\over\sin\theta$
	void ssint() ;		    
		
	/// $Id\over\cos\theta$
	void scost() ;		    

	/// $\cos\theta \, Id$
	void mult_ct() ;		    

	/// $\sin\theta \, Id$
	void mult_st() ;		    

	/// ${\partial \over \partial \phi}$ 
	void dsdp() ;		    

	/// ${\partial^2 \over \partial \phi^2}$ 
	void d2sdp2() ;		    

	/// $\cos\phi \, Id$
	void mult_cp() ;		    

	/// $\sin\phi \, Id$
	void mult_sp() ;		    

	/// Angular Laplacian
	void lapang() ;
	
	// PDE resolution
	//---------------
	public: 
	/** Resolution of an angular Poisson equation.
	 * The angular Poisson equation is $\Delta_{\theta\varphi} u = \sigma$,
	 * where $\Delta_{\theta\varphi} u := \frac{\partial^2 u}
	 *  {\partial \theta^2} + \frac{1}{\tan \theta} \frac{\partial u}
	 *  {\partial \theta} +\frac{1}{\sin^2 \theta}\frac{\partial^2 u}
	 *  {\partial \varphi^2}$.
	 * 
	 * Before the call to {\tt poisson\_angu()}, {\tt *this} contains the
	 * coefficients of the source $\sigma$; after the call, it contains the
	 * coefficients of the solution $u$.
	 */
	void poisson_angu() ; 
	
} ;
ostream& operator<<(ostream& , const Mtbl_cf& ) ;   

/**
 * @name Mtbl\_cf Mathematics
 */
    //@{
Mtbl_cf operator+(const Mtbl_cf& ) ;			/// + Mtbl\_cf
Mtbl_cf operator-(const Mtbl_cf& ) ;			/// - Mtbl\_cf
Mtbl_cf operator+(const Mtbl_cf&, const Mtbl_cf& ) ;	/// Mtbl\_cf + Mtbl\_cf
Mtbl_cf operator-(const Mtbl_cf&, const Mtbl_cf& ) ;	/// Mtbl\_cf - Mtbl\_cf
Mtbl_cf operator*(const Mtbl_cf&, double ) ;		/// Mtbl\_cf * double
Mtbl_cf operator*(double, const Mtbl_cf& ) ;		/// double * Mtbl\_cf
Mtbl_cf operator*(const Mtbl_cf&, int ) ;		/// Mtbl\_cf * int
Mtbl_cf operator*(int, const Mtbl_cf& ) ;		/// int * Mtbl\_cf
Mtbl_cf operator/(const Mtbl_cf&, double ) ;		/// Mtbl\_cf / double
Mtbl_cf operator/(const Mtbl_cf&, int ) ;		/// Mtbl\_cf / int

Mtbl_cf abs(const Mtbl_cf& ) ;	    /// Absolute value

/**
 * Maximum values of the {\tt Mtbl\_cf} elements in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Mtbl_cf& ) ;   

/**
 * Minimum values of the {\tt Mtbl\_cf} elements in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Mtbl_cf& ) ;   

/**
 * Sums of the absolute values of all the {\tt Mtbl\_cf} elements in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Mtbl_cf& ) ;   

/**
 * Relative difference between two {\tt Mtbl\_cf} (norme version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt norme[a(l)-b(l)]/norme[b(l)]} if {\tt b(l)!=0} and
 *	   {\tt norme[a(l)-b(l)]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrel(const Mtbl_cf& a, const Mtbl_cf& b) ; 

/**
 * Relative difference between two {\tt Mtbl\_cf} (max version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt max[abs(a(l)-b(l))]/max[abs(b(l))]} if {\tt b(l)!=0} and
 *	   {\tt max[abs(a(l)-b(l))]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrelmax(const Mtbl_cf& a, const Mtbl_cf& b) ; 

    //@}

#endif
