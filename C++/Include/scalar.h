/*
 *  Definition of Lorene class Scalar
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *   
 *   Copyright (c) 1999-2000 Jean-Alain Marck (for previous class Cmp)
 *   Copyright (c) 1999-2002 Eric Gourgoulhon (for previous class Cmp)
 *   Copyright (c) 1999-2001 Philippe Grandclement (for previous class Cmp)
 *   Copyright (c) 2000-2002 Jerome Novak (for previous class Cmp)
 *   Copyright (c) 2000-2001 Keisuke Taniguchi (for previous class Cmp)
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


#ifndef __SCALAR_H_ 
#define __SCALAR_H_ 


/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/09/24 15:10:54  j_novak
 * Suppression of the etat flag in class Tensor (still present in Scalar)
 *
 * Revision 1.3  2003/09/24 12:01:44  j_novak
 * Added friend functions for math.
 *
 * Revision 1.2  2003/09/24 10:22:01  e_gourgoulhon
 * still in progress...
 *
 * Revision 1.1  2003/09/22 12:50:47  e_gourgoulhon
 * First version: not ready yet!
 *
 *
 * $Header$
 *
 */

// Headers Lorene 
#include "valeur.h"

class Param ; 

/**
 * Tensor field of valence 0 (or component of a tensorial field).
 * 
 * @version #$Id$#
 * 
 */

class Scalar : public Tensor {
  
  // Data : 
  // -----
 protected:
  
  /** The logical state {\tt ETATNONDEF} (undefined), {\tt ETATZERO}(null)
   *  or {\tt ETATQCQ}(ordinary).
   */
  int etat ; 
  
	/**
	 * Power of {\it r} by which the quantity represented by {\tt this} 
	 * must be divided in the  compactified external domain (CED) in order 
	 * to get the correct physical values
	 */
	int dzpuis ;	

	Valeur va ;		/// The numerical value of the {\tt Scalar}    

    // Derived data : 
    // ------------
    protected:
	/// Pointer on $\partial/\partial r$ of {\tt *this}
	mutable Scalar* p_dsdr ;	
	/// Pointer on $1/r \partial/\partial \theta$ of {\tt *this}
	mutable Scalar* p_srdsdt ;	
	/// Pointer on $1/(r\sin\theta) \partial/\partial \phi$ of {\tt *this}
	mutable Scalar* p_srstdsdp ;
	
	/** Pointer on $\partial/\partial x$ of {\tt *this},
	 *  where $x=r\sin\theta \cos\phi$
	 */
	mutable Scalar* p_dsdx ;	

	/** Pointer on $\partial/\partial y$ of {\tt *this},
	 *  where $y=r\sin\theta \sin\phi$
	 */
	mutable Scalar* p_dsdy ;	

	/** Pointer on $\partial/\partial z$ of {\tt *this},
	 *  where $z=r\cos\theta$
	 */
	mutable Scalar* p_dsdz ;	

	/** Pointer on the Laplacian of {\tt *this}
	 */
	mutable Scalar* p_lap ;	

	/** Power of {\it r} by which the last computed Laplacian has been 
	 *  multiplied in the compactified external domain.  
	 */
	mutable int ind_lap ; 

	/** Pointer on the space integral of {\tt *this} (values in each 
	 *  domain)
	 */
	mutable Tbl* p_integ ; 

    // Constructors - Destructor
    // -------------------------
	
    public:

	explicit Scalar(const Map& ) ;	/// Constructor from mapping

	explicit Scalar(const Map* ) ;	/// Constructor from mapping
	
	Scalar(const Scalar& ) ;		/// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	Scalar(const Map&, const Mg3d&, FILE* ) ;    		

	~Scalar() ;			/// Destructor


    // Extraction of information
    // -------------------------
    public:
	/** Returns the logical state {\tt ETATNONDEF} (undefined), 
	 * {\tt ETATZERO}(null) or {\tt ETATQCQ}(ordinary).
	 */
	int get_etat() const {return etat;} ; 

	int get_dzpuis() const {return dzpuis;} ; /// Returns {\tt dzpuis}
	
	/// Returns {\tt va} (read only version)
	const Valeur& get_spectral_va() const {return va;} ; 

	/// Returns {\tt va} (read/write version)
	Valeur& set_spectral_va() {return va;} ; 
	
	/** Returns {\tt true} if the last domain is compactified and
	 *  {\tt *this} is not zero in this domain
	 */
	bool dz_nonzero() const ; 
	
	/** Returns {\tt false} if the last domain is compactified 
	 *  and {\tt *this} is not zero in this domain and {\tt dzpuis}
	 *  is not equal to {\tt dzi}, otherwise return true. 
	 */
	bool check_dzpuis(int dzi) const ; 
	
	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point $(r, \theta, \phi)$, by means of the spectral 
	*   expansion.
	*	 @param r [input] value of the coordinate {\it r}
	*	 @param theta [input] value of the coordinate $\theta$
	*	 @param phi [input] value of the coordinate $\phi$
	*	 @return value at the point $(r, \theta, \phi)$ 
	*		 of the field represented by {\tt *this}. 
	*/
	double val_point(double r, double theta, double phi) const ; 


    // Assignment
    // -----------
    public: 
	/// Assignment to another {\tt Scalar} defined on the same mapping
	void operator=(const Scalar&) ;	

	/// Assignment to a {\tt Tensor} (of valence 0)
	virtual void operator=(const Tensor&) ; 

	void operator=(const Valeur& ) ; /// Assignment to a {\tt Valeur}
	void operator=(const Mtbl& ) ;	 /// Assignment to a {\tt Mtbl}
	void operator=(double ) ;	 /// Assignment to a {\tt double}
	void operator=(int ) ;		 /// Assignment to an {\tt int}
	    	


    // Access to individual elements
    // -----------------------------
    public:

	/** Read/write of the value in a given domain.
	 * NB: to gain in efficiency, the method {\tt del\_deriv()} (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *
	 * @param l [input] domain index
	 * @return Tbl containing the value of the field in domain {\tt l}.
	 */ 
	Tbl& set_domain(int l) {
	    assert(etat == ETATQCQ) ;
	    return va.set(l) ;
	};
	
	/** Read-only of the value in a given domain.
	 * @param l [input] domain index
	 * @return Tbl containing the value of the field in domain {\tt l}.
	 */ 
	const Tbl& domain(int l) const {
	    assert(etat == ETATQCQ) ;
	    return va(l) ;
	};


	/** Read/write of a particular element.
	 * NB: to gain in efficiency, the method {\tt del\_deriv()} (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *     
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] {\it r} ($\xi$) index
	 */ 
	double& set_point(int l, int k, int j, int i) {
	    assert(etat == ETATQCQ) ;
	    return va.set(l, k, j, i) ;
	};
	
	
	/** Read-only of a particular element.
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] {\it r} ($\xi$) index
	 */ 
	double point(int l, int k, int j, int i) const {
	    assert(etat != ETATNONDEF) ;
	    if (etat == ETATZERO) {
			double zero = 0. ;
			return zero ; 
	    }
	    else{ 	    
			return va(l, k, j, i) ;
	    }
	};

    // Memory management
    // -----------------
    protected:
	void del_t() ;		    /// Logical destructor
	virtual void del_deriv() ;	    /// Logical destructor of the derivatives
	virtual void set_der_0x0() ;	    /// Sets the pointers for derivatives to 0x0

    public:

    /**
     * Sets the logical state to {\tt ETATNONDEF} (undefined). 
     * Calls the logical destructor of the {\tt Valeur va} and
     * deallocates the memory occupied by all the derivatives. 
     */
	virtual void set_etat_nondef() ;   

    /**
     * Sets the logical state to {\tt ETATZERO} (zero). 
     * Calls the logical destructor of the {\tt Valeur va} and
     * deallocates the memory occupied by all the derivatives. 
     */
	virtual void set_etat_zero() ;	    
	
    /**
     * Sets the logical state to {\tt ETATQCQ} (ordinary state).
     * If the state is already {\tt ETATQCQ}, this function does nothing.
     * Otherwise, it calls the logical destructor of the {\tt Valeur va} and
     * deallocates the memory occupied by all the derivatives.
     */
	virtual void set_etat_qcq() ;	    

    /**
     * Sets the logical state to {\tt ETATQCQ} (ordinary state)
     *  and performs the memory allocation of all the 
     *  elements, down to the {\tt double} arrays of the {\tt Tbl}s. 
     *  This function performs in fact recursive calls to {\tt set\_etat\_qcq()}
     *  on each element of the chain {\tt Scalar} ->
     *  {\tt Valeur} -> {\tt Mtbl} -> {\tt Tbl}. 
     */
	virtual void allocate_all() ; 

    /**
     * Sets the {\tt Scalar} to zero in a hard way. 
     * 1/ Sets the logical state to {\tt ETATQCQ}, i.e. to an ordinary state.
     * 2/ Fills the {\tt Valeur va} with zeros. 
     * NB: this function must be used for debugging purposes only.
     * For other operations, the functions {\tt set\_etat\_zero()}
     * or {\tt annule(int, int)} must be perferred. 
     */
	void annule_hard() ;

    /**
     * Sets the {\tt Scalar} to zero in several domains.
     *	@param l_min [input] The {\tt Scalar} will be set (logically) to zero
     *			     in the domains whose indices are in the range
     *			     {\tt [l\_min, l\_max]}.
     *	@param l_max [input] see the comments for {\tt l\_min}.
     * 
     * Note that {\tt annule(0, va.mg->get\_nzone()-1)} is equivalent to
     *	 {\tt set\_etat\_zero()}.
     */
	virtual void annule(int l_min, int l_max) ; 
    
    /**
     * Gives the spectrum in terms of multipolar modes {\it l}.
     *  @return a {\tt Tbl} of size (nzone, lmax), where lmax is the
     *  maximal multipolar momentum over all domains. The {\it l}-th
     *  element contains the L1 norm of the {\it l}-th multipole 
     *  ({\it i.e.} a sum over all {\it m} of the norms (coefficient space)
     *  of the component of a given $Y_l^m$.
     */
	Tbl multipole_spectrum () ;
	
    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file
    
	/** Prints the spectral coefficients.
	 *  Prints only the spectral coefficients greater than a given threshold.
	 *   @param ostr [input] Output stream used for the printing
	 *   @param threshold [input] Value above which an array element is printed
	 *    (default: 1.e-7)
	 *   @param type [input] Type of display : 0 = prints only the
	 *     coefficients,  1 = prints only the values in configuration 
	 *     space, 2 = prints both
	 *   @param precision [input] Number of printed digits (default: 4)
	 */
	void spectral_display(ostream& ostr, double threshold = 1.e-7,
			int type = 0, int precision = 4) const ;

	/// Display
	friend ostream& operator<<(ostream& , const Scalar & ) ;	


    // Member arithmetics
    // ------------------
    public:
	void operator+=(const Scalar &) ;		    /// += Scalar
	void operator-=(const Scalar &) ;		    /// -= Scalar
	void operator*=(const Scalar &) ;		    /// *= Scalar

    // Manipulation of spectral bases
    // ------------------------------    
    /** Sets the spectral bases of the {\tt Valeur va} to the standard ones 
     *  for a scalar field
     */
    void std_spectral_base_scal() ;	 

    /** Sets the spectral bases of the {\tt Valeur va} 
     */
    void set_spectral_base(const Base_val& ) ;	 

    /** Modifies the {\tt dzpuis} flag.
     *  NB: this method does not change the field values stored in
     *  the compactified external domain (use methods {\tt dec\_dzpuis()},
     *  etc... for this purpose).  
     */
    void set_dzpuis(int ) ; 

	


friend Scalar operator-(const Scalar& ) ;			
friend Scalar operator+(const Scalar&, const Scalar &) ;	
friend Scalar operator+(const Scalar&, double ) ;		
friend Scalar operator-(const Scalar &, const Scalar &) ;
friend Scalar operator-(const Scalar&, double ) ;		
friend Scalar operator*(const Scalar &, const Scalar &) ;
friend Scalar operator%(const Scalar &, const Scalar &) ;
friend Scalar operator*(double, const Scalar &) ;		
friend Scalar operator/(const Scalar &, const Scalar &) ;
friend Scalar operator/(const Scalar&, double ) ;	       
friend Scalar operator/(double, const Scalar &) ;

friend Scalar sin(const Scalar& ) ;
friend Scalar cos(const Scalar& ) ;
friend Scalar tan(const Scalar& ) ;
friend Scalar asin(const Scalar& ) ;
friend Scalar acos(const Scalar& ) ;
friend Scalar atan(const Scalar& ) ;
friend Scalar exp(const Scalar& ) ;	
friend Scalar log(const Scalar& ) ;	
friend Scalar log10(const Scalar& ) ;	
friend Scalar sqrt(const Scalar& ) ;	
friend Scalar racine_cubique (const Scalar& ) ;
friend Scalar pow(const Scalar& , int ) ;	
friend Scalar pow(const Scalar& , double ) ; 
friend Scalar abs(const Scalar& ) ;	


};
ostream& operator<<(ostream& , const Scalar & ) ;	

// Prototypage de l'arithmetique
/**
 * @name Scalar mathematics
 */
    //@{
Scalar operator+(const Scalar& ) ;			/// + Scalar
Scalar operator-(const Scalar& ) ;			/// - Scalar
Scalar operator+(const Scalar&, const Scalar &) ;	/// Scalar + Scalar
Scalar operator+(const Scalar&, double ) ;		/// Scalar + double
Scalar operator+(double, const Scalar& ) ;		/// double + Scalar 
Scalar operator+(const Scalar&, int ) ;		/// Scalar + int
Scalar operator+(int, const Scalar& ) ;		/// int + Scalar 
Scalar operator-(const Scalar &, const Scalar &) ;	/// Scalar - Scalar
Scalar operator-(const Scalar&, double ) ;		/// Scalar - double
Scalar operator-(double, const Scalar& ) ;		/// double - Scalar 
Scalar operator-(const Scalar&, int ) ;		/// Scalar - int
Scalar operator-(int, const Scalar& ) ;		/// int - Scalar 
Scalar operator*(const Scalar &, const Scalar &) ;	/// Scalar * Scalar
Scalar operator%(const Scalar &, const Scalar &) ;	/// Scalar * Scalar with desaliasing
Scalar operator*(const Scalar&, double ) ;		/// Scalar * double
Scalar operator*(double, const Scalar &) ;		/// double * Scalar
Scalar operator*(const Scalar&, int ) ;		/// Scalar * int
Scalar operator*(int, const Scalar& ) ;		/// int * Scalar 
Scalar operator/(const Scalar &, const Scalar &) ;	/// Scalar / Scalar
Scalar operator/(const Scalar&, double ) ;		/// Scalar / double
Scalar operator/(double, const Scalar &) ;		/// double / Scalar
Scalar operator/(const Scalar&, int ) ;		/// Scalar / int
Scalar operator/(int, const Scalar &) ;		/// int / Scalar

Scalar sin(const Scalar& ) ;		/// Sine
Scalar cos(const Scalar& ) ;		/// Cosine
Scalar tan(const Scalar& ) ;		/// Tangent
Scalar asin(const Scalar& ) ;		/// Arcsine
Scalar acos(const Scalar& ) ;		/// Arccosine
Scalar atan(const Scalar& ) ;		/// Arctangent
Scalar exp(const Scalar& ) ;		/// Exponential
Scalar log(const Scalar& ) ;		/// Neperian logarithm
Scalar log10(const Scalar& ) ;	/// Basis 10 logarithm
Scalar sqrt(const Scalar& ) ;		/// Square root
Scalar racine_cubique (const Scalar& ) ;		/// Cube root
Scalar pow(const Scalar& , int ) ;	/// Power ${\tt Scalar}^{\tt int}$
Scalar pow(const Scalar& , double ) ; /// Power ${\tt Scalar}^{\tt double}$
Scalar abs(const Scalar& ) ;		/// Absolute value

/**
 * Maximum values of a {\tt Scalar} in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Scalar& ) ;   

/**
 * Minimum values of a {\tt Scalar} in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Scalar& ) ;   

/**
 * Sums of the absolute values of all the values of the {\tt Scalar} 
 * in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Scalar& ) ;   

/**
 * Relative difference between two {\tt Scalar} (norme version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt norme[a(l)-b(l)]/norme[b(l)]} if {\tt b(l)!=0} and
 *	   {\tt norme[a(l)-b(l)]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrel(const Scalar& a, const Scalar& b) ; 

/**
 * Relative difference between two {\tt Scalar} (max version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt max[abs(a(l)-b(l))]/max[abs(b(l))]} if {\tt b(l)!=0} and
 *	   {\tt max[abs(a(l)-b(l))]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrelmax(const Scalar& a, const Scalar& b) ; 

    //@}
#endif
