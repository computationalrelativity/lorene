/*
 *  Definition of Lorene class Base_val
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


#ifndef __BASE_VAL_H_ 
#define __BASE_VAL_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2003/10/20 06:41:43  e_gourgoulhon
 * Corrected documentation.
 *
 * Revision 1.6  2003/10/19 19:42:50  e_gourgoulhon
 * -- member nzone now private
 * -- introduced new methods get_nzone() and get_b()
 * -- introduced new methods name_r, name_theta and name_phi.
 *
 * Revision 1.5  2003/09/16 08:53:05  j_novak
 * Addition of the T_LEG_II base (odd in theta, only for odd m) and the
 * transformation functions from and to T_SIN_P.
 *
 * Revision 1.4  2002/10/16 14:36:28  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/09/13 09:17:31  j_novak
 * Modif. commentaires
 *
 * Revision 1.2  2002/06/17 14:05:15  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.20  2001/10/12  14:56:09  novak
 * Ajout des methodes dsdx(), dsdt(), etc...
 *
 * Revision 2.19  2001/05/29 16:08:45  eric
 * Modif commentaires (mise en conformite Doc++ 3.4.7).
 *
 * Revision 2.18  2000/09/28  10:19:50  eric
 * Modif commentaires (nouvelles bases T_LEG_IP et T_LEG_PI).
 *
 * Revision 2.17  2000/09/08  11:42:55  eric
 * *** empty log message ***
 *
 * Revision 2.16  2000/09/08  10:14:11  eric
 * *** empty log message ***
 *
 * Revision 2.15  2000/09/08  10:11:49  eric
 * *** empty log message ***
 *
 * Revision 2.14  2000/09/08  10:10:24  eric
 * *** empty log message ***
 *
 * Revision 2.13  2000/09/08  10:06:19  eric
 * Ajout des methodes set_base_r, etc... et get_base_r, etc...
 *
 * Revision 2.12  1999/12/29  10:49:12  eric
 * theta_functions et phi_functions declarees const.
 *
 * Revision 2.11  1999/12/29  10:36:58  eric
 * Modif commentaires.
 *
 * Revision 2.10  1999/12/28  12:57:44  eric
 * Introduction des methodes theta_functions et phi_functions.
 *
 * Revision 2.9  1999/11/19  14:53:11  phil
 * *** empty log message ***
 *
 * Revision 2.8  1999/11/18  13:48:35  phil
 * *** empty log message ***
 *
 * Revision 2.7  1999/11/18  12:51:12  phil
 * commentaire de operator*
 *
 * Revision 2.6  1999/10/26  13:23:12  phil
 * *** empty log message ***
 *
 * Revision 2.5  1999/10/26  13:18:06  phil
 * ajout de l'operator*
 *
 * Revision 2.4  1999/10/13  15:49:12  eric
 * *** empty log message ***
 *
 * Revision 2.3  1999/10/12  10:02:17  eric
 * Amelioration des commentaires: description de toutes les bases.
 *
 * Revision 2.2  1999/10/01  15:55:58  eric
 * Depoussierage.
 * Documentation
 *
 * Revision 2.1  1999/09/13  14:38:08  phil
 * ajout de l'operateur ==
 *
 * Revision 2.0  1999/02/22  15:17:47  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

#include <stdio.h>
#include <assert.h>
#include "headcpp.h"


#include "type_parite.h"
class Tbl ;

/**
 * Bases of the spectral expansions.
 * 
 * The class {\tt Base\_val} describes, in each domain, on which basis the
 * spectral expansion of a given function is performed. The corresponding
 * coefficients will be stored in a {\tt Mtbl\_cf}. 
 * 
 * The various bases in each of the three dimensions {\it r}, $\theta$ and
 * $\phi$ are identified by an integer defined by a macro in the
 * file {\tt type\_parite.h}. These three integers are then merged (by means of 
 * the bitwise OR operator) to give a single integer, stored in 
 * {\tt Base\_val::b[l]},  
 * {\tt l} being the domain index. The bases in {\it r}, $\theta$ and $\phi$ can
 * be restored by applying the bitwise AND operator with the masks 
 * {\tt MSQ\_R}, {\tt MSQ\_T} and {\tt MSQ\_P} defined in {\tt type\_parite.h}.
 * 
 * The basis functions for expansion with respect to the radial coordinate $\xi$
 * are coded as follows, {\it m} being the order of the Fourier expansion in $\phi$:
 * \begin{itemize}
 *   \item {\tt R\_CHEB} (0x00000001) : Chebyshev polynomials 
 *					    ({\it r}-sampling: {\tt FIN});
 *   \item {\tt R\_CHEBP} (0x00000002) : Even Chebyshev polynomials 
 *					    ({\it r}-sampling: {\tt RARE});
 *   \item {\tt R\_CHEBI} (0x00000003) : Odd Chebyshev polynomials 
 *						({\it r}-sampling: {\tt RARE});
 *   \item {\tt R\_CHEBPIM\_P} (0x00000006) : Even (resp. odd) Chebyshev
 *	     polynomials for {\it m} even (resp. odd) ({\it r}-sampling: {\tt RARE});
 *   \item {\tt R\_CHEBPIM\_I} (0x00000007) : Odd (resp. even) Chebyshev
 *	     polynomials for {\it m} even (resp. odd) ({\it r}-sampling: {\tt RARE}) ;
 *   \item {\tt R\_CHEBU} (0x00000008) : Chebyshev polynomials 
 *					    ({\it r}-sampling: {\tt UNSURR}).
 * \end{itemize}
 * 
 * The basis functions for expansion with respect to the co-latitude coordinate 
 * $\theta$ are coded as follows, {\it m} being the order of the Fourier expansion 
 * in $\phi$:  
 * \begin{itemize}
 *   \item {\tt T\_COS\_P} (0x00000500) : $\cos(2j \theta)$
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_SIN\_P} (0x00000600) : $\sin(2j \theta)$
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_COS\_I} (0x00000700) : $\cos((2j+1) \theta)$
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_SIN\_I} (0x00000800) : $\sin((2j+1) \theta)$
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_COSSIN\_CP} (0x00000900) : $\cos(2j \theta)$ for {\it m} even,
 *					      $\sin((2j+1) \theta)$ for {\it m} odd
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_COSSIN\_SP} (0x00000a00) : $\sin(2j \theta)$ for {\it m} even,
 *					      $\cos((2j+1) \theta)$ for {\it m} odd
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_COSSIN\_CI} (0x00000b00) : $\cos((2j+1) \theta)$ for {\it m} even,
 *					      $\sin(2j \theta)$ for {\it m} odd
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_COSSIN\_SI} (0x00000c00) : $\sin((2j+1) \theta)$ for {\it m} even,
 *					      $\cos(2j \theta)$ for {\it m} odd
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_LEG\_P} (0x00000d00) : Associated Legendre functions
 *					  $P_{2j}^m(\cos\theta)$ for {\it m} even,
 *					  $P_{2j+1}^m(\cos\theta)$ for {\it m} odd
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_LEG\_PP} (0x00000e00) : Associated Legendre functions
 *					  $P_{2j}^m(\cos\theta)$, {\it m} being 
 *					   always even
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_LEG\_I} (0x00000f00) : Associated Legendre functions
 *					  $P_{2j+1}^m(\cos\theta)$ for {\it m} even,
 *					  $P_{2j}^m(\cos\theta)$ for {\it m} odd
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_LEG\_IP} (0x00001000) : Associated Legendre functions
 *					  $P_{2j+1}^m(\cos\theta)$, {\it m} being 
 *					   always even
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_LEG\_PI} (0x00001100) : Associated Legendre functions
 *					  $P_{2j+1}^m(\cos\theta)$, {\it m} being 
 *					   always odd
 *						($\theta$-sampling: {\tt SYM});
 *   \item {\tt T\_LEG\_II} (0x00001200) : Associated Legendre functions
 *					  $P_{2j}^m(\cos\theta)$, {\it m} being 
 *					   always odd
 *						($\theta$-sampling: {\tt SYM});
 * \end{itemize}
 *
 * The basis functions for expansion with respect to the azimuthal coordinate 
 * $\phi$ are coded as follows
 * \begin{itemize}
 *   \item {\tt P\_COSSIN} (0x00010000) : Fourrier series 
 *					   $(\cos(m\phi),\ \sin(m\phi))$
 *					   ($\phi$-sampling: {\tt NONSYM});
 *   \item {\tt P\_COSSIN\_P} (0x00020000) : Fourrier series 
 *					   $(\cos(m\phi),\ \sin(m\phi))$ with
 * 					   even harmonics only (i.e. {\it m} even)
 *					   ($\phi$-sampling: {\tt SYM});
 *   \item {\tt P\_COSSIN\_I} (0x00030000) : Fourrier series 
 *					   $(\cos(m\phi),\ \sin(m\phi))$ with
 * 					   odd harmonics only (i.e. {\it m} odd)
 *					   ($\phi$-sampling: {\tt SYM});
 * \end{itemize}
 *
 */

class Base_val {

    // Data 
    // ----
    private:
	int nzone ;	/// Number of domains (zones)

	public:	
	/// Array (size: {\tt nzone}) of the spectral basis in each domain
	int* b ;	

    // Constructors - Destructor
    // -------------------------
    public:
	explicit Base_val(int nb_of_domains) ;	    /// Standard constructor
	Base_val(const Base_val& ) ;		    /// Copy constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	explicit Base_val(FILE* ) ;    	   

	~Base_val() ;			    /// Destructor


    // Mutators / assignment
    // ---------------------
    public:
	void set_base_nondef() ;    /// Sets the spectral bases to {\tt NONDEF}

	/** Sets the expansion basis for {\it r} ($\xi$) functions in a 
	 *  given domain.
	 *  
	 *  @param l	    Domain index
	 *  @param base_r   type of basis functions in {\it r} ($\xi$)
	 *		    (e.g. {\tt R\_CHEB\_P}, etc..., 
	 *		     see general documentation of class {\tt Base\_val}
	 *		     for denomination of the various bases). 
	 */
	void set_base_r(int l, int base_r) ; 

	/** Sets the expansion basis for $\theta$ functions in all
	 *  domains.
	 *  
	 *  @param base_t   type of basis functions in $\theta$
	 *		    (e.g. {\tt T\_COS\_P}, etc..., 
	 *		     see general documentation of class {\tt Base\_val}
	 *		     for denomination of the various bases). 
	 */
	void set_base_t(int base_t) ; 

	/** Sets the expansion basis for $\phi$ functions in all
	 *  domains.
	 *  
	 *  @param base_p   type of basis functions in $\phi$
	 *		    (e.g. {\tt P\_COSSIN}, etc..., 
	 *		     see general documentation of class {\tt Base\_val}
	 *		     for denomination of the various bases). 
	 */
	void set_base_p(int base_p) ; 

	void operator=(const Base_val& ) ;	/// Assignment
    
    // Accessors
    // ---------
    public:
	bool operator==(const Base_val& ) const ;  /// Comparison operator

	/// Returns the code for the expansion basis in domain no. {\tt l}
	int get_b(int l) const {
	    assert( (l>=0) && (l<nzone) ) ;
	    return b[l] ; 
	}

	/** Returns the expansion basis for {\it r} ($\xi$) functions in the 
	 *  domain of index {\tt l} 
	 *  (e.g. {\tt R\_CHEB\_P}, etc..., 
	 *   see general documentation of class {\tt Base\_val}
	 *   for denomination of the various bases). 
	 */
	int get_base_r(int l) const {
	    assert( (l>=0) && (l<nzone) ) ;
	    return b[l] & MSQ_R ; 
	} ; 

	/** Returns the expansion basis for $\theta$ functions in the
	 *  domain of index {\tt l} 
	 *  (e.g. {\tt T\_COS\_P}, etc..., 
	 *   see general documentation of class {\tt Base\_val}
	 *   for denomination of the various bases). 
	 */
	int get_base_t(int l) const {
	    assert( (l>=0) && (l<nzone) ) ;
	    return b[l] & MSQ_T ; 
	} ; 

	/** Returns the expansion basis for $\phi$ functions in the 
	 *  domain of index {\tt l} 
	 *  (e.g. {\tt P\_COSSIN}, etc..., 
	 *   see general documentation of class {\tt Base\_val}
	 *   for denomination of the various bases). 
	 */
	int get_base_p(int l) const {
	    assert( (l>=0) && (l<nzone) ) ;
	    return b[l] & MSQ_P ; 
	} ; 

	/** Name of the basis function in {\it r} ($\xi$)
	 *
	 *	@param l [input] domain index
	 *	@param k [input] phi index (for the basis in {\it r} may depend upon
	 *		the phi index)
	 *	@param j [input] theta index (for the basis in {\it r} may depend upon
	 *		the theta index)
	 *	@param i [input] r index
	 *  @param basename [output] string containing the name of the basis function;
	 *		this {\tt char} array must have a size of (at least) 8 elements
	 *		and must have been allocated before the call
	 *		to {\tt name\_r}.
	 */
	void name_r(int l, int k, int j, int i, char* basename) const ; 

	/** Name of the basis function in $\theta$
	 *
	 *	@param l [input] domain index
	 *	@param k [input] phi index (for the basis in $\theta$ may depend upon
	 *		the phi index)
	 *	@param j [input] theta index
	 *  @param basename [output] string containing the name of the basis function;
	 *		this {\tt char} array must have a size of (at least) 8 elements
	 *		and must have been allocated before the call
	 *		to {\tt name\_theta}.
	 */
	void name_theta(int l, int k, int j, char* basename) const ; 

	/** Name of the basis function in $\varphi$
	 *
	 *	@param l [input] domain index
	 *	@param k [input] phi index
	 *  @param basename [output] string containing the name of the basis function;
	 *		this {\tt char} array must have a size of (at least) 8 elements
	 *		and must have been allocated before the call
	 *		to {\tt name\_phi}.
	 */
	void name_phi(int l, int k, char* basename) const ; 


	/** Values of the theta basis functions at the theta collocation points.
	 *  @param l [input] domain index
	 *  @param nt [input] number of theta collocation points 
	 *    (or equivalently number of theta basis functions) in the 
	 *     domain of index {\tt l}
	 *  @return resu : {\tt Tbl} 3-D containing the values 
	 *    $B_i(\theta_j)$ of 
	 *    the {\tt nt} theta basis functions $B_i$ at the 
	 *    {\tt nt} collocation points
	 *    $\theta_j$ in the domain of index {\tt l}. The storage convention
	 *    is the following one : \\
	 *    {\tt resu(ind\_phi, i, j)} = $B_i(\theta_j)$ with {\tt ind\_phi}
	 *    being a supplementary dimension in case of a dependence in 
	 *    phi of the theta basis : 
	 *    for example, if the theta basis is {\tt T\_COS\_P} (no dependence
	 *    in phi), {\tt resu.get\_dim(2)} = 1 and {\tt ind\_phi} can take
	 *    only the value 0; if the theta basis is
	 *    {\tt T\_COSSIN\_CP}, {\tt resu.get\_dim(2)} = 2, with 
	 *    {\tt ind\_phi} = 0 for {\it m} even and  
	 *    {\tt ind\_phi} = 1 for {\it m} odd. 
	 * 
	 */
	const Tbl& theta_functions(int l, int nt) const ; 
	
	/** Values of the phi basis functions at the phi collocation points.
	 *  @param l [input] domain index
	 *  @param np [input] number of phi collocation points 
	 *   (or equivalently number of phi basis functions) in the 
	 *     domain of index {\tt l}
	 *  @return resu : {\tt Tbl} 2-D containing the values 
	 *    $B_i(\phi_k)$ of 
	 *    the {\tt np} phi basis functions $B_i$ at the {\tt np} 
	 *    collocation points
	 *    $\phi_k$ in the domain of index {\tt l}. The storage convention
	 *    is the following one : \\
	 *    {\tt resu(i, k)} = $B_i(\phi_k)$. 
	 * 
	 */
	const Tbl& phi_functions(int l, int np) const ; 

	/// Returns the number of domains
	int get_nzone() const {return nzone; } ; 
	
     // Manipulation of basis
    // ----------------------
    public:
	/**
	 * The basis is transformed as with a $\frac{\partial}{\partial \xi}$
	 * operation.
	 */
	void dsdx() ; 
 
	/**
	 * The basis is transformed as with a $\frac{1}{\xi}$
	 * multiplication.
	 */
	void sx() { this->dsdx(); };

	/**
	 * The basis is transformed as with a 
	 * $\frac{\partial}{\partial \theta}$ operation.
	 */
	void dsdt() ;  

	/**
	 * The basis is transformed as with a 
	 * $\frac{1}{\sin \theta}$ multiplication.
	 */
	void ssint() ;  
  
	/**
	 * The basis is transformed as with a transformation to 
	 * $Y^l_m$ basis.
	 */
	void ylm() ;  

   // Outputs
    // -------
    public:	    
	void sauve(FILE *) const ;	    /// Save in a file
    
	friend ostream& operator<<(ostream& , const Base_val& ) ; /// Display	
	friend Base_val operator*(const Base_val&, const Base_val&) ;
};
ostream& operator<<(ostream& , const Base_val& ) ;
/**
 * @name Base\_val arithmetic.
 */
    //@{

/**
 * This operator is used when calling multiplication or division of {\tt Valeur}.
 * It returns the appropriate base, taking into account the symmetry of the result.
 * The calculation propreties are those of the multiplication of -1 for 
 * antisymmetry and 1 for symmetry.
 * 
 * Should the product of the {\tt Base\_val} not be possible, the result is set to
 * {\tt ETATNONDEF}, (state not defined). It would be the case,  for example,  if
 * one and only one of the {\tt Valeur} is given in spherical harmonics.
 * 
 */

Base_val operator*(const Base_val&, const Base_val&) ;

   //@}

#endif
