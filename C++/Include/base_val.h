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
 * Revision 1.14  2005/09/07 13:09:48  j_novak
 * New method for determining the highest multipole that can be described on a 3D
 * grid.
 *
 * Revision 1.13  2004/11/23 15:03:14  m_forot
 * Corrected error in comments.
 *
 * Revision 1.12  2004/11/04 15:39:45  e_gourgoulhon
 * Modif documentation: added description of bases without any symmetry
 *  in theta.
 *
 * Revision 1.11  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.10  2004/03/22 13:12:40  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.9  2004/01/27 14:31:25  j_novak
 * New method Base_val::mult_sint()
 *
 * Revision 1.8  2004/01/27 14:13:58  j_novak
 * Added method Base_val::mult_x()
 *
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
class Mg3d ;

/**
 * Bases of the spectral expansions. \ingroup (spec)
 * 
 * The class \c Base_val  describes, in each domain, on which basis the
 * spectral expansion of a given function is performed. The corresponding
 * coefficients will be stored in a \c Mtbl_cf. 
 * 
 * The various bases in each of the three dimensions \e r , \f$\theta\f$ and
 * \f$\phi\f$ are identified by an integer defined by a macro in the
 * file \c type_parite.h. These three integers are then merged (by means of 
 * the bitwise OR operator) to give a single integer, stored in 
 * \c Base_val::b[l],  
 * \c l being the domain index. The bases in \e r , \f$\theta\f$ and \f$\phi\f$ can
 * be restored by applying the bitwise AND operator with the masks 
 * \c MSQ_R, \c MSQ_T and \c MSQ_P defined in \c type_parite.h.
 * 
 * The basis functions for expansion with respect to the radial coordinate \f$\xi\f$
 * are coded as follows, \e m being the order of the Fourier expansion in \f$\phi\f$:
 *   \li \c R_CHEB (0x00000001) : Chebyshev polynomials 
 *					    (\e r -sampling: \c FIN);
 *   \li \c R_CHEBP (0x00000002) : Even Chebyshev polynomials 
 *					    (\e r -sampling: \c RARE);
 *   \li \c R_CHEBI (0x00000003) : Odd Chebyshev polynomials 
 *						(\e r -sampling: \c RARE );
 *   \li \c R_CHEBPI_P (0x00000004) : Even (resp. odd) Chebyshev
 *	     polynomials for \e l (spherical harmonic index) even (resp. odd) (\e r -sampling: \c RARE );
 *   \li \c R_CHEBPI_I (0x00000005) : Odd (resp. even) Chebyshev
 *	     polynomials for \e l (spherical harmonic index) even (resp. odd) (\e r -sampling: \c RARE ) ;
 *   \li \c R_CHEBPIM_P (0x00000006) : Even (resp. odd) Chebyshev
 *	     polynomials for \e m even (resp. odd) (\e r -sampling: \c RARE );
 *   \li \c R_CHEBPIM_I (0x00000007) : Odd (resp. even) Chebyshev
 *	     polynomials for \e m even (resp. odd) (\e r -sampling: \c RARE ) ;
 *   \li \c R_CHEBU  (0x00000008) : Chebyshev polynomials 
 *					    (\e r -sampling: \c UNSURR ).
 * 
 * The basis functions for expansion with respect to the co-latitude coordinate 
 * \f$\theta\f$ are coded as follows, \e m being the order of the Fourier expansion 
 * in \f$\phi\f$:  
 *   \li \c T_COSSIN_C  (0x00000100) : \f$\cos(j \theta)\f$ for \e m even,
 *					      \f$\sin(j \theta)\f$ for \e m odd
 *						(\f$\theta\f$-sampling: \c NONSYM);
 *   \li \c T_COSSIN_S  (0x00000200) : \f$\sin(j \theta)\f$ for \e m even,
 *					      \f$\cos(j \theta)\f$ for \e m odd
 *						(\f$\theta\f$-sampling: \c NONSYM);
 *   \li \c T_COS_P  (0x00000500) : \f$\cos(2j \theta)\f$
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_SIN_P  (0x00000600) : \f$\sin(2j \theta)\f$
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_COS_I  (0x00000700) : \f$\cos((2j+1) \theta)\f$
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_SIN_I  (0x00000800) : \f$\sin((2j+1) \theta)\f$
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_COSSIN_CP  (0x00000900) : \f$\cos(2j \theta)\f$ for \e m even,
 *					      \f$\sin((2j+1) \theta)\f$ for \e m odd
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_COSSIN_SP  (0x00000a00) : \f$\sin(2j \theta)\f$ for \e m even,
 *					      \f$\cos((2j+1) \theta)\f$ for \e m odd
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_COSSIN_CI  (0x00000b00) : \f$\cos((2j+1) \theta)\f$ for \e m even,
 *					      \f$\sin(2j \theta)\f$ for \e m odd
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_COSSIN_SI  (0x00000c00) : \f$\sin((2j+1) \theta)\f$ for \e m even,
 *					      \f$\cos(2j \theta)\f$ for \e m odd
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_LEG_P  (0x00000d00) : Associated Legendre functions
 *					  \f$P_{2j}^m(\cos\theta)\f$ for \e m even,
 *					  \f$P_{2j+1}^m(\cos\theta)\f$ for \e m odd
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_LEG_PP  (0x00000e00) : Associated Legendre functions
 *					  \f$P_{2j}^m(\cos\theta)\f$, \e m being 
 *					   always even
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_LEG_I  (0x00000f00) : Associated Legendre functions
 *					  \f$P_{2j+1}^m(\cos\theta)\f$ for \e m even,
 *					  \f$P_{2j}^m(\cos\theta)\f$ for \e m odd
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_LEG_IP  (0x00001000) : Associated Legendre functions
 *					  \f$P_{2j+1}^m(\cos\theta)\f$, \e m being 
 *					   always even
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_LEG_PI  (0x00001100) : Associated Legendre functions
 *					  \f$P_{2j+1}^m(\cos\theta)\f$, \e m being 
 *					   always odd
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_LEG_II  (0x00001200) : Associated Legendre functions
 *					  \f$P_{2j}^m(\cos\theta)\f$, \e m being 
 *					   always odd
 *						(\f$\theta\f$-sampling: \c SYM);
 *   \li \c T_LEG  (0x00001700) : Associated Legendre functions 
 *                              \f$P_l^m(\cos\theta)\f$ of all types
 *
 * The basis functions for expansion with respect to the azimuthal coordinate 
 * \f$\phi\f$ are coded as follows
 *   \li \c P_COSSIN  (0x00010000) : Fourrier series 
 *					   \f$(\cos(m\phi),\ \sin(m\phi))\f$
 *					   (\f$\phi\f$-sampling: \c NONSYM);
 *   \li \c P_COSSIN_P  (0x00020000) : Fourrier series 
 *					   \f$(\cos(m\phi),\ \sin(m\phi))\f$ with
 * 					   even harmonics only (i.e. \e m even)
 *					   (\f$\phi\f$-sampling: \c SYM);
 *   \li \c P_COSSIN_I  (0x00030000) : Fourrier series 
 *					   \f$(\cos(m\phi),\ \sin(m\phi))\f$ with
 * 					   odd harmonics only (i.e. \e m odd)
 *					   (\f$\phi\f$-sampling: \c SYM);
 *
 */

class Base_val {

    // Data 
    // ----
    private:
	int nzone ;	///< Number of domains (zones)

	public:	
	/// Array (size: \c nzone ) of the spectral basis in each domain
	int* b ;	

    // Constructors - Destructor
    // -------------------------
    public:
	explicit Base_val(int nb_of_domains) ;	    ///< Standard constructor
	Base_val(const Base_val& ) ;		    ///< Copy constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	explicit Base_val(FILE* ) ;    	   

	~Base_val() ;			    ///< Destructor


    // Mutators / assignment
    // ---------------------
    public:
	void set_base_nondef() ;    ///< Sets the spectral bases to \c NONDEF 

	/** Sets the expansion basis for \e r  (\f$\xi\f$) functions in a 
	 *  given domain.
	 *  
	 *  @param l	    Domain index
	 *  @param base_r   type of basis functions in \e r  (\f$\xi\f$)
	 *		    (e.g. \c R_CHEB_P , etc..., 
	 *		     see general documentation of class \c Base_val 
	 *		     for denomination of the various bases). 
	 */
	void set_base_r(int l, int base_r) ; 

	/** Sets the expansion basis for \f$\theta\f$ functions in all
	 *  domains.
	 *  
	 *  @param base_t   type of basis functions in \f$\theta\f$
	 *		    (e.g. \c T_COS_P , etc..., 
	 *		     see general documentation of class \c Base_val 
	 *		     for denomination of the various bases). 
	 */
	void set_base_t(int base_t) ; 

	/** Sets the expansion basis for \f$\phi\f$ functions in all
	 *  domains.
	 *  
	 *  @param base_p   type of basis functions in \f$\phi\f$
	 *		    (e.g. \c P_COSSIN , etc..., 
	 *		     see general documentation of class \c Base_val 
	 *		     for denomination of the various bases). 
	 */
	void set_base_p(int base_p) ; 

	void operator=(const Base_val& ) ;	///< Assignment
    
    // Accessors
    // ---------
    public:
	bool operator==(const Base_val& ) const ;  ///< Comparison operator

	/// Returns the code for the expansion basis in domain no. \c l 
	int get_b(int l) const {
	    assert( (l>=0) && (l<nzone) ) ;
	    return b[l] ; 
	}

	/** Returns the expansion basis for \e r  (\f$\xi\f$) functions in the 
	 *  domain of index \c l  
	 *  (e.g. \c R_CHEB_P , etc..., 
	 *   see general documentation of class \c Base_val 
	 *   for denomination of the various bases). 
	 */
	int get_base_r(int l) const {
	    assert( (l>=0) && (l<nzone) ) ;
	    return b[l] & MSQ_R ; 
	} ; 

	/** Returns the expansion basis for \f$\theta\f$ functions in the
	 *  domain of index \c l  
	 *  (e.g. \c T_COS_P , etc..., 
	 *   see general documentation of class \c Base_val 
	 *   for denomination of the various bases). 
	 */
	int get_base_t(int l) const {
	    assert( (l>=0) && (l<nzone) ) ;
	    return b[l] & MSQ_T ; 
	} ; 

	/** Returns the expansion basis for \f$\phi\f$ functions in the 
	 *  domain of index \c l  
	 *  (e.g. \c P_COSSIN , etc..., 
	 *   see general documentation of class \c Base_val 
	 *   for denomination of the various bases). 
	 */
	int get_base_p(int l) const {
	    assert( (l>=0) && (l<nzone) ) ;
	    return b[l] & MSQ_P ; 
	} ; 

	/** Name of the basis function in \e r  (\f$\xi\f$)
	 *
	 *	@param l [input] domain index
	 *	@param k [input] phi index (for the basis in \e r  may depend upon
	 *		the phi index)
	 *	@param j [input] theta index (for the basis in \e r  may depend upon
	 *		the theta index)
	 *	@param i [input] r index
	 *  @param basename [output] string containing the name of the basis function;
	 *		this \c char  array must have a size of (at least) 8 elements
	 *		and must have been allocated before the call
	 *		to \c name_r .
	 */
	void name_r(int l, int k, int j, int i, char* basename) const ; 

	/** Name of the basis function in \f$\theta\f$
	 *
	 *	@param l [input] domain index
	 *	@param k [input] phi index (for the basis in \f$\theta\f$ may depend upon
	 *		the phi index)
	 *	@param j [input] theta index
	 *  @param basename [output] string containing the name of the basis function;
	 *		this \c char  array must have a size of (at least) 8 elements
	 *		and must have been allocated before the call
	 *		to \c name_theta .
	 */
	void name_theta(int l, int k, int j, char* basename) const ; 

	/** Name of the basis function in \f$\varphi\f$
	 *
	 *	@param l [input] domain index
	 *	@param k [input] phi index
	 *  @param basename [output] string containing the name of the basis function;
	 *		this \c char  array must have a size of (at least) 8 elements
	 *		and must have been allocated before the call
	 *		to \c name_phi .
	 */
	void name_phi(int l, int k, char* basename) const ; 


	/** Values of the theta basis functions at the theta collocation points.
	 *  @param l [input] domain index
	 *  @param nt [input] number of theta collocation points 
	 *    (or equivalently number of theta basis functions) in the 
	 *     domain of index \c l 
	 *  @return resu : \c Tbl  3-D containing the values 
	 *    \f$B_i(\theta_j)\f$ of 
	 *    the \c nt  theta basis functions \f$B_i\f$ at the 
	 *    \c nt  collocation points
	 *    \f$\theta_j\f$ in the domain of index \c l . The storage convention
	 *    is the following one : \\
	 *    \c resu(ind_phi,i,j)  = \f$B_i(\theta_j)\f$ with \c ind_phi 
	 *    being a supplementary dimension in case of a dependence in 
	 *    phi of the theta basis : 
	 *    for example, if the theta basis is \c T_COS_P  (no dependence
	 *    in phi), \c resu.get_dim(2)  = 1 and \c ind_phi  can take
	 *    only the value 0; if the theta basis is
	 *    \c T_COSSIN_CP , \c resu.get_dim(2)  = 2, with 
	 *    \c ind_phi  = 0 for \e m even and  
	 *    \c ind_phi  = 1 for \e m odd. 
	 * 
	 */
	const Tbl& theta_functions(int l, int nt) const ; 
	
	/** Values of the phi basis functions at the phi collocation points.
	 *  @param l [input] domain index
	 *  @param np [input] number of phi collocation points 
	 *   (or equivalently number of phi basis functions) in the 
	 *     domain of index \c l 
	 *  @return resu : \c Tbl  2-D containing the values 
	 *    \f$B_i(\phi_k)\f$ of 
	 *    the \c np  phi basis functions \f$B_i\f$ at the \c np  
	 *    collocation points
	 *    \f$\phi_k\f$ in the domain of index \c l . The storage convention
	 *    is the following one : \\
	 *    \c resu(i,k)  = \f$B_i(\phi_k)\f$. 
	 * 
	 */
	const Tbl& phi_functions(int l, int np) const ; 

	/// Returns the number of domains
	int get_nzone() const {return nzone; } ; 
	
     // Manipulation of basis
    // ----------------------
    public:
	/**
	 * The basis is transformed as with a \f$\frac{\partial}{\partial \xi}\f$
	 * operation.
	 */
	void dsdx() ; 
 
	/**
	 * The basis is transformed as with a \f$\frac{1}{\xi}\f$
	 * multiplication.
	 */
	void sx() ;

	/**
	 * The basis is transformed as with a
	 * multiplication by \f$\xi\f$.
	 */
	void mult_x() ;

	/**
	 * The basis is transformed as with a 
	 * \f$\frac{\partial}{\partial \theta}\f$ operation.
	 */
	void dsdt() ;  

	/**
	 * The basis is transformed as with a 
	 * \f$\frac{1}{\sin \theta}\f$ multiplication.
	 */
	void ssint() ;  
  
	/**
	 * The basis is transformed as with a 
	 * \f$\sin \theta\f$ multiplication.
	 */
	void mult_sint() ;  
  
	/**
	 * The basis is transformed as with a transformation to 
	 * \f$Y^l_m\f$ basis.
	 */
	void ylm() ;  

	/**
	 * Computes the various quantum numbers and 1d radial base
	 **/
	void give_quant_numbers (int, int, int, 
				 int&, int&, int&) const ;

	/** 
	 * Returns the highest multipole for a given grid.
	 * @param mgrid : the \c Mg3d 
	 * @param lz : the domain to consider
	 * @return the highest multipolar momentum \e l that can be
	 * described by \c this base on the grid.
	 */
	int give_lmax(const Mg3d& mgrid, int lz) const ;

   // Outputs
    // -------
    public:	    
	void sauve(FILE *) const ;	    ///< Save in a file
    
	friend ostream& operator<<(ostream& , const Base_val& ) ; ///< Display	
	friend Base_val operator*(const Base_val&, const Base_val&) ;
};
ostream& operator<<(ostream& , const Base_val& ) ;
/**
 * \defgroup baseval_m Base_val Mathematics
 * \ingroup (spec)
 * @{
 */

/**
 * This operator is used when calling multiplication or division of \c Valeur .
 * It returns the appropriate base, taking into account the symmetry of the result.
 * The calculation propreties are those of the multiplication of -1 for 
 * antisymmetry and 1 for symmetry.
 * 
 * Should the product of the \c Base_val  not be possible, the result is set to
 * \c ETATNONDEF , (state not defined). It would be the case,  for example,  if
 * one and only one of the \c Valeur is given in spherical harmonics.
 * 
 */
Base_val operator*(const Base_val&, const Base_val&) ;

/** @}*/

#endif
