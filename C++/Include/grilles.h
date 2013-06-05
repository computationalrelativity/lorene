/*
 *  Definition of Lorene classes Grille3d and Mg3d
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


#ifndef __GRILLES_H_ 
#define __GRILLES_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.21  2013/06/05 15:00:26  j_novak
 * Suppression of all classes derived from Grille3d. Now Grille3d is no
 * longer an abstract class. r-samplings are only one of RARE, FIN or
 * UNSURR (FINJAC has been removed). Instead, Mg3d possesses a new member
 * colloc_r[nzone] defining the type of collocation points (spectral
 * bases) in each domain.
 *
 * Revision 1.20  2013/01/11 15:44:53  j_novak
 * Addition of Legendre bases (part 2).
 *
 * Revision 1.19  2012/01/17 10:10:15  j_penner
 * added a constructor in which the nucleus and outer domain are both of type FIN
 *
 * Revision 1.18  2008/10/29 08:17:51  jl_cornou
 * Standard spectral bases for pseudo vectors added
 *
 * Revision 1.17  2008/02/18 13:53:37  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.16  2007/12/11 15:28:05  jl_cornou
 * Jacobi(0,2) polynomials partially implemented
 *
 * Revision 1.15  2006/05/17 13:17:02  j_novak
 * New member g_angu_1dom, the one-domain angular grid associated with the
 * current grid.
 *
 * Revision 1.14  2005/10/25 08:56:34  p_grandclement
 * addition of std_spectral_base in the case of odd functions near the origin
 *
 * Revision 1.13  2005/10/07 08:47:20  j_novak
 * Addition of the pointer g_non_axi on a grid, with at least 5 points in the
 * theta direction and 4 in the phi one (for tensor rotations).
 *
 * Revision 1.12  2005/03/25 14:54:04  e_gourgoulhon
 * Corrected documentation.
 *
 * Revision 1.11  2004/07/06 13:36:27  j_novak
 * Added methods for desaliased product (operator |) only in r direction.
 *
 * Revision 1.10  2004/06/22 08:49:56  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.9  2004/03/22 13:12:41  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.8  2003/06/20 14:16:57  f_limousin
 * Add the operator== to compare two Mg3d
 *
 * Revision 1.7  2003/06/18 08:45:26  j_novak
 * In class Mg3d: added the member get_radial, returning only a radial grid
 * For dAlembert solver: the way the coefficients of the operator are defined has been changed.
 *
 * Revision 1.6  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.5  2002/08/13 08:02:45  j_novak
 * Handling of spherical vector/tensor components added in the classes
 * Mg3d and Tenseur. Minor corrections for the class Metconf.
 *
 * Revision 1.4  2002/06/17 14:05:16  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.3  2001/12/12 09:23:46  e_gourgoulhon
 * Parameter compact added to the simplified constructor of class Mg3d
 *
 * Revision 1.2  2001/12/11 06:47:42  e_gourgoulhon
 * Simplified constructor for class Mg3d
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.10  2001/06/13  14:23:40  eric
 * Les fonctions Mg3d::del_deriv() et Mg3d::set_deriv_0x0() ne sont plus
 * virtuelles puisque Mg3d n'a aucune classe derivee.
 *
 * Revision 2.9  2001/05/26  13:24:49  eric
 * Ajout du membre g_twice (grille double pour le desaliasing)
 * Modif de la declaration de g_angu (pointeur mutable)
 *   g_twice et g_angu ne sont calcules que si necessaire (cad si
 *   on appelle la fonction get_twice() ou get_angu()).
 *
 * Revision 2.8  1999/11/16  14:15:57  eric
 * Ajout de la fonction Mg3d::get_angu().
 *
 * Revision 2.7  1999/10/12  14:54:11  eric
 * Ajout du membre Base_val std_base_scal() const.
 *
 * Revision 2.6  1999/10/01  10:35:42  eric
 * Amelioration des commentaires.
 *
 * Revision 2.5  1999/09/30  14:58:00  eric
 * Operator!= declare const/
 *
 * Revision 2.4  1999/09/30  14:11:43  eric
 * sauve et std_base_vect_cart declarees const.
 *
 * Revision 2.3  1999/09/30  12:52:38  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.2  1999/09/24  14:23:24  eric
 * Declaration de methodes const.
 *
 * Revision 2.1  1999/09/14  15:24:04  phil
 * ajout de std_base_vect_cart
 *
 * Revision 2.0  1999/02/15  10:41:51  hyc
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/15  09:59:50  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1998/12/01  14:28:00  hyc
 * Version 2
 *
 *
 * $Header$
 *
 */

// Classes utilisees

// Fichiers includes
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "headcpp.h"

#include "type_parite.h"

class Base_val ; 

		    	//-------------//
		    	// Mono-grille //
		    	//-------------//

// Classe de base
/**
 * 3D grid class in one domain.\ingroup (spec)
 * 
 * Basic 3D spherical grid class in spherical coordinates \f$(r,\theta,\phi)\f$. 
 * The radial coordinate \f$\xi\f$ lies in the range [0, 1] or [-1, 1]
 * depending upon the sampling (\c RARE or \c FIN ). Its relation with the 
 * physical radial coordinate \e r is defined by the mapping (cf. class \c Map) 
 * and is described in Bonazzola, Gourgoulhon \& Marck, \a Phys. \a Rev. \a D
 * \b 58, 104020 (1998).
 * Note: this monogrid should not be used. Use instead \c Mg3d.
 *
 * @version #$Id$#
 */

class Grille3d {
    protected:
	const int nr ;	///< Number of points in \e r (\f$\xi\f$)
	const int nt ;	///< Number of points in \f$\theta\f$
	const int np ;	///< Number of points in \f$\phi\f$

	///Type of sampling in \e r (\f$\xi\f$) (\c RARE,\c FIN,\c UNSURR )
	int type_r ;	
	int type_t ;	///< Type of sampling in \f$\theta\f$ (\c SYM,\c NONSYM)
	int type_p ;	///< Type of sampling in \f$\phi\f$ (\c SYM,\c NONSYM)
	/// Type of radial spectral basis (\c BASE_CHEB, \c BASE_LEG, BASE_JAC02 )
	int base_r ; 
   
    public:
	/// Array of values of \f$\xi\f$ at the \c nr collocation points
	double* x ;	
	/// Array of values of \f$\theta\f$ at the \c nt collocation points
	double* tet ;	
	/// Array of values of \f$\phi\f$ at the \c np collocation points
	double* phi ;	

	/// Constructor 
	Grille3d(int n_r, int n_t, int n_p, int typer, int typet, 
		 int typep, int baser) ;
    
	/// Copy constructor 
	Grille3d(const Grille3d& ) ;
	
	/// Assignement operator 
	void operator=(const Grille3d& ) ;
	 	
    public:
	virtual ~Grille3d() ;		///< Destructor

    public:
	/// Returns \c nr
    	int get_nr() const {return nr ;} ;
	/// Returns \c nt	    
    	int get_nt() const {return nt ;} ;
	/// Returns \c np
    	int get_np() const {return np ;} ;

	/// Returns \c type_r
    	int get_type_r() const {return type_r ;} ; 
	/// Returns \c type_t
    	int get_type_t() const {return type_t ;} ; 
	/// Returns \c type_p
    	int get_type_p() const {return type_p ;} ; 
	/// Returns \c base_r
	int get_base_r() const {return base_r ;} ;

     protected:
	/// Computes the collocation point coordinates in the radial direction
	void compute_radial_grid() ;
};



		    	//---------------//
		    	// Multi-grilles //
		    	//---------------//

/**
 * Multi-domain grid. \ingroup (spec)
 *
 * A multi-domain grid is a set of 3D grids for the implementation of 
 * the multi-domain spectral method described in Bonazzola, Gourgoulhon 
 * \& Marck, \a Phys. \a Rev. \a D \b 58, 104020 (1998).
 * Each domain is represented by a 3D mono grid (\c Grille3d) and
 * is called a \e zone. 
 * For each direction, the Number of Degrees of Freedom (NDF) is 
 * \e a \e priori independent of the zone. However, some methods or 
 * routines may refuse to work if the NDF of some \f$(\theta, \phi)\f$ 
 * direction is not identical in all the zones. 
 * This holds for the type of sampling (symmetry) too.
 *
 * @version #$Id$#
 */
 
class Mg3d {

    // Data  
    // ----
    protected:
	int nzone ;	///< Number of domains (zones)
	
	int* nr ;	///< Array (size: \c nzone) of nb. of points in \e r (\f$\xi\f$)
	int* nt ;	///< Array (size: \c nzone) of nb. of points in \f$\theta\f$
	int* np ;	///< Array (size: \c nzone) of nb. of points in \f$\phi\f$
	
	/** Array (size: \c nzone) of type of sampling in \e r (\f$\xi\f$) 
	 *(\c RARE,\c FIN, c UNSURR)
	 */
	int* type_r ;	

	/// Type of sampling in \f$\theta\f$ (\c SYM, \c NONSYM)
	int type_t ;

	/// Type of sampling in \f$\phi\f$ (\c SYM, \c NONSYM)	
	int type_p ;

	/** Array (size: \c nzone)  of type of collocation points 
	 * in \e r (\f$\xi\f$) and related decompoisition bases 
	 * (\c BASE_CHEB, \c BASE_LEG, \c BASE_JAC02 ).
	 */
	int* colloc_r ;
	
	/// Array (size: \c nzone) of pointers on the \c Grille3d's
	Grille3d** g ;	

	mutable Mg3d* g_angu ;	///< Pointer on the associated angular grid
	///Pointer on the associated angular grid with only one domain
	mutable Mg3d* g_angu_1dom ; 
	mutable Mg3d* g_radial ; ///< Pointer on the associated radial grid
	
	/** Pointer on the grid which has twice the number of points in
	 *  each dimension (for desaliasing).
	 */
	mutable Mg3d* g_twice ; 

	/** Pointer on the grid which has 50% more points in
	 *  \e r dimension (for desaliasing).
	 */
	mutable Mg3d* g_plus_half ; 

	/** Pointer on the grid which has at least 4 points in
	 *  the \f$\phi\f$ direction and at least 5 in the \f$\theta\f$ 
	 *  direction (for tensor rotations).
	 */
	mutable Mg3d* g_non_axi ; 

    // Constructors - Destructor
    // -------------------------
	
    public:

/**
 * General constructor.
 *
 * @param   nz	    [input] Number of domains (zones).
 * @param   nbr[]   [input] Array (size: \c nz ) of number of degree of
 *		    freedom (NDF) in \e r-direction
 * @param   typr[]  [input] Array (size: \c nz ) of type of sampling in \e r -direction
 * @param   nbt[]   [input] Array (size: \c nz ) of NDF in \f$\theta\f$-direction
 * @param   typt    [input] Type of sampling in \f$\theta\f$-direction
 * @param   nbp[]   [input] Array (size: \c nz) of NDF in \f$\phi\f$-direction
 * @param   typp    [input] Type of sampling in \f$\phi\f$-direction
 * @param   base_r  [input] Types of \e r bases in each domain, to define the
 *                          collocation points: \c BASE_CHEB , \c BASE_LEG or 
 *                          \c BASE_JAC02. If the pointer is null, BASE_CHEB is set.
 *
 */
	Mg3d(int nz, int nbr[], int typr[], int nbt[], int typt, int nbp[],
	     int typp, int* base_r = 0x0) ;

/**
 * Simplified constructor for a standard multi-grid.
 * This provides a multi-grid with the same number of degrees of freedom
 * in all the domains. \n
 * The domain of index \c l = 0 will be a nucleus:
 * \f$\xi\in [0,1]\f$, rarefied sampling (type \c RARE) near the origin; \n
 * domains of indices \f$ 1 \le {\tt l} \le {\tt nz}-2\f$ will be shells:
 * \f$\xi\in [-1,1]\f$, dense sampling (type \c FIN) near -1 and 1; \n
 * if \c compact == true, 
 * the domains of index \c l = \c nz-1 will be the outermost compactified
 * shell:
 * \f$\xi\in [-1,1]\f$, dense sampling (type \c UNSURR) near -1 and 1
 * for a \e 1/r discretization.
 *
 *
 * @param   nz	    [input] Number of domains (zones).
 * @param   nbr     [input] Number of degree of freedom (NDF) in
 *				\e r -direction in each domain
 * @param   nbt     [input] Number of degree of freedom (NDF) in
 *				\f$\theta\f$-direction in each domain
 * @param   nbp     [input] Number of degree of freedom (NDF) in
 *				\f$\phi\f$-direction in each domain
 * @param   typt    [input] Type of sampling in \f$\theta\f$-direction:  \n
 *				\c SYM  for a sampling in \f$[0,\pi/2]\f$
 *			(symmetry with respect to the equatorial plane), \n
 *			\c NONSYM  for a sampling in \f$[0,\pi]\f$
 * @param   typp    [input] Type of sampling in \f$\phi\f$-direction: \n
 *			\c SYM  for a sampling in \f$[0,\pi[\f$
 *			(symmetry with respect to a \f$\pi\f$ translation
 *			 in \f$\phi\f$) \n
 *			\c NONSYM  for a sampling in \f$[0,2\pi[\f$
 * @param  compact [input] \c true  for the last domain to have 
 *			a \e 1/r  sampling (\c UNSURR ) instead of a
 *			\e r  sampling (\c FIN ).
 * @param  legendre [input] \c true if the nucleus and shells have Legendre-type
 *                          collocation points. If \c false, Chebyshev ones shall
 *                          be used.
 */
	Mg3d(int nz, int nbr, int nbt, int nbp, int typt, int typp, 
	     bool compact, bool legendre=false) ;

/**
 * Simplified constructor for a standard multi-grid only with shells and Chebyshev.
 * This provides a multi-grid with the same number of degrees of freedom
 * in all the domains. \n
 * ALL DOMAINS ARE TREATED AS SHELLS
 * \f$\xi\in [-1,1]\f$, dense sampling (type \c FIN) near -1 and 1; \n and only 
 * Chebyshev bases are considered.
 *
 * @param   nz	    [input] Number of domains (zones).
 * @param   nbr     [input] Number of degree of freedom (NDF) in
 *				\e r -direction in each domain
 * @param   nbt     [input] Number of degree of freedom (NDF) in
 *				\f$\theta\f$-direction in each domain
 * @param   nbp     [input] Number of degree of freedom (NDF) in
 *				\f$\phi\f$-direction in each domain
 * @param   typt    [input] Type of sampling in \f$\theta\f$-direction:  \n
 *				\c SYM  for a sampling in \f$[0,\pi/2]\f$
 *			(symmetry with respect to the equatorial plane), \n
 *			\c NONSYM  for a sampling in \f$[0,\pi]\f$
 * @param   typp    [input] Type of sampling in \f$\phi\f$-direction: \n
 *			\c SYM  for a sampling in \f$[0,\pi[\f$
 *			(symmetry with respect to a \f$\pi\f$ translation
 *			 in \f$\phi\f$) \n
 *			\c NONSYM  for a sampling in \f$[0,2\pi[\f$
 */
	Mg3d(int nz, int nbr, int nbt, int nbp, int typt, int typp) ;

	/**Constructor from a file (see \c sauve(FILE*)).
	 * If the boolean flag \c read_base is \c false (default) the 
	 * spectral basis is not saved, for compatibility with older 
	 * versions. If it is set to \c true , the basis is saved.
	 */
	Mg3d(FILE* fd, bool read_base=false) ;	 
	
   public:
	/**
	 * Copy constructor (private and not implemented to make \c Mg3d  a 
	 * non-copyable class)
	 */
	Mg3d(const Mg3d& ) ;
	
    public:
	
	~Mg3d() ;   ///< Destructor
		
    // Assignement
    // -----------
    private:
	/** 
	 * Assignement operator (private and not implemented to make \c Mg3d 
	 * a non-copyable class)
	 */
	void operator=(const Mg3d& ) ;
	 	
    // Extraction of information
    // -------------------------
    public:
   	/// Returns the number of domains
	int get_nzone() const { 	 
	    return nzone ;
	} ;
	/// Returns the number of points in the radial direction (\f$\xi\f$) in domain no. \e l 
	int get_nr(int l) const { 
	    assert(l>=0 && l<nzone) ;
	    return nr[l] ;
	} ;
	/// Returns the number of points in the co-latitude direction (\f$\theta\f$) in domain no. \e l 
	int get_nt(int l) const { 
	    assert(l>=0 && l<nzone) ;
	    return nt[l] ;
	} ;
	/// Returns the number of points in the azimuthal direction (\f$\phi\f$) in domain no. \e l 
	int get_np(int l) const { 
	    assert(l>=0 && l<nzone) ;
	    return np[l] ;
	} ;
    
	/** Returns the type of sampling in the radial direction
     * in domain no.\e l : \n
     *   \c RARE : \f$\xi\in[0,1]\f$ : rarefied at the origin  \n
     *   \c FIN : \f$\xi\in[-1,1]\f$ :  dense at the two extremities \n
     *   \c UNSURR : \f$\xi\in[-1,1]\f$ : dense at the two extremities, 
     *      in view of using \f$u=1/r\f$ as radial variable 
     */
	int get_type_r(int l) const {
	    assert(l>=0 && l<nzone) ;
	    return type_r[l] ;
	} ;

	/** Returns the type of sampling in the \f$\theta\f$ direction: \n
     *   \c SYM : \f$\theta\in[0,\pi/2]\f$ : symmetry with respect to 
     *      the equatorial plane \n
     *   \c NONSYM : \f$\theta\in[0,\pi]\f$ : no symmetry with respect 
     *              to the equatorial plane 
     */
	int get_type_t() const { 	
	    return type_t ;
	} ;

	/** Returns the type of sampling in the \f$\phi\f$ direction: \n
     *   \c SYM : \f$\phi\in[0,\pi[\f$ : symmetry with respect to the 
     *              transformation   \f$ \phi \mapsto \phi + \pi\f$   \n
     *   \c NONSYM : \f$\phi\in[0,2\pi[\f$ :  no symmetry with respect to 
     *              the transformation    \f$ \phi \mapsto \phi + \pi\f$  
     */
	int get_type_p() const {
	    return type_p ;
	} ;
	
	/// Returns a pointer on the 3D mono-grid for domain no. l
	const Grille3d* get_grille3d(int l) const { 
	    assert(l>=0 && l<nzone) ;
	    return g[l] ;
	} ;

	/** Returns the type of collocation points used in domain no. \e l.
	 * This type depends on the spectral decomposition basis:
	 * \li \c BASE_CHEB : Chebyshev basis and collocation points
	 * \li \c BASE_LEG : Legendre basis and collocation points
	 * \li \c BASE_JAC02 : Jacobi (0,2) polynomials
	 */
	int get_colloc_r(int l) const {
	  assert(l>=0 && l<nzone) ;
	  return colloc_r[l] ;
	}

	/// Returns the pointer on the associated angular grid
	const Mg3d* get_angu() const ;
	
	/** Returns the pointer on the associated mono-domain angular grid.
	 * This angular grid corresponds to the first domain (index 0).
	 */
	const Mg3d* get_angu_1dom() const ;
	
	/// Returns the pointer on the associated radial grid
	const Mg3d* get_radial() const ;
	
	/** Returns the pointer on the grid which has twice the number
	 *  of points in each dimension (for desaliasing).
	 */
	const Mg3d* get_twice() const ;

	/** Returns the pointer on the grid which has 50% more points in
	 *  \e r dimension (for desaliasing).
	 */
	const Mg3d* plus_half() const ;

	/** Returns the pointer on the grid which has at least 4 points in
	 *  the \f$\phi\f$ direction and at least 5 in the \f$\theta\f$ 
	 *  direction (for tensor rotations).
	 */
	const Mg3d* get_non_axi() const ;

	/// Comparison operator (egality)
	bool operator==(const Mg3d& ) const ;  


	
    // Outputs
    // -------
    public: 
	/** Saves into a file.
	 * By default, if \c save_base is false, the spectral decomposition
	 * basis \c colloc_r is not saved (for compatibility reasons).
	 * If \c save_base is \c true, then \c colloc_r is save, too.
	 * The same value for the flag \c read_base must be given in the
	 * contructor reading the file.
	 */
	void sauve(FILE* fd, bool save_base=false) const ; 
	
	friend ostream& operator<<(ostream& , const Mg3d & ) ;	///< Display
	
    // Management of derived quantities
    // --------------------------------
    protected:
	/** Deletes all the derived quantities 
	 *   (\c g_radial , \c g_angu, \c g_twice, ...)
	 */
	void del_deriv() const ; 
	
	/** Sets to \c 0x0  all the pointers on derived quantities
	 *   (\c g_radial , \c g_angu, \c g_twice, ... )
	 */
	void set_deriv_0x0() const ; 


    // Miscellaneous
    // -------------
    public:
	bool operator!=(const Mg3d & ) const ;  ///< Operator !=

	/// Returns the standard spectral bases for a scalar
	Base_val std_base_scal() const ; 	
	
	/// Returns the standard odd spectral bases for a scalar
	Base_val std_base_scal_odd() const ; 	
	
	/** Returns the standard spectral bases for the Cartesian components 
	 *  of a vector
	 */ 
	Base_val** std_base_vect_cart() const ;

	/** Returns the standard spectral bases for the spherical components 
	 *  of a vector
	 */ 
	Base_val** std_base_vect_spher() const ;

	/** Returns the standard spectral bases for the Cartesian components 
	 *  of a pseudo-vector
	 */ 
	Base_val** pseudo_base_vect_cart() const ;

	/** Returns the standard spectral bases for the spherical components 
	 *  of a pseudo-vector
	 */ 
	Base_val** pseudo_base_vect_spher() const ;

};
ostream& operator<<(ostream& , const Mg3d & ) ;


//======================================
// One domain standard bases definitions
//======================================
int std_base_scal_1z(int type_r, int type_t, int type_p) ; 
int std_base_scal_odd_1z(int type_r, int type_t, int type_p) ; 
int leg_base_scal_1z(int type_r, int type_t, int type_p) ; 
int leg_base_scal_odd_1z(int type_r, int type_t, int type_p) ; 
int jac02_base_scal_1z(int type_r, int type_t, int type_p) ; 
int jac02_base_scal_odd_1z(int type_r, int type_t, int type_p) ; 

#endif 
