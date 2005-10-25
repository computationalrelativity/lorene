/*
 *  Definition of Lorene classes Grille3d
 *				 Grille3d_*
 *				 Mg3d
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

#include "indent.h"
#include "type_parite.h"

class Base_val ; 

		    	//-------------//
		    	// Mono-grille //
		    	//-------------//

// Classe de base
/**
 * Base 3D grid class.\ingroup (spec)
 * 
 * This class is an abstract one. It can't be instanciated. 
 * There is no \c ostream 
 * operator. For historical reasons, the mnemonics refer to spherical 
 * coordinates \f$(r,\theta,\phi)\f$. However all this stuff can be used for 
 * any coordinate system if the appropriated constructors are created.
 *
 * The radial coordinate \f$\xi\f$ lies in the range [0, 1] or [-1, 1]
 * depending upon the sampling. Its relation with the physical radial 
 * coordinate \e r is defined by the mapping (cf. class \c Map) and 
 * is described in Bonazzola, Gourgoulhon \& Marck, \a Phys. \a Rev. \a D
 * \b 58, 104020 (1998).
 *
 * @version #$Id$#
 */

class Grille3d {
    protected:
	const int nr ;	///< Number of points in \e r (\f$\xi\f$)
	const int nt ;	///< Number of points in \f$\theta\f$
	const int np ;	///< Number of points in \f$\phi\f$

	int type_r ;	///< Type of sampling in \e r (\f$\xi\f$) (\c RARE, FIN, UNSURR)
	int type_t ;	///< Type of sampling in \f$\theta\f$ (\c SYM, NONSYM)
	int type_p ;	///< Type of sampling in \f$\phi\f$ (\c SYM, NONSYM)
    public:
	/// Array of values of \f$\xi\f$ at the \c nr collocation points
	double* x ;	
	/// Array of values of \f$\theta\f$ at the \c nt collocation points
	double* tet ;	
	/// Array of values of \f$\phi\f$ at the \c np collocation points
	double* phi ;	

    protected: 
	/// Constructor (protected to make \c Grille3d an abstract class) 
	Grille3d(int n_r, int n_t, int n_p) ;
    
    private:
	/** Copy constructor (private and not implemented to make \c Grille3d
	 * a non-copyable class)
	 */ 
	Grille3d(const Grille3d& ) ;
	
	/** Assignement operator (private and not implemented to make 
	 *   \c Grille3d a non-copyable class)
	 */
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

};


		    // -------------------- //
    	    	    // Les classes derivees //
		    // -------------------- //

/**@name Derived classes of class Grille3d.\ingroup (spec)
 * These derived classes differ only by their constructors 
 * (see the corresponding constructors for details).
 * 
 * Note: the monogrids should not be used. Use instead \c Mg3d.
 */
//@{

// Cas rare + sans symetrie
// ------------------------
/**
 * 3D grid for a spherical kernel without symmetry.
 * It contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [0,1]\f$, \f$\theta \in [0,\pi]\f$ and 
 * \f$\phi \in [0,2\pi[\f$.
 *
 */
class Grille3d_r : public Grille3d {
    public:
	Grille3d_r(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_r() ; ///< Destructor
};

// cas fin + sans symetrie
// -----------------------
/**
 * 3D grid for a spherical shell without symmetry.
 * It contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [-1,1]\f$, \f$\theta \in [0,\pi]\f$ and 
 * \f$\phi \in [0,2\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_f : public Grille3d {
    public:
	Grille3d_f(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_f() ; ///< Destructor
};

// cas echantillonnage (fin) en 1/r + sans symetrie
// ------------------------------------------------
/**
 * 3D grid for a compatified spherical shell without symmetry.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [-1,1]\f$, \f$\theta \in [0,\pi]\f$ and 
 * \f$\phi \in [0,2\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_i : public Grille3d {
    public:
	Grille3d_i(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_i() ; ///< Destructor
};

// Cas rare + symetrie equatoriale
// -------------------------------
/**
 * 3D grid for a spherical kernel with equatorial symmetry.
 * \f$z \rightarrow -z\f$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [0,1]\f$, \f$\theta \in [0,\pi/2]\f$ and 
 * \f$\phi \in [0,2\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_req : public Grille3d {
    public:
	Grille3d_req(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_req() ; ///< Destructor
};

// cas fin + symetrie equatoriale
// ------------------------------
/**
 * 3D grid for a spherical shell with equatorial symmetry 
 * \f$z \rightarrow -z\f$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [-1,1]\f$, \f$\theta \in [0,\pi/2]\f$ and 
 * \f$\phi \in [0,2\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_feq : public Grille3d {
    public:
	Grille3d_feq(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_feq() ; ///< Destructor
};

// cas echantillonnage (fin) en 1/r + symetrie equatoriale
// -------------------------------------------------------
/**
 * 3D grid for a compatified spherical shell with equatorial symmetry
 * \f$z \rightarrow -z\f$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [-1,1]\f$, \f$\theta \in [0,\pi/2]\f$ and 
 * \f$\phi \in [0,2\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_ieq : public Grille3d {
    public:
	Grille3d_ieq(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_ieq() ; ///< Destructor
};

// Cas rare supersymetrique
// ------------------------
/**
 * 3D grid for a spherical kernel with super-symmetry 
 * \f$(x,y) \rightarrow (-x,-y)\f$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [0,1]\f$, \f$\theta \in [0,\pi/2]\f$ and 
 * \f$\phi \in [0,\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_rs : public Grille3d {
    public:
	Grille3d_rs(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_rs() ; ///< Destructor
};

// cas fin supersymetrique
// -----------------------
/**
 * 3D grid for a spherical shell with super-symmetry 
 * \f$(x,y) \rightarrow (-x,-y)\f$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [-1,1]\f$, \f$\theta \in [0,\pi/2]\f$ and 
 * \f$\phi \in [0,\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_fs : public Grille3d {
    public:
	Grille3d_fs(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_fs() ; ///< Destructor
};

// cas echantillonnage (fin) en 1/r supersymetrique
// --------------------------
/**
 * 3D grid for a compactified spherical shell with super-symmetry 
 * \f$(x,y) \rightarrow (-x,-y)\f$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [-1,1]\f$, \f$\theta \in [0,\pi/2]\f$ and 
 * \f$\phi \in [0,\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_is : public Grille3d {
    public:
	Grille3d_is(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_is() ; ///< Destructor
};

// Cas rare + 2 phi
// ----------------
/**
 * 3D grid for a spherical kernel with only even harmonics in \f$\phi\f$. 
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [0,1]\f$, \f$\theta \in [0,\pi]\f$ and 
 * \f$\phi \in [0,\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_r2p : public Grille3d {
    public:
	Grille3d_r2p(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_r2p() ; ///< Destructor
};

// cas fin + 2 phi
// ---------------
/**
 * 3D grid for a spherical shell with only even harmonics in \f$\phi\f$. 
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [-1,1]\f$, \f$\theta \in [0,\pi]\f$ and 
 * \f$\phi \in [0,\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_f2p : public Grille3d {
    public:
	Grille3d_f2p(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_f2p() ; ///< Destructor
};

// cas echantillonnage (fin) en 1/r + 2 phi
// ----------------------------------------
/**
 * 3D grid for a compactified spherical kernel 
 * with only even harmonics in \f$\phi\f$. 
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: \f$\xi \in [-1,1]\f$, \f$\theta \in [0,\pi]\f$ and 
 * \f$\phi \in [0,\pi[\f$.
 *
 * @version #$Id$#
 */
class Grille3d_i2p : public Grille3d {
    public:
	Grille3d_i2p(int n_r, int n_t, int n_p) ; ///< Constructor
	~Grille3d_i2p() ; ///< Destructor
};

//@}

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
    private:
	int nzone ;	///< Number of domains (zones)
	
	int* nr ;	///< Array (size: \c nzone) of nb. of points in \e r (\f$\xi\f$)
	int* nt ;	///< Array (size: \c nzone) of nb. of points in \f$\theta\f$
	int* np ;	///< Array (size: \c nzone) of nb. of points in \f$\phi\f$
	
	/** Array (size: \c nzone) of type of sampling in \e r (\f$\xi\f$) 
     *(\c RARE,\c FIN, \c UNSURR)
     */
	int* type_r ;	
	/// Type of sampling in \f$\theta\f$ (\c SYM, \c NONSYM)
	int type_t ;
	/// Type of sampling in \f$\phi\f$ (\c SYM, \c NONSYM)	
	int type_p ;	
	
	/// Array (size: \c nzone) of pointers on the \c Grille3d's
	Grille3d** g ;	

	mutable Mg3d* g_angu ;	///< Pointer on the associated angular grid
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
 *
 */
	Mg3d(int nz, int nbr[], int typr[], int nbt[], int typt, int nbp[],
	     int typp) ;

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
 */
	Mg3d(int nz, int nbr, int nbt, int nbp, int typt, int typp, 
	     bool compact) ;

	Mg3d(FILE* ) ;	 ///< Constructor from a file (see \c sauve(FILE*))
	
   private:
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

	/// Returns the pointer on the associated angular grid
	const Mg3d* get_angu() const ;
	
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
	void sauve(FILE* ) const ; ///< Save in a file
	
	friend ostream& operator<<(ostream& , const Mg3d & ) ;	///< Display
	
    // Management of derived quantities
    // --------------------------------
    protected:
	/** Deletes all the derived quantities 
	 *   (\c g_radial , \c g_angu  and \c g_twice )
	 */
	void del_deriv() const ; 
	
	/** Sets to \c 0x0  all the pointers on derived quantities
	 *   (\c g_radial , \c g_angu  and \c g_twice )
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

};
ostream& operator<<(ostream& , const Mg3d & ) ;

#endif 
