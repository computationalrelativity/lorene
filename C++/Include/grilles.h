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

class Base_val ; 

		    	//-------------//
		    	// Mono-grille //
		    	//-------------//

// Classe de base
/**
 * Base 3D grid class.
 * 
 * This class is an abstract one. It can't be instanciated. 
 * There is no {\tt ostream} 
 * operator. For historical reasons, the mnemonics refer to spherical 
 * coordinates $(r,\theta,\phi)$. However all this stuff can be used for 
 * any coordinate system if the appropriated constructors are created.
 *
 * The radial coordinate $\xi$ lies in the range [0, 1] or [-1, 1]
 * depending upon the sampling. Its relation with the physical radial 
 * coordinate {\it r} is defined by the mapping (cf. class {\tt Map}) and 
 * is described in Bonazzola, Gourgoulhon \& Marck, {\sl Phys. Rev. D}
 * {\bf 58}, 104020 (1998).
 *
 * @version #$Id$#
 */

class Grille3d {
    protected:
	const int nr ;	/// Number of points in {\it r} ($\xi$)
	const int nt ;	/// Number of points in $\theta$
	const int np ;	/// Number of points in $\phi$

	int type_r ;	/// Type of sampling in {\it r} ($\xi$) ({\tt RARE, FIN, UNSURR})
	int type_t ;	/// Type of sampling in $\theta$ ({\tt SYM, NONSYM})
	int type_p ;	/// Type of sampling in $\phi$ ({\tt SYM, NONSYM})
    public:
	/// Array of values of $\xi$ at the {\tt nr} collocation points
	double* x ;	
	/// Array of values of $\theta$ at the {\tt nt} collocation points
	double* tet ;	
	/// Array of values of $\phi$ at the {\tt np} collocation points
	double* phi ;	

    protected: 
	/// Constructor (protected to make {\tt Grille3d} an abstract class) 
	Grille3d(int n_r, int n_t, int n_p) ;
    
    private:
	/** Copy constructor (private and not implemented to make {\tt Grille3d}
	 * a non-copyable class)
	 */ 
	Grille3d(const Grille3d& ) ;
	
	/** Assignement operator (private and not implemented to make 
	 *   {\tt Grille3d} a non-copyable class)
	 */
	void operator=(const Grille3d& ) ;
	 	
    public:
	virtual ~Grille3d() ;		/// Destructor

    public:
    	int get_nr() const {return nr ;} ;	    /// Returns {\tt nr}
    	int get_nt() const {return nt ;} ;	    /// Returns {\tt nt}
    	int get_np() const {return np ;} ;	    /// Returns {\tt np}

    	int get_type_r() const {return type_r ;} ;    /// Returns {\tt type\_r}
    	int get_type_t() const {return type_t ;} ;    /// Returns {\tt type\_t}
    	int get_type_p() const {return type_p ;} ;    /// Returns {\tt type\_p}

};


		    // -------------------- //
    	    	    // Les classes derivees //
		    // -------------------- //

/**@name Derived classes of class Grille3d.
 * These derived classes differ only by their constructors 
 * (see the corresponding constructors for details).
 * 
 * Note: the monogrids should not be used. Use instead {\tt Mg3d}.
 */
//@{

// Cas rare + sans symetrie
// ------------------------
/**
 * 3D grid for a spherical kernel without symmetry.
 * It contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [0,1]$, $\theta \in [0,\pi]$ and 
 * $\phi \in [0,2\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_r : public Grille3d {
    public:
	Grille3d_r(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_r() ; /// Destructor
};

// cas fin + sans symetrie
// -----------------------
/**
 * 3D grid for a spherical shell without symmetry.
 * It contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [-1,1]$, $\theta \in [0,\pi]$ and 
 * $\phi \in [0,2\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_f : public Grille3d {
    public:
	Grille3d_f(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_f() ; /// Destructor
};

// cas echantillonnage (fin) en 1/r + sans symetrie
// ------------------------------------------------
/**
 * 3D grid for a compatified spherical shell without symmetry.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [-1,1]$, $\theta \in [0,\pi]$ and 
 * $\phi \in [0,2\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_i : public Grille3d {
    public:
	Grille3d_i(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_i() ; /// Destructor
};

// Cas rare + symetrie equatoriale
// -------------------------------
/**
 * 3D grid for a spherical kernel with equatorial symmetry.
 * $z \rightarrow -z$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [0,1]$, $\theta \in [0,\pi/2]$ and 
 * $\phi \in [0,2\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_req : public Grille3d {
    public:
	Grille3d_req(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_req() ; /// Destructor
};

// cas fin + symetrie equatoriale
// ------------------------------
/**
 * 3D grid for a spherical shell with equatorial symmetry 
 * $z \rightarrow -z$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [-1,1]$, $\theta \in [0,\pi/2]$ and 
 * $\phi \in [0,2\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_feq : public Grille3d {
    public:
	Grille3d_feq(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_feq() ; /// Destructor
};

// cas echantillonnage (fin) en 1/r + symetrie equatoriale
// -------------------------------------------------------
/**
 * 3D grid for a compatified spherical shell with equatorial symmetry
 * $z \rightarrow -z$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [-1,1]$, $\theta \in [0,\pi/2]$ and 
 * $\phi \in [0,2\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_ieq : public Grille3d {
    public:
	Grille3d_ieq(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_ieq() ; /// Destructor
};

// Cas rare supersymetrique
// ------------------------
/**
 * 3D grid for a spherical kernel with super-symmetry 
 * $(x,y) \rightarrow (-x,-y)$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [0,1]$, $\theta \in [0,\pi/2]$ and 
 * $\phi \in [0,\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_rs : public Grille3d {
    public:
	Grille3d_rs(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_rs() ; /// Destructor
};

// cas fin supersymetrique
// -----------------------
/**
 * 3D grid for a spherical shell with super-symmetry 
 * $(x,y) \rightarrow (-x,-y)$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [-1,1]$, $\theta \in [0,\pi/2]$ and 
 * $\phi \in [0,\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_fs : public Grille3d {
    public:
	Grille3d_fs(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_fs() ; /// Destructor
};

// cas echantillonnage (fin) en 1/r supersymetrique
// --------------------------
/**
 * 3D grid for a compactified spherical shell with super-symmetry 
 * $(x,y) \rightarrow (-x,-y)$.
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [-1,1]$, $\theta \in [0,\pi/2]$ and 
 * $\phi \in [0,\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_is : public Grille3d {
    public:
	Grille3d_is(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_is() ; /// Destructor
};

// Cas rare + 2 phi
// ----------------
/**
 * 3D grid for a spherical kernel with only even harmonics in $\phi$. 
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [0,1]$, $\theta \in [0,\pi]$ and 
 * $\phi \in [0,\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_r2p : public Grille3d {
    public:
	Grille3d_r2p(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_r2p() ; /// Destructor
};

// cas fin + 2 phi
// ---------------
/**
 * 3D grid for a spherical shell with only even harmonics in $\phi$. 
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [-1,1]$, $\theta \in [0,\pi]$ and 
 * $\phi \in [0,\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_f2p : public Grille3d {
    public:
	Grille3d_f2p(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_f2p() ; /// Destructor
};

// cas echantillonnage (fin) en 1/r + 2 phi
// ----------------------------------------
/**
 * 3D grid for a compactified spherical kernel 
 * with only even harmonics in $\phi$. 
 * This class contains only a constructor and a destructor. 
 * Coordinates ranges are: $\xi \in [-1,1]$, $\theta \in [0,\pi]$ and 
 * $\phi \in [0,\pi[$.
 *
 * @version #$Id$#
 */
class Grille3d_i2p : public Grille3d {
    public:
	Grille3d_i2p(int n_r, int n_t, int n_p) ; /// Constructor
	~Grille3d_i2p() ; /// Destructor
};

//@}

		    	//---------------//
		    	// Multi-grilles //
		    	//---------------//

/**
 * Multi-domain grid.
 *
 * A multi-domain grid is a set of 3D grids for the implementation of 
 * the multi-domain spectral method described in Bonazzola, Gourgoulhon 
 * \& Marck, {\sl Phys. Rev. D} {\bf 58}, 104020 (1998).
 * Each domain is represented by a 3D mono grid ({\tt Grille3d}) and
 * is called a {\em zone}. 
 * For each direction, the Number of Degrees of Freedom (NDF) is 
 * {\em a priori} independent of the zone. However, some methods or 
 * routines may refuse to work if the NDF of some $(\theta, \phi)$ 
 * direction is not identical in all the zones. 
 * This holds for the type of sampling (symmetry) too.
 *
 * @version #$Id$#
 */
 
class Mg3d {

    // Data  
    // ----
    private:
	int nzone ;	/// Number of zones
	
	int* nr ;	/// Array (size: {\tt nzone}) of nb. of points in {\it r} ($\xi$)
	int* nt ;	/// Array (size: {\tt nzone}) of nb. of points in $\theta$
	int* np ;	/// Array (size: {\tt nzone}) of nb. of points in $\phi$
	
	/// Array (size: {\tt nzone}) of type of sampling in {\it r} ($\xi$) ({\tt RARE, FIN, UNSURR})
	int* type_r ;	
	/// Type of sampling in $\theta$ ({\tt SYM, NONSYM})
	int type_t ;
	/// Type of sampling in $\phi$ ({\tt SYM, NONSYM})	
	int type_p ;	
	
	/// Array (size: {\tt nzone}) of pointers on the {\tt Grille3d}'s
	Grille3d** g ;	

	mutable Mg3d* g_angu ;	/// Pointer on the associated angular grid
	mutable Mg3d* g_radial ; /// Pointer on the associated radial grid
	
	/** Pointer on the grid which has twice the number of points in
	 *  each dimension (for desaliasing).
	 */
	mutable Mg3d* g_twice ; 

    // Constructors - Destructor
    // -------------------------
	
    public:

/**
 * General constructor.
 *
 * @param   nz	    [input] Number of domains (zones).
 * @param   nbr[]   [input] Array (size: {\tt nz}) of number of degree of
 *		    freedom (NDF) in {\it r}-direction
 * @param   typr[]  [input] Array (size: {\tt nz}) of type of sampling in {\it r}-direction
 * @param   nbt[]   [input] Array (size: {\tt nz}) of NDF in $\theta$-direction
 * @param   typt    [input] Type of sampling in $\theta$-direction
 * @param   nbp[]   [input] Array (size: {\tt nz}) of NDF in $\phi$-direction
 * @param   typp    [input] Type of sampling in $\phi$-direction
 *
 */
	Mg3d(int nz, int nbr[], int typr[], int nbt[], int typt, int nbp[],
	     int typp) ;

/**
 * Simplified constructor for a standard multi-grid.
 * This provides a multi-grid with the same number of degrees of freedom
 * in all the domains. \\
 * The domain of index {\tt l = 0} will be a nucleus:
 * $\xi\in [0,1]$, rarefied sampling (type {\tt RARE}) near the origin; \\
 * domains of indices $ 1 \le {\tt l} \le {\tt nz}-2$ will be shells:
 * $\xi\in [-1,1]$, dense sampling (type {\tt FIN}) near -1 and 1; \\
 * if {\tt compact == true}, 
 * the domains of index {\tt l = nz-1} will be the outermost compactified
 * shell:
 * $\xi\in [-1,1]$, dense sampling (type {\tt UNSURR}) near -1 and 1
 * for a {\it 1/r} discretization.
 *
 *
 * @param   nz	    [input] Number of domains (zones).
 * @param   nbr     [input] Number of degree of freedom (NDF) in
 *				{\it r}-direction in each domain
 * @param   nbt     [input] Number of degree of freedom (NDF) in
 *				$\theta$-direction in each domain
 * @param   nbp     [input] Number of degree of freedom (NDF) in
 *				$\phi$-direction in each domain
 * @param   typt    [input] Type of sampling in $\theta$-direction:  \\
 *				{\tt SYM} for a sampling in $[0,\pi/2]$
 *			(symmetry with respect to the equatorial plane), \\
 *			{\tt NONSYM} for a sampling in $[0,\pi]$
 * @param   typp    [input] Type of sampling in $\phi$-direction: \\
 *			{\tt SYM} for a sampling in $[0,\pi[$
 *			(symmetry with respect to a $\pi$ translation
 *			 in $\phi$) \\
 *			{\tt NONSYM} for a sampling in $[0,2\pi[$
 * @param  compact [input] {\tt true} for the last domain to have 
 *			a {\it 1/r} sampling ({\tt UNSURR}) instead of a
 *			{\it r} sampling ({\tt FIN}). 
 */
	Mg3d(int nz, int nbr, int nbt, int nbp, int typt, int typp, 
	     bool compact) ;

	Mg3d(FILE* ) ;	 /// Constructor from a file (see {\tt sauve(FILE* )})
	
   private:
	/**
	 * Copy constructor (private and not implemented to make {\tt Mg3d} a 
	 * non-copyable class)
	 */
	Mg3d(const Mg3d& ) ;
	
    public:
	
	~Mg3d() ;   /// Destructor
		
    // Assignement
    // -----------
    private:
	/** 
	 * Assignement operator (private and not implemented to make {\tt Mg3d}
	 * a non-copyable class)
	 */
	void operator=(const Mg3d& ) ;
	 	
    // Extraction of information
    // -------------------------
    public:
	int get_nzone() const { 	    	/// Returns {\tt nzone}
	    return nzone ;
	} ;

	int get_nr(int l) const { 	    	/// Returns {\tt nr[l]}
	    assert(l>=0 && l<nzone) ;
	    return nr[l] ;
	} ;
	
	int get_nt(int l) const { 	    	/// Returns {\tt nt[l]}
	    assert(l>=0 && l<nzone) ;
	    return nt[l] ;
	} ;
	
	int get_np(int l) const { 	    	/// Returns {\tt np[l]}
	    assert(l>=0 && l<nzone) ;
	    return np[l] ;
	} ;
	
	int get_type_r(int l) const { 	/// Returns {\tt type\_r[l]}
	    assert(l>=0 && l<nzone) ;
	    return type_r[l] ;
	} ;
	
	int get_type_t() const { 	    	/// Returns {\tt type\_t}
	    return type_t ;
	} ;
	
	int get_type_p() const {		/// Returns {\tt type\_p}
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

	/// Comparison operator (egality)
	bool operator==(const Mg3d& ) const ;  


	
    // Outputs
    // -------
    public: 
	void sauve(FILE* ) const ; /// Save in a file
	
	friend ostream& operator<<(ostream& , const Mg3d & ) ;	/// Display
	
    // Management of derived quantities
    // --------------------------------
    protected:
	/** Deletes all the derived quantities 
	 *   ({\tt g\_radial}, {\tt g\_angu} and {\tt g\_twice})
	 */
	void del_deriv() const ; 
	
	/** Sets to {\tt 0x0} all the pointers on derived quantities
	 *   ({\tt g\_radial}, {\tt g\_angu} and {\tt g\_twice})
	 */
	void set_deriv_0x0() const ; 


    // Miscellaneous
    // -------------
    public:
	bool operator!=(const Mg3d & ) const ;  /// Operator !=

	/// Returns the standard spectral bases for a scalar
	Base_val std_base_scal() const ; 	

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
