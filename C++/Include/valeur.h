/*
 *  Definition of Lorene class Valeur
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 1999-2001 Jerome Novak
 *   Copyright (c) 2001 Keisuke Taniguchi
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


#ifndef __VALEUR_H_ 
#define __VALEUR_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.12  2004/11/23 12:46:04  f_limousin
 * Add functiun filtre_tp(int nn, int nz1, int nz2).
 *
 * Revision 1.11  2004/08/24 09:14:40  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.10  2004/07/06 13:36:27  j_novak
 * Added methods for desaliased product (operator |) only in r direction.
 *
 * Revision 1.9  2004/03/22 13:12:44  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.8  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.7  2003/10/19 19:48:31  e_gourgoulhon
 * Introduced new method display_coef.
 *
 * Revision 1.6  2003/10/13 20:48:52  e_gourgoulhon
 * Added new method get_base() (to prepare the encapsulation of the member
 * base).
 *
 * Revision 1.5  2003/09/23 08:52:53  e_gourgoulhon
 * Added Scalar as a friend class.
 *
 * Revision 1.4  2002/10/16 14:36:30  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
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
 * Revision 2.51  2001/05/29  16:11:58  eric
 * Modif commentaires (mise en conformite Doc++ 3.4.7).
 *
 * Revision 2.50  2001/05/26  14:49:37  eric
 * Ajout de l'operator% : produit de deux Valeur avec desaliasage.
 *
 * Revision 2.49  2001/01/16  15:09:13  keisuke
 * Change the argument of the function smooth.
 *
 * Revision 2.48  2001/01/16  14:54:01  keisuke
 * Ajout de la fonction smooth.
 *
 * Revision 2.47  2000/11/10  13:31:53  eric
 * Ajout de la fonction equipot_outward.
 *
 * Revision 2.46  2000/09/11  13:52:39  eric
 * Ajout des methodes mult_cp() et mult_sp() ainsi que des membres associes
 *  p_mult_cp et p_mult_sp
 *
 * Revision 2.45  2000/09/08  11:43:26  eric
 * Modif commentaires.
 *
 * Revision 2.44  2000/09/08  10:07:02  eric
 * Ajout des methodes set_base_r, etc...
 *
 * Revision 2.43  2000/08/04  11:53:59  eric
 * Ajout de l'operateur (int l) et de la fonction set(int l) pour l'acces
 * individuel aux Tbl.
 *
 * Revision 2.42  1999/12/29  13:19:51  eric
 * Modif commentaires.
 *
 * Revision 2.41  1999/12/29  13:11:08  eric
 * Ajout de la fonction val_point_jk.
 *
 * Revision 2.40  1999/12/20  16:35:47  eric
 * Ajout de la fonction set_base.
 *
 * Revision 2.39  1999/12/10  16:19:40  eric
 * Modif commentaires.
 *
 * Revision 2.38  1999/12/10  16:11:13  eric
 * Fonction set: suppression de l'appel a set_etat_c_qcq() pour
 *  augmenter l'efficacite.
 *
 * Revision 2.37  1999/12/07  14:52:43  eric
 * Changement ordre des arguments (phi,theta,xi) --> (xi,theta,phi)
 *  dans la routine val_point.
 *
 * Revision 2.36  1999/12/06  16:46:46  eric
 * Ajout de la fonction val_point.
 *
 * Revision 2.35  1999/11/30  12:41:21  eric
 * Le membre base est desormais un objet de type Base_val et non plus
 *  un pointeur vers une Base_val.
 *
 * Revision 2.34  1999/11/29  10:25:32  eric
 * Ajout de Valeur/Mtbl et Mtbl / Valeur dans l'arithmetique.
 *
 * Revision 2.33  1999/11/29  10:05:47  eric
 * Ajout de Valeur*Mtbl dans l'arithmetique.
 *
 * Revision 2.32  1999/11/23  16:15:24  eric
 * Suppression du membre statique Valeur_Zero.
 * Suppression du constructeur par defaut.
 *
 * Revision 2.31  1999/11/23  14:30:47  novak
 * Ajout des membres mult_ct et mult_st
 *
 * Revision 2.30  1999/11/22  15:40:48  eric
 * Ajout des operateurs set(l,k,j,i) et (l,k,j,i).
 * Ajout de la fonction annule(int l).
 *
 * Revision 2.29  1999/11/19  11:21:48  eric
 * Ajout du membre p_stdsdp et de la fonction correspondante stdsdp().
 *
 * Revision 2.28  1999/11/19  09:28:38  eric
 * Les valeurs de retour des operateurs differentiels sont desormais
 *   const Valeur &
 * et non plus Valeur.
 * Le laplacien angulaire (lapang) a desormais le meme statut que les
 * autres operateurs differentiels.
 *
 * Revision 2.27  1999/11/16  13:09:48  novak
 * Ajout de mult_x et scost
 *
 * Revision 2.26  1999/11/09  15:24:08  phil
 * ajout de la fonction mathematique calculant la racine cubique
 *
 * Revision 2.25  1999/10/29  15:14:33  eric
 * Ajout de fonctions mathematiques (abs, norme, etc...).
 *
 * Revision 2.24  1999/10/27  08:48:46  eric
 * La classe Cmp est desormais amie (pour que Cmp::del_t() puisse appeler
 * Valeur::del_t()).
 *
 * Revision 2.23  1999/10/21  14:21:06  eric
 * Constructeur par lecture de fichier.
 *
 * Revision 2.22  1999/10/20  15:38:52  eric
 * *** empty log message ***
 *
 * Revision 2.21  1999/10/20  15:31:04  eric
 * Ajout de l'arithmetique.
 *
 * Revision 2.20  1999/10/19  15:30:24  eric
 * Ajout de la fonction affiche_seuil.
 *
 * Revision 2.19  1999/10/18  15:07:02  eric
 *  La fonction membre annule() est rebaptisee annule_hard().
 * Introduction de la fonction membre annule(int, int).
 *
 * Revision 2.18  1999/10/18  13:39:38  eric
 * Routines de derivation --> const
 * Suppression de sxdsdx (non implemente).
 *
 * Revision 2.17  1999/10/13  15:49:57  eric
 * Depoussierage.
 * Documentation.
 *
 * Revision 2.16  1999/09/14  17:17:47  phil
 * *** empty log message ***
 *
 * Revision 2.15  1999/09/14  17:15:38  phil
 * ajout de Valeur operator* (double, const Valeur&)
 *
 * Revision 2.14  1999/09/13  14:53:26  phil
 * *** empty log message ***
 *
 * Revision 2.13  1999/09/13  14:17:52  phil
 * ajout de Valeur friend operator+ (Valeur, Valeur)
 *
 * Revision 2.12  1999/04/26  16:24:23  phil
 * ajout de mult2_xm1_zec()
 *
 * Revision 2.11  1999/04/26  16:12:45  phil
 * ajout de mult_xm1_zec()
 *
 * Revision 2.10  1999/04/26  15:48:11  phil
 * ajout de sxm1_zec()
 *
 * Revision 2.9  1999/04/26  12:57:24  phil
 * ajout de lapang()
 *
 * Revision 2.8  1999/04/13  16:44:55  phil
 * ajout de ylm_i()
 *
 * Revision 2.7  1999/04/13  16:31:46  phil
 * *** empty log message ***
 *
 * Revision 2.6  1999/04/13  16:26:08  phil
 * ajout ylm
 *
 * Revision 2.5  1999/02/24  15:24:34  hyc
 * *** empty log message ***
 *
 * Revision 2.4  1999/02/23  15:55:46  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Fichier includes
#include <stdio.h>

#include "mtbl.h"
#include "mtbl_cf.h"

class Coord ;
class Itbl ; 

/**
 * Values and coefficients of a (real-value) function.
 * \ingroup (spec)
 *
 */

class Valeur {

    // Data : 
    // -----
    private:
	const Mg3d* mg ;  ///< Multi-grid \c Mgd3  on which \c this  is defined

	/// Logical state (\c ETATNONDEF , \c ETATQCQ  or \c ETATZERO ).
	int etat ;	

    public:
	/// Values of the function at the points of the multi-grid	
	mutable Mtbl* c ;	    

	/// Coefficients of the spectral expansion of the function
	mutable Mtbl_cf* c_cf ;	    

	/// Bases on which the spectral expansion is performed
	Base_val base ;	   

    // Derived data : 
    // ------------
    private:
	mutable Valeur* p_dsdx ;    ///< Pointer on \f$\partial / \partial \xi\f$ 
	mutable Valeur* p_d2sdx2 ;  ///< Pointer on \f$\partial^2 / \partial \xi^2\f$
	mutable Valeur* p_sx ;	    ///< Pointer on \f$1 / \xi\f$
	mutable Valeur* p_sx2 ;	    ///< Pointer on \f$1 / \xi^2\f$
	mutable Valeur* p_mult_x ;  ///< Pointer on \f$\xi \, Id\f$

	mutable Valeur* p_dsdt ;    ///< Pointer on \f$\partial / \partial \theta\f$
	mutable Valeur* p_d2sdt2 ;  ///< Pointer on \f$\partial^2 / \partial \theta^2\f$
	mutable Valeur* p_ssint ;   ///< Pointer on \f$1 / \sin(\theta)\f$
        mutable Valeur* p_scost ;   ///< Pointer on \f$1 / \cos(\theta)\f$	  
	mutable Valeur* p_mult_ct ; ///< Pointer on \f$\cos(\theta) \, Id\f$
	mutable Valeur* p_mult_st ; ///< Pointer on \f$\sin(\theta) \, Id\f$

	mutable Valeur* p_dsdp ;    ///< Pointer on \f$\partial / \partial \phi\f$
	mutable Valeur* p_stdsdp ;  ///< Pointer on \f$1/\sin\theta \partial / \partial \phi\f$
	mutable Valeur* p_d2sdp2 ;  ///< Pointer on \f$\partial^2 / \partial \phi^2\f$
	mutable Valeur* p_mult_cp ; ///< Pointer on \f$\cos(\phi) \, Id\f$
	mutable Valeur* p_mult_sp ; ///< Pointer on \f$\sin(\phi) \, Id\f$

	mutable Valeur* p_lapang ;  ///< Pointer on the angular Laplacian

    // Constructors - Destructor
    // -------------------------
	
    public:
	explicit Valeur(const Mg3d& mgrid) ;		///< Constructor
	explicit Valeur(const Mg3d* p_mgrid) ;		///< Constructor

	/// Constructor from a file (see \c sauve(FILE*) )
	Valeur(const Mg3d&, FILE* ) ;    		

	Valeur(const Valeur& ) ;	///< Copy constructor
	~Valeur() ;			///< Destructor

    // Assignement
    // -----------
    public: 
	void operator=(const Valeur& a) ; ///< Assignement to another \c Valeur 
	void operator=(const Mtbl& mt) ;	 ///< Assignement to a \c Mtbl 
	void operator=(const Mtbl_cf& mtcf) ; ///< Assignement to a \c Mtbl_cf 
	void operator=(double ) ;	///< Assignement to a \c double 	    

    // Access to individual elements
    // -----------------------------
    public:
	/** Read/write of the value in a given domain (configuration space).
	 * NB: to gain in efficiency, the method \c del_deriv()  (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *
	 * @param l [input] domain index
	 * @return  Tbl containing the value of the field in domain \c l .
	 */ 
	Tbl& set(int l) {
	    assert(l < mg->get_nzone()) ;
	    assert(etat == ETATQCQ) ;
	    if (c == 0x0) {
		coef_i() ;
	    }
	    if (c_cf != 0x0) {
		delete c_cf ; 
		c_cf = 0 ; 
	    }
	    return c->set(l) ;
	};
	
	
	/** Read-only of the value in a given domain (configuration space).
	 * @param l [input] domain index
	 * @return  Tbl containing the value of the field in domain \c l .
	 */ 
	const Tbl& operator()(int l) const {
	    assert(l < mg->get_nzone()) ;
	    assert(etat == ETATQCQ) ;
	    if (c == 0x0) {
		coef_i() ;
	    }
	    return (*c)(l) ;
	};


	/** Read/write of a particular element (configuration space).
	 * NB: to gain in efficiency, the method \c del_deriv()  (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *
	 * @param l [input] domain index
	 * @param k [input] \f$\phi\f$ index
	 * @param j [input] \f$\theta\f$ index
	 * @param i [input] \e r  (\f$\xi\f$) index
	 */ 
	double& set(int l, int k, int j, int i) {
	    assert(l < mg->get_nzone()) ;
	    assert(etat == ETATQCQ) ;
	    if (c == 0x0) {
		coef_i() ;
	    }
	    if (c_cf != 0x0) {
		delete c_cf ; 
		c_cf = 0 ; 
	    }
	    return c->set(l, k, j, i) ;
	};
	
	
	/** Read-only of a particular element (configuration space).
	 * @param l [input] domain index
	 * @param k [input] \f$\phi\f$ index
	 * @param j [input] \f$\theta\f$ index
	 * @param i [input] \e r  (\f$\xi\f$) index
	 */ 
	double operator()(int l, int k, int j, int i) const {
	    assert(etat != ETATNONDEF) ;
	    assert(l < mg->get_nzone()) ;
	    if (etat == ETATZERO) {
		double zero = 0. ;
		return zero ; 
	    }
	    else{ 	    
		if (c == 0x0) {
		    coef_i() ;
		}
	    	return (*c)(l, k, j, i) ;
	    }
	};

	/** Computes the value of the field represented by \c *this  at an
	*   arbitrary point, by means of the spectral expansion.
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate \f$\xi\f$
	*	 @param theta [input] value of the coordinate \f$\theta'\f$
	*	 @param phi [input] value of the coordinate \f$\phi'\f$
	*	 @return value at the point \f$(\xi, \theta', \phi')\f$ in
	*	    the domain no. \e l  of the field represented by \c *this . 
	*/
	double val_point(int l, double x, double theta, double phi) const ; 

	/** Computes the value of the field represented by \c *this  at an
	*   arbitrary point in \f$\xi\f$, but collocation point in 
	*   \f$(\theta', \phi')\f$, by means of the spectral expansion.
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate \f$\xi\f$
	*	 @param j [input] index of the collocation point in \f$\theta'\f$
	*	 @param k [input] index of the collocation point in \f$\phi'\f$
	*	 @return value at the point 
	*		    \f$(\xi, {\theta'}_j, {\phi'}_k)\f$ in
	*	    the domain no. \e l  of the field represented by \c *this . 
	*/
	double val_point_jk(int l, double x, int j, int k) const ; 

   
    // Operations on coefficients
    // --------------------------
    public:
	void coef() const ;	///< Computes the coeffcients of \c *this 
	void coef_i() const ;	///< Computes the physical value of \c *this 
	void ylm() ;	    ///< Computes the coefficients \f$Y_l^m\f$ of \c *this 
	void ylm_i() ;	    ///< Inverse of \c ylm()  

	/**
	 * Set the basis to the eigenvalues of 
	 * \f$ \partial^2_theta - \cos\theta/\sin\theta \partial_\theta\f$.
	 **/
	void val_propre_1d() ;
	/**
	 * Inverse transformation of \c val_propre_1d.
	 **/
	void val_propre_1d_i() ;
	
	/// Return the bases for spectral expansions (member \c base )
	const Base_val& get_base() const {return base; } ; 

	/// Sets the bases for spectral expansions (member \c base ) 
	void set_base(const Base_val& ) ; 

	/** Sets the bases for spectral expansions (member \c base ) 
	 *  to the standard ones for a scalar
	 */
	void std_base_scal() ;	 

	/** Sets the expansion basis for \e r  (\f$\xi\f$) functions in a 
	 *  given domain.
	 *  
	 *  @param l	    Domain index
	 *  @param base_r   type of basis functions in \e r  (\f$\xi\f$)
	 *		    (e.g. \c R_CHEBP , etc..., 
	 *		     see documentation of class \c Base_val 
	 *		     for denomination of the various bases). 
	 */
	void set_base_r(int l, int base_r) ; 

	/** Sets the expansion basis for \f$\theta\f$ functions in all
	 *  domains.
	 *  
	 *  @param base_t   type of basis functions in \f$\theta\f$
	 *		    (e.g. \c T_COS_P , etc..., 
	 *		     see documentation of class \c Base_val 
	 *		     for denomination of the various bases). 
	 */
	void set_base_t(int base_t) ; 

	/** Sets the expansion basis for \f$\phi\f$ functions in all
	 *  domains.
	 *  
	 *  @param base_p   type of basis functions in \f$\phi\f$
	 *		    (e.g. \c P_COSSIN , etc..., 
	 *		     see documentation of class \c Base_val 
	 *		     for denomination of the various bases). 
	 */
	void set_base_p(int base_p) ; 

	/**
	 * Sets the \c n  lasts coefficients in \f$\theta\f$ to 0 in the 
	 * domain \c nz1 to \c nz2 when expressed in spherical harmonics.
	 */
	void filtre_tp(int nn, int nz1, int nz2) ;



    // Differential operators
    // ----------------------
    public:
	/// Returns \f$\partial / \partial \xi\f$ of \c *this 
	const Valeur& dsdx() const ;	    
	/// Returns \f$\partial^2 / \partial \xi^2\f$ of \c *this 
	const Valeur& d2sdx2() const ;     

	/// Returns \f$\partial / \partial \theta\f$ of \c *this 
	const Valeur& dsdt() const ;	
	/// Returns \f$\partial^2 / \partial \theta^2\f$ of \c *this 
	const Valeur& d2sdt2() const ;	
	/// Returns \f$1 / \sin(\theta)\f$ of \c *this 
	const Valeur& ssint() const ;	
	/// Returns \f$1 / \cos(\theta)\f$ of \c *this 
        const Valeur& scost() const ;  
	/// Returns \f$\cos(\theta) \, Id\f$ applied to \c *this 
	const Valeur& mult_ct() const ;
	/// Returns \f$\sin(\theta) \, Id\f$ applied to \c *this 
	const Valeur& mult_st() const ;

	/// Returns \f$\partial / \partial \phi\f$ of \c *this 
	const Valeur& dsdp() const ;	
	/// Returns \f$1/\sin(\theta) \, \partial / \partial \phi\f$ of \c *this 
	const Valeur& stdsdp() const ;	
	/// Returns \f$\partial^2 / \partial \phi^2\f$ of \c *this 
	const Valeur& d2sdp2() const ;	
	/// Returns \f$\cos(\phi) \, Id\f$ applied to \c *this 
	const Valeur& mult_cp() const ;
	/// Returns \f$\sin(\phi) \, Id\f$ applied to \c *this 
	const Valeur& mult_sp() const ;
	
	/// Returns the angular Laplacian of \c *this 
	const Valeur& lapang() const ;		

	/** Returns \f${1 \over \xi}\f$ (\e r -sampling = \c RARE ) \\
	 * Id (\e r  sampling = \c FIN ) \\
	 * \f${1 \over \xi-1}\f$ (\e r -sampling = \c UNSURR )
	 */
	const Valeur& sx() const ;
		    
	/** Returns \f${1 \over \xi^2}\f$ (\e r -sampling = \c RARE ) \\
	 * Id (\e r  sampling = \c FIN ) \\
	 * \f${1 \over (\xi-1)^2}\f$ (\e r -sampling = \c UNSURR )
	 */
	const Valeur& sx2() const ;	    
	
	/** Returns \f$\xi \, Id\f$ (\e r -sampling = \c RARE ) \\
	 * Id (\e r  sampling = \c FIN ) \\
	 * \f$(\xi-1) \, Id \f$ (\e r -sampling = \c UNSURR )
	 */
	const Valeur& mult_x() const ;
	
	/** Applies the following operator to \c *this : \\
	 * Id (\e r  sampling = \c RARE, FIN ) \\
	 * \f${1 \over (\xi-1)}\f$ (\e r -sampling = \c UNSURR )
	 */
	void sxm1_zec() ;

	/** Applies the following operator to \c *this : \\
	 * Id (\e r  sampling = \c RARE, FIN ) \\
	 * \f$(\xi-1) \, Id\f$ (\e r -sampling = \c UNSURR )
	 */
	void mult_xm1_zec() ;	

	/** Applies the following operator to \c *this : \\
	 * Id (\e r  sampling = \c RARE, FIN ) \\
	 * \f$(\xi-1)^2 \, Id\f$ (\e r -sampling = \c UNSURR )
	 */
	void mult2_xm1_zec() ;	
	
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    ///< Save in a file
    
	/** Displays the spectral coefficients and the associated
	 *  basis functions. This function shows only the values greater than a 
	 *  given threshold.
	 *   @param threshold [input] Value above which a coefficient is printed
	 *    (default: 1.e-7)
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param ostr [input] Output stream used for the printing (default: cout)
	 */
	void display_coef(double threshold = 1.e-7, int precision = 4, 
			   ostream& ostr = cout) const ;

	/** Prints only the values greater than a given threshold.
	 *   @param ostr [input] Output stream used for the printing
	 *   @param type [input] Type of display : 0 = prints only the
	 *     coefficients,  1 = prints only the values in configuration 
	 *     space, 2 = prints both
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param threshold [input] Value above which an array element is printed
	 *    (default: 1.e-7)
	 */
	void affiche_seuil(ostream& ostr, int type = 0, int precision = 4, 
			   double threshold = 1.e-7) const ;

	/// Display
	friend ostream& operator<<(ostream& , const Valeur& ) ;   

    // Memory management
    // -----------------
    private:
	void nouveau() ;	    ///< Memory allocation
	void del_t() ;		    ///< Logical destructor
	void del_deriv() ;	    ///< Logical destructor of the derivatives
	void set_der_0x0() ;	    ///< Sets the pointers for derivatives to 0x0

    // State manipulations
    public:

    /**
     * Sets the logical state to \c ETATNONDEF  (undefined). 
     * Deallocates the memory occupied by the \c Mtbl  \c c  and
     * the \c Mtbl_cf  \c c_cf ,  as well as by all the derivatives. 
     */
	void set_etat_nondef() ;  

    /**
     * Sets the logical state to \c ETATZERO  (zero). 
     * Deallocates the memory occupied by the \c Mtbl  \c c  and
     * the \c Mtbl_cf  \c c_cf ,  as well as by all the derivatives. 
     */
	void set_etat_zero() ;	    

    /**
     * Sets the logical state to \c ETATQCQ  (ordinary state) for
     * values in the configuration space (\c Mtbl  \c c ). 
     * If \c c  is not 0x0, this 
     * function does nothing on \c c . Otherwise, it performs the memory 
     * allocation for \c c .
     * In all cases, this function deallocates the memory occupied by the 
     * \c Mtbl_cf  \c c_cf ,  as well as by all the derivatives. 
     *
     */
	void set_etat_c_qcq() ;	    

    /**
     * Sets the logical state to \c ETATQCQ  (ordinary state) for
     * values in the configuration space (\c Mtbl_cf  \c c_cf ). 
     * If \c c_cf  is not 0x0, this 
     * function does nothing on \c c_cf . Otherwise, it performs the memory 
     * allocation for \c c_cf .
     * In all cases, this function deallocates the memory occupied by the 
     * \c Mtbl  \c c ,  as well as by all the derivatives. 
     *
     */
	void set_etat_cf_qcq() ;    

    /**
     * Sets the \c Valeur  to zero in a hard way. 
     * 1/ Sets the logical state to \c ETATQCQ , i.e. to an ordinary state.
     * 2/ Allocates the memory  for \c c  and \c c_cf , and fills it
     * with zeros. NB: this function must be used for debugging purposes only.
     * For other operations, the functions \c set_etat_zero() 
     * or \c annule(int,int)  must be perferred. 
     */
	void annule_hard() ;		    

    /**
     * Sets the \c Valeur  to zero in a given domain.
     *	@param l [input]  Index of the domain in which the \c Valeur 
     *			  will be set (logically) to zero.
     */
	void annule(int l) ; 

    /**
     * Sets the \c Valeur  to zero in several domains.
     *	@param l_min [input] The \c Valeur  will be set (logically) to zero
     *			     in the domains whose indices are in the range
     *			     \c [l_min,l_max] .
     *	@param l_max [input] see the comments for \c l_min .
     * 
     * Note that \c annule(0,mg->get_nzone()-1)  is equivalent to
     *	 \c set_etat_zero() .
     */
	void annule(int l_min, int l_max) ; 


    // Extraction of information
    // -------------------------
    public:
	/// Returns the logical state
	int get_etat() const {return etat ; };   

	/// Returns the \c Mg3d  on which the \c this  is defined
	const Mg3d* get_mg() const { return mg ; }; 
	
    // Member arithmetics
    // ------------------
    public:
	void operator+=(const Valeur& ) ;	///< += Valeur
	void operator-=(const Valeur& ) ;	///< -= Valeur
	void operator*=(const Valeur& ) ;	///< *= Valeur

    // Miscellaneous
    // -------------
    public:
	/** Determines an equipotential surface of the field represented 
	 *   by \c *this  (inward search).
	 *  The equipotential is supposed to have the form \\ 
	 *   \f$l=L(\theta', \phi') \qquad (1) \f$    \\
	 *   \f$\xi = X(\theta', \phi') \qquad (2)\f$ \\
	 *  where \e l  is the domain index and \f$\xi\f$ the radial variable in
	 *  each domain. 
	 * @param uu0 [input] Value defining the equipotential by
	 *			\e u  = const = \c uu0  where \e u  is 
	 *			the field represented by \c *this .
	 * @param nz_search [input] Number of domains where the equipotential is
	 *			searched : the routine scans inward 
	 *			the \c nz_search 
	 *			innermost domains, starting from the domain
	 *			of index \c nz_search-1 
	 * @param precis [input] Required absolute precision in the 
	 *			  determination of zeros by the secant 
	 *			  method (standard value: 1.e-14).
	 * @param nitermax [input] Maximum number of iterations in the secant 
	 *			    method (standard value: 100).
	 * @param niter [output] Number of iterations effectively used 
	 *			    in the secant method
	 * @param l_iso [output] 2-D \c Itbl  containing the values
	 *	of \e l  on the equipotential surface at the collocation points
	 *	 in \f$(\theta', \phi')\f$ [Eq. (1)], with the following storage
	 *	 convention \\
	 *	 \c l_iso(k,j)  = \f$L({\theta'}_j, {\phi'}_k)\f$ 
	 * @param xi_iso [output] 2-D \c Tbl  containing the values
	 *	of \f$\xi\f$ on the equipotential surface at the collocation points
	 *	 in \f$(\theta', \phi')\f$ [Eq. (2)], with the following storage
	 *	 convention \\
	 *	 \c xi_iso(k,j)  = \f$X({\theta'}_j, {\phi'}_k)\f$ 
	 * 
	 */
	void equipot(double uu0, int nz_search, double precis, int nitermax,
		     int& niter, Itbl& l_iso, Tbl& xi_iso) const ; 
    
	/** Determines an equipotential surface of the field represented 
	 *   by \c *this  (outward search).
	 *  The equipotential is supposed to have the form \\ 
	 *   \f$l=L(\theta', \phi') \qquad (1) \f$    \\
	 *   \f$\xi = X(\theta', \phi') \qquad (2)\f$ \\
	 *  where \e l  is the domain index and \f$\xi\f$ the radial variable in
	 *  each domain. 
	 * @param uu0 [input] Value defining the equipotential by
	 *			\e u  = const = \c uu0  where \e u  is 
	 *			the field represented by \c *this .
	 * @param nz_search [input] Number of domains where the equipotential is
	 *			searched : the routine scans outward the 
	 *			\c nz_search  innermost domains, starting 
	 *			from the domain of index 0
	 * @param precis [input] Required absolute precision in the 
	 *			  determination of zeros by the secant 
	 *			  method (standard value: 1.e-14).
	 * @param nitermax [input] Maximum number of iterations in the secant 
	 *			    method (standard value: 100).
	 * @param niter [output] Number of iterations effectively used 
	 *			    in the secant method
	 * @param l_iso [output] 2-D \c Itbl  containing the values
	 *	of \e l  on the equipotential surface at the collocation points
	 *	 in \f$(\theta', \phi')\f$ [Eq. (1)], with the following storage
	 *	 convention \\
	 *	 \c l_iso(k,j)  = \f$L({\theta'}_j, {\phi'}_k)\f$ 
	 * @param xi_iso [output] 2-D \c Tbl  containing the values
	 *	of \f$\xi\f$ on the equipotential surface at the collocation points
	 *	 in \f$(\theta', \phi')\f$ [Eq. (2)], with the following storage
	 *	 convention \\
	 *	 \c xi_iso(k,j)  = \f$X({\theta'}_j, {\phi'}_k)\f$ 
	 * 
	 */
	void equipot_outward(double uu0, int nz_search, double precis, 
			     int nitermax, int& niter, Itbl& l_iso, 
			     Tbl& xi_iso) const ; 

	/** Changes the function \c *this  as a smooth one
	 *   when there exists a discontinuity
	 *   between the nucleus and the first shell.
	 *
	 * @param nzet [input] Number of domains covering a star.
	 * @param uuva [output] Smoothed function.
	 * 
	 */
	void smooth(int nzet, Valeur& uuva) const ;

    friend class Cmp ;	    ///< Friend class
    friend class Scalar ;	    ///< Friend class
    friend void rotate_propre_pair (Valeur&, bool) ; ///< Friend fonction.
    friend void rotate_propre_impair (Valeur&, bool) ; ///< Friend fonction.
};
ostream& operator<<(ostream& , const Valeur& ) ;   

/**
 * \defgroup val_m Valeur Mathematics
 * \ingroup (spec)
 * @{
 */
Valeur operator+(const Valeur& ) ;	///< + Valeur
Valeur operator-(const Valeur& ) ;	///< \c - Valeur
Valeur operator+(const Valeur&, const Valeur& ) ; ///< Valeur + Valeur
Valeur operator+(const Valeur&, double ) ;	  ///< Valeur + double
Valeur operator+(double, const Valeur& ) ;	  ///< double + Valeur 
Valeur operator+(const Valeur&, int ) ;		  ///< Valeur + int
Valeur operator+(int, const Valeur& ) ;		  ///< int + Valeur 
Valeur operator-(const Valeur&, const Valeur& ) ; ///< Valeur - Valeur
Valeur operator-(const Valeur&, double ) ;	  ///< Valeur - double
Valeur operator-(double, const Valeur& ) ;	  ///< double - Valeur 
Valeur operator-(const Valeur&, int ) ;		  ///< Valeur - int
Valeur operator-(int, const Valeur& ) ;		  ///< int - Valeur 
Valeur operator*(const Valeur&, const Valeur& ) ; ///< Valeur * Valeur

/// Valeur * Valeur with desaliasing
Valeur operator%(const Valeur&, const Valeur& ) ; 

/// Valeur * Valeur with desaliasing only in \e r direction
Valeur operator|(const Valeur&, const Valeur& ) ; 

Valeur operator*(const Valeur&, double ) ;	  ///< Valeur * double
Valeur operator*(double, const Valeur& ) ;	  ///< double * Valeur 
Valeur operator*(const Valeur&, int ) ;		  ///< Valeur * int
Valeur operator*(int, const Valeur& ) ;		  ///< int * Valeur 
Valeur operator*(const Valeur& a, const Mtbl& b) ;	  ///< Valeur * Mtbl
Valeur operator*(const Mtbl& b, const Valeur& a) ;	  ///< Mtbl * Valeur
Valeur operator*(const Valeur& a, const Coord& c) ;  ///< Valeur * Coord
Valeur operator*(const Coord& c, const Valeur& a) ;  ///< Coord * Coord
Valeur operator/(const Valeur& a, const Valeur& b) ; ///< Valeur / Valeur
Valeur operator/(const Valeur&, double ) ;	  ///< Valeur / double
Valeur operator/(double, const Valeur& ) ;	  ///< double / Valeur 
Valeur operator/(const Valeur&, int ) ;		  ///< Valeur / int
Valeur operator/(int, const Valeur& ) ;		  ///< int / Valeur 
Valeur operator/(const Valeur& a, const Mtbl& b) ;	  ///< Valeur / Mtbl
Valeur operator/(const Mtbl& b, const Valeur& a) ;	  ///< Mtbl / Valeur

Valeur sin(const Valeur& ) ;	    ///< Sine
Valeur cos(const Valeur& ) ;	    ///< Cosine
Valeur tan(const Valeur& ) ;	    ///< Tangent
Valeur asin(const Valeur& ) ;	    ///< Arcsine
Valeur acos(const Valeur& ) ;	    ///< Arccosine
Valeur atan(const Valeur& ) ;	    ///< Arctangent
Valeur exp(const Valeur& ) ;	    ///< Exponential
Valeur log(const Valeur& ) ;	    ///< Neperian logarithm
Valeur log10(const Valeur& ) ;      ///< Basis 10 logarithm
Valeur sqrt(const Valeur& ) ;	    ///< Square root
Valeur pow(const Valeur& , int ) ;  ///< Power \f${\tt Valeur}^{\tt int}\f$
Valeur pow(const Valeur& , double ) ; ///< Power \f${\tt Valeur}^{\tt double}\f$
Valeur abs(const Valeur& ) ;	    ///< Absolute value
Valeur racine_cubique (const Valeur&) ; ///< Cube root

/**
 * Maximum values of the \c Valeur  (configuration space)
 * in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Valeur& ) ;   

/**
 * Minimum values of the \c Valeur  (configuration space)
 * in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Valeur& ) ;   

/**
 * Sums of the absolute values of all the \c Valeur  (configuration space)
 * in each domain.
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Valeur& ) ;   

/**
 * Relative difference between two \c Valeur  (configuration space)
 * (norme version).
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are \c norme[a(l)-b(l)]/norme[b(l)]  if \c b(l)!=0  and
 *	   \c norme[a(l)-b(l)]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c a  and \c b  
 *	   in domain no. \c l . 
 */
Tbl diffrel(const Valeur& a, const Valeur& b) ; 

/**
 * Relative difference between two \c Valeur  (configuration space)
 * (max version).
 * @return 1-D \c Tbl  of size the number of domains, the elements of which 
 *	   are \c max[abs(a(l)-b(l))]/max[abs(b(l))]  if \c b(l)!=0  and
 *	   \c max[abs(a(l)-b(l))]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c a  and \c b  
 *	   in domain no. \c l . 
 */
Tbl diffrelmax(const Valeur& a, const Valeur& b) ; 


/** @} */
	
#endif
