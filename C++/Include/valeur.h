/*
 *  Definition of Lorene class Valeur
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
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
#include <iostream.h>


#include "mtbl.h"
#include "mtbl_cf.h"

class Coord ;
class Itbl ; 

/**
 * Values and coefficients of a (real-value) function.
 * 
 * @version #$Id$#
 *
 */

class Valeur {

    // Data : 
    // -----
    private:
	const Mg3d* mg ;  /// Multi-grid {\tt Mgd3} on which {\tt this} is defined

	/// Logical state ({\tt ETATNONDEF}, {\tt ETATQCQ} or {\tt ETATZERO}).
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
	mutable Valeur* p_dsdx ;    /// Pointer on $\partial / \partial \xi$ 
	mutable Valeur* p_d2sdx2 ;  /// Pointer on $\partial^2 / \partial \xi^2$
	mutable Valeur* p_sx ;	    /// Pointer on $1 / \xi$
	mutable Valeur* p_sx2 ;	    /// Pointer on $1 / \xi^2$
	mutable Valeur* p_mult_x ;  /// Pointer on $\xi \, Id$

	mutable Valeur* p_dsdt ;    /// Pointer on $\partial / \partial \theta$
	mutable Valeur* p_d2sdt2 ;  /// Pointer on $\partial^2 / \partial \theta^2$
	mutable Valeur* p_ssint ;   /// Pointer on $1 / \sin(\theta)$
        mutable Valeur* p_scost ;   /// Pointer on $1 / \cos(\theta)$	  
	mutable Valeur* p_mult_ct ; /// Pointer on $\cos(\theta) \, Id$
	mutable Valeur* p_mult_st ; /// Pointer on $\sin(\theta) \, Id$

	mutable Valeur* p_dsdp ;    /// Pointer on $\partial / \partial \phi$
	mutable Valeur* p_stdsdp ;  /// Pointer on $1/\sin\theta \partial / \partial \phi$
	mutable Valeur* p_d2sdp2 ;  /// Pointer on $\partial^2 / \partial \phi^2$
	mutable Valeur* p_mult_cp ; /// Pointer on $\cos(\phi) \, Id$
	mutable Valeur* p_mult_sp ; /// Pointer on $\sin(\phi) \, Id$

	mutable Valeur* p_lapang ;  /// Pointer on the angular Laplacian

    // Constructors - Destructor
    // -------------------------
	
    public:
	explicit Valeur(const Mg3d& ) ;		/// Constructor
	explicit Valeur(const Mg3d* ) ;		/// Constructor

	/// Constructor from a file (see {\tt sauve(FILE* )})
	Valeur(const Mg3d&, FILE* ) ;    		

	Valeur(const Valeur& ) ;	/// Copy constructor
	~Valeur() ;			/// Destructor

    // Assignement
    // -----------
    public: 
	void operator=(const Valeur& ) ; /// Assignement to another {\tt Valeur}
	void operator=(const Mtbl& ) ;	 /// Assignement to a {\tt Mtbl}
	void operator=(const Mtbl_cf& ) ; /// Assignement to a {\tt Mtbl\_cf}
	void operator=(double ) ;	/// Assignement to a {\tt double}	    

    // Access to individual elements
    // -----------------------------
    public:
	/** Read/write of the value in a given domain (configuration space).
	 * NB: to gain in efficiency, the method {\tt del\_deriv()} (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *
	 * @param l [input] domain index
	 * @return  Tbl containing the value of the field in domain {\tt l}.
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
	 * @return  Tbl containing the value of the field in domain {\tt l}.
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
	 * NB: to gain in efficiency, the method {\tt del\_deriv()} (to delete
	 *     the derived members) is not called by this function. It must
	 *     thus be invoqued by the user.  
	 *
	 * @param l [input] domain index
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] $r$ ($\xi$) index
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
	 * @param k [input] $\phi$ index
	 * @param j [input] $\theta$ index
	 * @param i [input] $r$ ($\xi$) index
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

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point, by means of the spectral expansion.
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate $\xi$
	*	 @param theta [input] value of the coordinate $\theta'$
	*	 @param phi [input] value of the coordinate $\phi'$
	*	 @return value at the point $(\xi, \theta', \phi')$ in
	*	    the domain no. $l$ of the field represented by {\tt *this}. 
	*/
	double val_point(int l, double x, double theta, double phi) const ; 

	/** Computes the value of the field represented by {\tt *this} at an
	*   arbitrary point in $\xi$, but collocation point in 
	*   $(\theta', \phi')$, by means of the spectral expansion.
	*	 @param l [input] index of the domain
	*	 @param x [input] value of the coordinate $\xi$
	*	 @param j [input] index of the collocation point in $\theta'$
	*	 @param k [input] index of the collocation point in $\phi'$
	*	 @return value at the point 
	*		    $(\xi, {\theta'}_j, {\phi'}_k)$ in
	*	    the domain no. $l$ of the field represented by {\tt *this}. 
	*/
	double val_point_jk(int l, double x, int j, int k) const ; 


    // Operations on coefficients
    // --------------------------
    public:
	void coef() const ;	/// Computes the coeffcients of {\tt *this}
	void coef_i() const ;	/// Computes the physical value of {\tt *this}
	void ylm() ;	    /// Computes the coefficients $Y_l^m$ of {\tt *this}
	void ylm_i() ;	    /// Inverse of {\tt ylm()} 
	
	/// Sets the bases for spectral expansions (member {\tt base}) 
	void set_base(const Base_val& ) ; 

	/** Sets the bases for spectral expansions (member {\tt base}) 
	 *  to the standard ones for a scalar
	 */
	void std_base_scal() ;	 

	/** Sets the expansion basis for $r$ ($\xi$) functions in a 
	 *  given domain.
	 *  
	 *  @param l	    Domain index
	 *  @param base_r   type of basis functions in $r$ ($\xi$)
	 *		    (e.g. {\tt R\_CHEB\_P}, etc..., 
	 *		     see documentation of class {\tt Base\_val}
	 *		     for denomination of the various bases). 
	 */
	void set_base_r(int l, int base_r) ; 

	/** Sets the expansion basis for $\theta$ functions in all
	 *  domains.
	 *  
	 *  @param base_t   type of basis functions in $\theta$
	 *		    (e.g. {\tt T\_COS\_P}, etc..., 
	 *		     see documentation of class {\tt Base\_val}
	 *		     for denomination of the various bases). 
	 */
	void set_base_t(int base_t) ; 

	/** Sets the expansion basis for $\phi$ functions in all
	 *  domains.
	 *  
	 *  @param base_p   type of basis functions in $\phi$
	 *		    (e.g. {\tt P\_COSSIN}, etc..., 
	 *		     see documentation of class {\tt Base\_val}
	 *		     for denomination of the various bases). 
	 */
	void set_base_p(int base_p) ; 

    // Differential operators
    // ----------------------
    public:
	/// Returns $\partial / \partial \xi$ of {\tt *this}
	const Valeur& dsdx() const ;	    
	/// Returns $\partial^2 / \partial \xi^2$ of {\tt *this}
	const Valeur& d2sdx2() const ;     

	/// Returns $\partial / \partial \theta$ of {\tt *this}
	const Valeur& dsdt() const ;	
	/// Returns $\partial^2 / \partial \theta^2$ of {\tt *this}
	const Valeur& d2sdt2() const ;	
	/// Returns $1 / \sin(\theta)$ of {\tt *this}
	const Valeur& ssint() const ;	
	/// Returns $1 / \cos(\theta)$ of {\tt *this}
        const Valeur& scost() const ;  
	/// Returns $\cos(\theta) \, Id$ applied to {\tt *this}
	const Valeur& mult_ct() const ;
	/// Returns $\sin(\theta) \, Id$ applied to {\tt *this}
	const Valeur& mult_st() const ;

	/// Returns $\partial / \partial \phi$ of {\tt *this}
	const Valeur& dsdp() const ;	
	/// Returns $1/\sin(\theta) \, \partial / \partial \phi$ of {\tt *this}
	const Valeur& stdsdp() const ;	
	/// Returns $\partial^2 / \partial \phi^2$ of {\tt *this}
	const Valeur& d2sdp2() const ;	
	/// Returns $\cos(\phi) \, Id$ applied to {\tt *this}
	const Valeur& mult_cp() const ;
	/// Returns $\sin(\phi) \, Id$ applied to {\tt *this}
	const Valeur& mult_sp() const ;
	
	/// Returns the angular Laplacian of {\tt *this}
	const Valeur& lapang() const ;		

	/** Returns ${1 \over \xi}$ ($r$-sampling = {\tt RARE}) \\
	 * Id ($r$ sampling = {\tt FIN}) \\
	 * ${1 \over \xi-1}$ ($r$-sampling = {\tt UNSURR})
	 */
	const Valeur& sx() const ;
		    
	/** Returns ${1 \over \xi^2}$ ($r$-sampling = {\tt RARE}) \\
	 * Id ($r$ sampling = {\tt FIN}) \\
	 * ${1 \over (\xi-1)^2}$ ($r$-sampling = {\tt UNSURR})
	 */
	const Valeur& sx2() const ;	    
	
	/** Returns $\xi \, Id$ ($r$-sampling = {\tt RARE}) \\
	 * Id ($r$ sampling = {\tt FIN}) \\
	 * $(\xi-1) \, Id $ ($r$-sampling = {\tt UNSURR})
	 */
	const Valeur& mult_x() const ;
	
	/** Applies the following operator to {\tt *this}: \\
	 * Id ($r$ sampling = {\tt RARE, FIN}) \\
	 * ${1 \over (\xi-1)}$ ($r$-sampling = {\tt UNSURR})
	 */
	void sxm1_zec() ;

	/** Applies the following operator to {\tt *this}: \\
	 * Id ($r$ sampling = {\tt RARE, FIN}) \\
	 * $(\xi-1) \, Id$ ($r$-sampling = {\tt UNSURR})
	 */
	void mult_xm1_zec() ;	

	/** Applies the following operator to {\tt *this}: \\
	 * Id ($r$ sampling = {\tt RARE, FIN}) \\
	 * $(\xi-1)^2 \, Id$ ($r$-sampling = {\tt UNSURR})
	 */
	void mult2_xm1_zec() ;	
	
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file
    
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
	void nouveau() ;	    /// Memory allocation
	void del_t() ;		    /// Logical destructor
	void del_deriv() ;	    /// Logical destructor of the derivatives
	void set_der_0x0() ;	    /// Sets the pointers for derivatives to 0x0

    // State manipulations
    public:

    /**
     * Sets the logical state to {\tt ETATNONDEF} (undefined). 
     * Deallocates the memory occupied by the {\tt Mtbl} {\tt c} and
     * the {\tt Mtbl\_cf} {\tt c\_cf},  as well as by all the derivatives. 
     */
	void set_etat_nondef() ;  

    /**
     * Sets the logical state to {\tt ETATZERO} (zero). 
     * Deallocates the memory occupied by the {\tt Mtbl} {\tt c} and
     * the {\tt Mtbl\_cf} {\tt c\_cf},  as well as by all the derivatives. 
     */
	void set_etat_zero() ;	    

    /**
     * Sets the logical state to {\tt ETATQCQ} (ordinary state) for
     * values in the configuration space ({\tt Mtbl} {\tt c}). 
     * If {\tt c} is not 0x0, this 
     * function does nothing on {\tt c}. Otherwise, it performs the memory 
     * allocation for {\tt c}.
     * In all cases, this function deallocates the memory occupied by the 
     * {\tt Mtbl\_cf} {\tt c\_cf},  as well as by all the derivatives. 
     *
     */
	void set_etat_c_qcq() ;	    

    /**
     * Sets the logical state to {\tt ETATQCQ} (ordinary state) for
     * values in the configuration space ({\tt Mtbl\_cf} {\tt c\_cf}). 
     * If {\tt c\_cf} is not 0x0, this 
     * function does nothing on {\tt c\_cf}. Otherwise, it performs the memory 
     * allocation for {\tt c\_cf}.
     * In all cases, this function deallocates the memory occupied by the 
     * {\tt Mtbl} {\tt c},  as well as by all the derivatives. 
     *
     */
	void set_etat_cf_qcq() ;    

    /**
     * Sets the {\tt Valeur} to zero in a hard way. 
     * 1/ Sets the logical state to {\tt ETATQCQ}, i.e. to an ordinary state.
     * 2/ Allocates the memory  for {\tt c} and {\tt c\_cf}, and fills it
     * with zeros. NB: this function must be used for debugging purposes only.
     * For other operations, the functions {\tt set\_etat\_zero()}
     * or {\tt annule(int, int)} must be perferred. 
     */
	void annule_hard() ;		    

    /**
     * Sets the {\tt Valeur} to zero in a given domain.
     *	@param l [input]  Index of the domain in which the {\tt Valeur}
     *			  will be set (logically) to zero.
     */
	void annule(int l) ; 

    /**
     * Sets the {\tt Valeur} to zero in several domains.
     *	@param l_min [input] The {\tt Valeur} will be set (logically) to zero
     *			     in the domains whose indices are in the range
     *			     {\tt [l\_min, l\_max]}.
     *	@param l_max [input] see the comments for {\tt l\_min}.
     * 
     * Note that {\tt annule(0, mg->get\_nzone()-1)} is equivalent to
     *	 {\tt set\_etat\_zero()}.
     */
	void annule(int l_min, int l_max) ; 


    // Extraction of information
    // -------------------------
    public:
	int get_etat() const {return etat ; };   /// Returns the logical state

	/// Returns the {\tt Mg3d} on which the {\tt this} is defined
	const Mg3d* get_mg() const { return mg ; }; 
	
    // Member arithmetics
    // ------------------
    public:
	void operator+=(const Valeur& ) ;	/// += Valeur
	void operator-=(const Valeur& ) ;	/// -= Valeur
	void operator*=(const Valeur& ) ;	/// *= Valeur

    // Miscellaneous
    // -------------
    public:
	/** Determines an equipotential surface of the field represented 
	 *   by {\tt *this} (inward search).
	 *  The equipotential is supposed to have the form \\ 
	 *   $l=L(\theta', \phi') \qquad (1) $    \\
	 *   $\xi = X(\theta', \phi') \qquad (2)$ \\
	 *  where $l$ is the domain index and $\xi$ the radial variable in
	 *  each domain. 
	 * @param uu0 [input] Value defining the equipotential by
	 *			$u = {\rm const} = {\tt uu0}$ where $u$ is 
	 *			the field represented by {\tt *this}.
	 * @param nz_search [input] Number of domains where the equipotential is
	 *			searched : the routine scans inward 
	 *			the {\tt nz\_search}
	 *			innermost domains, starting from the domain
	 *			of index {\tt nz\_search-1}
	 * @param precis [input] Required absolute precision in the 
	 *			  determination of zeros by the secant 
	 *			  method (standard value: 1.e-14).
	 * @param nitermax [input] Maximum number of iterations in the secant 
	 *			    method (standard value: 100).
	 * @param niter [output] Number of iterations effectively used 
	 *			    in the secant method
	 * @param l_iso [output] 2-D {\tt Itbl} containing the values
	 *	of $l$ on the equipotential surface at the collocation points
	 *	 in $(\theta', \phi')$ [Eq. (1)], with the following storage
	 *	 convention \\
	 *	 {\tt l\_iso(k, j)} = $L({\theta'}_j, {\phi'}_k)$ 
	 * @param xi_iso [output] 2-D {\tt Tbl} containing the values
	 *	of $\xi$ on the equipotential surface at the collocation points
	 *	 in $(\theta', \phi')$ [Eq. (2)], with the following storage
	 *	 convention \\
	 *	 {\tt xi\_iso(k, j)} = $X({\theta'}_j, {\phi'}_k)$ 
	 * 
	 */
	void equipot(double uu0, int nz_search, double precis, int nitermax,
		     int& niter, Itbl& l_iso, Tbl& xi_iso) const ; 
    
	/** Determines an equipotential surface of the field represented 
	 *   by {\tt *this} (outward search).
	 *  The equipotential is supposed to have the form \\ 
	 *   $l=L(\theta', \phi') \qquad (1) $    \\
	 *   $\xi = X(\theta', \phi') \qquad (2)$ \\
	 *  where $l$ is the domain index and $\xi$ the radial variable in
	 *  each domain. 
	 * @param uu0 [input] Value defining the equipotential by
	 *			$u = {\rm const} = {\tt uu0}$ where $u$ is 
	 *			the field represented by {\tt *this}.
	 * @param nz_search [input] Number of domains where the equipotential is
	 *			searched : the routine scans outward the 
	 *			{\tt nz\_search} innermost domains, starting 
	 *			from the domain of index 0
	 * @param precis [input] Required absolute precision in the 
	 *			  determination of zeros by the secant 
	 *			  method (standard value: 1.e-14).
	 * @param nitermax [input] Maximum number of iterations in the secant 
	 *			    method (standard value: 100).
	 * @param niter [output] Number of iterations effectively used 
	 *			    in the secant method
	 * @param l_iso [output] 2-D {\tt Itbl} containing the values
	 *	of $l$ on the equipotential surface at the collocation points
	 *	 in $(\theta', \phi')$ [Eq. (1)], with the following storage
	 *	 convention \\
	 *	 {\tt l\_iso(k, j)} = $L({\theta'}_j, {\phi'}_k)$ 
	 * @param xi_iso [output] 2-D {\tt Tbl} containing the values
	 *	of $\xi$ on the equipotential surface at the collocation points
	 *	 in $(\theta', \phi')$ [Eq. (2)], with the following storage
	 *	 convention \\
	 *	 {\tt xi\_iso(k, j)} = $X({\theta'}_j, {\phi'}_k)$ 
	 * 
	 */
	void equipot_outward(double uu0, int nz_search, double precis, 
			     int nitermax, int& niter, Itbl& l_iso, 
			     Tbl& xi_iso) const ; 

	/** Changes the function {\tt *this} as a smooth one
	 *   when there exists a discontinuity
	 *   between the nucleus and the first shell.
	 *
	 * @param nzet [input] Number of domains covering a star.
	 * @param uuva [output] Smoothed function.
	 * 
	 */
	void smooth(int nzet, Valeur& uuva) const ;

    friend class Cmp ;	    /// Friend class
};

/**
 * @name Valeur mathematics
 */
    //@{
Valeur operator+(const Valeur& ) ;	/// + Valeur
Valeur operator-(const Valeur& ) ;	/// - Valeur
Valeur operator+(const Valeur&, const Valeur& ) ; /// Valeur + Valeur
Valeur operator+(const Valeur&, double ) ;	  /// Valeur + double
Valeur operator+(double, const Valeur& ) ;	  /// double + Valeur 
Valeur operator+(const Valeur&, int ) ;		  /// Valeur + int
Valeur operator+(int, const Valeur& ) ;		  /// int + Valeur 
Valeur operator-(const Valeur&, const Valeur& ) ; /// Valeur - Valeur
Valeur operator-(const Valeur&, double ) ;	  /// Valeur - double
Valeur operator-(double, const Valeur& ) ;	  /// double - Valeur 
Valeur operator-(const Valeur&, int ) ;		  /// Valeur - int
Valeur operator-(int, const Valeur& ) ;		  /// int - Valeur 
Valeur operator*(const Valeur&, const Valeur& ) ; /// Valeur * Valeur

/// Valeur * Valeur with desaliasing
Valeur operator%(const Valeur&, const Valeur& ) ; 

Valeur operator*(const Valeur&, double ) ;	  /// Valeur * double
Valeur operator*(double, const Valeur& ) ;	  /// double * Valeur 
Valeur operator*(const Valeur&, int ) ;		  /// Valeur * int
Valeur operator*(int, const Valeur& ) ;		  /// int * Valeur 
Valeur operator*(const Valeur&, const Mtbl& ) ;	  /// Valeur * Mtbl
Valeur operator*(const Mtbl&, const Valeur& ) ;	  /// Mtbl * Valeur
Valeur operator*(const Valeur&, const Coord& ) ;  /// Valeur * Coord
Valeur operator*(const Coord&, const Valeur& ) ;  /// Coord * Coord
Valeur operator/(const Valeur&, const Valeur& ) ; /// Valeur / Valeur
Valeur operator/(const Valeur&, double ) ;	  /// Valeur / double
Valeur operator/(double, const Valeur& ) ;	  /// double / Valeur 
Valeur operator/(const Valeur&, int ) ;		  /// Valeur / int
Valeur operator/(int, const Valeur& ) ;		  /// int / Valeur 
Valeur operator/(const Valeur&, const Mtbl& ) ;	  /// Valeur / Mtbl
Valeur operator/(const Mtbl&, const Valeur& ) ;	  /// Mtbl / Valeur

Valeur sin(const Valeur& ) ;	    /// Sine
Valeur cos(const Valeur& ) ;	    /// Cosine
Valeur tan(const Valeur& ) ;	    /// Tangent
Valeur asin(const Valeur& ) ;	    /// Arcsine
Valeur acos(const Valeur& ) ;	    /// Arccosine
Valeur atan(const Valeur& ) ;	    /// Arctangent
Valeur exp(const Valeur& ) ;	    /// Exponential
Valeur log(const Valeur& ) ;	    /// Neperian logarithm
Valeur log10(const Valeur& ) ;      /// Basis 10 logarithm
Valeur sqrt(const Valeur& ) ;	    /// Square root
Valeur pow(const Valeur& , int ) ;  /// Power ${\tt Valeur}^{\tt int}$
Valeur pow(const Valeur& , double ) ; /// Power ${\tt Valeur}^{\tt double}$
Valeur abs(const Valeur& ) ;	    /// Absolute value
Valeur racine_cubique (const Valeur&) ; /// Cube root

/**
 * Maximum values of the {\tt Valeur} (configuration space)
 * in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the maximum values in each domain.  
 */
Tbl max(const Valeur& ) ;   

/**
 * Minimum values of the {\tt Valeur} (configuration space)
 * in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the minimum values in each domain.  
 */
Tbl min(const Valeur& ) ;   

/**
 * Sums of the absolute values of all the {\tt Valeur} (configuration space)
 * in each domain.
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are the set of the sums of the absolute values in each domain.  
 */
Tbl norme(const Valeur& ) ;   

/**
 * Relative difference between two {\tt Valeur} (configuration space)
 * (norme version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt norme[a(l)-b(l)]/norme[b(l)]} if {\tt b(l)!=0} and
 *	   {\tt norme[a(l)-b(l)]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrel(const Valeur& a, const Valeur& b) ; 

/**
 * Relative difference between two {\tt Valeur} (configuration space)
 * (max version).
 * @return 1-D {\tt Tbl} of size the number of domains, the elements of which 
 *	   are {\tt max[abs(a(l)-b(l))]/max[abs(b(l))]} if {\tt b(l)!=0} and
 *	   {\tt max[abs(a(l)-b(l))]} if  {\tt b(l)=0},  where {\tt a(l)} and 
 *	   {\tt b(l)} denote symbolically the values of {\tt a} and {\tt b} 
 *	   in domain no. {\tt l}. 
 */
Tbl diffrelmax(const Valeur& a, const Valeur& b) ; 


    //@}
	
#endif
