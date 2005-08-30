/*
 *  Definition of Lorene classes Tenseur
 *				 Tenseur_sym
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
 *   Copyright (c) 2002 Jerome Novak
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


#ifndef __TENSEUR_H_
#define __TENSEUR_H_


/*
 * $Id$
 * $Log$
 * Revision 1.16  2005/08/30 08:35:10  p_grandclement
 * Addition of the Tau version of the vectorial Poisson equation for the Tensors
 *
 * Revision 1.15  2004/12/29 16:24:03  k_taniguchi
 * Addition of the method for vectorial Poisson equations with a multipole
 * falloff condition at the outer boundary.
 *
 * Revision 1.14  2004/12/22 18:25:12  k_taniguchi
 * Change an argument of poisson_vect_falloff
 *
 * Revision 1.13  2004/11/30 20:43:32  k_taniguchi
 * Addition of the method for vectorial Poisson equations with falloff
 * condition at the outer boundary.
 *
 * Revision 1.12  2004/03/22 13:12:43  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.11  2003/11/06 12:17:31  r_prix
 * fixed mini-bug in documentation: without explicit argument in function-prototype,
 * doc++ seemed to merge docu of Tensor::operator=(Cmp&) and operator=(Tenseur&)
 *
 * Revision 1.10  2003/06/20 14:23:38  f_limousin
 * Add the functions compare().
 *
 * Revision 1.9  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.8  2002/09/19 14:12:37  e_gourgoulhon
 * Modif documentation for LaTeX compliance.
 *
 * Revision 1.7  2002/09/10 13:44:17  j_novak
 * The method "manipule" of one indice has been removed for Tenseur_sym objects
 * (the result cannot be a Tenseur_sym).
 * The method "sans_trace" now computes the traceless part of a Tenseur (or
 * Tenseur_sym) of valence 2.
 *
 * Revision 1.6  2002/09/06 14:49:25  j_novak
 * Added method lie_derive for Tenseur and Tenseur_sym.
 * Corrected various errors for derive_cov and arithmetic.
 *
 * Revision 1.5  2002/08/14 13:46:14  j_novak
 * Derived quantities of a Tenseur can now depend on several Metrique's
 *
 * Revision 1.4  2002/08/08 15:10:44  j_novak
 * The flag "plat" has been added to the class Metrique to show flat metrics.
 *
 * Revision 1.3  2002/08/07 16:14:11  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.2  2002/06/17 14:05:17  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.46  2001/08/27  10:03:56  eric
 * Ajout de l'operator% (produit tensoriel avec desaliasing)
 *
 * Revision 2.45  2001/06/19  15:38:00  eric
 * Modif commentaires: mise en conformite Doc++ 3.4.8.
 *
 * Revision 2.44  2001/06/18  13:56:07  novak
 * Ajout de la fonction abs()
 *
 * Revision 2.43  2001/05/29 16:10:43  eric
 * Modif commentaires (mise en conformite Doc++ 3.4.7).
 *
 * Revision 2.42  2001/05/26  15:48:17  eric
 * *** empty log message ***
 *
 * Revision 2.41  2001/05/26  15:42:54  eric
 * Ajout de la fonction flat_scalar_prod_desal (desaliasage)
 *
 * Revision 2.40  2000/10/19  10:37:51  phil
 * *** empty log message ***
 *
 * Revision 2.39  2000/10/19  09:47:15  phil
 * *** empty log message ***
 *
 * Revision 2.38  2000/10/19  09:41:48  phil
 * ajout de inverse_poisson_vect
 *
 * Revision 2.37  2000/10/06  12:37:09  keisuke
 * Add a vectorial Poisson equation Tenseur::poisson_vect_regu.
 *
 * Revision 2.36  2000/09/27  08:52:34  eric
 * Modif commentaires.
 *
 * Revision 2.35  2000/09/13  12:21:36  eric
 * Modif commentaires.
 *
 * Revision 2.34  2000/09/13  12:11:26  eric
 * Ajout de la fonction allocate_all().
 *
 * Revision 2.33  2000/05/22  14:39:11  phil
 * ajout de inc_dzpuis et dec_dzpuis
 *
 * Revision 2.32  2000/04/03  15:18:54  phil
 * suppression de poisson_vect_dirichlet
 *
 * Revision 2.31  2000/03/31  13:29:43  phil
 * poisson_vect_neumann devient poisson_vect_dirichlet
 *
 * Revision 2.30  2000/03/30  16:10:46  phil
 * *** empty log message ***
 *
 * Revision 2.29  2000/02/18  10:45:05  eric
 * Modif commentaires.
 *
 * Revision 2.28  2000/02/11  19:11:55  phil
 * commentaires
 *
 * Revision 2.27  2000/02/10  16:26:48  eric
 * Modif commentaires.
 *
 * Revision 2.26  2000/02/10  16:10:49  eric
 * Ajout de la fonction change_triad.
 *
 * Revision 2.25  2000/02/09  19:29:25  eric
 * MODIF IMPORTANTE: la triade de decomposition est desormais passee en
 * argument des constructeurs.
 *
 * Revision 2.24  2000/02/09  09:53:56  phil
 * ajout de poisson_vect_oohara
 * ,
 *
 * Revision 2.23  2000/02/08  19:04:25  eric
 * Les fonctions arithmetiques ne sont plus amies.
 * Ajout de nouvelles fonctions arithmetiques.
 *
 * Revision 2.22  2000/02/01  15:40:19  eric
 * Ajout de la fonction sqrt
 *
 * Revision 2.21  2000/02/01  14:13:56  eric
 * Modif commentaires.
 * Ajout de la fonction amie flat_scalar_prod.
 *
 * Revision 2.20  2000/01/21  12:48:55  phil
 * changement prototypage de Tenseur::poisson_vect
 *
 * Revision 2.19  2000/01/20  16:02:20  eric
 * Ajout des operator=(double ) et operator=(int ).
 *
 * Revision 2.18  2000/01/20  14:52:08  phil
 * *** empty log message ***
 *
 * Revision 2.17  2000/01/20  14:50:56  phil
 * *** empty log message ***
 *
 * Revision 2.16  2000/01/20  14:39:31  phil
 * *** empty log message ***
 *
 * Revision 2.15  2000/01/20  13:10:56  phil
 * Ajout de Tenseur::poisson_vect (double)
 *
 * Revision 2.14  2000/01/20  11:20:17  phil
 * changement prototypage
 *
 * Revision 2.13  2000/01/20  10:31:30  phil
 * ajout de xksk
 *
 * Revision 2.12  2000/01/14  14:17:47  eric
 * Modif commentaires.
 *
 * Revision 2.11  2000/01/14  14:04:07  eric
 * Ajout de la fonction annule.
 * classe Tenseur: constructeurs pour les classes derivees et fonctions
 *  de gestion memoire (del_t(), ...) declarees protected et non plus
 *  private.
 *
 * Revision 2.10  2000/01/13  14:15:30  eric
 * Modif commentaires.
 *
 * Revision 2.9  2000/01/13  14:10:18  eric
 * Ajout du constructeur par copie d'un Cmp (pour un scalaire)
 * ainsi que l'affectation a un Cmp.
 *
 * Revision 2.8  2000/01/13  13:45:56  eric
 * Ajout du membre p_gradient_spher et des fonctions fait_gradient_spher(),
 *  gradient_spher() pour le calcul du gradient d'un scalaire en
 *  coordonnees spheriques sur la triade spherique associee.
 *
 * Revision 2.7  2000/01/12  13:17:49  eric
 * Les operator::(...) renvoie desormais une reference const sur le c[...]
 * correspondant et non plus un Cmp copie de c[...].
 * (ceci grace a Map::cmp_zero()).
 *
 * Revision 2.6  2000/01/11  11:13:16  eric
 * Changement de nom pour la base vectorielle : base --> triad
 *
 * Revision 2.5  2000/01/10  17:22:26  eric
 *  Modif des #include
 *
 * Revision 2.4  2000/01/10  15:14:58  eric
 * Ajout du membre base (base vectorielle sur laquelle sont definies
 *   les composantes).
 *
 * Revision 2.3  1999/12/09  12:39:39  phil
 * changement prototypage des derivees
 *
 * Revision 2.2  1999/12/07  15:24:17  phil
 * ajout include
 *
 * Revision 2.1  1999/12/03  09:37:11  phil
 * *** empty log message ***
 *
 * Revision 2.0  1999/12/02  17:15:32  phil
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/02  17:13:29  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */

#define COV -1
#define CON +1

#define N_MET_MAX 5

// Headers Lorene 
#include "cmp.h"
#include "itbl.h"
#include "base_vect.h"

class Metrique ;
class Tenseur_sym ;

			//---------------------------------//
			//	class Tenseur		   //
			//---------------------------------//
			

/**
 * Tensor handling. \ingroup (otens)
 * 
 * This class is intended to store the components of a tensorial field in 
 * a specific basis. \a Tensor \a densities can also be stored. A tensor
 * density \f$\tau^{i_1\ldots i_p}_{j_1\ldots j_q}\f$ is defined by:
 * \f$ \tau^{i_1\ldots i_p}_{j_1\ldots j_q} = \gamma^{\frac{n}{2}} 
 * T^{i_1\ldots i_p}_{j_1\ldots j_q}\f$ where \e T is a \e q -covariant
 * \e p -contravariant tensor and \f$\gamma\f$ is the determinant of the 
 * used 3-metric. \e n is called the weight of the tensor density.
 * 
 * All this is \b 3D meaning that the indices go from 0 to 2. Moreover,
 * the components are described in orthonormal bases.
 * 
 * When first constructed, the memory for each component is not allocated.
 * 
 */
class Tenseur { 

    // Data : 
    // -----
    protected:
	const Map* const mp ;	///< Reference mapping
	int valence ;		///< Valence
	
	/** Vectorial basis (triad) with respect to which the tensor
	 *  components are defined. 
	 */
	const Base_vect* triad ; 

	/** Array of size \c valence  contening the type of each index, 
	 * \c COV  for a covariant one and \c CON  for a contravariant one.
	 * 
	 */	
	Itbl type_indice ;	
	
	int n_comp ;	///< Number of components, depending on the symmetry.
	int etat ;  ///< Logical state \c ETATZERO , \c ETATQCQ or \c ETATNONDEF
	Cmp** c ;   ///< The components.
	double poids ; ///< For tensor densities: the weight
	/// For tensor densities: the metric defining the conformal factor
	const Metrique* metric ;      
	
    // Derived data : 
    // ------------
    protected:
	/** Array of pointers on the \c Metrique 's used to calculate 
	 * derivatives members. If no such member has been calculated 
	 * the pointer is set to zero. The array has the size 
	 * \c N_MET_MAX  abd the i-th corresponds to the i-th ones
	 * in \c p_derive_cov and p_derive_con.
	 * 
	 */
	const Metrique** met_depend ;
	
	/** Pointer on the gradient of \c *this .
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	mutable Tenseur* p_gradient ;
	
	/** Pointer on the gradient of \c *this  in a spherical orthonormal
	 * basis (scalar field only).
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	mutable Tenseur* p_gradient_spher ;
	
	/** Array of pointers on the covariant derivatives of \c *this  
	 * with respect to the corresponding metric in \c *met_depend . 
	 * It is set to zero if it has not been calculated yet, or 
	 * for a scalar field.
	 * 
	 */
	Tenseur** p_derive_cov ;

	/** Array of pointers on the contravariant derivatives of \c *this  
	 * with respect to the corresponding metric in \c *met_depend . 
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	Tenseur** p_derive_con ;

	/** Array of pointers on the scalar squares of \c *this  
	 * with respect to the corresponding metric in \c *met_depend . 
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	Tenseur** p_carre_scal ;
   
    // Constructors - Destructor :
    // -------------------------
    protected:
	/// Returns false for a tensor density without a defined metric
	bool verif() const ; 
		
	/** Builds the arrays \c met_depend , \c p_derive_cov , 
	 *  \c p_derive_con and  \c p_carre_scal and fills them with
	 *  null pointers.
	 *
	 */
	void new_der_met() ;
	
    public:
	explicit Tenseur (const Map& map, const Metrique* met = 0x0, 
		     double weight = 0) ; ///< Constructor for a scalar field. 

	/// Constructor for a scalar field and from a \c Cmp . 
	explicit Tenseur (const Cmp& cmp, const Metrique* met = 0x0, 
		     double weight = 0) ; 

	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor
	 * @param tipe  1-D \c Itbl  of size \c valence  containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index 
	 *			\li and so on... 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined 
	 * @param met   for tensor densities only: a pointer on the metric
	 *              defining the conformal factor
	 * @param weight for tensor densities: the weight
	 */
	Tenseur (const Map& map, int val, const Itbl& tipe, 
		 const Base_vect& triad_i, const Metrique* met = 0x0, 
		     double weight = 0) ;

	/** Standard constructor with the triad passed as a pointer.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor
	 * @param tipe  1-D \c Itbl  of size \c valence  containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index
	 *			\li and so on... 
	 * @param triad_i  pointer on the vectorial basis (triad) with respect 
	 *		    to which the tensor components are defined 
	 *		    (can be set to 0x0 for a scalar field)
	 * @param met   for tensor densities only: a pointer on the metric
	 *              defining the conformal factor
	 * @param weight for tensor densities: the weight
	 */
	Tenseur (const Map& map, int val, const Itbl& tipe, 
		 const Base_vect* triad_i, const Metrique* met = 0x0, 
		     double weight = 0) ;

	/** Standard constructor when all the indices are of 
	 *  the same type.
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  the type of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined.
	 * @param met   for tensor densities only: a pointer on the metric
	 *              defining the conformal factor
	 * @param weight for tensor densities: the weight
	 */
	Tenseur (const Map& map, int val, int tipe, const 
		 Base_vect& triad_i, const Metrique* met = 0x0, 
		     double weight = 0) ;

	Tenseur (const Tenseur&) ;  ///< Copy constructor

	/// Constructor from a symmetric tensor.
	explicit Tenseur (const Tenseur_sym&) ;
	
	/** Constructor from a file (see \c  sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been created by 
	 *			    the function \c sauve(FILE*) .
	 * @param met   for tensor densities only: a pointer on the metric
	 *              defining the conformal factor
	 */
	Tenseur (const Map& map, const Base_vect& triad_i, FILE* fich, 
		 const Metrique* met = 0x0) ;

	/** Constructor from a file for a scalar field
	 *  (see \c sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param fich  file which has been created by 
	 *			    the function \c sauve(FILE*) .
	 * @param met   for tensor densities only: a pointer on the metric
	 *              defining the conformal factor
	 */
	Tenseur (const Map& map, FILE* fich, const Metrique* met = 0x0) ;

	
    protected:
	/**
	 * Constructor used by the derived classes.
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  1-D \c Itbl  of size \c valence  containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index 
	 *			\li and so on... 
	 * @param n_comp  the number of components.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 * @param met   for tensor densities only: a pointer on the metric
	 *              defining the conformal factor
	 * @param weight for tensor densities: the weight
	 */	 
	Tenseur (const Map& map, int val, const Itbl& tipe, int n_comp,
		 const Base_vect& triad_i, const Metrique* met = 0x0, 
		     double weight = 0) ;

	/**
	 * Constructor used by the derived classes when all the indices are of 
	 * the same type.
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  the type of the indices.
	 * @param n_comp  the number of components.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 * @param met   for tensor densities only: a pointer on the metric
	 *              defining the conformal factor
	 * @param weight for tensor densities: the weight
	 */
	Tenseur (const Map&, int val, int tipe, int n_comp, 
		 const Base_vect& triad_i, const Metrique* met = 0x0, 
		     double weight = 0) ;


    public: 

	virtual ~Tenseur() ;	///< Destructor
	
    // Memory management
    // -----------------
    protected:
	void del_t() ;	///< Logical destructor

	/**
	 * Logical destructor of the derivatives depending on the i-th
	 * element of \c *met_depend .
	 */	
	void del_derive_met(int i) const ;

	/**
	 * Logical destructor of all the derivatives.
	 */
	void del_derive() const ;
	
	/**
	 * Sets the pointers of the derivatives depending on the i-th
	 * element of \c *met_depend  to zero (as well as that i-th 
	 * element).
	 */
	void set_der_met_0x0(int i) const ;

	/**
	 * Sets the pointers of all the derivatives
	 * to zero.
	 */
	void set_der_0x0() const ;

    // Mutators / assignment
    // ---------------------
    public:
	/**
	 * Sets the logical state to \c ETATNONDEF  (undefined state).
	 * The components are not allocated.
	 */
	void set_etat_nondef() ;
	
	/**
	 * Sets the logical state to \c ETATZERO  (zero state).
	 * The components are not allocated.
	 */
	void set_etat_zero() ;

	/**
	 * Sets the logical state to \c ETATQCQ  (ordinary state).
	 * The components are now allocated and set to \c ETATNONDEF .
	 */
	void set_etat_qcq() ;
	
	/**
	 * Sets the logical state to \c ETATQCQ  (ordinary state)
	 *  and performs the memory allocation of all the 
	 *  elements, down to the \c double  arrays of the \c Tbl s. 
	 *  This function performs in fact recursive calls to 
	 *  \c set_etat_qcq() 
	 *  on each element of the chain \c Tenseur  -> \c Cmp  ->
	 *  \c Valeur  -> \c Mtbl  -> \c Tbl . 
	 */
	void allocate_all() ; 

	/** Sets a new vectorial basis (triad) of decomposition and modifies
	 *  the components accordingly. 
	 */
	void change_triad(const Base_vect& new_triad) ; 
    
	/** Assigns a new vectorial basis (triad) of decomposition. 
	 *  NB: this function modifies only the member \c triad  and
	 *  leave unchanged the components (member \c c ). In order to 
	 *  change them coherently with the new basis, the function 
	 *  \c change_triad(const Base_vect\&)  must be called instead. 
	 */
	void set_triad(const Base_vect& new_triad) ; 
	void set_poids(double weight) ; ///<Sets the weight for a tensor density
	/// Sets the pointer on the metric for a tensor density
	void set_metric(const Metrique& met) ;
    
	/// Assignment to another \c Tenseur 
	virtual void operator=(const Tenseur& tens) ; 
	
	/// Assignment to a \c Cmp  (scalar field only)
	void operator=(const Cmp& field) ; 
	
	 /// Assignment to a \c double  (scalar field only, except for zero)
	void operator=(double ) ;	

	 /// Assignment to a \c int  (scalar field only, except for zero)
	void operator=(int ) ;	

	/// Read/write for a scalar (see also \c operator=(const \c Cmp\&) ). 
	Cmp& set () ;  
	Cmp& set (int) ; ///< Read/write for a vector.
	Cmp& set (int, int) ; ///< Read/write for a tensor of valence 2.
	Cmp& set (int, int, int) ; ///< Read/write for a tensor of valence 3.
	Cmp& set (const Itbl&) ; ///< Read/write in the general case.
	    
	/**
	 * Sets the \c Tenseur  to zero in a given domain.
	 *	@param l [input]  Index of the domain in which the \c Tenseur 
	 *			  will be set (logically) to zero.
	 */
	void annule(int l) ; 

	/**
	 * Sets the \c Tenseur  to zero in several domains.
	 *	@param l_min [input] The \c Tenseur  will be set (logically) 
	 *			     to zero
	 *			     in the domains whose indices are in the range
	 *			     \c [l_min,l_max] .
	 *	@param l_max [input] see the comments for \c l_min .
	 * 
         * Note that \c annule(0,nz-1), where \c nz  is the total number
	 * of domains, is equivalent to \c set_etat_zero() .
         */
	void annule(int l_min, int l_max) ; 

	/**
	 * Set the standard spectal basis of decomposition for each component.
	 * To be used only with \c valence  strictly lower than 3.
	 * 
	 */
	void set_std_base() ; 
	
	void dec_dzpuis() ;	///< dzpuis -= 1 ;
	void inc_dzpuis() ;	///< dzpuis += 1 ;
	void dec2_dzpuis() ;	///< dzpuis -= 2 ;
	void inc2_dzpuis() ;	///< dzpuis += 2 ;
	void mult_r_zec() ; ///< Multiplication by \e r  in the external zone.
	
	/**
	 * Compute \f$\Delta + \lambda \nabla\nabla\f$ of \c *this , \c *this 
	 * being of valence 1.
	 */
	Tenseur inverse_poisson_vect (double lambda) const ;
	
    // Accessors
    // ---------
    public:
	/**
	 * Returns the position in the \c Cmp  1-D array \c c  of a 
	 * component given by its indices.  
	 *
	 * @return position in the \c Cmp  1-D array \c c   
	 * corresponding to the indices given in \c idx . \c idx 
	 * must be a 1-D \c Itbl  of size \c valence , 
	 * each element of which must be 0, 1 or 2, 
	 * corresponding to spatial indices 1, 2 or 3 respectively. 
	 */
	virtual int donne_place (const Itbl& idx) const ;

	/**
	 * Returns the indices of a component given by its position in the 
	 * \c Cmp  1-D array \c c . 
	 *
	 * @return 1-D array of integers (\c Itbl ) of
	 *         size \c valence  giving the value of each index 
	 *	   for the component located at the position \c place 
	 *	   in the \c Cmp  1-D array \c c . 
	 *	   Each element of this \c Itbl  is 0, 1 or 2, which 
	 *	   corresponds to spatial indices 1, 2 or 3 respectively. 
	 * If \c (*this)  is a scalar the function returns an undefined 
	 * \c Itbl .
	 */
	virtual Itbl donne_indices (int place) const ;
	
	/// Returns pointer on the mapping.
	const Map* get_mp() const {return mp ;} ; 

	/** Returns the vectorial basis (triad) on which the components
	 *  are defined.  
	 */
	const Base_vect* get_triad() const {return triad;} ; 
    
	/// Returns the logical state.
	int get_etat() const {return etat ;} ;
 
	///Returns the valence.
	int get_valence() const {return valence ; } ;

	///Returns the number of components. 
	int get_n_comp() const {return n_comp ;} ; 
	
	/**
	 *  Returns the type of the index number \c i . \c i  must be
	 *  strictly lower than \c valence  and obey the following
	 *		      convention: 
	 *			\li \c i  = 0 : first index 
	 *			\li \c i  = 1 : second index
	 *			\li and so on... 
	 * 
	 *  @return COV for a covariant index, CON for a
	 *	    contravariant one. 
	 */
	int get_type_indice (int i) const {return type_indice(i) ;};

	/**
	 * Returns the types of all the indices.
	 * 
	 *  @return 1-D \c Itbl  of size \c valence  containing the type 
	 *  of each index, \c COV  for a covariant one and \c CON  
	 *  for a contravariant one.
	 */
	Itbl get_type_indice () const {return type_indice ; } ;

	///Returns the weight
	double get_poids() const {return poids ; } ; 

	/**
	 * Returns a pointer on the metric defining the conformal factor 
	 * for tensor densities. Otherwise (case of a pure tensor), it
	 * returns 0x0.
	 */
	const Metrique* get_metric() const {return metric ; } ; 
	
	const Cmp& operator()() const ; ///< Read only for a scalar.
	const Cmp& operator()(int) const ; ///< Read only for a vector.
	const Cmp& operator()(int, int) const ; ///< Read only for a tensor of valence 2.
	const Cmp& operator()(int, int, int) const ; ///< Read only for a tensor of valence 3.
	const Cmp& operator()(const Itbl&) const ; ///< Read only in the general case.
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    ///< Save in a file
	friend ostream& operator<<(ostream& , const Tenseur & ) ;
	
    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the gradient of \c *this .
	 * The result is in \c *p_gradient 
	 */
	virtual void fait_gradient () const ;

	/**
	 * Calculates, if needed, the gradient of \c *this  in a 
	 * spherical orthonormal basis (scalar field only). 
	 * The result is in \c *p_gradient_spher 
	 */
	void fait_gradient_spher () const ;

	/**
	 * Calculates, if needed, the covariant derivative of \c *this , 
	 * with respect to \c met .
	 * The result is in \c *p_derive_cov[i] 
	 */
	virtual void fait_derive_cov (const Metrique& met, int i) const ;

	/**
	 * Calculates, if needed, the contravariant derivative of \c *this ,
	 * with respect to \c met .
	 * The result is in \c *p_derive_con[i] 
	 */
	virtual void fait_derive_con (const Metrique&, int i) const ;

	/**
	 * Calculates, if needed, the scalar square of \c *this ,
	 * with respect to \c met .
	 * The result is in \c *p_carre_scal[i] 
	 */
	void fait_carre_scal (const Metrique&, int i) const ;
	
	/**
	 * To be used to describe the fact that the derivatives members have
	 * been calculated with \c met .
	 * 
	 * First it sets a null element of \c met_depend  to 
	 * \c \&met  and puts \c this  in 
	 * the list of the dependancies of \c met .
	 * 
	 */
	void set_dependance (const Metrique& met) const ;

	/**
	 * Returns the position of the pointer on \c metre  in 
	 * the array \c met_depend .
	 *
	 */
	int get_place_met(const Metrique& metre) const ;
	
    // Differential operators
    // ----------------------
    public:
	/// Returns the gradient of \c *this  (Cartesian coordinates)
	const Tenseur& gradient() const ; 
	
	/** Returns the gradient of \c *this  (Spherical coordinates)
	 *	(scalar field only). 
	 */
	const Tenseur& gradient_spher() const ; 
	
	/**
	 * Returns the covariant derivative of \c *this , with respect to 
	 * \c met .
	 */
	const Tenseur& derive_cov (const Metrique& met) const ;

	/**
	 * Returns the contravariant derivative of \c *this , with respect to 
	 * \c met .
	 */
	const Tenseur& derive_con (const Metrique&) const ;	

	/**
	 * Returns the scalar square of \c *this , with respect to 
	 * \c met .
	 */
	const Tenseur& carre_scal (const Metrique&) const ;
    
    // Resolution d'EDP :
    /**
     * Solves the vectorial Poisson equation :
     * \f$\Delta N^i +\lambda \nabla^i \nabla_k N^k = S^i\f$.
     * with \f$\lambda \not= 1\f$.
     * 
     * \c *this  must be given with \c dzpuis  = 4.
     * 
     * It uses the Shibata scheme,  where \f$N^i\f$ is given by :
     * \f[
     *  N^i = \frac{1}{2}\frac{\lambda+2}{\lambda+1}W^i-\frac{1}{2}
     * \frac{\lambda}{\lambda+1}\left(\nabla^i\chi+\nabla^iW^kx_k\right)
     * \f]
     * with \f$\Delta W^i = S^i\f$ and \f$\Delta \chi = -x_kS^k\f$.
     * 
     * @param lambda [input] \f$\lambda\f$.
     * @param par [input/output] see Map::donne_para_poisson_vect.
     * @param shift [input] solution \f$N^i\f$ at the previous step.
     *			Zero if nothing is known.
     * @param shift [output] solution at this step.
     * @param vect [input/output] the same thing than for \c shift  but for 
     * \f$W^i\f$.
     * @param scal [input/output] the same thing than for \c shift  but for 
     * \f$\chi\f$.
     */

    void poisson_vect(double lambda, Param& par, Tenseur& shift, Tenseur& vect
			, Tenseur& scal) const ;
    
	/*
	* Same as poisson_vect with a Tau method
	**/
    void poisson_vect_tau(double lambda, Param& par, Tenseur& shift, Tenseur& vect
			, Tenseur& scal) const ;
			
    void poisson_vect_falloff(double lambda, Param& par, Tenseur& shift, 
			      Tenseur& vect, Tenseur& scal, int* k_falloff) const ;

    void poisson_vect_ylm(double lambda, Param& para, Tenseur& shift,
			  Tenseur& vecteur, Tenseur& scalaire, int nylm,
			  double* intvec) const ;

    /**
     * Solves the vectorial Poisson equation \f$\Delta N^i +\lambda \nabla^i
     * \nabla_k N^k = S^i\f$.
     * with \f$\lambda \not= 1\f$.
     *  
     * \c *this  must be given with \c dzpuis  = 4.
     * 
     * It uses the Shibata scheme,  where \f$N^i\f$ is given by :
     * \f[
     * N^i = \frac{1}{2}\frac{\lambda+2}{\lambda+1}W^i-\frac{1}{2}
     * \frac{\lambda}{\lambda+1}\left(\nabla^i\chi+\nabla^iW^kx_k\right)
     * \f]
     * with \f$\Delta W^i = S^i\f$ and \f$\Delta \chi = -x_kS^k\f$.
     * 
     * This version is to be used only with an affine mapping.
     * 
     * @param lambda [input] \f$\lambda\f$.
     * @param vect [input] \f$W^i\f$ at the previous step.
     *			Zero if nothing is known.
     * @param vect [output] \f$W^i\f$ at this step.
     * @param scal [input/output] the same thing than for \c shift  but for 
     * \f$\chi\f$.
     *
     * @return the solution \f$N^i\f$.
     */
    

    Tenseur poisson_vect(double lambda, Tenseur& vect , Tenseur& scal ) const ;
    
	/*
	* Same as poisson_vect with a Tau method
	**/
    Tenseur poisson_vect_tau(double lambda, Tenseur& vect , Tenseur& scal ) const ;
			
    Tenseur poisson_vect_falloff(double lambda, Tenseur& vect , 
				 Tenseur& scal, int* k_falloff ) const ;

    Tenseur poisson_vect_ylm(double lambda, Tenseur& vecteur, 
			     Tenseur& scalaire, int nylm, double* intvec) const ;

    /**
     * Solves the vectorial Poisson equation \f$\Delta N^i +\lambda \nabla^i
     * \nabla_k N^k = S^i\f$.
     * with \f$\lambda \not= 1\f$.
     *
     * \c *this  must be given with \c dzpuis  = 3 or 4 and be continuous.
     * 
     * It uses the Oohara scheme, where \f$N^i\f$ is given by 
     * \f[
     *  \Delta N^i = S^i-\lambda \nabla^i \chi 
     * \f]
     * with \f$\chi\f$ solution of :
     * \f[
     *   \Delta \chi = \frac{1}{\lambda+1}\nabla_k S^k
     * \f]
     *
     * @param lambda [input] \f$\lambda\f$.
     * @param par [input/output] see Map::donne_para_poisson_vect.
     * @param shift [input] solution \f$N^i\f$ at the previous step.
     *			Zero if nothing is known.
     * @param shift [output] solution at this step.
     * @param scal [input/output] the same thing than for \c shift  but for 
     * \f$\chi\f$.
     */

    void poisson_vect_oohara(double lambda, Param& par, Tenseur& shift, 
				    Tenseur& scal) const ;
       
	/*
	* Same as poisson_vect_oohara with a Tau method
	**/
    void poisson_vect_oohara_tau(double lambda, Param& par, Tenseur& shift, 
				    Tenseur& scal) const ;
			
  /**
     * Solves the vectorial Poisson equation \f$\Delta N^i +\lambda \nabla^i
     * \nabla_k N^k = S^i\f$.
     * with \f$\lambda \not= 1\f$.
     * 
     * \c *this  must be given with \c dzpuis  = 3 or 4 and be continuous.
     *
     * This version is to be used only with an affine mapping.
     *
     * It uses the Oohara scheme,  where \f$N^i\f$ is given by :
     * \f[
     *  \Delta N^i = S^i-\lambda \nabla^i \chi 
     * \f]
     * 
     * with \f$\chi\f$ solution of :
     * \f[
     *  \Delta \chi = \frac{1}{\lambda+1}\nabla_k S^k
     * \f]
     * 
     * This version is to be used only with an affine mapping.
     *
     * @param lambda [input] \f$\lambda\f$.
     * @param scal [input]\f$\chi\f$ at the previous step.
     *			Zero if nothing is known.
     * @param scal [output] \f$\chi\f$ at this step.
     * @return the solution \f$N^i\f$.
     */
			    
    Tenseur poisson_vect_oohara(double lambda, Tenseur& scal) const ;
   /*
	* Same as poisson_vect_oohara with a Tau method
	**/
    Tenseur poisson_vect_oohara_tau(double lambda, Tenseur& scal) const ;
			
    /**
     * Solves the vectorial Poisson equation :
     * \f$\Delta N^i +\lambda \nabla^i \nabla_k N^k = S^i\f$.
     * with \f$\lambda \not= 1\f$ by regularizing the source term.
     * 
     * \c *this  must be given with \c dzpuis  = 4.
     * 
     * It uses the Shibata scheme,  where \f$N^i\f$ is given by :
     * \f[
     *  N^i = \frac{1}{2}\frac{\lambda+2}{\lambda+1}W^i-\frac{1}{2}
     * \frac{\lambda}{\lambda+1}\left(\nabla^i\chi+\nabla^iW^kx_k\right)
     * \f]
     * with \f$\Delta W^i = S^i\f$ and \f$\Delta \chi = -x_kS^k\f$.
     * 
     * @param k_div [input] regularization degree.
     * @param nzet [input] number of domains covering a star.
     * @param unsgam1 [input] \f$1/(\gamma - 1)\f$.
     * @param lambda [input] \f$\lambda\f$.
     * @param par [input/output] see Map::donne_para_poisson_vect.
     * @param shift [input] solution \f$N^i\f$ at the previous step.
     *			Zero if nothing is known.
     * @param shift [output] solution at this step.
     * @param vect [input/output] the same thing than for \c shift but for 
     * \f$W^i\f$.
     * @param scal [input/output] the same thing than for \c shift  but for 
     * \f$\chi\f$.
     */

    void poisson_vect_regu(int k_div, int nzet, double unsgam1,
                           double lambda, Param& par, Tenseur& shift,
                           Tenseur& vect, Tenseur& scal) const ;
   
    
    /// Functions to compare the values of tensors.

    void compare( const Tenseur& tens, const char* name ) ;
    void compare( FILE* fich, const char* name_i ) ;
    


    // Friend classes 
    // ---------------
    friend class Tenseur_sym ;
    friend class Metrique ;
    
    // Mathematical operators
    // ----------------------
    
    friend Tenseur operator* (const Tenseur&, const Tenseur&) ; 
    friend Tenseur operator% (const Tenseur&, const Tenseur&) ; 
    friend Tenseur contract(const Tenseur&, int id1, int id2) ;
    friend Tenseur contract(const Tenseur&, int id1, const Tenseur&, 
			    int id2) ;
    friend Tenseur contract_desal(const Tenseur&, int id1, const Tenseur&, 
			    int id2) ;
    friend Tenseur flat_scalar_prod(const Tenseur& t1, const Tenseur& t2) ;
    friend Tenseur flat_scalar_prod_desal(const Tenseur& t1, 
					  const Tenseur& t2) ;
    friend Tenseur manipule(const Tenseur&, const Metrique&, int idx) ;
    friend Tenseur manipule(const Tenseur&, const Metrique&) ;
    friend Tenseur skxk (const Tenseur&) ;
    friend Tenseur lie_derive(const Tenseur& , const Tenseur& , 
			      const Metrique* ) ;

};


/**
 * \defgroup tens_cal Tenseur calculus 
 *  \ingroup (otens)
 *
 *@{
 */
/// Tensorial product.
Tenseur operator*(const Tenseur&, const Tenseur&) ; 

/// Tensorial product with desaliasing.
Tenseur operator%(const Tenseur&, const Tenseur&) ; 

/**
 * Self contraction of two indices of a \c Tenseur .
 * 
 * The two indices must be of different type, i.e. covariant and
 * contravariant, or contravariant and covariant.
 *
 * @param id1 [input] number of the first index for the contraction;
 *		      \c id1  must be strictly lower than the
 *		      valence of the tensor and obeys the following
 *		      convention: 
 *			\li \c id1  = 0 : first index 
 *			\li \c id1  = 1 : second index 
 *			\li and so on...
 * @param id2 [input] number of the second index for the contraction;
 *		      \c id2  must be strictly lower than the
 *		      valence of the tensor and obeys the following
 *		      convention: 
 *			\li \c id2  = 0 : first index 
 *			\li \c id2  = 1 : second index 
 *			\li and so on...
 * 
 */
Tenseur contract(const Tenseur&, int id1, int id2) ;

/**
 * Contraction of two \c Tenseur .
 * 
 * The two indices must be of different type, i.e. covariant and
 * contravariant, or contravariant and covariant.
 *
 * @param id1 [input] number of the index of contraction for
 *		      the first \c Tenseur ;
 *		      \c id1  must be strictly lower than the
 *		      valence of the tensor and obeys the following
 *		      convention: 
 *			\li \c id1  = 0 : first index 
 *			\li \c id1  = 1 : second index 
 *			\li and so on...
 * @param id2 [input] number of index of contraction for the second one;
 *		      \c id2  must be strictly lower than the
 *		      valence of the tensor and obeys the following
 *		      convention: 
 *			\li \c id2  = 0 : first index 
 *			\li \c id2  = 1 : second index 
 *			\li and so on...
 */	
Tenseur contract(const Tenseur&, int id1, const Tenseur&, int id2) ;

/**
 *  Scalar product of two \c Tenseur  when the metric is
 *  \f$\delta_{ij}\f$: performs the contraction of the 
 *  last index of \c t1  with the first one of \c t2 , irrespective
 *  of the type of these indices. 
 */
Tenseur flat_scalar_prod(const Tenseur& t1, const Tenseur& t2) ;
	
/**
 *  Same as \c flat_scalar_prod  but with desaliasing. 
 */
Tenseur flat_scalar_prod_desal(const Tenseur& t1, const Tenseur& t2) ;
	
/**
 * Raise or lower the index \c idx  depending on its type, using the
 * given \c Metrique .
 */
Tenseur manipule(const Tenseur&, const Metrique&, int idx) ;

/**
 * Raise or lower all the indices, depending on their type,  using the given
 * \c Metrique .
 */
Tenseur manipule(const Tenseur&, const Metrique&) ;
	
/**
 * Contraction of the last index of (*this) with \f$x^k\f$ or \f$x_k\f$, depending
 * on the type of \e S . 
 * 
 * The calculation is performed to avoid singularities in the external 
 * zone. This is done only for a flat metric.
 */
Tenseur skxk (const Tenseur&) ;

/**
 * Lie Derivative of \c t  with respect to \c x . If no other argument
 * is given, it uses partial derivatives with respect to cartesian coordinates
 * to calculate the result (this is the default). Otherwise, it uses the 
 * covariant derivative associated to the metric given as last argument.
 */
Tenseur lie_derive (const Tenseur& t, const Tenseur& x, const Metrique* = 0x0);

/**
 *  Computes the traceless part of a \c Tenseur  of valence 2.
 * 
 * @param tens [input] the \c Tenseur  of valence 2
 * @param metre [input] the metric used to raise or lower the indices
 * 
 * @return The traceless part of the input \c Tenseur 
 */
Tenseur sans_trace(const Tenseur& tens, const Metrique& metre) ;


/** @} */

/**
 * \defgroup tens_ma Tenseur mathematics 
 * \ingroup (otens)
 *
 *  @{
 */
Tenseur operator+(const Tenseur& ) ;			///< + Tenseur
Tenseur operator-(const Tenseur& ) ;			///< \c - Tenseur
Tenseur operator+(const Tenseur&, const Tenseur &) ;	///< Tenseur + Tenseur

/// Tenseur + double (the \c Tenseur  must be a scalar)
Tenseur operator+(const Tenseur&, double ) ;		

/// double + Tenseur (the \c Tenseur  must be a scalar)
Tenseur operator+(double, const Tenseur& ) ;		

/// Tenseur + int (the \c Tenseur  must be a scalar)
Tenseur operator+(const Tenseur&, int ) ;		

/// int + Tenseur (the \c Tenseur  must be a scalar)
Tenseur operator+(int, const Tenseur& ) ;		

Tenseur operator-(const Tenseur &, const Tenseur &) ;	///< Tenseur - Tenseur

/// Tenseur - double (the \c Tenseur  must be a scalar)
Tenseur operator-(const Tenseur&, double ) ;		

/// double - Tenseur (the \c Tenseur  must be a scalar)
Tenseur operator-(double, const Tenseur& ) ;		

/// Tenseur - int (the \c Tenseur  must be a scalar)
Tenseur operator-(const Tenseur&, int ) ;		

/// int - Tenseur (the \c Tenseur  must be a scalar)
Tenseur operator-(int, const Tenseur& ) ;		

/// Tenseur * double 
Tenseur operator*(const Tenseur&, double ) ;		

/// double * Tenseur 
Tenseur operator*(double, const Tenseur& ) ;		

/// Tenseur * int 
Tenseur operator*(const Tenseur&, int ) ;		

/// int * Tenseur 
Tenseur operator*(int, const Tenseur& ) ;		

/// Tenseur / Tenseur (\c b  must be a scalar)
Tenseur operator/(const Tenseur& a, const Tenseur& b) ;	

Tenseur operator/(const Tenseur&, double ) ;	///< Tenseur / double

/// double / Tenseur  (the \c Tenseur  must be a scalar)
Tenseur operator/(double, const Tenseur &) ;		

Tenseur operator/(const Tenseur&, int ) ;		///< Tenseur / int

/// int / Tenseur  (the \c Tenseur  must be a scalar)
Tenseur operator/(int, const Tenseur &) ;		

Tenseur exp(const Tenseur& ) ;		///< Exponential (for a scalar only)
Tenseur log(const Tenseur& ) ;		///< Neperian logarithm (for a scalar only)
Tenseur sqrt(const Tenseur& ) ;		///< Square root (for a scalar only)
Tenseur abs(const Tenseur& ) ;		///< Absolute value (for a scalar only)
Tenseur pow(const Tenseur&, int ) ;	///< Power (for a scalar only)
Tenseur pow(const Tenseur&, double ) ;	///< Power (for a scalar only)

/** @} */





			//---------------------------------//
			//	class Tenseur_sym	   //
			//---------------------------------//
			
/**
 * Class intended to describe tensors with a symmetry on the two last indices.
 * The storage and the calculations are different and quicker than with an 
 * usual \c Tenseur . \ingroup (otens)
 * 
 * The valence must be >1.
 */
class Tenseur_sym : public Tenseur {

    // Constructors - Destructor :
    // -------------------------
	
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor; must be greater or equal to 2.
	 * @param tipe  1-D \c Itbl  of size \c valence  containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index 
	 *			\li and so on... 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */
	Tenseur_sym (const Map& map, int val, const Itbl& tipe, 
		     const Base_vect& triad_i, const Metrique* met = 0x0,
		     double weight = 0) ;

	/** Standard constructor when all the indices are of the same type.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor; must be greater or equal to 2.
	 * @param tipe  the type of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 * 
	 */
	Tenseur_sym (const Map& map, int val, int tipe, 
		     const Base_vect& triad_i, const Metrique* met = 0x0,
		     double weight = 0) ;

	Tenseur_sym (const Tenseur_sym&) ; ///< Copy constructor

	/** Constructor from a \c Tenseur .
	 *  The symmetry is assumed to be true but not checked.
	 */
	explicit Tenseur_sym (const Tenseur&) ;
	
	/** Constructor from a file (see \c sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been created by 
	 *			    the function \c sauve(FILE*) .
	 */
	Tenseur_sym (const Map& map, const Base_vect& triad_i, FILE* fich,
		     const Metrique* met = 0x0) ;

	virtual ~Tenseur_sym() ;    ///< Destructor
	
    // Mutators / assignment
    // ---------------------
    public:
	/**
	 * Assignment from a \c Tenseur .
	 * 
	 * The symmetry is assumed but not checked.
	 */
	virtual void operator= (const Tenseur&) ;
    

    // Accessors
    // ---------
    public:
	/**
	 * Returns the position in the \c Cmp  1-D array \c c  of a 
	 * component given by its indices.  
	 *
	 * @return position in the \c Cmp  1-D array \c c   
	 * corresponding to the indices given in \c idx . \c idx 
	 * must be a 1-D \c Itbl  of size \c valence , 
	 * each element of which must be 0, 1 or 2, 
	 * corresponding to spatial indices 1, 2 or 3 respectively. 
	 */
	virtual int donne_place (const Itbl& idx) const ;

	/**
	 * Returns the indices of a component given by its position in the 
	 * \c Cmp  1-D array \c c . 
	 *
	 * @return 1-D array of integers (\c Itbl ) of
	 *         size \c valence  giving the value of each index 
	 *	   for the component located at the position \c place 
	 *	   in the \c Cmp  1-D array \c c . 
	 *	   Each element of this \c Itbl  is 0, 1 or 2, which 
	 *	   corresponds to spatial indices 1, 2 or 3 respectively. 
	 */
	virtual Itbl donne_indices (int place) const ;
		
    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the gradient of \c *this .
	 * The result is in \c *p_gradient 
	 */
	virtual void fait_gradient () const ;

	/**
	 * Calculates, if needed, the covariant derivative of \c *this , with 
	 * respect to \c met .
	 * The result is in \c *p_derive_cov[i] 
	 */
	virtual void fait_derive_cov (const Metrique& met, int i) const ;

	/**
	 * Calculates, if needed, the contravariant derivative of \c *this ,
	 * with respect to \c met .
	 * The result is in \c *p_derive_con[i] 
	 */
	virtual void fait_derive_con (const Metrique&, int i) const ;
	
    // Mathematical operators
    // ----------------------
	friend Tenseur_sym operator* (const Tenseur&, const Tenseur_sym&) ; 
	friend Tenseur_sym manipule(const Tenseur_sym&, const Metrique&) ;
	friend Tenseur lie_derive (const Tenseur& , const Tenseur& , 
			    const Metrique* );
 
} ;
/**
 * \defgroup tsym_cal Tenseur_sym calculus 
 * \ingroup (otens)
 *
 * @{
 */
/// Tensorial product.
Tenseur_sym operator* (const Tenseur&, const Tenseur_sym&) ; 

/**
 * Raise or lower all the indices, depending on their type,  using the given
 * \c Metrique .
 */
Tenseur_sym manipule(const Tenseur_sym&, const Metrique&) ;

/**
 * Lie Derivative of \c t  with respect to \c x . If no other 
 * argument is given, it uses partial derivatives with respect to 
 * cartesian coordinates to calculate the result (this is the 
 * default). Otherwise, it uses the covariant derivative associated 
 * to the metric given as last argument.
 */
Tenseur_sym lie_derive (const Tenseur_sym& t, const Tenseur& x, 
			    const Metrique* = 0x0);

/**
 *  Computes the traceless part of a \c Tenseur_sym  of valence 2.
 * 
 * @param tens [input] the \c Tenseur_sym  of valence 2
 * @param metre [input] the metric used to raise or lower the indices
 * 
 * @return The traceless part of the input \c Tenseur_sym 
 */
Tenseur_sym sans_trace(const Tenseur_sym& tens, const Metrique& metre) ;

/** @} */



/**
 * \defgroup tsym_mat Tenseur_sym mathematics 
 *  \ingroup (otens)
 *
 * @{
 */
Tenseur_sym operator+(const Tenseur_sym& ) ;	///< + Tenseur_sym
Tenseur_sym operator-(const Tenseur_sym& ) ;	///< \c - Tenseur_sym

/// Tenseur_sym + Tenseur_sym
Tenseur_sym operator+(const Tenseur_sym&, const Tenseur_sym &) ;	

/// Tenseur_sym - Tenseur_sym
Tenseur_sym operator-(const Tenseur_sym &, const Tenseur_sym &) ;	

/// Tenseur_sym * double 
Tenseur_sym operator*(const Tenseur_sym&, double ) ;		

/// double * Tenseur_sym 
Tenseur_sym operator*(double, const Tenseur_sym& ) ;		

/// Tenseur_sym * int 
Tenseur_sym operator*(const Tenseur_sym&, int ) ;		

/// int * Tenseur_sym 
Tenseur_sym operator*(int, const Tenseur_sym& ) ;		

/// Tenseur_sym / Tenseur (\c b  must be a scalar)
Tenseur_sym operator/(const Tenseur_sym& a, const Tenseur& b) ;	

Tenseur_sym operator/(const Tenseur_sym&, double ) ;  ///< Tenseur_sym / double

Tenseur_sym operator/(const Tenseur_sym&, int ) ;     ///< Tenseur_sym / int

/** @} */





#endif
