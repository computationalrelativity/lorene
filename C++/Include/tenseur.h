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

// Headers C++
#include <iostream.h>

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
 * Tensor handling.
 * 
 * This class is intended to store the components of a tensorial field in 
 * a specific basis. {\em Tensor densities} can also be stored. A tensor
 * density $\tau^{i_1\ldots i_p}_{j_1\ldots j_q}$ is defined by:
 * $ \tau^{i_1\ldots i_p}_{j_1\ldots j_q} = \gamma^{\frac{n}{2}} 
 * T^{i_1\ldots i_p}_{j_1\ldots j_q}$ where {\it T} is a {\it q}-covariant
 * {\it p}-contravariant tensor and $\gamma$ is the determinant of the 
 * used 3-metric. {\it n} is called the weight of the tensor density.
 * 
 * All this is {\it 3D} meaning that the indices go from 0 to 2. Moreover,
 * the components are described in orthonormal bases.
 * 
 * When first constructed, the memory for each component is not allocated.
 * 
 * @version #$Id$#
 */
class Tenseur { 

    // Data : 
    // -----
    protected:
	const Map* const mp ;	/// Reference mapping
	int valence ;		/// Valence
	
	/** Vectorial basis (triad) with respect to which the tensor
	 *  components are defined. 
	 */
	const Base_vect* triad ; 

	/** Array of size {\tt valence} contening the type of each index, 
	 * {\tt COV} for a covariant one and {\tt CON} for a contravariant one.
	 * 
	 */	
	Itbl type_indice ;	
	
	int n_comp ;	/// Number of components, depending on the symmetry.
	int etat ;  /// Logical state {\tt (ETATZERO, ETATQCQ or ETATNONDEF)}
	Cmp** c ;   /// The components.
	double poids ; /// For tensor densities: the weight
	/// For tensor densities: the metric defining the conformal factor
	const Metrique* metric ;      
	
    // Derived data : 
    // ------------
    protected:
	/** Array of pointers on the {\tt Metrique}'s used to calculate 
	 * derivatives members. If no such member has been calculated 
	 * the pointer is set to zero. The array has the size 
	 * {\tt N\_MET\_MAX} abd the i-th corresponds to the i-th ones
	 * in {\tt p\_derive\_cov} and {\tt p\_derive\_con}.
	 * 
	 */
	const Metrique** met_depend ;
	
	/** Pointer on the gradient of {\tt *this}.
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	mutable Tenseur* p_gradient ;
	
	/** Pointer on the gradient of {\tt *this} in a spherical orthonormal
	 * basis (scalar field only).
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	mutable Tenseur* p_gradient_spher ;
	
	/** Array of pointers on the covariant derivatives of {\tt *this} 
	 * with respect to the corresponding metric in {\tt *met\_depend}. 
	 * It is set to zero if it has not been calculated yet, or 
	 * for a scalar field.
	 * 
	 */
	Tenseur** p_derive_cov ;

	/** Array of pointers on the contravariant derivatives of {\tt *this} 
	 * with respect to the corresponding metric in {\tt *met\_depend}. 
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	Tenseur** p_derive_con ;

	/** Array of pointers on the scalar squares of {\tt *this} 
	 * with respect to the corresponding metric in {\tt *met\_depend}. 
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	Tenseur** p_carre_scal ;
   
    // Constructors - Destructor :
    // -------------------------
    protected:
	/// Returns false for a tensor density without a defined {\tt metric} 
	bool verif() const ; 
		
	/** Builds the arrays {\tt met\_depend}, {\tt p\_derive\_cov}, 
	 *  {\tt p\_derive\con} and {\tt p\_carre\_scal} and fills them with
	 *  null pointers.
	 *
	 */
	void new_der_met() ;
	
    public:
	explicit Tenseur (const Map& map, const Metrique* met = 0x0, 
		     double weight = 0) ; /// Constructor for a scalar field. 

	/// Constructor for a scalar field and from a {\tt Cmp}. 
	explicit Tenseur (const Cmp& cmp, const Metrique* met = 0x0, 
		     double weight = 0) ; 

	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor
	 * @param tipe  1-D {\tt Itbl} of size {\tt valence} containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index \\
	 *			and so on... 
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
	 * @param tipe  1-D {\tt Itbl} of size {\tt valence} containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index \\
	 *			and so on... 
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

	Tenseur (const Tenseur&) ;  /// Copy constructor

	/// Constructor from a symmetric tensor.
	explicit Tenseur (const Tenseur_sym&) ;
	
	/** Constructor from a file (see {\tt sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been created by 
	 *			    the function {\tt sauve(FILE* )}.
	 * @param met   for tensor densities only: a pointer on the metric
	 *              defining the conformal factor
	 */
	Tenseur (const Map& map, const Base_vect& triad_i, FILE* fich, 
		 const Metrique* met = 0x0) ;

	/** Constructor from a file for a scalar field
	 *  (see {\tt sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param fich  file which has been created by 
	 *			    the function {\tt sauve(FILE* )}.
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
	 * @param tipe  1-D {\tt Itbl} of size {\tt valence} containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index \\
	 *			and so on... 
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

	virtual ~Tenseur() ;	/// Destructor
	
    // Memory management
    // -----------------
    protected:
	void del_t() ;	/// Logical destructor

	/**
	 * Logical destructor of the derivatives depending on the i-th
	 * element of {\tt *met\_depend}.
	 */	
	void del_derive_met(int i) const ;

	/**
	 * Logical destructor of all the derivatives.
	 */
	void del_derive() const ;
	
	/**
	 * Sets the pointers of the derivatives depending on the i-th
	 * element of {\tt *met\_depend} to zero (as well as that i-th 
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
	 * Sets the logical state to {\tt ETATNONDEF} (undefined state).
	 * The components are not allocated.
	 */
	void set_etat_nondef() ;
	
	/**
	 * Sets the logical state to {\tt ETATZERO} (zero state).
	 * The components are not allocated.
	 */
	void set_etat_zero() ;

	/**
	 * Sets the logical state to {\tt ETATQCQ} (ordinary state).
	 * The components are now allocated and set to {\tt ETATNONDEF}.
	 */
	void set_etat_qcq() ;
	
	/**
	 * Sets the logical state to {\tt ETATQCQ} (ordinary state)
	 *  and performs the memory allocation of all the 
	 *  elements, down to the {\tt double} arrays of the {\tt Tbl}s. 
	 *  This function performs in fact recursive calls to 
	 *  {\tt set\_etat\_qcq()}
	 *  on each element of the chain {\tt Tenseur} -> {\tt Cmp} ->
	 *  {\tt Valeur} -> {\tt Mtbl} -> {\tt Tbl}. 
	 */
	void allocate_all() ; 

	/** Sets a new vectorial basis (triad) of decomposition and modifies
	 *  the components accordingly. 
	 */
	void change_triad(const Base_vect& new_triad) ; 
    
	/** Assigns a new vectorial basis (triad) of decomposition. 
	 *  NB: this function modifies only the member {\tt triad} and
	 *  leave unchanged the components (member {\tt c}). In order to 
	 *  change them coherently with the new basis, the function 
	 *  {\tt change\_triad(const Base\_vect\& )} must be called instead. 
	 */
	void set_triad(const Base_vect& new_triad) ; 
	void set_poids(double weight) ; ///Sets the weight for a tensor density
	/// Sets the pointer on the metric for a tensor density
	void set_metric(const Metrique& met) ;
    
	/// Assignment to another {\tt Tenseur}
	virtual void operator=(const Tenseur&) ; 
	
	/// Assignment to a {\tt Cmp} (scalar field only)
	void operator=(const Cmp&) ; 
	
	 /// Assignment to a {\tt double} (scalar field only, except for zero)
	void operator=(double ) ;	

	 /// Assignment to a {\tt int} (scalar field only, except for zero)
	void operator=(int ) ;	

	/// Read/write for a scalar (see also {\tt operator=(const Cmp\&)}). 
	Cmp& set () ;  
	Cmp& set (int) ; /// Read/write for a vector.
	Cmp& set (int, int) ; /// Read/write for a tensor of valence 2.
	Cmp& set (int, int, int) ; /// Read/write for a tensor of valence 3.
	Cmp& set (const Itbl&) ; /// Read/write in the general case.
	    
	/**
	 * Sets the {\tt Tenseur} to zero in a given domain.
	 *	@param l [input]  Index of the domain in which the {\tt Tenseur}
	 *			  will be set (logically) to zero.
	 */
	void annule(int l) ; 

	/**
	 * Sets the {\tt Tenseur} to zero in several domains.
	 *	@param l_min [input] The {\tt Tenseur} will be set (logically) 
	 *			     to zero
	 *			     in the domains whose indices are in the range
	 *			     {\tt [l\_min, l\_max]}.
	 *	@param l_max [input] see the comments for {\tt l\_min}.
	 * 
         * Note that {\tt annule(0, nz-1)}, where {\tt nz} is the total number
	 * of domains, is equivalent to {\tt set\_etat\_zero()}.
         */
	void annule(int l_min, int l_max) ; 

	/**
	 * Set the standard spectal basis of decomposition for each component.
	 * To be used only with {\tt valence} strictly lower than 3.
	 * 
	 */
	void set_std_base() ; 
	
	void dec_dzpuis() ;	/// dzpuis -= 1 ;
	void inc_dzpuis() ;	/// dzpuis += 1 ;
	void dec2_dzpuis() ;	/// dzpuis -= 2 ;
	void inc2_dzpuis() ;	/// dzpuis += 2 ;
	void mult_r_zec() ; /// Multiplication by {\it r} in the external zone.
	
	/**
	 * Compute $\Delta + \lambda \nabla\nabla$ of {\tt *this}, {\tt *this}
	 * being of valence 1.
	 */
	Tenseur inverse_poisson_vect (double lambda) const ;
	
    // Accessors
    // ---------
    public:
	/**
	 * Returns the position in the {\tt Cmp} 1-D array {\tt c} of a 
	 * component given by its indices.  
	 *
	 * @return position in the {\tt Cmp} 1-D array {\tt c}  
	 * corresponding to the indices given in {\tt idx}. {\tt idx}
	 * must be a 1-D {\tt Itbl} of size {\tt valence}, 
	 * each element of which must be 0, 1 or 2, 
	 * corresponding to spatial indices 1, 2 or 3 respectively. 
	 */
	virtual int donne_place (const Itbl& idx) const ;

	/**
	 * Returns the indices of a component given by its position in the 
	 * {\tt Cmp} 1-D array {\tt c}. 
	 *
	 * @return 1-D array of integers ({\tt Itbl}) of
	 *         size {\tt valence} giving the value of each index 
	 *	   for the component located at the position {\tt place}
	 *	   in the {\tt Cmp} 1-D array {\tt c}. 
	 *	   Each element of this {\tt Itbl} is 0, 1 or 2, which 
	 *	   corresponds to spatial indices 1, 2 or 3 respectively. 
	 * If {\tt (*this)} is a scalar the function returns an undefined 
	 * {\tt Itbl}.
	 */
	virtual Itbl donne_indices (int place) const ;
	
	const Map* get_mp() const {return mp ;} ; /// Returns pointer on the mapping.

	/** Returns the vectorial basis (triad) on which the components
	 *  are defined.  
	 */
	const Base_vect* get_triad() const {return triad;} ; 
    
	int get_etat() const {return etat ;} ; /// Returns the logical state.
	int get_valence() const {return valence ; } ; ///Returns the valence.
	int get_n_comp() const {return n_comp ;} ; ///Returns the number of components.
	
	/**
	 *  Returns the type of the index number {\tt i}. {\tt i} must be
	 *  strictly lower than {\tt valence} and obey the following
	 *		      convention: \\
	 *			{\tt i} = 0 : first index \\
	 *			{\tt i} = 1 : second index \\
	 *			and so on... 
	 * 
	 *  @return COV for a covariant index, CON for a
	 *	    contravariant one. 
	 */
	int get_type_indice (int i) const {return type_indice(i) ;};

	/**
	 * Returns the types of all the indices.
	 * 
	 *  @return 1-D {\tt Itbl} of size {\tt valence} containing the type 
	 *  of each index, {\tt COV} for a covariant one and {\tt CON} 
	 *  for a contravariant one.
	 */
	Itbl get_type_indice () const {return type_indice ; } ;

	double get_poids() const {return poids ; } ; ///Returns the weight

	/**
	 * Returns a pointer on the metric defining the conformal factor 
	 * for tensor densities. Otherwise (case of a pure tensor), it
	 * returns 0x0.
	 */
	const Metrique* get_metric() const {return metric ; } ; 
	
	const Cmp& operator()() const ; /// Read only for a scalar.
	const Cmp& operator()(int) const ; /// Read only for a vector.
	const Cmp& operator()(int, int) const ; /// Read only for a tensor of valence 2.
	const Cmp& operator()(int, int, int) const ; /// Read only for a tensor of valence 3.
	const Cmp& operator()(const Itbl&) const ; /// Read only in the general case.
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file

	friend ostream& operator<<(ostream& , const Tenseur & ) ;
	
    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the gradient of {\tt *this}.
	 * The result is in {\tt *p\_gradient}
	 */
	virtual void fait_gradient () const ;

	/**
	 * Calculates, if needed, the gradient of {\tt *this} in a 
	 * spherical orthonormal basis (scalar field only). 
	 * The result is in {\tt *p\_gradient\_spher}
	 */
	void fait_gradient_spher () const ;

	/**
	 * Calculates, if needed, the covariant derivative of {\tt *this}, 
	 * with respect to {\tt met}.
	 * The result is in {\tt *p\_derive\_cov[i]}
	 */
	virtual void fait_derive_cov (const Metrique& met, int i) const ;

	/**
	 * Calculates, if needed, the contravariant derivative of {\tt *this},
	 * with respect to {\tt met}.
	 * The result is in {\tt *p\_derive\_con[i]}
	 */
	virtual void fait_derive_con (const Metrique&, int i) const ;

	/**
	 * Calculates, if needed, the scalar square of {\tt *this},
	 * with respect to {\tt met}.
	 * The result is in {\tt *p\_carre\_scal[i]}
	 */
	void fait_carre_scal (const Metrique&, int i) const ;
	
	/**
	 * To be used to describe the fact that the derivatives members have
	 * been calculated with {\tt met}.
	 * 
	 * First it sets a null element of {\tt met\_depend} to 
	 * {\tt \&met} and puts {\tt this} in 
	 * the list of the dependancies of {\tt met}.
	 * 
	 */
	void set_dependance (const Metrique& met) const ;

	/**
	 * Returns the position of the pointer on {\tt metre} in 
	 * the array {\tt met\_depend}.
	 *
	 */
	int get_place_met(const Metrique& metre) const ;
	
    // Differential operators
    // ----------------------
    public:
	/// Returns the gradient of {\tt *this} (Cartesian coordinates)
	const Tenseur& gradient() const ; 
	
	/** Returns the gradient of {\tt *this} (Spherical coordinates)
	 *	(scalar field only). 
	 */
	const Tenseur& gradient_spher() const ; 
	
	/**
	 * @return the covariant derivative of {\tt *this}, with respect to 
	 * {\tt met}.
	 */
	const Tenseur& derive_cov (const Metrique& met) const ;

	/**
	 * @return the contravariant derivative of {\tt *this}, with respect to 
	 * {\tt met}.
	 */
	const Tenseur& derive_con (const Metrique&) const ;	

	/**
	 * @return the scalar square of {\tt *this}, with respect to 
	 * {\tt met}.
	 */
	const Tenseur& carre_scal (const Metrique&) const ;
    
    // Resolution d'EDP :
    /**
     * Solves the vectorial Poisson equation :
     * $\Delta N^i +\lambda \nabla^i \nabla_k N^k = S^i$.
     * with $\lambda \not= 1$.
     * 
     * {\tt *this} must be given with {\tt dzpuis} = 4.
     * 
     * It uses the Shibata scheme,  where $N^i$ is given by :
     * \begin{equation}
     *  N^i = \frac{1}{2}\frac{\lambda+2}{\lambda+1}W^i-\frac{1}{2}
     * \frac{\lambda}{\lambda+1}\left(\nabla^i\chi+\nabla^iW^kx_k\right)
     * \end{equation}
     * with $\Delta W^i = S^i$ and $\Delta \chi = -x_kS^k$.
     * 
     * @param lambda [input] $\lambda$.
     * @param par [input/output] see Map::donne\_para\_poisson\_vect.
     * @param shift [input] solution $N^i$ at the previous step.
     *			Zero if nothing is known.
     * @param shift [output] solution at this step.
     * @param vect [input/output] the same thing than for {\tt shift} but for 
     * $W^i$.
     * @param scal [input/output] the same thing than for {\tt shift} but for 
     * $\chi$.
     */

    void poisson_vect(double lambda, Param& par, Tenseur& shift, Tenseur& vect
			, Tenseur& scal) const ;
    
    /**
     * Solves the vectorial Poisson equation $\Delta N^i +\lambda \nabla^i
     * \nabla_k N^k = S^i$.
     * with $\lambda \not= 1$.
     *  
     * {\tt *this} must be given with {\tt dzpuis} = 4.
     * 
     * It uses the Shibata scheme,  where $N^i$ is given by :
     * \begin{equation}
     * N^i = \frac{1}{2}\frac{\lambda+2}{\lambda+1}W^i-\frac{1}{2}
     * \frac{\lambda}{\lambda+1}\left(\nabla^i\chi+\nabla^iW^kx_k\right)
     * \end{equation}
     * with $\Delta W^i = S^i$ and $\Delta \chi = -x_kS^k$.
     * 
     * This version is to be used only with an affine mapping.
     * 
     * @param lambda [input] $\lambda$.
     * @param vect [input] $W^i$ at the previous step.
     *			Zero if nothing is known.
     * @param vect [output] $W^i$ at this step.
     * @param scal [input/output] the same thing than for {\tt shift} but for 
     * $\chi$.
     *
     * @return the solution $N^i$.
     */
    

    Tenseur poisson_vect(double lambda, Tenseur& vect , Tenseur& scal ) const ;
    
    /**
     * Solves the vectorial Poisson equation $\Delta N^i +\lambda \nabla^i
     * \nabla_k N^k = S^i$.
     * with $\lambda \not= 1$.
     *
     * {\tt *this} must be given with {\tt dzpuis} = 3 or 4 and be continuous.
     * 
     * It uses the Oohara scheme, where $N^i$ is given by 
     * \begin{equation}
     *  \Delta N^i = S^i-\lambda \nabla^i \chi 
     * \end{equation}
     * with $\chi$ solution of :
     * \begin{equation}
     *   \Delta \chi = \frac{1}{\lambda+1}\nabla_k S^k
     * \end{equation}
     *
     * @param lambda [input] $\lambda$.
     * @param par [input/output] see Map::donne\_para\_poisson\_vect.
     * @param shift [input] solution $N^i$ at the previous step.
     *			Zero if nothing is known.
     * @param shift [output] solution at this step.
     * @param scal [input/output] the same thing than for {\tt shift} but for 
     * $\chi$.
     */

    void poisson_vect_oohara(double lambda, Param& par, Tenseur& shift, 
				    Tenseur& scal) const ;
    
  /**
     * Solves the vectorial Poisson equation $\Delta N^i +\lambda \nabla^i
     * \nabla_k N^k = S^i$.
     * with $\lambda \not= 1$.
     * 
     * {\tt *this} must be given with {\tt dzpuis} = 3 or 4 and be continuous.
     *
     * This version is to be used only with an affine mapping.
     *
     * It uses the Oohara scheme,  where $N^i$ is given by :
     * \begin{equation}
     *  \Delta N^i = S^i-\lambda \nabla^i \chi 
     * \end{equation}
     * 
     * with $\chi$ solution of :
     * \begin{equation}
     *  \Delta \chi = \frac{1}{\lambda+1}\nabla_k S^k
     * \end{equation}
     * 
     * This version is to be used only with an affine mapping.
     *
     * @param lambda [input] $\lambda$.
     * @param scal [input]$\chi$ at the previous step.
     *			Zero if nothing is known.
     * @param scal [output] $\chi$ at this step.
     * @return the solution $N^i$.
     */
			    
    Tenseur poisson_vect_oohara(double lambda, Tenseur& scal) const ;
   
    /**
     * Solves the vectorial Poisson equation :
     * $\Delta N^i +\lambda \nabla^i \nabla_k N^k = S^i$.
     * with $\lambda \not= 1$ by regularizing the source term.
     * 
     * {\tt *this} must be given with {\tt dzpuis} = 4.
     * 
     * It uses the Shibata scheme,  where $N^i$ is given by :
     * \begin{equation}
     *  N^i = \frac{1}{2}\frac{\lambda+2}{\lambda+1}W^i-\frac{1}{2}
     * \frac{\lambda}{\lambda+1}\left(\nabla^i\chi+\nabla^iW^kx_k\right)
     * \end{equation}
     * with $\Delta W^i = S^i$ and $\Delta \chi = -x_kS^k$.
     * 
     * @param k_div [input] regularization degree.
     * @param nzet [input] number of domains covering a star.
     * @param unsgam1 [input] $1/(\gamma - 1)$.
     * @param lambda [input] $\lambda$.
     * @param par [input/output] see Map::donne\_para\_poisson\_vect.
     * @param shift [input] solution $N^i$ at the previous step.
     *			Zero if nothing is known.
     * @param shift [output] solution at this step.
     * @param vect [input/output] the same thing than for {\tt shift} but for 
     * $W^i$.
     * @param scal [input/output] the same thing than for {\tt shift} but for 
     * $\chi$.
     */

    void poisson_vect_regu(int k_div, int nzet, double unsgam1,
                           double lambda, Param& par, Tenseur& shift,
                           Tenseur& vect, Tenseur& scal) const ;

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
 * @name Tenseur calculus
 */
//@{
/// Tensorial product.
Tenseur operator*(const Tenseur&, const Tenseur&) ; 

/// Tensorial product with desaliasing.
Tenseur operator%(const Tenseur&, const Tenseur&) ; 

/**
 * Self contraction of two indices of a {\tt Tenseur}.
 * 
 * The two indices must be of different type, i.e. covariant and
 * contravariant, or contravariant and covariant.
 *
 * @param id1 [input] number of the first index for the contraction;
 *		      {\tt id1} must be strictly lower than the
 *		      valence of the tensor and obeys the following
 *		      convention: \\
 *			{\tt id1} = 0 : first index \\
 *			{\tt id1} = 1 : second index \\
 *			and so on...
 * @param id2 [input] number of the second index for the contraction;
 *		      {\tt id2} must be strictly lower than the
 *		      valence of the tensor and obeys the following
 *		      convention: \\
 *			{\tt id2} = 0 : first index \\
 *			{\tt id2} = 1 : second index \\
 *			and so on...
 * 
 */
Tenseur contract(const Tenseur&, int id1, int id2) ;

/**
 * Contraction of two {\tt Tenseur}.
 * 
 * The two indices must be of different type, i.e. covariant and
 * contravariant, or contravariant and covariant.
 *
 * @param id1 [input] number of the index of contraction for
 *		      the first {\tt Tenseur};
 *		      {\tt id1} must be strictly lower than the
 *		      valence of the tensor and obeys the following
 *		      convention: \\
 *			{\tt id1} = 0 : first index \\
 *			{\tt id1} = 1 : second index \\
 *			and so on...
 * @param id2 [input] number of index of contraction for the second one;
 *		      {\tt id2} must be strictly lower than the
 *		      valence of the tensor and obeys the following
 *		      convention: \\
 *			{\tt id2} = 0 : first index \\
 *			{\tt id2} = 1 : second index \\
 *			and so on...
 */	
Tenseur contract(const Tenseur&, int id1, const Tenseur&, int id2) ;

/**
 *  Scalar product of two {\tt Tenseur} when the metric is
 *  $\delta_{ij}$: performs the contraction of the 
 *  last index of {\tt t1} with the first one of {\tt t2}, irrespective
 *  of the type of these indices. 
 */
Tenseur flat_scalar_prod(const Tenseur& t1, const Tenseur& t2) ;
	
/**
 *  Same as {\tt flat\_scalar\_prod} but with desaliasing. 
 */
Tenseur flat_scalar_prod_desal(const Tenseur& t1, const Tenseur& t2) ;
	
/**
 * Raise or lower the index {\tt idx} depending on its type, using the
 * given {\tt Metrique}.
 */
Tenseur manipule(const Tenseur&, const Metrique&, int idx) ;

/**
 * Raise or lower all the indices, depending on their type,  using the given
 * {\tt Metrique}.
 */
Tenseur manipule(const Tenseur&, const Metrique&) ;
	
/**
 * Contraction of the last index of (*this) with $x^k$ or $x_k$, depending
 * on the type of {\it S}. 
 * 
 * The calculation is performed to avoid singularities in the external 
 * zone. This is done only for a flat metric.
 */
Tenseur skxk (const Tenseur&) ;

/**
 * Lie Derivative of {\tt t} with respect to {\tt x}. If no other argument
 * is given, it uses partial derivatives with respect to cartesian coordinates
 * to calculate the result (this is the default). Otherwise, it uses the 
 * covariant derivative associated to the metric given as last argument.
 */
Tenseur lie_derive (const Tenseur& t, const Tenseur& x, const Metrique* = 0x0);


//@}

/**
 * @name Tenseur mathematics
 */
    //@{
Tenseur operator+(const Tenseur& ) ;			/// + Tenseur
Tenseur operator-(const Tenseur& ) ;			/// - Tenseur
Tenseur operator+(const Tenseur&, const Tenseur &) ;	/// Tenseur + Tenseur

/// Tenseur + double (the {\tt Tenseur} must be a scalar)
Tenseur operator+(const Tenseur&, double ) ;		

/// double + Tenseur (the {\tt Tenseur} must be a scalar)
Tenseur operator+(double, const Tenseur& ) ;		

/// Tenseur + int (the {\tt Tenseur} must be a scalar)
Tenseur operator+(const Tenseur&, int ) ;		

/// int + Tenseur (the {\tt Tenseur} must be a scalar)
Tenseur operator+(int, const Tenseur& ) ;		

Tenseur operator-(const Tenseur &, const Tenseur &) ;	/// Tenseur - Tenseur

/// Tenseur - double (the {\tt Tenseur} must be a scalar)
Tenseur operator-(const Tenseur&, double ) ;		

/// double - Tenseur (the {\tt Tenseur} must be a scalar)
Tenseur operator-(double, const Tenseur& ) ;		

/// Tenseur - int (the {\tt Tenseur} must be a scalar)
Tenseur operator-(const Tenseur&, int ) ;		

/// int - Tenseur (the {\tt Tenseur} must be a scalar)
Tenseur operator-(int, const Tenseur& ) ;		

/// Tenseur * double 
Tenseur operator*(const Tenseur&, double ) ;		

/// double * Tenseur 
Tenseur operator*(double, const Tenseur& ) ;		

/// Tenseur * int 
Tenseur operator*(const Tenseur&, int ) ;		

/// int * Tenseur 
Tenseur operator*(int, const Tenseur& ) ;		

/// Tenseur / Tenseur ({\tt b} must be a scalar)
Tenseur operator/(const Tenseur& a, const Tenseur& b) ;	

Tenseur operator/(const Tenseur&, double ) ;	/// Tenseur / double

/// double / Tenseur  (the {\tt Tenseur} must be a scalar)
Tenseur operator/(double, const Tenseur &) ;		

Tenseur operator/(const Tenseur&, int ) ;		/// Tenseur / int

/// int / Tenseur  (the {\tt Tenseur} must be a scalar)
Tenseur operator/(int, const Tenseur &) ;		

Tenseur exp(const Tenseur& ) ;		/// Exponential (for a scalar only)
Tenseur log(const Tenseur& ) ;		/// Neperian logarithm (for a scalar only)
Tenseur sqrt(const Tenseur& ) ;		/// Square root (for a scalar only)
Tenseur abs(const Tenseur& ) ;		/// Absolute value (for a scalar only)
Tenseur pow(const Tenseur&, int ) ;	/// Power (for a scalar only)
Tenseur pow(const Tenseur&, double ) ;	/// Power (for a scalar only)

    //@}





			//---------------------------------//
			//	class Tenseur_sym	   //
			//---------------------------------//
			
/**
 * Class intended to describe tensors with a symmetry on the two last indices.
 * The storage and the calculations are different and quicker than with an 
 * usual {\tt Tenseur}.
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
	 * @param tipe  1-D {\tt Itbl} of size {\tt valence} containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index \\
	 *			and so on... 
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

	Tenseur_sym (const Tenseur_sym&) ; /// Copy constructor

	/** Constructor from a {\tt Tenseur}.
	 *  The symmetry is assumed to be true but not checked.
	 */
	explicit Tenseur_sym (const Tenseur&) ;
	
	/** Constructor from a file (see {\tt sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been created by 
	 *			    the function {\tt sauve(FILE* )}.
	 */
	Tenseur_sym (const Map& map, const Base_vect& triad_i, FILE* fich,
		     const Metrique* met = 0x0) ;

	virtual ~Tenseur_sym() ;    /// Destructor
	
    // Mutators / assignment
    // ---------------------
    public:
	/**
	 * Assignment from a {\tt Tenseur}.
	 * 
	 * The symmetry is assumed but not checked.
	 */
	virtual void operator= (const Tenseur&) ;
    

    // Accessors
    // ---------
    public:
	/**
	 * Returns the position in the {\tt Cmp} 1-D array {\tt c} of a 
	 * component given by its indices.  
	 *
	 * @return position in the {\tt Cmp} 1-D array {\tt c}  
	 * corresponding to the indices given in {\tt idx}. {\tt idx}
	 * must be a 1-D {\tt Itbl} of size {\tt valence}, 
	 * each element of which must be 0, 1 or 2, 
	 * corresponding to spatial indices 1, 2 or 3 respectively. 
	 */
	virtual int donne_place (const Itbl& idx) const ;

	/**
	 * Returns the indices of a component given by its position in the 
	 * {\tt Cmp} 1-D array {\tt c}. 
	 *
	 * @return 1-D array of integers ({\tt Itbl}) of
	 *         size {\tt valence} giving the value of each index 
	 *	   for the component located at the position {\tt place}
	 *	   in the {\tt Cmp} 1-D array {\tt c}. 
	 *	   Each element of this {\tt Itbl} is 0, 1 or 2, which 
	 *	   corresponds to spatial indices 1, 2 or 3 respectively. 
	 */
	virtual Itbl donne_indices (int place) const ;
		
    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the gradient of {\tt *this}.
	 * The result is in {\tt *p\_gradient}
	 */
	virtual void fait_gradient () const ;

	/**
	 * Calculates, if needed, the covariant derivative of {\tt *this}, with 
	 * respect to {\tt met}.
	 * The result is in {\tt *p\_derive\_cov[i]}
	 */
	virtual void fait_derive_cov (const Metrique& met, int i) const ;

	/**
	 * Calculates, if needed, the contravariant derivative of {\tt *this},
	 * with respect to {\tt met}.
	 * The result is in {\tt *p\_derive\_con[i]}
	 */
	virtual void fait_derive_con (const Metrique&, int i) const ;
	
    // Mathematical operators
    // ----------------------
	friend Tenseur_sym operator* (const Tenseur&, const Tenseur_sym&) ; 
	friend Tenseur_sym manipule(const Tenseur_sym&, const Metrique&, 
				    int idx) ;   
	friend Tenseur_sym manipule(const Tenseur_sym&, const Metrique&) ;
	friend Tenseur lie_derive (const Tenseur& , const Tenseur& , 
			    const Metrique* );
 
} ;
/**
 * @name Tenseur\_sym calculus
 */
//@{
/// Tensorial product.
Tenseur_sym operator* (const Tenseur&, const Tenseur_sym&) ; 

/**
 * Raise or lower the index {\tt idx} depending on its type, using the
 * given {\tt Metrique}.
 */
Tenseur_sym manipule(const Tenseur_sym&, const Metrique&, int idx) ;

/**
 * Raise or lower all the indices, depending on their type,  using the given
 * {\tt Metrique}.
 */
Tenseur_sym manipule(const Tenseur_sym&, const Metrique&) ;

/**
 * Lie Derivative of {\tt t} with respect to {\tt x}. If no other 
 * argument is given, it uses partial derivatives with respect to 
 * cartesian coordinates to calculate the result (this is the 
 * default). Otherwise, it uses the covariant derivative associated 
 * to the metric given as last argument.
 */
Tenseur_sym lie_derive (const Tenseur_sym& t, const Tenseur& x, 
			    const Metrique* = 0x0);

//@}



/**
 * @name Tenseur\_sym mathematics
 */
    //@{
Tenseur_sym operator+(const Tenseur_sym& ) ;	/// + Tenseur\_sym
Tenseur_sym operator-(const Tenseur_sym& ) ;	/// - Tenseur\_sym

/// Tenseur\_sym + Tenseur\_sym
Tenseur_sym operator+(const Tenseur_sym&, const Tenseur_sym &) ;	

/// Tenseur\_sym - Tenseur\_sym
Tenseur_sym operator-(const Tenseur_sym &, const Tenseur_sym &) ;	

/// Tenseur\_sym * double 
Tenseur_sym operator*(const Tenseur_sym&, double ) ;		

/// double * Tenseur\_sym 
Tenseur_sym operator*(double, const Tenseur_sym& ) ;		

/// Tenseur\_sym * int 
Tenseur_sym operator*(const Tenseur_sym&, int ) ;		

/// int * Tenseur\_sym 
Tenseur_sym operator*(int, const Tenseur_sym& ) ;		

/// Tenseur\_sym / Tenseur ({\tt b} must be a scalar)
Tenseur_sym operator/(const Tenseur_sym& a, const Tenseur& b) ;	

Tenseur_sym operator/(const Tenseur_sym&, double ) ;  /// Tenseur\_sym / double

Tenseur_sym operator/(const Tenseur_sym&, int ) ;     /// Tenseur\_sym / int

    //@}





#endif
