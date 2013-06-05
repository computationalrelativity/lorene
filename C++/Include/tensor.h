/*
 *  Definition of Lorene classes Tensor and Sym_tensor
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (for preceding class Tenseur)
 *   Copyright (c) 2000-2001 Eric Gourgoulhon      (for preceding class Tenseur)
 *   Copyright (c) 2002 Jerome Novak               (for preceding class Tenseur)
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


#ifndef __TENSOR_H_ 
#define __TENSOR_H_ 


/*
 * $Id$
 * $Log$
 * Revision 1.60  2013/06/05 15:43:49  j_novak
 * Suppression of leg_spectral_base()
 *
 * Revision 1.59  2013/01/11 15:44:54  j_novak
 * Addition of Legendre bases (part 2).
 *
 * Revision 1.58  2008/12/05 08:44:02  j_novak
 * New flag to control the "verbosity" of maxabs.
 *
 * Revision 1.57  2007/12/21 16:06:16  j_novak
 * Methods to filter Tensor, Vector and Sym_tensor objects.
 *
 * Revision 1.56  2006/06/07 14:08:58  j_novak
 * New methods set_index_type( / int).
 *
 * Revision 1.55  2005/10/25 08:56:34  p_grandclement
 * addition of std_spectral_base in the case of odd functions near the origin
 *
 * Revision 1.54  2004/07/08 12:21:51  j_novak
 * Replaced tensor::annule_extern_c2 with tensor::annule_extern_cn for a
 * more general transition.
 *
 * Revision 1.53  2004/06/17 06:54:23  e_gourgoulhon
 * Added method annule_extern_c2.
 *
 * Revision 1.52  2004/05/13 21:29:27  e_gourgoulhon
 * Added (external) functions central_value, max_all_domains,
 * min_all_domains and maxabs_all_domains.
 *
 * Revision 1.51  2004/03/24 14:53:39  j_novak
 * Double declarations suppressed
 *
 * Revision 1.50  2004/03/22 13:12:43  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.49  2004/02/27 21:12:44  e_gourgoulhon
 * Suppressed function contract_desal (since contract has now the
 * boolean argument "desaliasing").
 *
 * Revision 1.48  2004/02/26 22:44:37  e_gourgoulhon
 * -- constructor of Tensor from Map is now declared explicit.
 * -- class Tensor: added methods compute_derive_lie and derive_lie
 * -- class Tensor_sym: added methods derive_cov, derive_con and derive_lie.
 *
 * Revision 1.47  2004/02/19 22:08:51  e_gourgoulhon
 * Added argument "comment" in method spectral_display,
 * as well as in external functions min, max, maxabs, etc...
 *
 * Revision 1.46  2004/02/18 18:42:41  e_gourgoulhon
 * -- Added methods trace.
 * -- Method scontract suppressed ( since it is the same as trace(int, int) ).
 *
 * Revision 1.45  2004/02/18 15:52:39  e_gourgoulhon
 * -- Added optional argument desaliasing in function contract.
 * -- Added new function contract for double contraction.
 *
 * Revision 1.44  2004/02/16 10:48:06  e_gourgoulhon
 * Added "class Tensor_sym;" at the beginning.
 *
 * Revision 1.43  2004/02/15 21:53:48  e_gourgoulhon
 * Modif. comments: suppressed the mention *** under development ***.
 *
 * Revision 1.42  2004/01/30 12:44:17  e_gourgoulhon
 * Added Tensor_sym operator*(const Tensor_sym&, const Tensor_sym& ).
 *
 * Revision 1.41  2004/01/27 13:05:10  j_novak
 * Removed the method Tensor::mult_r_ced()
 *
 * Revision 1.40  2004/01/19 16:31:40  e_gourgoulhon
 * Added operator()(int, int, int, int) and set(int, int, int, int)
 * for direct access to components of valence 4 tensors.
 *
 * Revision 1.39  2004/01/15 11:09:27  f_limousin
 * Modif in method contract_desal
 *
 * Revision 1.38  2004/01/15 11:00:44  f_limousin
 * Added method contract_desal for the contraction of two tensors with desaliasing
 *
 * Revision 1.37  2004/01/14 11:39:00  f_limousin
 * Added method contract for one tensor
 *
 * Revision 1.36  2004/01/08 09:21:39  e_gourgoulhon
 * Added arithmetics of Tensor_sym.
 * Added arithmetics with Scalar (to solve some ambiguities with respect
 * to the Scalar arithmetics).
 * Added Tensor_sym tensorial product.
 *
 * Revision 1.35  2004/01/04 20:47:37  e_gourgoulhon
 * -- Introduction of new derived class Tensor_sym to store tensor with
 *    two symmetric indices
 * -- Suppression of class Tensor_delta (now a special case of Tensor_sym).
 *
 * Revision 1.34  2003/12/27 14:58:01  e_gourgoulhon
 * Improved documentation. In particular, better description of methods
 * derive_cov(), derive_con() and divergence(), taking into account the
 * new index convention for covariant derivatives.
 *
 * Revision 1.33  2003/12/05 16:41:05  f_limousin
 * Added method operator*
 *
 * Revision 1.32  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.31  2003/11/05 15:25:57  e_gourgoulhon
 * Added declaration of external functions:
 * max, min, maxabs, diffrel and diffrelmax.
 *
 * Revision 1.30  2003/11/03 10:58:00  j_novak
 * Suppressed the constructor from a Sym_tensor.
 *
 * Revision 1.29  2003/10/29 11:00:42  e_gourgoulhon
 * Virtual functions dec_dzpuis and inc_dzpuis have now an integer argument to
 *  specify by which amount dzpuis is to be increased.
 * Accordingly virtual methods dec2_dzpuis and inc2_dzpuis have been suppressed.
 *
 * Revision 1.28  2003/10/28 21:21:50  e_gourgoulhon
 * Member function Tensor::contract(int, int) renamed
 *  Tensor::scontract(int, int) in order not to mask
 * the non-member function contract.
 *
 * Revision 1.27  2003/10/27 10:44:00  e_gourgoulhon
 * Declaration of class Sym_tensor is now in file sym_tensor.h.
 *
 * Revision 1.26  2003/10/24 15:00:19  j_novak
 * Forgotten Class declaration... thanks IBM aix!
 *
 * Revision 1.25  2003/10/20 14:26:02  j_novak
 * New assignement operators.
 *
 * Revision 1.24  2003/10/20 09:32:10  j_novak
 * Members p_potential and p_div_free of the Helmholtz decomposition
 * + the method decompose_div(Metric).
 *
 * Revision 1.23  2003/10/19 19:47:31  e_gourgoulhon
 * Introduced new virtual method spectral_display.
 *
 * Revision 1.22  2003/10/16 15:24:30  e_gourgoulhon
 * Name of method annule(int ) changed to annule_domain(int ).
 *
 * Revision 1.21  2003/10/16 14:21:33  j_novak
 * The calculation of the divergence of a Tensor is now possible.
 *
 * Revision 1.20  2003/10/13 13:52:39  j_novak
 * Better managment of derived quantities.
 *
 * Revision 1.19  2003/10/08 14:24:08  j_novak
 * replaced mult_r_zec with mult_r_ced
 *
 * Revision 1.18  2003/10/06 20:48:23  e_gourgoulhon
 * Added methods down and up_down.
 *
 * Revision 1.17  2003/10/06 16:17:29  j_novak
 * Calculation of contravariant derivative and Ricci scalar.
 *
 * Revision 1.16  2003/10/06 15:12:56  e_gourgoulhon
 * Added tensor contraction and raising of index.
 *
 * Revision 1.15  2003/10/06 13:58:45  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.14  2003/10/05 21:07:27  e_gourgoulhon
 * Method std_spectral_base() is now virtual.
 *
 * Revision 1.13  2003/10/03 11:21:45  j_novak
 * More methods for the class Metric
 *
 * Revision 1.12  2003/10/02 15:45:48  j_novak
 * New class Metric
 *
 * Revision 1.11  2003/10/01 15:41:14  e_gourgoulhon
 * class name Delta changed to Tensor_delta.
 *
 * Revision 1.10  2003/10/01 13:03:52  e_gourgoulhon
 * The method get_mp() returns now a reference (and not a pointer)
 * onto a mapping.
 *
 * Revision 1.9  2003/09/29 13:48:17  j_novak
 * New class Delta.
 *
 * Revision 1.8  2003/09/26 14:33:51  j_novak
 * Arithmetic functions for the class Tensor
 *
 * Revision 1.7  2003/09/26 08:05:29  j_novak
 * New class Vector.
 *
 * Revision 1.6  2003/09/25 21:01:50  e_gourgoulhon
 * Improved comments.
 *
 * Revision 1.5  2003/09/25 13:37:38  j_novak
 * Symmetric tensors of valence 2 are now implemented (not tested yet).
 *
 * Revision 1.4  2003/09/24 15:10:54  j_novak
 * Suppression of the etat flag in class Tensor (still present in Scalar)
 *
 * Revision 1.3  2003/09/24 08:46:31  j_novak
 * Added tensor.h and scalar.h to the documentation
 *
 * Revision 1.2  2003/09/23 08:53:11  e_gourgoulhon
 * not ready yet
 *
 * Revision 1.1  2003/09/22 12:50:47  e_gourgoulhon
 * First version: not ready yet!
 *
 *
 * $Header$
 *
 */

#define COV -1
#define CON +1

#define N_MET_MAX 5

// Headers Lorene 
#include "itbl.h"
#include "base_vect.h"
#include "map.h"

class Scalar ;
class Vector ; 
class Tensor_sym ;
class Sym_tensor ;
class Metric ;

			//-------------------------//
			//       class Tensor      //
			//-------------------------//
			

/**
 * Tensor handling. \ingroup (tensor)
 *
 * This class has been devised to replace \c Tenseur  and \c Cmp  (the
 *  latter via the derived class \c Scalar ).
 * 
 * The \c Tensor  class is intended to store the components of a tensorial 
 * field with respect to a specific basis (triad).  
 * 
 * All this is \e 3D meaning that the indices go from 1 to 3. 
 * 
 * 
 */
class Tensor { 

    // Data : 
    // -----
    protected:
	
        /// Mapping on which the numerical values at the grid points are defined
	const Map* const mp ;	

    /// Valence of the tensor (0 = scalar, 1 = vector, etc...)
	int valence ;	
        	
	/** Vectorial basis (triad) with respect to which the tensor
	 *  components are defined. 
	 */
	const Base_vect* triad ; 

	/** 1D array of integers (class \c Itbl ) of size \c valence  
	 *  containing the type of each index: 
	 *  \c COV  for a covariant one and \c CON  for a contravariant one.
	 * 
	 */	
	Itbl type_indice ;	
	
	int n_comp ;	///< Number of stored components, depending on the symmetry.

	/// Array of size \c n_comp  of pointers onto the components.
	Scalar** cmp ;   


    // Derived data : 
    // ------------
     protected:
	/**
	 * Array on the \c Metric 's which were used to compute derived
	 * quantities, like \c p_derive_cov , etc... 
	 * The i-th element of this array is the \c Metric  used to
	 * compute the i-th element of \c p_derive_cov , etc..
	 */
	mutable const Metric* met_depend[N_MET_MAX] ; 

	/** Array of pointers on the covariant derivatives of \c this 
         * with respect to various metrics.
	 * See the comments of \c met_depend . See also the comments
         * of method \c derive_cov()  for the index convention of the
         * covariant derivation.  
	 */
	mutable Tensor* p_derive_cov[N_MET_MAX];
	
	/** Array of pointers on the contravariant derivatives of \c this 
         * with respect to various metrics.
	 * See the comments of \c met_depend . See also the comments
         * of method \c derive_con()  for a precise definition of a
         * "contravariant" derivative.
 	 */
	mutable Tensor* p_derive_con[N_MET_MAX];

	/** Array of pointers on the divergence of \c this 
         * with respect to various metrics.
	 * See the comments of \c met_depend . See also the comments
         * of method \c divergence()  for a precise definition of a
         * the divergence with respect to a given metric.
	 */
	mutable Tensor* p_divergence[N_MET_MAX];

   
    // Constructors - Destructor :
    // -------------------------
	
	public: 

	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor
	 * @param tipe  1-D array of integers (class \c Itbl ) 
	 *		of size \c valence  containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index 
	 *			\li and so on... 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined 
	 */
	Tensor(const Map& map, int val, const Itbl& tipe, 
		 	const Base_vect& triad_i) ;

	/** Standard constructor with the triad passed as a pointer.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor
	 * @param tipe  1-D array of integers (class \c Itbl ) 
	 *		of size \c valence  containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index 
	 *			\li and so on... 
	 * @param triad_i  pointer on the vectorial basis (triad) with respect 
	 *		    to which the tensor components are defined 
	 *		    (can be set to 0x0 for a scalar field)
	 */
	Tensor(const Map& map, int val, const Itbl& tipe, 
		 	const Base_vect* triad_i) ;

	/** Standard constructor when all the indices are of 
	 *  the same type.
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  the type (\c COV  or \c CON ) of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined.
	 */
	Tensor(const Map& map, int val, int tipe, 
			const Base_vect& triad_i) ;

	Tensor(const Tensor&) ;  ///< Copy constructor

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
	Tensor(const Map& map, const Base_vect& triad_i, FILE* fich) ;

    protected:
	/**
	 *  Constructor for a scalar field: to be used only by the derived
	 *  class \c Scalar .
	 *
	 */
	 explicit Tensor(const Map& map) ;

	/**
	 * Constructor to be used by derived classes, with symmetries among
	 *  the components. The number of independent components is
	 *  given as an argument (\c n_comp_i ), and not computed
	 *	from the valence, as in the standard constructor.  
	 *
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  1-D array of integers (class \c Itbl ) 
	 *		of size \c valence  containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index 
	 *			\li and so on... 
	 * @param n_comp_i number of components to be stored
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */	 
	Tensor(const Map& map, int val, const Itbl& tipe, int n_comp_i,
		    const Base_vect& triad_i) ;

	/**
	 * Constructor used by derived classes, with symmetries among
	 *  the components, when all the indices are of 
	 *  the same type. The number of independent components is
	 *  given as a argument (\c n_comp_i ), and not computed
	 *	from the valence, as in the standard constructor.  
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  the type of the indices.
	 * @param n_comp_i  number of components to be stored
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */
	Tensor(const Map& map, int val, int tipe, int n_comp_i, 
		 const Base_vect& triad_i) ;


    public: 

	virtual ~Tensor() ;	///< Destructor
	
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const ;	///< Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 

	/**
	 * Logical destructor of the derivatives depending on the i-th
	 * element of \c met_depend .
	 */	
	virtual void del_derive_met(int) const ;

	/**
	 * Sets all the i-th components of \c met_depend , 
	 * \c p_derive_cov , etc... to 0x0.
	 */
	void set_der_met_0x0(int) const ;

	/**
	 * To be used to describe the fact that the derivatives members have
	 * been calculated with \c met .
	 * 
	 * First it sets a null element of \c met_depend  to 
	 * \c \&met  and puts \c this  in 
	 * the list of the dependancies of \c met .
	 * 
	 */
	void set_dependance (const Metric&) const ;

	/**
	 * Returns the position of the pointer on \c metre  in 
	 * the array \c met_depend .
	 *
	 */
	int get_place_met(const Metric&) const ;
	
    // Mutators / assignment
    // ---------------------
    public:
	/**
	 * Sets the logical state of all components to \c ETATNONDEF  
	 * (undefined state).
	 */
	virtual void set_etat_nondef() ;
	
	/**
	 * Sets the logical state of all components to \c ETATZERO  
	 *(zero state).
	 */
	virtual void set_etat_zero() ;

	/**
	 * Sets the logical state of all components to \c ETATQCQ  
	 * (ordinary state).
	 */
	virtual void set_etat_qcq() ;
	
	/**
	 *  Performs the memory allocation of all the 
	 *  elements, down to the \c double  arrays of the \c Tbl s. 
	 *  This function performs in fact recursive calls to 
	 *  \c set_etat_qcq() 
	 *  on each element of the chain \c Scalar  ->
	 *  \c Valeur  -> \c Mtbl  -> \c Tbl . 
	 */
	virtual void allocate_all() ; 

	/** Sets a new vectorial basis (triad) of decomposition and modifies
	 *  the components accordingly. 
	 */
	virtual void change_triad(const Base_vect& new_triad) ; 
    
	/** Assigns a new vectorial basis (triad) of decomposition. 
	 *  NB: this function modifies only the member \c triad  and
	 *  leave unchanged the components (member \c cmp ). In order to 
	 *  change them coherently with the new basis, the function 
	 *  \c change_triad(const Base_vect\& ) must be called instead. 
	 */
	void set_triad(const Base_vect& new_triad) ; 

	
	virtual void operator=(const Tensor&) ;///< Assignment to a \c Tensor 
	
	/** Returns the value of a component (read/write version).
	 *
	 * @param ind  1-D \c Itbl  of size \c valence  containing the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: 
	 *			\li \c ind(0)  : value of the first index (1, 2 or 3) 
	 *			\li \c ind(1)  : value of the second index (1, 2 or 3) 
	 *			\li and so on... 
	 * @return modifiable reference on the component specified by \c ind 
	 *
	 */
	Scalar& set(const Itbl& ind) ; 
	
	/** Returns the value of a component for a tensor of valence 2
	 *  (read/write version).
	 *
	 * @param i1  value of the first index (1, 2 or 3)
	 * @param i2  value of the second index (1, 2 or 3)
	 *
	 * @return modifiable reference on the component specified by \c (i1,i2) 
	 *
	 */
	Scalar& set(int i1, int i2) ; 
	
	
	/** Returns the value of a component for a tensor of valence 3
	 *  (read/write version).
	 *
	 * @param i1  value of the first index (1, 2 or 3)
	 * @param i2  value of the second index (1, 2 or 3)
	 * @param i3  value of the third index (1, 2 or 3)
	 *
	 * @return modifiable reference on the component specified by \c (i1,i2,i3) 
	 *
	 */
	Scalar& set(int i1, int i2, int i3) ; 
	
	/** Returns the value of a component for a tensor of valence 4
	 *  (read/write version).
	 *
	 * @param i1  value of the first index (1, 2 or 3)
	 * @param i2  value of the second index (1, 2 or 3)
	 * @param i3  value of the third index (1, 2 or 3)
	 * @param i4  value of the fourth index (1, 2 or 3)
	 *
	 * @return modifiable reference on the component specified by 
     *   \c (i1,i2,i3,i4) 
	 *
	 */
	Scalar& set(int i1, int i2, int i3, int i4) ; 
	
	/**
	 * Sets the \c Tensor  to zero in a given domain.
	 *	@param l [input]  Index of the domain in which the \c Tensor 
	 *			  will be set (logically) to zero.
	 */
	void annule_domain(int l) ; 

	/**
	 * Sets the \c Tensor  to zero in several domains.
	 *	@param l_min [input] The \c Tensor  will be set (logically) 
	 *			     to zero
	 *			     in the domains whose indices are in the range
	 *			     \c [l_min,l_max] .
	 *	@param l_max [input] see the comments for \c l_min .
	 * 
	 * Note that \c annule(0,nz-1) , where \c nz  is the total number
	 * of domains, is equivalent to \c set_etat_zero() .
	 */
	virtual void annule(int l_min, int l_max) ; 

        /** Performs a smooth (C^n) transition in a given domain 
         *  to zero.
         *  @param l_0 [input] in the domain of index l0 the tensor is
         *  multiplied by the right polynomial (of degree \e 2n+1), to 
         *  ensure continuty of the function and its \e n first 
         *  derivative at both ends of this domain. The tensor is unchanged
         *  in the domains \e l < \e l_0 and set to zero in domains \e l > 
	 *  \e l_0.
	 *  @param deg [input] the degree \e n of smoothness of the 
	 *  transition.
         */
         void annule_extern_cn(int l_0, int deg) ;

	/**
	 * Sets the standard spectal bases of decomposition for each component.
	 * To be used only with \c valence  lower than or equal 2.
	 */
	virtual void std_spectral_base() ; 
	
	/**
	 * Sets the standard odd spectal bases of decomposition for each component.
	 * Currently only implemented for a scalar.
	 */
	virtual void std_spectral_base_odd() ; 
	
	/** Decreases by \c dec  units the value of \c dzpuis  and 
	 *  changes accordingly the values in the 
	 *  compactified external domain (CED).
	 */
	virtual void dec_dzpuis(int dec = 1) ; 

	/** Increases by \c inc  units the value of \c dzpuis  and 
	 *  changes accordingly the values in the 
	 *  compactified external domain (CED).
	 */
	virtual void inc_dzpuis(int inc = 1) ; 
	
	/** Applies exponential filters to all components 
	 * (see \c Scalar::exponential_filter_r ). Works only for Cartesian 
	 * components.
	 */
	virtual void  exponential_filter_r(int lzmin, int lzmax, int p, 
			    double alpha= -16.) ;

	/** Applies exponential filters to all components 
	 * (see \c Scalar::exponential_filter_ylm ). Works only for Cartesian 
	 * components.
	 */
	virtual void exponential_filter_ylm(int lzmin, int lzmax, int p, 
			    double alpha= -16.) ;

    // Computational methods
    // ---------------------
    
    protected: 
        /** Computes the Lie derivative of \c this  with respect to some
         *  vector field \c v  (protected method; the public interface
         *  is method \c derive_lie ).
         */
        void compute_derive_lie(const Vector& v, Tensor& resu) const ; 

    
    public:
	/** Returns the covariant derivative of \c this  with respect to some 
         * metric \f$\gamma\f$.
         * \f$T\f$ denoting the tensor represented by \c this  and
         * \f$\nabla T\f$ its covariant derivative with respect to 
         * the metric \f$\gamma\f$, 
         * the extra index (with respect to the indices of \f$T\f$)
         * of \f$\nabla T\f$ is chosen to be the \b last  one.
         * This convention agrees with that of MTW (see Eq. (10.17) of MTW).
         * For instance, if \f$T\f$ is a 1-form, whose components
         * w.r.t. the triad \f$e^i\f$ are \f$T_i\f$: \f$T=T_i \; e^i\f$,
         * then the covariant derivative of \f$T\f$ is the bilinear form
         * \f$\nabla T\f$ whose components \f$\nabla_j T_i\f$ are
         * such that 
         * \f[
         *  \nabla T = \nabla_j T_i \; e^i \otimes e^j
         * \f]
         *
         * @param gam metric \f$\gamma\f$
         * @return covariant derivative \f$\nabla T\f$ of \c this  with 
         *  respect to the connection \f$\nabla\f$ associated with the
         *  metric \f$\gamma\f$ 
	 */
	const Tensor& derive_cov(const Metric& gam) const ; 

	/** Returns the "contravariant" derivative of \c this  with respect 
	 * to some metric \f$\gamma\f$, by raising the last index of the
         * covariant derivative (cf. method \c derive_cov() ) with 
         * \f$\gamma\f$.
	 */
	const Tensor& derive_con(const Metric& gam) const ; 

	/** Computes the divergence of \c this  with respect to 
         * some metric \f$\gamma\f$. 
         * The divergence is taken with respect of the last index of \c this 
         * which thus must be contravariant.
         * For instance if the tensor \f$T\f$ represented by \c this 
         * is a twice contravariant tensor, whose 
         * components w.r.t. the
         * triad \f$e_i\f$ are \f$T^{ij}\f$: \f$T = T^{ij} \; e_i \otimes e_j\f$,
         * the divergence of \f$T\f$ w.r.t. \f$\gamma\f$ is the vector 
         * \f[
         *   {\rm div}\,  T = \nabla_k T^{ik} \; e_i
         * \f]
         * where \f$\nabla\f$ denotes the connection associated with the metric
         * \f$\gamma\f$. 
         * @param gam metric \f$\gamma\f$
         * @return divergence of \c this  with respect to \f$\gamma\f$. 
	 */
	const Tensor& divergence(const Metric& gam) const ; 


        /** Computes the Lie derivative of \c this  with respect to some
         *  vector field \c v 
         */
        Tensor derive_lie(const Vector& v) const ; 

	/** Computes a new tensor by raising an index of \c *this 
	 *
	 *  @param ind index to be raised, with the 
 	 *   following convention : 
 	 *    \li \c ind1  = 0 : first index of the tensor 
 	 *    \li \c ind1  = 1 : second index of the tensor 
	 *    \li and so on... 
	 *   (\c ind  must be of covariant type (\c COV )).
	 *  @param gam metric used to raise the index (contraction with the
	 *    twice contravariant form of the metric on the index \c ind ). 
	 * 
	 */
	Tensor up(int ind, const Metric& gam) const ; 

	/** Computes a new tensor by lowering an index of \c *this 
	 *
	 *  @param ind index to be lowered, with the 
 	 *   following convention : 
 	 *    \li \c ind1  = 0 : first index of the tensor 
 	 *    \li \c ind1  = 1 : second index of the tensor 
	 *    \li and so on... 
	 *   (\c ind  must be of covariant type (\c CON )).
	 *  @param gam metric used to lower the index (contraction with the
	 *    twice covariant form of the metric on the index \c ind ). 
	 * 
	 */
	Tensor down(int ind, const Metric& gam) const ; 

	/** Computes a new tensor by raising or lowering all the indices 
	 *  of \c *this .
	 *
	 *  @param gam metric used to lower the contravariant indices
	 *    and raising the covariant ones. 
	 * 
	 */
	Tensor up_down(const Metric& gam) const ; 

    /** Trace on two different type indices.
     *  @param ind1 first index for the contraction, with the 
 	 *   following convention : 
 	 *    \li \c ind1  = 0 : first index of the tensor 
 	 *    \li \c ind1  = 1 : second index of the tensor 
	 *    \li and so on... 
     *  @param ind2 second index for the contraction 
     */
    Tensor trace(int ind1, int ind2) const ; 

    /** Trace with respect to a given metric.
     *  @param ind1 first index for the contraction, with the 
 	 *   following convention : 
 	 *    \li \c ind1  = 0 : first index of the tensor 
 	 *    \li \c ind1  = 1 : second index of the tensor 
	 *    \li and so on... 
     *  @param ind2 second index for the contraction 
     *  @param gam metric used to raise or lower ind1 in order that it
     *      has a opposite type than ind2
     */
    Tensor trace(int ind1, int ind2, const Metric& gam) const ; 

    /** Trace on two different type indices for a valence 2 tensor.
     */
    Scalar trace() const ; 

    /** Trace with respect to a given metric for a valence 2 tensor.
     *  @param gam metric used to raise or lower one of the indices, 
     *      in order to take the trace 
     */
    Scalar trace(const Metric& gam) const ; 


    // Accessors
    // ---------
        public:
	/**
	 * Returns the position in the array \c cmp  of a 
	 * component given by its indices.  
	 *
	 * @param ind [input] 1-D array of integers (class \c Itbl )
	 *		 of size \c valence  giving the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: 
	 *			\li \c ind(0)  : value of the first index (1, 2 or 3) 
	 *			\li \c ind(1)  : value of the second index (1, 2 or 3) 
	 *			\li and so on... 
	 *
	 * @return position in the array \c cmp  of the pointer to the
	 *  	\c Scalar  containing the component specified by \c ind 
	 */
	virtual int position(const Itbl& ind) const ;

	/**
	 * Returns the indices of a component given by its position in the 
	 * array \c cmp . 
	 *
	 * @param pos [input] position in the array \c cmp 
	 *		of the pointer to the \c Scalar  representing a component
	 *
	 * @return 1-D array of integers (class \c Itbl ) of
	 *         size \c valence  giving the value of each index 
	 *	   for the component located at the position \c pos  in
	 *		the array \c cmp, with the 
	 *		following storage convention: 
	 *			\li \c Itbl(0)  : value of the first index (1, 2 or 3) 
	 *			\li \c Itbl(1)  : value of the second index (1, 2 or 3) 
	 *			\li and so on... 
	 */
	virtual Itbl indices(int pos) const ;
	
	public:
	/// Returns the mapping.
	const Map& get_mp() const {return *mp ;} ; 

	/** Returns the vectorial basis (triad) on which the components
	 *  are defined.  
	 */
	const Base_vect* get_triad() const {return triad;} ; 
    
	/// Returns the valence.
	int get_valence() const {return valence ; } ; 

	/// Returns the number of stored components.
	int get_n_comp() const {return n_comp ;} ; 
	
	/**
	 *  Gives the type (covariant or contravariant)
	 *  of the index number \c i . \c i  must be
	 *  strictly lower than \c valence  and obey the following
	 *		      convention: 
	 *			\li \c i  = 0 : first index 
	 *			\li \c i  = 1 : second index 
	 *			\li and so on... 
	 * 
	 *  @return COV for a covariant index, CON for a
	 *	    contravariant one. 
	 */
	int get_index_type(int i) const {return type_indice(i) ;};

	/**
	 * Returns the types of all the indices.
	 * 
	 *  @return 1-D array of integers (class \c Itbl ) of size \c valence  
	 *  containing the type of each index, 
	 *  \c COV  for a covariant one and \c CON  
	 *  for a contravariant one.
	 */
	Itbl get_index_type() const {return type_indice ; } ;

	/**
	 *  Sets the type of the index number \c i . \c i  must be
	 *  strictly lower than \c valence  and obey the following
	 *		      convention: 
	 *			\li \c i  = 0 : first index 
	 *			\li \c i  = 1 : second index 
	 *			\li and so on... 
	 * 
	 *  @return reference on the type that can be modified 
	 *  (\c COV for a covariant index, \c CON for a contravariant one)
	 */
	int& set_index_type(int i) {return type_indice.set(i) ;};

	/**
	 * Sets the types of all the indices.
	 * 
	 *  @return a reference on the 1-D array of integers (class \c Itbl ) 
	 *  of size \c valence that can be modified
	 *  (\c COV  for a covariant one and \c CON  for a contravariant one)
	 */
	Itbl& set_index_type() {return type_indice ; } ;

	
	/** Returns the value of a component (read-only version).
	 *
	 * @param ind  1-D \c Itbl  of size \c valence  containing the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: 
	 *			\li \c ind(0)  : value of the first index (1, 2 or 3) 
	 *			\li \c ind(1)  : value of the second index (1, 2 or 3) 
	 *			\li and so on... 
	 * @return reference on the component specified by \c ind 
	 *
	 */
	const Scalar& operator()(const Itbl& ind) const ; 

	/** Returns the value of a component for a tensor of valence 2
	 *  (read-only version).
	 *
	 * @param i1  value of the first index (1, 2 or 3)
	 * @param i2  value of the second index (1, 2 or 3)
	 *
	 * @return reference on the component specified by \c (i1,i2) 
	 *
	 */
	const Scalar& operator()(int i1, int i2) const ; 

	/** Returns the value of a component for a tensor of valence 3
	 *  (read-only version).
	 *
	 * @param i1  value of the first index (1, 2 or 3)
	 * @param i2  value of the second index (1, 2 or 3)
	 * @param i3  value of the third index (1, 2 or 3)
	 *
	 * @return reference on the component specified by \c (i1,i2,i3) 
	 *
	 */
	const Scalar& operator()(int i1, int i2, int i3) const ; 
	
	/** Returns the value of a component for a tensor of valence 4
	 *  (read-only version).
	 *
	 * @param i1  value of the first index (1, 2 or 3)
	 * @param i2  value of the second index (1, 2 or 3)
	 * @param i3  value of the third index (1, 2 or 3)
	 * @param i4  value of the fourth index (1, 2 or 3)
	 *
	 * @return reference on the component specified by \c (i1,i2,i3,i4) 
	 *
	 */
	const Scalar& operator()(int i1, int i2, int i3, int i4) const ; 
	
    // Member arithmetics
    // ------------------
    public:
	void operator+=(const Tensor &) ;		    ///< += Tensor
	void operator-=(const Tensor &) ;		    ///< -= Tensor

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    ///< Save in a binary file

	/** Displays the spectral coefficients and the associated
	 *  basis functions of each component. This function shows 
	 *  only the values greater than a given threshold.
         *   @param comment comment to be printed at top of the display
         *      (default: 0x0 = nothing printed)
	 *   @param threshold [input] Value above which a coefficient is printed
	 *    (default: 1.e-7)
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param ostr [input] Output stream used for the printing (default: cout)
	 */
	virtual void spectral_display(const char* comment = 0x0, 
                            double threshold = 1.e-7, int precision = 4, 
			    ostream& ostr = cout) const ;

	friend ostream& operator<<(ostream& , const Tensor & ) ;
	

    // Friend classes 
    // ---------------
	friend class Scalar ;
	friend class Vector ;
	friend class Sym_tensor ;
    friend class Tensor_sym ; 
	friend class Metric ;
  
    // Mathematical operators
    // ----------------------
    
    friend Scalar operator+(const Tensor&, const Scalar&) ;
    friend Scalar operator+(const Scalar&, const Tensor&) ;
    friend Scalar operator-(const Tensor&, const Scalar&) ;
    friend Scalar operator-(const Scalar&, const Tensor&) ;
    friend Tensor operator*(const Tensor&, const Tensor&) ; 
    friend Tensor_sym operator*(const Tensor&, const Tensor_sym&) ; 
    friend Tensor_sym operator*(const Tensor_sym&, const Tensor&) ; 	
    friend Tensor_sym operator*(const Tensor_sym&, const Tensor_sym&) ; 	
   
};



			//-------------------------//
			//    class Tensor_sym     //
			//-------------------------//

/**
 * Symmetric tensors (with respect to two of their arguments).
 *
 * This subclass of \c Tensor  is intended to store the components of a 
 * tensorial field with respect to a specific basis (triad), in the case
 * the tensor has a valence at least 2 and is symmetric with respect 
 * to two of its arguments (or in other words, the components are
 * symmetric with respect to two of their indices).   
 * \ingroup (tensor)
 * 
 */
class Tensor_sym : public Tensor { 

    // Data : 
    // -----
    protected:
	
        /// Number of the first symmetric index (\c 0<= \c id_sym1 < \c valence )
        int id_sym1 ;
        
        /** Number of the second symmetric index 
         * (\c id_sym1 < \c id_sym2 < \c valence )
         */
        int id_sym2 ;
        
         
    // Constructors - Destructor :
    // -------------------------
	
	public: 

	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor (must be at least 2)
	 * @param tipe  1-D array of integers (class \c Itbl ) 
	 *		of size \c valence  containing the type 
	 *		of each index, \c COV  for a covariant one 
	 *		and \c CON  for a contravariant one,  with the 
	 *		following storage convention: 
	 *			\li \c tipe(0)  : type of the first index 
	 *			\li \c tipe(1)  : type of the second index 
	 *			\li and so on... 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined 
         * @param index_sym1 number of the first symmetric index 
         *                      (\c 0<= \c index_sym1 < \c valence )
         * @param index_sym2 number of the second symmetric index 
         *                      (\c index_sym1 < \c index_sym2 < \c valence )
	 */
	Tensor_sym(const Map& map, int val, const Itbl& tipe, 
		 	const Base_vect& triad_i, int index_sym1, 
                        int index_sym2) ;

	/** Standard constructor when all the indices are of 
	 *  the same type.
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor (must be at least 2)
	 * @param tipe  the type (\c COV  or \c CON ) of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined.
         * @param index_sym1 number of the first symmetric index 
         *                      (\c 0<= \c index_sym1 < \c valence )
         * @param index_sym2 number of the second symmetric index 
         *                      (\c index_sym1 < \c index_sym2 < \c valence )
	 */
	Tensor_sym(const Map& map, int val, int tipe, const Base_vect& triad_i,
                    int index_sym1, int index_sym2) ;

	/** Constructor for a valence 3 symmetric tensor.
	 * 
	 * @param map  the mapping
	 * @param tipe0 type (\c COV  or \c CON ) of the first index.
	 * @param tipe1 type (\c COV  or \c CON ) of the second index.
	 * @param tipe2 type (\c COV  or \c CON ) of the third index.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined.
         * @param index_sym1 number of the first symmetric index 
         *                      (\c 0<= \c index_sym1 \c <=2 )
         * @param index_sym2 number of the second symmetric index 
         *                      (\c index_sym1 < \c index_sym2 \c <=2 )
	 */
	Tensor_sym(const Map& map, int tipe0, int tipe1, int tipe2, 
                   const Base_vect& triad_i,
                   int index_sym1, int index_sym2) ;

	Tensor_sym(const Tensor_sym& a) ;  ///< Copy constructor

	/** Constructor from a file (see \c sauve(FILE*) ).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been created by 
	 *			    the function \c sauve(FILE*).
	 */
	Tensor_sym(const Map& map, const Base_vect& triad_i, FILE* fich) ;

    public: 

	virtual ~Tensor_sym() ;	///< Destructor
	
    // Mutators / assignment
    // ---------------------
    public:
	
        /// Assignment to another \c Tensor_sym 
	virtual void operator=(const Tensor_sym& a) ;
	
        /** Assignment to a \c Tensor 
         * NB: the symmetry about the indices \c id_sym1  and 
         * \c id_sym2  of the input tensor is assumed but is not checked
         */
	virtual void operator=(const Tensor& a) ; 
	

    // Accessors
    // ---------
        public:
        /// Number of the first symmetric index (\c 0<= \c id_sym1 < \c valence )
        int sym_index1() const {return id_sym1;} ;
        
        /** Number of the second symmetric index 
         * (\c  id_sym1 < \c id_sym2 < \c valence )
         */
        int sym_index2() const {return id_sym2;} ;
        
	/**
	 * Returns the position in the array \c cmp  of a 
	 * component given by its indices.  
	 *
	 * @param ind [input] 1-D array of integers (class \c Itbl )
	 *		 of size \c valence  giving the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: 
	 *			\li \c ind(0)  : value of the first index (1, 2 or 3) 
	 *			\li \c ind(1)  : value of the second index (1, 2 or 3) 
	 *			\li and so on... 
	 *
	 * @return position in the array \c cmp  of the pointer to the
	 *  	\c Scalar  containing the component specified by \c ind 
	 */
	virtual int position(const Itbl& ind) const ;

	/**
	 * Returns the indices of a component given by its position in the 
	 * array \c cmp . 
	 *
	 * @param pos [input] position in the array \c cmp 
	 *		of the pointer to the \c Scalar  representing a component
	 *
	 * @return 1-D array of integers (class \c Itbl ) of
	 *         size \c valence  giving the value of each index 
	 *	   for the component located at the position \c pos  in
	 *		the array \c cmp, with the 
	 *		following storage convention: 
	 *			\li \c Itbl(0)  : value of the first index (1, 2 or 3) 
	 *			\li \c Itbl(1)  : value of the second index (1, 2 or 3) 
	 *			\li and so on... 
	 */
	virtual Itbl indices(int pos) const ;
        
            
    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;      ///< Save in a binary file
	

    // Tensor calculus
    // ---------------
    public:
    
	/** Returns the covariant derivative of \c this  with respect to some 
         * metric \f$\gamma\f$.
         * \f$T\f$ denoting the tensor represented by \c this  and
         * \f$\nabla T\f$ its covariant derivative with respect to 
         * the metric \f$\gamma\f$, 
         * the extra index (with respect to the indices of \f$T\f$)
         * of \f$\nabla T\f$ is chosen to be the \b last  one.
         * This convention agrees with that of MTW (see Eq. (10.17) of MTW).
         *
         * @param gam metric \f$\gamma\f$
         * @return covariant derivative \f$\nabla T\f$ of \c this  with 
         *  respect to the connection \f$\nabla\f$ associated with the
         *  metric \f$\gamma\f$ 
	 */
	const Tensor_sym& derive_cov(const Metric& gam) const ; 

	/** Returns the "contravariant" derivative of \c this  with respect 
	 * to some metric \f$\gamma\f$, by raising the last index of the
         * covariant derivative (cf. method \c derive_cov() ) with 
         * \f$\gamma\f$.
	 */
	const Tensor_sym& derive_con(const Metric& gam) const ; 

        /** Computes the Lie derivative of \c this  with respect to some
         *  vector field \c v 
         */
        Tensor_sym derive_lie(const Vector& v) const ; 


    // Mathematical operators
    // ----------------------
    
    friend Tensor_sym operator*(const Tensor&, const Tensor_sym&) ; 
    friend Tensor_sym operator*(const Tensor_sym&, const Tensor&) ; 
 	
}; 



/**
 * \defgroup tenso_cal Tensor calculus
 * \ingroup (tensor)
 * @{
 */

/// Tensorial product
Tensor operator*(const Tensor& a, const Tensor& b) ; 

/// Tensorial product with symmetries
Tensor_sym operator*(const Tensor& a, const Tensor_sym& b) ; 

/// Tensorial product with symmetries
Tensor_sym operator*(const Tensor_sym& a, const Tensor& b) ; 

/** Tensorial product of two symmetric tensors.
 * NB: the output is an object of class \c Tensor_sym , with
 * the two symmetric indices corresponding to the symmetric indices
 * of tensor \c a . This means that the symmetries of tensor
 * \c b  indices are not used in the storage, since 
 * there is currently no class in Lorene to manage
 * tensors with more than two symmetric indices. 
 */
Tensor_sym operator*(const Tensor_sym& a, const Tensor_sym& b) ; 


/** Contraction of two tensors. 
 *
 * @param t1 [input] first tensor 
 * @param ind1 [input] index of the first tensor for the contraction, 
 *    obeying to the following convention : 
 *    \li \c ind1  = 0 : first index of the tensor 
 *    \li \c ind1  = 1 : second index of the tensor 
 *    \li and so on... 
 *  (\c ind1  must thus be in the range 0...t1.valence-1)  
 * @param t2 [input] second tensor 
 * @param ind2 [input] index of the second tensor for the contraction, with 
 *   the same convention as \c ind1 
 * @param desaliasing [input] determines whether the products are performed
 *      with desaliasing or not  
 * @return tensor resulting of the contraction of the index \c ind1  of
 *  \c t1  with the index \c ind2  of \c t2 .
 * NB: the types (\c COV  or \c CON ) of the indices \c ind1  and
 * \c ind2  must be different. 
 */
Tensor contract(const Tensor& t1, int ind1, const Tensor& t2, int ind2, 
        bool desaliasing = false) ;

/** Double contraction of two tensors. 
 *
 * @param t1 [input] first tensor 
 * @param ind_i1 [input] position of the first index \e i1 in the first 
 *      tensor for the contraction, 
 *    obeying to the following convention : 
 *    \li \c ind_i1  = 0 : first index of the tensor 
 *    \li \c ind_i1  = 1 : second index of the tensor 
 *    \li and so on... 
 *  (\c ind_i1  must thus be in the range 0...t1.valence-1)  
 * @param ind_j1 [input] position of the second index \e j1 in the first 
 *      tensor for the contraction; one must have \c ind_i1 < \c ind_ j1
 * @param t2 [input] second tensor 
 * @param ind_i2 [input] position of the first index \e i2 in the second 
 *      tensor for the contraction 
 * @param ind_j2 [input] position of the second index \e j2 in the second 
 *      tensor for the contraction; one must have \c ind_i2 < \c ind_ j2
 * @param desaliasing [input] determines whether the products are performed
 *      with desaliasing or not  
 * @return tensor resulting of the contraction of the index \c ind_i1  
 *  of \c t1  with the index \c ind_i2  of \c t2  and of the 
 * contraction of the index \c ind_j1  
 *  of \c t1  with the index \c ind_j2  of \c t2 
 * NB: the types (\c COV  or \c CON ) of the indices \c ind_i1  and
 * \c ind_i2  (resp. \c ind_j1  and \c ind_j2 ) must be different. 
 */
Tensor contract(const Tensor& t1, int ind_i1, int ind_j1, 
                const Tensor& t2, int ind_i2, int ind_j2,
                bool desaliasing = false) ;


/** Contraction on two indices of a single tensor (trace). 
 *
 * @param t1 [input] tensor 
 * @param ind1 [input] first index of the tensor for the contraction, 
 *    obeying to the following convention : 
 *    \li \c ind1  = 0 : first index of the tensor 
 *    \li \c ind1  = 1 : second index of the tensor 
 *    \li and so on... 
 *  (\c ind1  must thus be in the range 0...t1.valence-1)  
 * @param ind2 [input] second index of the tensor for the contraction, with 
 *   the same convention as \c ind1  
 * @return tensor resulting of the contraction 
 * NB: the types (\c COV  or \c CON ) of the indices \c ind1  and
 * \c ind2  must be different. 
 */
Tensor contract(const Tensor& t1, int ind1, int ind2) ;


/** Maxima in each domain of the values of the tensor components
 * @param aa tensor
 * @param comment comment to be printed on \c ost  before the result
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result
 * @return 2-D \c Tbl  of size the number of independent components
 *	times the number of domains, the elements \c (i,l) 
 *     of which are \c max(a(l)) , where \c a(l)  
 *      denotes symbolically the values of \c aa  
 *	   in domain no. \c l  and for component no.\c i . 
 */
Tbl max(const Tensor& aa, const char* comment = 0x0, ostream& ost = cout) ; 


/** Minima in each domain of the values of the tensor components
 * @param aa tensor
 * @param comment comment to be printed on \c ost  before the result
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result
 * @return 2-D \c Tbl  of size the number of independent components
 *	times the number of domains, the elements \c (i,l) 
 *     of which are \c min(a(l)), where \c a(l)  
 *      denotes symbolically the values of \c aa  
 *	   in domain no. \c l  and for component no.\c i . 
 */
Tbl min(const Tensor& aa, const char* comment = 0x0, ostream& ost = cout) ; 

/** Maxima in each domain of the absolute values of the tensor components
 * @param aa tensor
 * @param comment comment to be printed on \c ost  before the result
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result
 * @return 2-D \c Tbl  of size the number of independent components
 *	times the number of domains, the elements \c (i,l) 
 *     of which are \c max[abs(a(l))] , where \c a(l)  
 *      denotes symbolically the values of \c aa  
 *	   in domain no. \c l  and for component no.\c i . 
 */
Tbl maxabs(const Tensor& aa, const char* comment = 0x0, ostream& ost = cout, 
	   bool verb = true) ; 


/** Relative difference between two \c Tensor  (\f$L^1\f$ version).
 * @param aa first tensor
 * @param bb second tensor
 * @param comment comment to be printed on \c ost  before the result
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result
 * @return 2-D \c Tbl  of size the number of independent components
 *	times the number of domains, the elements \c (i,l) 
 *     of which 
 *	   are \c norme[a(l)-b(l)]/norme[b(l)]  if \c b(l)!=0  and
 *	   \c norme[a(l)-b(l)]  if  \c b(l)=0 ,  where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c aa  and \c bb  
 *	   in domain no. \c l  and for component no.\c i . 
 */
Tbl diffrel(const Tensor& aa, const Tensor& bb, const char* comment = 0x0,
            ostream& ost = cout) ; 

/** Relative difference between two \c Tensor  (max version).
 * @param aa first tensor
 * @param bb second tensor
 * @param comment comment to be printed on \c ost  before the result
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result
 * @return 2-D \c Tbl  of size the number of independent components
 *	times the number of domains, the elements \c (i,l) 
 *     of which 
 *	   are \c max[abs(a(l)-b(l))]/max[abs(b(l))]  if \c b(l)!=0  and
 *	   \c max[abs(a(l)-b(l))]  if  \c b(l)=0 , where \c a(l)  and 
 *	   \c b(l)  denote symbolically the values of \c aa  and \c bb  
 *	   in domain no. \c l  and for component no.\c i . 
 */
Tbl diffrelmax(const Tensor& aa, const Tensor& bb, const char* comment = 0x0,
               ostream& ost = cout) ; 

/** Central value of each component of a tensor. 
 * @param aa tensor
 * @param comment comment to be printed on \c ost
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result; used
 *  only if \c comment != \c 0x0. 
 * @return 1-D \c Tbl  of size the number of independent components
 *  (\c aa.get_ncomp()), the elements of which are the central values 
 *  of the various components. 
 */
Tbl central_value(const Tensor& aa, const char* comment = 0x0, ostream& ost = cout) ; 

/** Maximum value of each component of a tensor over all the domains. 
 * @param aa tensor
 * @param l_excluded index of domain to be excluded from the computation:
 *   the default = -1 corresponds to no excluded domain
 * @param comment comment to be printed on \c ost
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result; used
 *  only if \c comment != \c 0x0. 
 * @return 1-D \c Tbl  of size the number of independent components
 *  (\c aa.get_ncomp()), the elements of which are the maximum values 
 *  of the various components. 
 */
Tbl max_all_domains(const Tensor& aa, int l_excluded = -1, const char* comment = 0x0, 
    ostream& ost = cout) ; 


/** Minimum value of each component of a tensor over all the domains. 
 * @param aa tensor
 * @param l_excluded index of domain to be excluded from the computation:
 *   the default = -1 corresponds to no excluded domain
 * @param comment comment to be printed on \c ost
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result; used
 *  only if \c comment != \c 0x0. 
 * @return 1-D \c Tbl  of size the number of independent components
 *  (\c aa.get_ncomp()), the elements of which are the minimum values 
 *  of the various components. 
 */
Tbl min_all_domains(const Tensor& aa, int l_excluded = -1, const char* comment = 0x0, 
    ostream& ost = cout) ; 

/** Maximum of the absolute value of each component of a tensor over all the domains. 
 * @param aa tensor
 * @param l_excluded index of domain to be excluded from the computation:
 *   the default = -1 corresponds to no excluded domain
 * @param comment comment to be printed on \c ost
 *    (default: 0x0 = nothing printed)
 * @param ost output stream for a formatted output of the result; used
 *  only if \c comment != \c 0x0. 
 * @return 1-D \c Tbl  of size the number of independent components
 *  (\c aa.get_ncomp()), the elements of which are the maximum of of the 
 *  absolute value of the various components. 
 */
Tbl maxabs_all_domains(const Tensor& aa, int l_excluded = -1, const char* comment = 0x0, 
    ostream& ost = cout, bool verb = true) ; 



/** @} */


/**
 * \defgroup tens_ari Tensor arithmetics
 * \ingroup (tensor)
 * @{
 */

Tensor operator+(const Tensor& ) ;			///< + Tensor
Tensor operator-(const Tensor& ) ;			///< \c - Tensor
Tensor operator+(const Tensor& a, const Tensor& b) ;	///< Tensor + Tensor

/// Tensor + Scalar. The \c Tensor  must be of valence 0.
Scalar operator+(const Tensor& a, const Scalar& b) ;	

/// Scalar + Tensor. The \c Tensor  must be of valence 0.
Scalar operator+(const Scalar& a, const Tensor& b) ;	

Tensor operator-(const Tensor& a, const Tensor& b) ;    ///< Tensor - Tensor

/// Tensor - Scalar. The \c Tensor  must be of valence 0.
Scalar operator-(const Tensor& a, const Scalar& b) ;	

/// Scalar - Tensor. The \c Tensor must be of valence 0.
Scalar operator-(const Scalar& a, const Tensor& b) ;	

Tensor operator*(const Scalar& a , const Tensor& b) ;   ///< Scalar * Tensor
Tensor operator*(const Tensor& a, const Scalar& b) ;    ///< Tensor * Scalar 
Tensor operator*(double , const Tensor&) ;              ///< double * Tensor
Tensor operator* (const Tensor&, double) ;              ///< Tensor * double
Tensor operator*(int, const Tensor &) ;                 ///< int* Tensor
Tensor operator*(const Tensor&, int) ;                 ///< Tensor * int
Tensor operator/(const Tensor&, const Scalar&) ;       ///< Tensor / Scalar
Tensor operator/(const Tensor&, double) ;              ///< Tensor / double
Tensor operator/(const Tensor&, int) ;                 ///< Tensor / int

    //@}

/**
 * @name Tensor_sym arithmetics
 */
    //@{
/** + Tensor_sym.  For efficiency reasons this function is 
 *  distinct from \c Tensor \c operator+(const Tensor\& ) .
 */
Tensor_sym operator+(const Tensor_sym&) ;  

/** - Tensor_sym.  For efficiency reasons this function is 
 *  distinct from \c Tensor \c operator+(const Tensor\& ).
 */
Tensor_sym operator-(const Tensor_sym&) ;  

/** Tensor_sym + Tensor_sym.  For efficiency reasons this function is 
 *  distinct
 *  from \c Tensor \c operator+(const Tensor\&, const Tensor\& ).
 */
Tensor_sym operator+(const Tensor_sym&, const Tensor_sym&) ;  

/** Tensor_sym - Tensor_sym. For efficiency reasons this function is 
 *  distinct
 *  from \c Tensor \c operator-(const Tensor\&, const Tensor\&) .
 */
Tensor_sym operator-(const Tensor_sym&, const Tensor_sym&) ;  

/** Scalar * Tensor_sym. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator*(const Scalar\&, const Tensor\&) .
 */
Tensor_sym operator*(const Scalar& a, const Tensor_sym& b) ;   

/** Tensor_sym * Scalar. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator*(const Tensor\&, const Scalar\&) .
 */
Tensor_sym operator*(const Tensor_sym& a, const Scalar& b) ;  

/** double * Tensor_sym. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator*(double, const Tensor\&) .
 */
Tensor_sym operator*(double, const Tensor_sym&) ;  

/**  Tensor_sym * double. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator*(const Tensor\&, double) .
 */
Tensor_sym operator*(const Tensor_sym&, double) ;  

/** int * Tensor_sym. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator*(int, const Tensor\&) .
 */
Tensor_sym operator*(int, const Tensor_sym&) ;  

/**  Tensor_sym * int. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator*(const Tensor\&, int) .
 */
Tensor_sym operator*(const Tensor_sym&, int) ;  

/** Tensor_sym / Scalar. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator/(const Tensor\&, const Scalar\&) .
 */
Tensor_sym operator/(const Tensor_sym&, const Scalar&) ;  

/**  Tensor_sym / double. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator/(const Tensor\&, double) .
 */
Tensor_sym operator/(const Tensor_sym&, double) ;  

/**  Tensor_sym / int. For efficiency reasons this function is distinct
 *  from \c Tensor \c operator/(const Tensor\&, int) .
 */
Tensor_sym operator/(const Tensor_sym&, int) ;  

/** @} */


#include "scalar.h"

#include "vector.h"

#include "sym_tensor.h"


#endif
