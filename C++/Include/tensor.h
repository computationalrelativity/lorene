/*
 *  Definition of Lorene classes Tensor and Sym_tensor
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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
class Sym_tensor ;
class Metric ;

			//-------------------------//
			//       class Tensor      //
			//-------------------------//
			

/**
 * Tensor handling *** UNDER DEVELOPMENT  ***.
 *
 * This class is intended to replace {\tt Tenseur} and {\tt Cmp} (the
 *  latter via the derived class {\tt Scalar}).
 * 
 * The {\tt Tensor} class is intended to store the components of a tensorial 
 * field in a specific basis (triad).  
 * 
 * All this is {\it 3D} meaning that the indices go from 1 to 3. 
 * 
 * 
 * @version #$Id$#
 */
class Tensor { 

    // Data : 
    // -----
    protected:
	
	const Map* const mp ;	/// Reference mapping

	int valence ;		/// Valence
	
	/** Vectorial basis (triad) with respect to which the tensor
	 *  components are defined. 
	 */
	const Base_vect* triad ; 

	/** 1D array of integers (class {\tt Itbl}) of size {\tt valence} 
	 *  containing the type of each index: 
	 *  {\tt COV} for a covariant one and {\tt CON} for a contravariant one.
	 * 
	 */	
	Itbl type_indice ;	
	
	int n_comp ;	/// Number of stored components, depending on the symmetry.

	/// Array of size {\tt n\_comp} of pointers onto the components.
	Scalar** cmp ;   


    // Derived data : 
    // ------------
     protected:
	/**
	 * Array on the {\tt Metric} which were used to compute derived
	 * quantities, like {\tt p\_derive\_cov}, etc... 
	 * The i-th element of this array is the {\tt Metric} used to
	 * compute the i-th element of {\tt p\_derive\_cov}, etc..
	 */
	mutable const Metric* met_depend[N_MET_MAX] ; 

	/**Array of pointers on the covariant derivatives of {\tt this}
	 * (see the comments of {\tt met\_depend}).
	 */
	mutable Tensor* p_derive_cov[N_MET_MAX];
	
	/**Array of pointers on the contravariant derivatives of {\tt this}
	 * (see the comments of {\tt met\_depend}).
	 */
	mutable Tensor* p_derive_con[N_MET_MAX];

	/**Array of pointers on the divergence of {\tt this}
	 * (see the comments of {\tt met\_depend}). It is assumed that the first 
	 * index of {\tt this} is contravariant.
	 */
	mutable Tensor* p_divergence[N_MET_MAX];

   
    // Constructors - Destructor :
    // -------------------------
	
	public: 

	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor
	 * @param tipe  1-D array of integers (class {\tt Itbl}) 
	 *		of size {\tt valence} containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index \\
	 *			and so on... 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined 
	 */
	Tensor(const Map& map, int val, const Itbl& tipe, 
		 	const Base_vect& triad_i) ;

	/** Standard constructor with the triad passed as a pointer.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor
	 * @param tipe  1-D array of integers (class {\tt Itbl}) 
	 *		of size {\tt valence} containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index \\
	 *			and so on... 
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
	 * @param tipe  the type ({\tt COV} or {\tt CON}) of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined.
	 */
	Tensor(const Map& map, int val, int tipe, 
			const Base_vect& triad_i) ;

	Tensor(const Tensor&) ;  /// Copy constructor

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
	Tensor(const Map& map, const Base_vect& triad_i, FILE* fich) ;

    protected:
	/**
	 *  Constructor for a scalar field: to be used only by the derived
	 *  class {\tt Scalar}.
	 *
	 */
	 Tensor(const Map& map) ;

	/**
	 * Constructor to be used by derived classes, with symmetries among
	 *  the components. The number of independent components is
	 *  given as an argument ({\tt n\_comp\_i}), and not computed
	 *	from the valence, as in the standard constructor.  
	 *
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  1-D array of integers (class {\tt Itbl}) 
	 *		of size {\tt valence} containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index \\
	 *			and so on... 
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
	 *  given as a argument ({\tt n\_comp\_i}), and not computed
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

	virtual ~Tensor() ;	/// Destructor
	
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const ;	/// Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 

	/**
	 * Logical destructor of the derivatives depending on the i-th
	 * element of {\tt met\_depend}.
	 */	
	virtual void del_derive_met(int) const ;

	/**
	 * Sets all the i-th components of {\tt met\_depend}, 
	 * {\tt p\_derive\_cov}, etc... to 0x0.
	 */
	void set_der_met_0x0(int) const ;

	/**
	 * To be used to describe the fact that the derivatives members have
	 * been calculated with {\tt met}.
	 * 
	 * First it sets a null element of {\tt met\_depend} to 
	 * {\tt \&met} and puts {\tt this} in 
	 * the list of the dependancies of {\tt met}.
	 * 
	 */
	void set_dependance (const Metric&) const ;

	/**
	 * Returns the position of the pointer on {\tt metre} in 
	 * the array {\tt met\_depend}.
	 *
	 */
	int get_place_met(const Metric&) const ;
	
    // Mutators / assignment
    // ---------------------
    public:
	/**
	 * Sets the logical state of all components to {\tt ETATNONDEF} 
	 * (undefined state).
	 */
	virtual void set_etat_nondef() ;
	
	/**
	 * Sets the logical state of all components to {\tt ETATZERO} 
	 *(zero state).
	 */
	virtual void set_etat_zero() ;

	/**
	 * Sets the logical state of all components to {\tt ETATQCQ} 
	 * (ordinary state).
	 */
	virtual void set_etat_qcq() ;
	
	/**
	 *  Performs the memory allocation of all the 
	 *  elements, down to the {\tt double} arrays of the {\tt Tbl}s. 
	 *  This function performs in fact recursive calls to 
	 *  {\tt set\_etat\_qcq()}
	 *  on each element of the chain {\tt Scalar} ->
	 *  {\tt Valeur} -> {\tt Mtbl} -> {\tt Tbl}. 
	 */
	virtual void allocate_all() ; 

	/** Sets a new vectorial basis (triad) of decomposition and modifies
	 *  the components accordingly. 
	 */
	virtual void change_triad(const Base_vect& new_triad) ; 
    
	/** Assigns a new vectorial basis (triad) of decomposition. 
	 *  NB: this function modifies only the member {\tt triad} and
	 *  leave unchanged the components (member {\tt cmp}). In order to 
	 *  change them coherently with the new basis, the function 
	 *  {\tt change\_triad(const Base\_vect\& )} must be called instead. 
	 */
	void set_triad(const Base_vect& new_triad) ; 

	
	virtual void operator=(const Tensor&) ;/// Assignment to a {\tt Tensor}
	
	/** Returns the value of a component (read/write version).
	 *
	 * @param ind  1-D {\tt Itbl} of size {\tt valence} containing the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: \\
	 *			{\tt ind(0)} : value of the first index (1, 2 or 3) \\
	 *			{\tt ind(1)} : value of the second index (1, 2 or 3) \\
	 *			and so on... 
	 * @return modifiable reference on the component specified by {\tt ind}
	 *
	 */
	Scalar& set(const Itbl& ind) ; 
	
	/** Returns the value of a component for a tensor of valence 2
	 *  (read/write version).
	 *
	 * @param i1  value of the first index (1, 2 or 3)
	 * @param i2  value of the second index (1, 2 or 3)
	 *
	 * @return modifiable reference on the component specified by {\tt (i1,i2)}
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
	 * @return modifiable reference on the component specified by {\tt (i1,i2,i3)}
	 *
	 */
	Scalar& set(int i1, int i2, int i3) ; 
	
	/**
	 * Sets the {\tt Tensor} to zero in a given domain.
	 *	@param l [input]  Index of the domain in which the {\tt Tensor}
	 *			  will be set (logically) to zero.
	 */
	void annule_domain(int l) ; 

	/**
	 * Sets the {\tt Tensor} to zero in several domains.
	 *	@param l_min [input] The {\tt Tensor} will be set (logically) 
	 *			     to zero
	 *			     in the domains whose indices are in the range
	 *			     {\tt [l\_min, l\_max]}.
	 *	@param l_max [input] see the comments for {\tt l\_min}.
	 * 
	 * Note that {\tt annule(0, nz-1)}, where {\tt nz} is the total number
	 * of domains, is equivalent to {\tt set\_etat\_zero()}.
	 */
	virtual void annule(int l_min, int l_max) ; 

	/**
	 * Sets the standard spectal bases of decomposition for each component.
	 * To be used only with {\tt valence} lower than or equal 2.
	 */
	virtual void std_spectral_base() ; 
	
	/** Decreases by {\tt dec} units the value of {\tt dzpuis} and 
	 *  changes accordingly the values in the 
	 *  compactified external domain (CED).
	 */
	virtual void dec_dzpuis(int dec = 1) ; 

	/** Increases by {\tt inc} units the value of {\tt dzpuis} and 
	 *  changes accordingly the values in the 
	 *  compactified external domain (CED).
	 */
	virtual void inc_dzpuis(int inc = 1) ; 
	
	/// Multiplication by {\it r} in the external domain.
	virtual void mult_r_ced() ; 

	///The covariant derivative of {\tt this} with respect to a {\tt Metric}
	const Tensor& derive_cov(const Metric&) const ; 

	/**The contravariant derivative of {\tt this} with respect 
	 * to a {\tt Metric}. It uses the covariant derivative and then
	 * raises the first index.
	 */
	const Tensor& derive_con(const Metric&) const ; 

	/**The divergence of {\tt this} with respect to a {\tt Metric}.
	 * The first index of {\tt this} is assumed to be contravariant.
	 */
	const Tensor& divergence(const Metric&) const ; 

	/** Index contraction on two indices of different type (contravariant
 	 *  or covariant) of the tensor. 
 	 * 
	 * @param ind1 [input] first index for the contraction, obeying to the
 	 *   following convention : \\
 	 *    {\tt ind1} = 0 : first index of the tensor \\
 	 *    {\tt ind1} = 1 : second index of the tensor \\
	 *    and so on... \\
 	 *  ({\tt ind1} must thus be in the range 0...valence-1)  
 	 * @param ind2 [input] second index for the contraction, with the same convention 
 	 *   as {\tt ind1} 
 	 * @return tensor resulting of the contraction of the index {\tt ind1} with
 	 *   the index {\tt ind2}.
 	 * NB: the types ({\tt COV} or {\tt CON}) of the indices {\tt ind1} and
 	 * {\tt ind2} must be different. 
 	 */
	Tensor scontract(int ind1, int ind2) const ; 
	
	/** Computes a new tensor by raising an index of {\tt *this}
	 *
	 *  @param ind index to be raised, with the 
 	 *   following convention : \\
 	 *    {\tt ind1} = 0 : first index of the tensor \\
 	 *    {\tt ind1} = 1 : second index of the tensor \\
	 *    and so on... \\
	 *   ({\tt ind} must be of covariant type ({\tt COV})).
	 *  @param met metric used to raise the index (contraction with the
	 *    twice contravariant form of the metric on the index {\tt ind}). 
	 * 
	 */
	Tensor up(int ind, const Metric& met) const ; 

	/** Computes a new tensor by lowering an index of {\tt *this}
	 *
	 *  @param ind index to be lowered, with the 
 	 *   following convention : \\
 	 *    {\tt ind1} = 0 : first index of the tensor \\
 	 *    {\tt ind1} = 1 : second index of the tensor \\
	 *    and so on... \\
	 *   ({\tt ind} must be of covariant type ({\tt CON})).
	 *  @param met metric used to lower the index (contraction with the
	 *    twice covariant form of the metric on the index {\tt ind}). 
	 * 
	 */
	Tensor down(int ind, const Metric& met) const ; 

	/** Computes a new tensor by raising or lowering all the indices 
	 *  of {\tt *this}.
	 *
	 *  @param met metric used to lower the contravariant indices
	 *    and raising the covariant ones. 
	 * 
	 */
	Tensor up_down(const Metric& met) const ; 

    // Accessors
    // ---------
        public:
	/**
	 * Returns the position in the array {\tt cmp} of a 
	 * component given by its indices.  
	 *
	 * @param ind [input] 1-D array of integers (class {\tt Itbl})
	 *		 of size {\tt valence} giving the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: \\
	 *			{\tt ind(0)} : value of the first index (1, 2 or 3) \\
	 *			{\tt ind(1)} : value of the second index (1, 2 or 3) \\
	 *			and so on... 
	 *
	 * @return position in the array {\tt cmp} of the pointer to the
	 *  	{\tt Scalar} containing the component specified by {\tt ind}
	 */
	virtual int position(const Itbl& ind) const ;

	/**
	 * Returns the indices of a component given by its position in the 
	 * array {\tt cmp}. 
	 *
	 * @param pos [input] position in the array {\tt cmp}
	 *		of the pointer to the {\tt Scalar} representing a component
	 *
	 * @return 1-D array of integers (class {\tt Itbl}) of
	 *         size {\tt valence} giving the value of each index 
	 *	   for the component located at the position {\tt pos} in
	 *		the array [\tt cmp}, with the 
	 *		following storage convention: \\
	 *			{\tt Itbl(0)} : value of the first index (1, 2 or 3) \\
	 *			{\tt Itbl(1)} : value of the second index (1, 2 or 3) \\
	 *			and so on... 
	 */
	virtual Itbl indices(int pos) const ;
	
	public:
	const Map& get_mp() const {return *mp ;} ; /// Returns the mapping.

	/** Returns the vectorial basis (triad) on which the components
	 *  are defined.  
	 */
	const Base_vect* get_triad() const {return triad;} ; 
    
	int get_valence() const {return valence ; } ; /// Returns the valence.

	/// Returns the number of stored components.
	int get_n_comp() const {return n_comp ;} ; 
	
	/**
	 *  Gives the type (covariant or contravariant)
	 *  of the index number {\tt i}. {\tt i} must be
	 *  strictly lower than {\tt valence} and obey the following
	 *		      convention: \\
	 *			{\tt i} = 0 : first index \\
	 *			{\tt i} = 1 : second index \\
	 *			and so on... 
	 * 
	 *  @return COV for a covariant index, CON for a
	 *	    contravariant one. 
	 */
	int get_index_type(int i) const {return type_indice(i) ;};

	/**
	 * Returns the types of all the indices.
	 * 
	 *  @return 1-D array of integers (class {\tt Itbl}) of size {\tt valence} 
	 *  containing the type of each index, 
	 *  {\tt COV} for a covariant one and {\tt CON} 
	 *  for a contravariant one.
	 */
	Itbl get_index_type() const {return type_indice ; } ;

	
	/** Returns the value of a component (read-only version).
	 *
	 * @param ind  1-D {\tt Itbl} of size {\tt valence} containing the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: \\
	 *			{\tt ind(0)} : value of the first index (1, 2 or 3) \\
	 *			{\tt ind(1)} : value of the second index (1, 2 or 3) \\
	 *			and so on... 
	 * @return reference on the component specified by {\tt ind}
	 *
	 */
	const Scalar& operator()(const Itbl& ind) const ; 

	/** Returns the value of a component for a tensor of valence 2
	 *  (read-only version).
	 *
	 * @param i1  value of the first index (1, 2 or 3)
	 * @param i2  value of the second index (1, 2 or 3)
	 *
	 * @return reference on the component specified by {\tt (i1,i2)}
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
	 * @return reference on the component specified by {\tt (i1,i2,i3)}
	 *
	 */
	const Scalar& operator()(int i1, int i2, int i3) const ; 
	
    // Member arithmetics
    // ------------------
    public:
	void operator+=(const Tensor &) ;		    /// += Tensor
	void operator-=(const Tensor &) ;		    /// -= Tensor

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file

	/** Displays the spectral coefficients and the associated
	 *  basis functions of each component. This function shows 
	 *  only the values greater than a given threshold.
	 *   @param threshold [input] Value above which a coefficient is printed
	 *    (default: 1.e-7)
	 *   @param precision [input] Number of printed digits (default: 4)
	 *   @param ostr [input] Output stream used for the printing (default: cout)
	 */
	virtual void spectral_display(double threshold = 1.e-7, int precision = 4, 
			   ostream& ostr = cout) const ;

	friend ostream& operator<<(ostream& , const Tensor & ) ;
	

    // Friend classes 
    // ---------------
	friend class Scalar ;
	friend class Vector ;
	friend class Sym_tensor ;
	friend class Tensor_delta ;
	friend class Metric ;
    
};


/**
 * @name Tensor calculus
 */
//@{
/** Contraction of two tensors. 
 *
 * @param t1 [input] first tensor 
 * @param ind1 [input] index of the first tensor for the contraction, 
 *    obeying to the following convention : \\
 *    {\tt ind1} = 0 : first index of the tensor \\
 *    {\tt ind1} = 1 : second index of the tensor \\
 *    and so on... \\
 *  ({\tt ind1} must thus be in the range 0...t1.valence-1)  
 * @param t2 [input] second tensor 
 * @param ind2 [input] index of the second tensor for the contraction, with 
 *   the same convention as {\tt ind1} 
 * @return tensor resulting of the contraction of the index {\tt ind1} of
 *  [\tt t1} with the index {\tt ind2} of {\tt t2}.
 * NB: the types ({\tt COV} or {\tt CON}) of the indices {\tt ind1} and
 * {\tt ind2} must be different. 
 */
Tensor contract(const Tensor& t1, int ind1, const Tensor& t2, int ind2) ;


//@}

/**
 * @name Tensor mathematics
 */
    //@{
Tensor operator+(const Tensor& ) ;			/// + Tensor
Tensor operator-(const Tensor& ) ;			/// - Tensor
Tensor operator+(const Tensor&, const Tensor &) ;	/// Tensor + Tensor
Tensor operator-(const Tensor&, const Tensor &) ;       /// Tensor - Tensor
Tensor operator*(double , const Tensor&) ;              /// double * Tensor
Tensor operator* (const Tensor&, double) ;              /// Tensor * double
Tensor operator*(int, const Tensor &) ;                 /// int* Tensor
Tensor operator* (const Tensor&, int) ;                 /// Tensor * int
Tensor operator/ (const Tensor&, const Scalar&) ;       /// Tensor / Scalar
Tensor operator/ (const Tensor&, double) ;              /// Tensor / double
Tensor operator/ (const Tensor&, int) ;                 /// Tensor / int

    //@}




			//---------------------------------//
			//        class Tensor_delta       //
			//---------------------------------//
			
/**
 * Class describing valence-3 tensors, symmetric on the last two indices.
 * These object typically represent the difference between two connections
 * (see the class {\tt Connection}). The storage and the calculations are 
 * different and quicker than with a usual {\tt Tensor}.
 * 
 * The valence must be 3.
 *
 * @version #$Id$#
 */
class Tensor_delta : public Tensor {

    // Constructors - Destructor :
    // -------------------------
	
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param tipe  1-D array of integers (class {\tt Itbl}) of size 3 
	 *		containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index, etc..
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */
	Tensor_delta(const Map& map, const Itbl& tipe, const Base_vect& triad_i) ;

	/** Standard constructor with indices types specified individually
	 * 
	 * @param map   the mapping 
	 * @param tipe0 : type of the first index \\
	 * @param tipe1 : type of the second index \\
	 * @param tipe2 : type of the third index \\
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */
	Tensor_delta(const Map& map, int tipe0, int tipe1, int tipe2,
				 const Base_vect& triad_i) ;

	/** Standard constructor when both indices are of the same type.
	 * 
	 * @param map   the mapping 
	 * @param tipe  the type of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 * 
	 */
	Tensor_delta(const Map& map, int tipe, const Base_vect& triad_i) ;

	Tensor_delta(const Tensor_delta&) ; /// Copy constructor

	/** Constructor from a {\tt Tensor}.
	 *  The symmetry of the input tensor is assumed to be true but not checked.
	 */
	Tensor_delta(const Tensor&) ;
	
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
	Tensor_delta(const Map& map, const Base_vect& triad_i, FILE* fich) ;

	virtual ~Tensor_delta() ;    /// Destructor
	
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	/// Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to a {\tt Tensor\_delta}
	void operator=(const Tensor_delta&) ;

	/**
	 * Assignment from a {\tt Tensor}.
	 * 
	 * The symmetry is assumed but not checked.
	 */
	virtual void operator= (const Tensor&) ;
    

    // Accessors
    // ---------
    public:
	/**
	 * Returns the position in the array {\tt cmp} of a 
	 * component given by its indices.  
	 *
	 * @param ind [input] 1-D array of integers (class {\tt Itbl})
	 *		 of size 3 giving the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: \\
	 *		   {\tt ind(0)} : value of the first index (1, 2 or 3) \\
	 *		   {\tt ind(1)} : value of the second index (1, 2 or 3)\\
	 *                 {\it etc ...}
	 *
	 * @return position in the array {\tt cmp} of the pointer to the
	 *  	{\tt Scalar} containing the component specified by {\tt ind}
	 */
	virtual int position(const Itbl& ind) const ;


	/**
	 * Returns the indices of a component given by its position in the 
	 * array {\tt cmp}. 
	 *
	 * @param pos [input] position in the array {\tt cmp}
	 *		of the pointer to the {\tt Scalar} representing a component
	 *
	 * @return 1-D array of integers (class {\tt Itbl}) of
	 *         size 3 giving the value of each index 
	 *	   for the component located at the position {\tt pos} in
	 *		the array [\tt cmp}, with the 
	 *		following storage convention: \\
	 *		     {\tt Itbl(0)} : value of the first index (1, 2 or 3) \\
	 *		     {\tt Itbl(1)} : value of the second index (1, 2 or 3)\\
	 *                   {\it etc...}
	 */
	virtual Itbl indices(int pos) const ;
	
} ;

#include "scalar.h"

#include "vector.h"

#include "sym_tensor.h"


#endif
