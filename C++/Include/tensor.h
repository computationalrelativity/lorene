/*
 *  Definition of Lorene class Tensor
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

			//-------------------------//
			//	class Tensor		   //
			//-------------------------//
			

/**
 * Tensor handling *** UNDER DEVELOPMENT  ***.
 * This class is intended to replace {\tt Tenseur} and {\tt Cmp} (via
 *  its derived class {\tt Scalar}).
 * 
 * This class is intended to store the components of a tensorial field in 
 * a specific basis. 
 * 
 * All this is {\it 3D} meaning that the indices go from 1 to 3. Moreover,
 * the components are described in orthonormal bases.
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

	/** Array of size {\tt valence} contening the type of each index, 
	 * {\tt COV} for a covariant one and {\tt CON} for a contravariant one.
	 * 
	 */	
	Itbl type_indice ;	
	
	int n_comp ;	/// Number of stored components, depending on the symmetry.

	Scalar** cmp ;   /// The components.


    // Derived data : 
    // ------------
    // protected:
   
    // Constructors - Destructor :
    // -------------------------
	
	public: 

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
	 */
	Tensor (const Map& map, int val, const Itbl& tipe, 
		 	const Base_vect& triad_i) ;

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
	 */
	Tensor (const Map& map, int val, const Itbl& tipe, 
		 	const Base_vect* triad_i) ;

	/** Standard constructor when all the indices are of 
	 *  the same type.
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  the type of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined.
	 */
	Tensor (const Map& map, int val, int tipe, 
			const Base_vect& triad_i) ;

	Tensor (const Tensor&) ;  /// Copy constructor

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
	Tensor (const Map& map, const Base_vect& triad_i, FILE* fich) ;

	
    protected:
	/**
	 *  Constructor for a scalar field: to be used by the derived
	 *  class {\tt Scalar}
	 *
	 */
	 Tensor(const Map& map) ;

	/**
	 * Constructor used by derived classes, with symmetries among
	 *  the components.
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
	 * @param n_comp  number of stored components.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */	 
	Tensor (const Map& map, int val, const Itbl& tipe, int n_comp,
		    const Base_vect& triad_i) ;

	/**
	 * Constructor used by derived classes, with symmetries among
	 *  the components, when all the indices are of 
	 *  the same type.
	 * 
	 * @param map  the mapping
	 * @param val   valence of the tensor
	 * @param tipe  the type of the indices.
	 * @param n_comp  the number of stored components.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */
	Tensor (const Map&, int val, int tipe, int n_comp, 
		 const Base_vect& triad_i) ;


    public: 

	virtual ~Tensor() ;	/// Destructor
	
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() ;	/// Deletes derived quantities

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
	void change_triad(const Base_vect& new_triad) ; 
    
	/** Assigns a new vectorial basis (triad) of decomposition. 
	 *  NB: this function modifies only the member {\tt triad} and
	 *  leave unchanged the components (member {\tt cmp}). In order to 
	 *  change them coherently with the new basis, the function 
	 *  {\tt change\_triad(const Base\_vect\& )} must be called instead. 
	 */
	void set_triad(const Base_vect& new_triad) ; 

    
	/// Assignment to another {\tt Tensor}
	virtual void operator=(const Tensor&) ; 
	
	
	Scalar& set(int, int) ; /// Read/write for a tensor of valence 2.
	Scalar& set(int, int, int) ; /// Read/write for a tensor of valence 3.
	Scalar& set(const Itbl&) ; /// Read/write in the general case.
	    
	/**
	 * Sets the {\tt Tensor} to zero in a given domain.
	 *	@param l [input]  Index of the domain in which the {\tt Tensor}
	 *			  will be set (logically) to zero.
	 */
	void annule(int l) ; 

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
	 * Set the standard spectal basis of decomposition for each component.
	 * To be used only with {\tt valence} strictly lower than 3.
	 * 
	 */
	void std_spectral_base() ; 
	
	virtual void dec_dzpuis() ;	/// dzpuis -= 1 ;
	virtual void inc_dzpuis() ;	/// dzpuis += 1 ;
	virtual void dec2_dzpuis() ;	/// dzpuis -= 2 ;
	virtual void inc2_dzpuis() ;	/// dzpuis += 2 ;
	virtual void mult_r_zec() ; /// Multiplication by {\it r} in the external zone.
	
    // Accessors
    // ---------
    protected:
	/**
	 * Returns the position in the {\tt Scalar} array {\tt cmp} of a 
	 * component given by its indices.  
	 *
	 * @return position in the {\tt Scalar} array {\tt cmp}  
	 * corresponding to the indices given in {\tt idx}. {\tt idx}
	 * must be a 1-D {\tt Itbl} of size {\tt valence}, 
	 * each element of which must be one of the spatial indices 
	 * 1, 2 or 3. 
	 */
	virtual int position(const Itbl& idx) const ;

	/**
	 * Returns the indices of a component given by its position in the 
	 * {\tt Scalar} array {\tt cmp}. 
	 *
	 * @return 1-D array of integers ({\tt Itbl}) of
	 *         size {\tt valence} giving the value of each index 
	 *	   for the component located at the position {\tt place}
	 *	   in the {\tt Scalar} array {\tt cmp}. 
	 *	   Each element of this {\tt Itbl} 
	 *	   contains the spatial index 1, 2 or 3. 
	 */
	virtual Itbl indices(int place) const ;
	
	public:
	const Map* get_mp() const {return mp ;} ; /// Returns pointer on the mapping.

	/** Returns the vectorial basis (triad) on which the components
	 *  are defined.  
	 */
	const Base_vect* get_triad() const {return triad;} ; 
    
	int get_valence() const {return valence ; } ; ///Returns the valence.

	/// Returns the number of stored components.
	int get_n_comp() const {return n_comp ;} ; 
	
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
	int get_index_type(int i) const {return type_indice(i) ;};

	/**
	 * Returns the types of all the indices.
	 * 
	 *  @return 1-D {\tt Itbl} of size {\tt valence} containing the type 
	 *  of each index, {\tt COV} for a covariant one and {\tt CON} 
	 *  for a contravariant one.
	 */
	Itbl get_index_type () const {return type_indice ; } ;

	
	const Scalar& operator()(int, int) const ; /// Read only for a tensor of valence 2.
	const Scalar& operator()(int, int, int) const ; /// Read only for a tensor of valence 3.
	const Scalar& operator()(const Itbl&) const ; /// Read only in the general case.
	
    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file

	friend ostream& operator<<(ostream& , const Tensor & ) ;
	
    // Computation of derived members
    // ------------------------------
    // protected:

    // Friend classes 
    // ---------------
	friend class Scalar ;
    
    // Mathematical operators
    // ----------------------
    

};


/**
 * @name Tensor calculus
 */
//@{


//@}

/**
 * @name Tensor mathematics
 */
    //@{
Tensor operator+(const Tensor& ) ;			/// + Tensor
Tensor operator-(const Tensor& ) ;			/// - Tensor
Tensor operator+(const Tensor&, const Tensor &) ;	/// Tensor + Tensor


    //@}






#include "scalar.h"


#endif
