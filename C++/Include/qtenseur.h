/*
 *  Definition of Lorene classes Qtenseur
 *				 Qtenseur_sym
 *
 */

/*
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


#ifndef __QTENSEUR_H_
#define __QTENSEUR_H_

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/11/06 14:43:37  e_gourgoulhon
 * Gave a name to const arguments in certain method prototypes (e.g.
 * constructors) to correct a bug of DOC++.
 *
 * Revision 1.2  2002/09/19 10:11:17  j_novak
 * Modif. commentaires
 *
 * Revision 1.1  2002/09/19 09:52:42  j_novak
 * Added objects Qtenseur and Qmetrique for 4D tensor and metric handling.
 *
 *
 * $Header$
 *
 */

// Headers Lorene 
#include "tenseur.h"
#include "param.h"

class Qmetrique ;
class Qtenseur_sym ;

			//---------------------------------//
			//	class Qtenseur		   //
			//---------------------------------//
			

/**
 * Four dimensional tensor handling.
 * 
 * This class is intended to store the components of a tensorial field in 
 * a specific basis. It is very similar to the class {\tt Tenseur}.
 * 
 * All this is {\it 4D} meaning that the indices go from 0 to 3. Moreover,
 * the components are described in orthonormal bases. Time derivatives
 * rely on a constant time-step (e.g. in the method {\tt derive\_cov}).
 * 
 * When first constructed, the memory for each component is not allocated.
 * 
 * @version #$Id$#
 */
class Qtenseur { 

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
	
    // Derived data : 
    // ------------
    protected:
	/** Pointer on the {\tt Metrique} used to calculate derivatives 
	 *  members.
	 * If no such member has been calculated the pointer is set to zero.
	 * 
	 */
	mutable const Qmetrique* met_depend ;
	
	/** Pointer on the gradient of {\tt *this}.
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	mutable Qtenseur* p_gradient ;
	
	/** Pointer on the gradient of {\tt *this} in a spherical orthonormal
	 * basis (scalar field only).
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	mutable Qtenseur* p_gradient_spher ;
	
	/** Pointer on the covariant derivative of {\tt *this} with respect to 
	 * {\tt *met\_depend}. 
	 * It is set to zero if it has not been calculated yet, or for a scalar 
	 * field.
	 * 
	 */
	mutable Qtenseur* p_derive_cov ;

	/** Pointer on the contravariant derivative of {\tt *this} with 
	 * respect to {\tt *met\_depend}. 
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	mutable Qtenseur* p_derive_con ;

	/** Pointer on the scalar square of {\tt *this} with respect to 
	 * {\tt *met\_depend}. 
	 * It is set to zero if it has not been calculated yet.
	 * 
	 */
	mutable Qtenseur* p_carre_scal ;
   
    // Constructors - Destructor :
    // -------------------------

    public:
	/// Constructor for a scalar field. 
	explicit Qtenseur (const Map& map) ; 

	/// Constructor for a scalar field and from a {\tt Cmp}. 
	explicit Qtenseur (const Cmp& cmp) ; 

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
	Qtenseur (const Map& map, int val, const Itbl& tipe, 
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
	Qtenseur (const Map& map, int val, const Itbl& tipe, 
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
	Qtenseur (const Map& map, int val, int tipe, const 
		 Base_vect& triad_i) ;

	Qtenseur (const Qtenseur&) ;  /// Copy constructor

	/// Constructor from a symmetric tensor.
	explicit Qtenseur (const Qtenseur_sym&) ;
	
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
	Qtenseur (const Map& map, const Base_vect& triad_i, FILE* fich) ;

	/** Constructor from a file for a scalar field
	 *  (see {\tt sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param fich  file which has been created by 
	 *			    the function {\tt sauve(FILE* )}.
	 */
	Qtenseur (const Map& map, FILE* fich) ;

	
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
	 */	 
	Qtenseur (const Map& map, int val, const Itbl& tipe, int n_comp,
		 const Base_vect& triad_i) ;

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
	 */
	Qtenseur (const Map&, int val, int tipe, int n_comp, 
		 const Base_vect& triad_i) ;


    public: 

	virtual ~Qtenseur() ;	/// Destructor
	
    // Memory management
    // -----------------
    protected:
	void del_t() ;	/// Logical destructor

	/**
	 * Logical destructor of the derivatives depending on 
	 * {\tt *met\_depend}.
	 */	
	void del_derive_met() const ;

	/**
	 * Logical destructor of all the derivatives.
	 */
	void del_derive() const ;
	
	/**
	 * Sets the pointers of the derivatives depending on 
	 * {\tt *met\_depend} to zero.
	 */
	void set_der_met_0x0() const ;

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
	 *  on each element of the chain {\tt Qtenseur} -> {\tt Cmp} ->
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
    
	/// Assignment to another {\tt Qtenseur}
	virtual void operator=(const Qtenseur& a) ; 
	
	/// Assignment to a {\tt Cmp} (scalar field only)
	void operator=(const Cmp& a) ; 
	
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
	 * Sets the {\tt Qtenseur} to zero in a given domain.
	 *	@param l [input]  Index of the domain in which the {\tt Qtenseur}
	 *			  will be set (logically) to zero.
	 */
	void annule(int l) ; 

	/**
	 * Sets the {\tt Qtenseur} to zero in several domains.
	 *	@param l_min [input] The {\tt Qtenseur} will be set (logically) 
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
	 * each element of which must be 0, 1, 2, 3 
	 * corresponding to respective indices.
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
	 *	   Each element of this {\tt Itbl} is 0, 1, 2, or 3 which 
	 *	   corresponds to respective indices.
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

	const Cmp& operator()() const ; /// Read only for a scalar.
	const Cmp& operator()(int) const ; /// Read only for a vector.
	const Cmp& operator()(int, int) const ; /// Read only for a tensor of valence 2.
	const Cmp& operator()(int, int, int) const ; /// Read only for a tensor of valence 3.
	const Cmp& operator()(const Itbl&) const ; /// Read only in the general case.
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file

	friend ostream& operator<<(ostream& , const Qtenseur & ) ;
	
    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the gradient of {\tt *this}.
	 * The result is in {\tt *p\_gradient}
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	virtual void fait_gradient(const Param& para) const ;

	/**
	 * Calculates, if needed, the gradient of {\tt *this} in a 
	 * spherical orthonormal basis (scalar field only). 
	 * The result is in {\tt *p\_gradient\_spher}
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	void fait_gradient_spher (const Param& para) const ;

	/**
	 * Calculates, if needed, the covariant derivative of {\tt *this}, 
	 * with respect to {\tt met}.
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 * The result is in {\tt *p\_derive\_cov}
	 */
	virtual void fait_derive_cov (const Param& para, const 
				      Qmetrique& met) const ;

	/**
	 * Calculates, if needed, the contravariant derivative of {\tt *this},
	 * with respect to {\tt met}.
	 * The result is in {\tt *p\_derive\_con}
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	virtual void fait_derive_con (const Param& para, const 
				      Qmetrique& met) const ;

	/**
	 * Calculates, if needed, the scalar square of {\tt *this},
	 * with respect to {\tt met}.
	 * The result is in {\tt *p\_carre\_scal}
	 */
	void fait_carre_scal (const Qmetrique& met) const ;
	
	/**
	 * To be used to describe the fact that the derivatives members have
	 * been calculated with {\tt met}.
	 * 
	 * First it sets a null element of {\tt met\_depend} to 
	 * {\tt \&met} and puts {\tt this} in 
	 * the list of the dependancies of {\tt met}.
	 * 
	 */
	void set_dependance (const Qmetrique& met) const ;

    // Differential operators
    // ----------------------
    public:
	/**
	 * Returns the gradient of {\tt *this} (Cartesian coordinates)
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constany.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	const Qtenseur& gradient(const Param& para) const ; 
	
	/** Returns the gradient of {\tt *this} (Spherical coordinates)
	 *	(scalar field only). 
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	const Qtenseur& gradient_spher(const Param& para) const ; 
	
	/**
	 * @return the covariant derivative of {\tt *this}, with respect to 
	 * {\tt met}.
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	const Qtenseur& derive_cov (const Param& para, const 
				    Qmetrique& met) const ;

	/**
	 * @return the contravariant derivative of {\tt *this}, with respect to 
	 * {\tt met}.
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	const Qtenseur& derive_con (const Param& para, const Qmetrique& met) const ;	

	/**
	 * @return the scalar square of {\tt *this}, with respect to 
	 * {\tt met}.
	 */
	const Qtenseur& carre_scal (const Qmetrique& met) const ;

    // Friend classes 
    // ---------------
    friend class Qtenseur_sym ;
    friend class Qmetrique ;
    
    // Mathematical operators
    // ----------------------
    
    friend Qtenseur operator* (const Qtenseur&, const Qtenseur&) ; 
    friend Qtenseur operator% (const Qtenseur&, const Qtenseur&) ; 
    friend Qtenseur contract(const Qtenseur&, int id1, int id2) ;
    friend Qtenseur contract(const Qtenseur&, int id1, const Qtenseur&, 
			    int id2) ;
    friend Qtenseur flat_scalar_prod(const Qtenseur& t1, const Qtenseur& t2) ;
    friend Qtenseur flat_scalar_prod_desal(const Qtenseur& t1, 
					  const Qtenseur& t2) ;
    friend Qtenseur manipule(const Qtenseur&, const Qmetrique&, int idx) ;
    friend Qtenseur manipule(const Qtenseur&, const Qmetrique&) ;
 
};


/**
 * @name Qtenseur calculus
 */
//@{
/// Tensorial product.
Qtenseur operator*(const Qtenseur&, const Qtenseur&) ; 

/// Tensorial product with desaliasing.
Qtenseur operator%(const Qtenseur&, const Qtenseur&) ; 

/**
 * Self contraction of two indices of a {\tt Qtenseur}.
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
Qtenseur contract(const Qtenseur&, int id1, int id2) ;

/**
 * Contraction of two {\tt Qtenseur}.
 * 
 * The two indices must be of different type, i.e. covariant and
 * contravariant, or contravariant and covariant.
 *
 * @param id1 [input] number of the index of contraction for
 *		      the first {\tt Qtenseur};
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
Qtenseur contract(const Qtenseur&, int id1, const Qtenseur&, int id2) ;

/**
 *  Scalar product of two {\tt Qtenseur} when the metric is
 *  $\delta_{ij}$: performs the contraction of the 
 *  last index of {\tt t1} with the first one of {\tt t2}, irrespective
 *  of the type of these indices. 
 */
Qtenseur flat_scalar_prod(const Qtenseur& t1, const Qtenseur& t2) ;
	
/**
 *  Same as {\tt flat\_scalar\_prod} but with desaliasing. 
 */
Qtenseur flat_scalar_prod_desal(const Qtenseur& t1, const Qtenseur& t2) ;
	
/**
 * Raise or lower the index {\tt idx} depending on its type, using the
 * given {\tt Qmetrique}.
 */
Qtenseur manipule(const Qtenseur&, const Qmetrique&, int idx) ;

/**
 * Raise or lower all the indices, depending on their type,  using the given
 * {\tt Qmetrique}.
 */
Qtenseur manipule(const Qtenseur&, const Qmetrique&) ;
	
//@}

/**
 * @name Qtenseur mathematics
 */
    //@{
Qtenseur operator+(const Qtenseur& ) ;			/// + Qtenseur
Qtenseur operator-(const Qtenseur& ) ;			/// - Qtenseur
Qtenseur operator+(const Qtenseur&, const Qtenseur &) ;	/// Qtenseur + Qtenseur

/// Qtenseur + double (the {\tt Qtenseur} must be a scalar)
Qtenseur operator+(const Qtenseur&, double ) ;		

/// double + Qtenseur (the {\tt Qtenseur} must be a scalar)
Qtenseur operator+(double, const Qtenseur& ) ;		

/// Qtenseur + int (the {\tt Qtenseur} must be a scalar)
Qtenseur operator+(const Qtenseur&, int ) ;		

/// int + Qtenseur (the {\tt Qtenseur} must be a scalar)
Qtenseur operator+(int, const Qtenseur& ) ;		

Qtenseur operator-(const Qtenseur &, const Qtenseur &) ;	/// Qtenseur - Qtenseur

/// Qtenseur - double (the {\tt Qtenseur} must be a scalar)
Qtenseur operator-(const Qtenseur&, double ) ;		

/// double - Qtenseur (the {\tt Qtenseur} must be a scalar)
Qtenseur operator-(double, const Qtenseur& ) ;		

/// Qtenseur - int (the {\tt Qtenseur} must be a scalar)
Qtenseur operator-(const Qtenseur&, int ) ;		

/// int - Qtenseur (the {\tt Qtenseur} must be a scalar)
Qtenseur operator-(int, const Qtenseur& ) ;		

/// Qtenseur * double 
Qtenseur operator*(const Qtenseur&, double ) ;		

/// double * Qtenseur 
Qtenseur operator*(double, const Qtenseur& ) ;		

/// Qtenseur * int 
Qtenseur operator*(const Qtenseur&, int ) ;		

/// int * Qtenseur 
Qtenseur operator*(int, const Qtenseur& ) ;		

/// Qtenseur / Qtenseur ({\tt b} must be a scalar)
Qtenseur operator/(const Qtenseur& a, const Qtenseur& b) ;	

Qtenseur operator/(const Qtenseur&, double ) ;	/// Qtenseur / double

/// double / Qtenseur  (the {\tt Qtenseur} must be a scalar)
Qtenseur operator/(double, const Qtenseur &) ;		

Qtenseur operator/(const Qtenseur&, int ) ;		/// Qtenseur / int

/// int / Qtenseur  (the {\tt Qtenseur} must be a scalar)
Qtenseur operator/(int, const Qtenseur &) ;		

Qtenseur exp(const Qtenseur& ) ;		/// Exponential (for a scalar only)
Qtenseur log(const Qtenseur& ) ;		/// Neperian logarithm (for a scalar only)
Qtenseur sqrt(const Qtenseur& ) ;		/// Square root (for a scalar only)
Qtenseur abs(const Qtenseur& ) ;		/// Absolute value (for a scalar only)
Qtenseur pow(const Qtenseur&, int ) ;	/// Power (for a scalar only)
Qtenseur pow(const Qtenseur&, double ) ;	/// Power (for a scalar only)

    //@}





			//---------------------------------//
			//	class Qtenseur_sym	   //
			//---------------------------------//
			
/**
 * Class intended to describe 4D tensors with a symmetry on the two last 
 * indices. The storage and the calculations are different and quicker than 
 * with an usual {\tt Qtenseur}.
 * 
 * The valence must be >1.
 */
class Qtenseur_sym : public Qtenseur {

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
	Qtenseur_sym (const Map& map, int val, const Itbl& tipe, 
		     const Base_vect& triad_i) ;

	/** Standard constructor when all the indices are of the same type.
	 * 
	 * @param map   the mapping 
	 * @param val   valence of the tensor; must be greater or equal to 2.
	 * @param tipe  the type of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 * 
	 */
	Qtenseur_sym (const Map& map, int val, int tipe, 
		     const Base_vect& triad_i) ;

	Qtenseur_sym (const Qtenseur_sym&) ; /// Copy constructor

	/** Constructor from a {\tt Qtenseur}.
	 *  The symmetry is assumed to be true but not checked.
	 */
	explicit Qtenseur_sym (const Qtenseur&) ;
	
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
	Qtenseur_sym (const Map& map, const Base_vect& triad_i, FILE* fich) ;

	virtual ~Qtenseur_sym() ;    /// Destructor
	
    // Mutators / assignment
    // ---------------------
    public:
	/**
	 * Assignment from a {\tt Qtenseur}.
	 * 
	 * The symmetry is assumed but not checked.
	 */
	virtual void operator= (const Qtenseur&) ;
    

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
	 * each element of which must be 0, 1, 2 or 3 
	 * corresponding to respective indices. 
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
	 *	   Each element of this {\tt Itbl} is 0, 1, 2 or 3,  which 
	 *	   corresponds to respective  indices. 
	 */
	virtual Itbl donne_indices (int place) const ;
		
    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the gradient of {\tt *this}.
	 * The result is in {\tt *p\_gradient}
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	virtual void fait_gradient (const Param& para) const ;

	/**
	 * Calculates, if needed, the covariant derivative of {\tt *this}, with 
	 * respect to {\tt met}.
	 * The result is in {\tt *p\_derive\_cov}
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	virtual void fait_derive_cov (const Param& para, const 
				      Qmetrique& met) const ;

	/**
	 * Calculates, if needed, the contravariant derivative of {\tt *this},
	 * with respect to {\tt met}.
	 * The result is in {\tt *p\_derive\_con}
	 *
	 * @param para [input] parameters for calculating the time 
	 * derivative: \\
	 *  {\tt par.get\_double(0)} : [input] the time step, supposed
	 *  to be constant.\\
	 *  {\tt par.get\_qtenseur\_mod(0)} : [input] the value of the
	 *  field at previous time-step. If set to {\tt 0x0}, the
	 *  method will suppose the field had the same value.\\
	 *  {\tt par.get\_qtenseur\_mod(1)} : [input] the value of the
	 *  field two time-steps before. If set to {\tt 0x0}, the
	 *  method will use a first ordre scheme. It is supposed to be
	 *  {\tt 0x0} if the previous one is null.
	 * 
	 */
	virtual void fait_derive_con (const Param& para, const 
				      Qmetrique&) const ;
	
    // Mathematical operators
    // ----------------------

	/// Tensorial product.
	friend Qtenseur_sym operator* (const Qtenseur&, const Qtenseur_sym&) ; 

};



/**
 * @name Qtenseur\_sym mathematics
 */
    //@{
Qtenseur_sym operator+(const Qtenseur_sym& ) ;	/// + Qtenseur\_sym
Qtenseur_sym operator-(const Qtenseur_sym& ) ;	/// - Qtenseur\_sym

/// Qtenseur\_sym + Qtenseur\_sym
Qtenseur_sym operator+(const Qtenseur_sym&, const Qtenseur_sym &) ;	

/// Qtenseur\_sym - Qtenseur\_sym
Qtenseur_sym operator-(const Qtenseur_sym &, const Qtenseur_sym &) ;	

/// Qtenseur\_sym * double 
Qtenseur_sym operator*(const Qtenseur_sym&, double ) ;		

/// double * Qtenseur\_sym 
Qtenseur_sym operator*(double, const Qtenseur_sym& ) ;		

/// Qtenseur\_sym * int 
Qtenseur_sym operator*(const Qtenseur_sym&, int ) ;		

/// int * Qtenseur\_sym 
Qtenseur_sym operator*(int, const Qtenseur_sym& ) ;		

/// Qtenseur\_sym / Qtenseur ({\tt b} must be a scalar)
Qtenseur_sym operator/(const Qtenseur_sym& a, const Qtenseur& b) ;	

Qtenseur_sym operator/(const Qtenseur_sym&, double ) ;  /// Qtenseur\_sym / double

Qtenseur_sym operator/(const Qtenseur_sym&, int ) ;     /// Qtenseur\_sym / int

    //@}


#endif
