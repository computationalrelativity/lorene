/*
 *  Definition of Lorene class Sym_tensor, 
 *  as well as derived classes Sym_tensor_trans and Sym_tensor_tt
 *
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon & Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

#ifndef __SYM_TENSOR_H_ 
#define __SYM_TEBSOR_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/11/03 22:29:54  e_gourgoulhon
 * Class Sym_tensor_tt: added functions set_eta_mu and update_tp.
 *
 * Revision 1.3  2003/11/03 17:09:30  e_gourgoulhon
 * Class Sym_tensor_tt: added the methods eta() and mu().
 *
 * Revision 1.2  2003/10/28 21:22:51  e_gourgoulhon
 * Class Sym_tensor_trans: added methods trace() and tt_part().
 *
 * Revision 1.1  2003/10/27 10:45:19  e_gourgoulhon
 * New derived classes Sym_tensor_trans and Sym_tensor_tt.
 *
 *
 * $Header$
 *
 */

class Sym_tensor_trans ;
class Sym_tensor_tt ;




			//---------------------------------//
			//        class Sym_tensor         //
			//---------------------------------//
			
/**
 * Class intended to describe valence-2 symmetric tensors.
 * The storage and the calculations are different and quicker than with an 
 * usual {\tt Tensor}.
 * 
 * The valence must be 2.
 *
 * @version #$Id$#
 */
class Sym_tensor : public Tensor {

    // Constructors - Destructor :
    // -------------------------
	
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param tipe  1-D array of integers (class {\tt Itbl}) of size 2 
	 *		containing the type 
	 *		of each index, {\tt COV} for a covariant one 
	 *		and {\tt CON} for a contravariant one,  with the 
	 *		following storage convention: \\
	 *			{\tt tipe(0)} : type of the first index \\
	 *			{\tt tipe(1)} : type of the second index 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 */
	Sym_tensor(const Map& map, const Itbl& tipe,const Base_vect& triad_i) ;

	/** Standard constructor when both indices are of the same type.
	 * 
	 * @param map   the mapping 
	 * @param tipe  the type of the indices.
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined
	 * 
	 */
	Sym_tensor(const Map& map, int tipe, const Base_vect& triad_i) ;

	Sym_tensor(const Sym_tensor&) ; /// Copy constructor

	/** Constructor from a {\tt Tensor}.
	 *  The symmetry of the input tensor is assumed to be true but not checked.
	 */
	Sym_tensor(const Tensor&) ;
	
	/** Constructor from a file (see {\tt sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param fich  file which has been used by 
	 *			    the function {\tt sauve(FILE* )}.
	 */
	Sym_tensor(const Map& map, const Base_vect& triad_i, FILE* fich) ;

	virtual ~Sym_tensor() ;    /// Destructor
	
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	/// Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to a {\tt Sym\_tensor}
	virtual void operator=(const Sym_tensor&) ;

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
	 *		 of size 2 giving the 
	 *		values of each index specifing the component,  with the 
	 *		following storage convention: \\
	 *			{\tt ind(0)} : value of the first index (1, 2 or 3) \\
	 *			{\tt ind(1)} : value of the second index (1, 2 or 3) 
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
	 *         size 2 giving the value of each index 
	 *	   for the component located at the position {\tt pos} in
	 *		the array [\tt cmp}, with the 
	 *		following storage convention: \\
	 *			{\tt Itbl(0)} : value of the first index (1, 2 or 3) \\
	 *			{\tt Itbl(1)} : value of the second index (1, 2 or 3) 
	 */
	virtual Itbl indices(int pos) const ;
	
		
	/**Returns the divergence of {\tt this} with respect to a {\tt Metric}.
	 * The indices are assumed to be contravariant.
	 */
	const Vector& divergence(const Metric&) const ; 

    // Computation of derived members
    // ------------------------------
	//    protected:

    // Mathematical operators
    // ----------------------
 protected:
	/**
	 * Returns a pointer on the inverse of the {\tt Sym\_tensor} 
	 * (seen as a matrix).
	 */
	Sym_tensor* inverse() const ;

    // Friend classes
    //-----------------
	friend class Metric ;
 
} ;


			//---------------------------------//
			//    class Sym_tensor_trans       //
			//---------------------------------//
			

/**
 * Transverse symmetric tensors of rank 2.
 *
 * This class is designed to store transverse (divergence-free) 
 * symmetric contravariant tensors of rank 2,
 * with the component expressed in an orthonormal spherical basis
 * $(e_r,e_\theta,e_\varphi)$.
 *
 * 
 * @version #$Id$#
 */
class Sym_tensor_trans: public Sym_tensor {

    // Data : 
    // -----
    protected:
	/// Metric with respect to which the divergence and the trace are defined
	const Metric* const met_div ; 
	
	/// Trace with respect to the metric {\tt *met\_div} 
	mutable Scalar* p_trace ; 
	
	/// Traceless part with respect to the metric {\tt *met\_div} 
	mutable Sym_tensor_tt* p_tt ;
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined 
	 * @param met the metric with respect to which the divergence is defined
	 */
	Sym_tensor_trans(const Map& map, const Base_vect& triad_i, 
		const Metric& met) ;

	Sym_tensor_trans(const Sym_tensor_trans& ) ;       /// Copy constructor

	/** Constructor from a file (see {\tt Tensor::sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param met the metric with respect to which the divergence is defined
	 * @param fich  file which has been used by 
	 *			    the function {\tt sauve(FILE* )}.
	 */
	Sym_tensor_trans(const Map& map, const Base_vect& triad_i, 
		const Metric& met, FILE* fich) ;

	virtual ~Sym_tensor_trans() ;			/// Destructor

 
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	/// Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------

	public:
	/// Assignment from another {\tt Sym\_tensor\_trans}
	virtual void operator=(const Sym_tensor_trans&) ;	
	
	/// Assignment from a {\tt Sym_tensor}
	virtual void operator=(const Sym_tensor&) ;	
	
	/// Assignment from a {\tt Tensor}
	virtual void operator=(const Tensor&) ;	
	
	// Computational methods
	// ---------------------
	/// Returns the trace of the tensor with respect to metric {\tt *met\_div}
	const Scalar& trace() const ; 
	
	/** Returns the transverse traceless part of the tensor, the trace being defined
	 * with respect to metric {\tt *met\_div}
	 */
	const Sym_tensor_tt& tt_part() const ; 
	
} ; 
	

			//------------------------------//
			//    class Sym_tensor_tt       //
			//------------------------------//
			

/**
 * Transverse and traceless symmetric tensors of rank 2.
 *
 * This class is designed to store transverse (divergence-free) 
 * and transverse symmetric contravariant tensors of rank 2,
 * with the component expressed in an orthonormal spherical basis
 * $(e_r,e_\theta,e_\varphi)$.
 *
 * 
 * @version #$Id$#
 */
class Sym_tensor_tt: public Sym_tensor_trans {

    // Data : 
    // -----

    protected:
	/** Field $\eta$ such that the components $(h^{r\theta}, h^{r\varphi})$
	 * of the tensor are written:
	 * \begin{equation}
	 *	h^{r\theta} =  {1\over r} \left( {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi} \right) 
	 * \end{equation} 
	 * \begin{equation}
	 *	h^{r\varphi} =  {1\over r} \left( {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} \right)
	 * \end{equation} 
	 */
	mutable Scalar* p_eta ;
	
	/** Field $\mu$ such that the components $(h^{r\theta}, h^{r\varphi})$
	 * of the tensor are written:
	 * \begin{equation}
	 *	h^{r\theta} =  {1\over r} \left( {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi} \right) 
	 * \end{equation} 
	 * \begin{equation}
	 *	h^{r\varphi} =  {1\over r} \left( {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} \right)
	 * \end{equation} 
	 */
	mutable Scalar* p_mu ;
	
	
    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor.
	 * 
	 * @param map   the mapping 
	 * @param triad_i  vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined 
	 * @param met the metric with respect to which the divergence is defined
	 */
	Sym_tensor_tt(const Map& map, const Base_vect& triad_i, 
		const Metric& met) ;

	Sym_tensor_tt(const Sym_tensor_tt& ) ;       /// Copy constructor

	/** Constructor from a file (see {\tt Tensor::sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *			  the tensor components are defined. It will
	 *			  be checked that it coincides with the basis
	 *			  saved in the file.
	 * @param met the metric with respect to which the divergence is defined
	 * @param fich  file which has been used by 
	 *			    the function {\tt sauve(FILE* )}.
	 */
	Sym_tensor_tt(const Map& map, const Base_vect& triad_i, 
		const Metric& met, FILE* fich) ;

	virtual ~Sym_tensor_tt() ;			/// Destructor

 
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() const;	/// Deletes the derived quantities

	/// Sets the pointers on derived quantities to 0x0
	void set_der_0x0() const ; 


    // Mutators / assignment
    // ---------------------

	public:
	/// Assignment from another {\tt Sym\_tensor\_tt}
	virtual void operator=(const Sym_tensor_tt&) ;	
	
	/// Assignment from a {\tt Sym_tensor\_trans}
	virtual void operator=(const Sym_tensor_trans&) ;	
	
	/// Assignment from a {\tt Sym_tensor}
	virtual void operator=(const Sym_tensor&) ;	
	
	/// Assignment from a {\tt Tensor}
	virtual void operator=(const Tensor&) ;	
	
	/** Sets the angular potentials $\eta$ and $\mu$ (see members
	 *  {\tt p\_eta} and {\tt p\_mu}). 
	 *  The components $h^{r\theta}$ and $h^{r\varphi}$ are updated consistently
	 *  by a call to the method {\tt update\_tp()}.
	 *
	 *	@param eta_i [input] angular potential $\eta$
	 *	@param mu_i [input] angular potential $\mu$
	 *
	 */
	void set_eta_mu(const Scalar& eta_i, const Scalar& mu_i) ; 
	
	
	// Computational methods
	// ---------------------
	
	/** Gives the field $\eta$ such that the components $(h^{r\theta}, h^{r\varphi})$
	 * of the tensor are written:
	 * \begin{equation}
	 *	h^{r\theta} =  {1\over r} \left( {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi} \right) 
	 * \end{equation} 
	 * \begin{equation}
	 *	h^{r\varphi} =  {1\over r} \left( {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} \right)
	 * \end{equation} 
	 */
	const Scalar& eta() const ;
	
	/** Gives the field $\mu$ such that the components $(h^{r\theta}, h^{r\varphi})$
	 * of the tensor are written:
	 * \begin{equation}
	 *	h^{r\theta} =  {1\over r} \left( {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi} \right) 
	 * \end{equation} 
	 * \begin{equation}
	 *	h^{r\varphi} =  {1\over r} \left( {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} \right)
	 * \end{equation} 
	 */
	const Scalar& mu() const ;
	

	/** Computes the components $(h^{r\theta}, h^{r\varphi})$ from the
	 *  potential $\eta$ and  $\mu$, according to:
	 * \begin{equation}
	 *	h^{r\theta} =  {1\over r} \left( {\partial \eta \over \partial\theta}
	 *		- {1\over\sin\theta} {\partial \mu \over \partial\varphi} \right) 
	 * \end{equation} 
	 * \begin{equation}
	 *	h^{r\varphi} =  {1\over r} \left( {1\over\sin\theta} 
	 *				{\partial \eta \over \partial\varphi}
	 *				+ {\partial \mu \over \partial\theta} \right)
	 * \end{equation} 
	 */
	void update_tp() ;
	

} ; 
	



#endif
