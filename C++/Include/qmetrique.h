/*
 *  Definition of Lorene class Qmetrique
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


#ifndef __QMETRIQUE_H_
#define __QMETRIQUE_H_


/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/09/25 12:08:02  j_novak
 * Tensors can be stored in Param objects
 *
 * Revision 1.3  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
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
#include "qtenseur.h"
#include "metrique.h"

/**
 * 4D-metric handling.
 * 
 * This class is used to described 4-dimensionnal metrics. 
 * 
 * This class includes a dynamical calculation of the Christoffel symbols, 
 * the Ricci-curvature and the Ricci-scalar. Although this is 4D, it is 
 * based on a 3+1 decomposition of the metric in terms of lapse, shift 
 * and 3-metric. (see also the documentation for {\tt Qtenseur}).
 * 
 * @version #$Id$#
 */
class Qmetrique {

    // Data : 
    // -----
    protected:
	const Map* const mp ;	/// Reference mapping.
	int etat ;  ///Logical state {\tt (ETATZERO, ETATQCQ or ETATNONDEF)} 
	Tenseur* lapse ; /// A pointer on the lapse.
	Tenseur* shift ; /// A pointer on the shift.
	Metrique* gamij ; ///A pointer on the spatial metric.
	bool plat ; ///Flag for a flat metric
    	
	/**
	 * Pointer on the contravariant representation.
	 */
	mutable Qtenseur_sym* p_met_con ;

	/**
	 * Pointer on the covariant representation.
	 */
	mutable Qtenseur_sym* p_met_cov ;

	/**
	 * Pointer on the covariant representation at previous time-step.
	 */
	mutable Qtenseur_sym* p_met_cov_jm1 ;

	/**
	 * Pointer on the covariant representation 2 time-steps before.
	 */
	mutable Qtenseur_sym* p_met_cov_jm2 ;

    // Derived data : 
    // ------------
	/**
	 * Pointer on the Christoffel symbols.
	 */
	mutable Qtenseur_sym* p_gamma ;

	/**
	 * Pointer on the Christoffel symbols at previous time-step.
	 */
	mutable Qtenseur_sym* p_gamma_jm1 ;

	/**
	 * Pointer on the Christoffel symbols 2 time-steps before.
	 */
	mutable Qtenseur_sym* p_gamma_jm2 ;

	/**
	 * Pointer on the Ricci curvature.
	 */
	mutable Qtenseur_sym* p_ricci ;

	/**
	 * Pointer on the Ricci scalar.
	 */
	mutable Qtenseur* p_ricci_scal ;
	
	/**
	 * Pointer on the determinant.
	 */
	mutable Qtenseur* p_determinant ;
	
	/**
	 * Pointer on the dependancies, that means the array contains pointers
	 * on all the {\tt Qtenseur} whose derivative members have been calculated
	 * using {\tt *this}.
	 */
	const Qtenseur** dependances ;
	
    // Constructors - Destructor :
    // -------------------------
    public:
	/** Standard constructor.
	 *  Nothing is allocated but {\tt dependances}.
	 *  By default, the {\tt Qmetrique} is not flat.
	 */
	explicit Qmetrique (const Map&, bool plate = false) ;

	Qmetrique (const Qmetrique&) ;    /// Constructor by copy.

	/**
	 * Constructor from lapse, shift and 3-metric.
	 * This should be used in most cases.
	 *
	 * @param alpha [input] the lapse (valence 0).
	 * @param beta [input] the shift (valence 1, contravariant).
	 * @param gamma [input] the 3-metric
	 */
	Qmetrique(const Tenseur& alpha, const Tenseur& beta,
		  const Metrique& gamma, bool plate = false) ;

	/** Constructor from a {\tt Qtenseur\_sym} of {\tt valence} = 2.
	 *  One representation is allocated depending on the 
	 *  type of {\tt source}. By default, the {\tt Qmetrique} is not flat.
	 */
	explicit Qmetrique (const Qtenseur_sym& source, bool plate = false) ;

	/** Constructor from a file (see {\tt sauve(FILE* )}).
	 * 
	 * @param map  the mapping
	 * @param triad_i   vectorial basis (triad) with respect to which 
	 *		    the tensor components are defined. It will
	 *		    be checked that it coincides with the basis
	 *		    saved in the file.
	 * @param fich  file which has been created by 
	 *			    the function {\tt sauve(FILE* )}.
	 */
	Qmetrique(const Map& map, const Base_vect& triad_i, FILE* fich) ;

	~Qmetrique() ;		    /// Destructor
	
    // Memory management
    // -----------------
    protected:
	void del_t() ;		    /// Logical destructor
	void del_deriv() ; /// Logical destructor of the derivative members.

	/**
	 * Sets all the pointer on the derivative members to zero.
	 */
	void set_der_0x0() ;

	/**
	 * Delete all the derivative members of the {\tt Qtenseur} contained in
	 * {\tt dependances}. Those quantities had been previously 
	 * calculated using {\tt *this}.
	 */
	void del_dependances() ;
		
    // Mutators / assignment
    // ---------------------
    public:	
	/**
	 * Sets the logical state to {\tt ETATNONDEF} (undefined state).
	 * Everything is deallocated.
	 */
	void set_etat_nondef() ;

	/**
	 * Sets the logical state to {\tt ETATZERO} (zero state).
	 * Everything is deallocated.
	 */
	void set_etat_zero() ;

	/**
	 * Sets the logical state to {\tt ETATQCQ} (ordinary state).
	 * The contravariant representation is allocated.
	 */
	void set_etat_con_qcq() ;

	/**
	 * Sets the logical state to {\tt ETATQCQ} (ordinary state).
	 * The covariant representation is allocated.
	 */
	void set_etat_cov_qcq() ;
	
	///Assignment from another {\tt Qmetrique}.
	void operator= (const Qmetrique&) ; 

	/**
	 * Assignment from a {\tt Qtenseur\_sym} of {\tt valence =2}.
	 * The allocated representation depends on the type of {\tt t}.
	 * All the other members are deleted.
	 */
	void operator= (const Qtenseur_sym& t) ;
	
	/**
	 * Set the standard spectral basis on the allocated representations.
	 */
	void set_std_base() ;

	/**
	 * Advance of one time-step.
	 *
	 * The values of the covariant representation are stored for
	 * the evaluation of the time derivatives. The same for the
	 * Christoffel symbols.
	 *
	 * @param alpha [input] the lapse (valence 0).
	 * @param beta [input] the shift (valence 1, contravariant).
	 * @param gamma [input] the 3-metric
	 * @param dt [input] the time-step (supposed to be constant).
	 */
	void avance_temps(const Tenseur& alpha, const Tenseur& beta,
		  const Metrique& gamma, double dt) ;

    // Accessors
    // ---------
    public:
	/// Returns the contravariant representation.
	const Qtenseur_sym& con() const ; 
	/// Returns the covariant representation.
	const Qtenseur_sym& cov() const ; 	
	/// Returns the Christoffel symbols. {\tt dt} is the time step.
	const Qtenseur_sym& gamma(double dt) const ; 
	/// Returns the Ricci-curvature. {\tt dt} is the time step.
	const Qtenseur_sym& ricci(double dt) const ; 
	/// Returns the Ricci-scalar. {\tt dt} is the time step.
	const Qtenseur& ricci_scal(double dt) const ; 
	/// Returns the determinant.
	const Qtenseur& determinant() const ; 
	
	/// Returns a pointer on the mapping.
	const Map* get_mp() const{return mp ; } ; 
	int get_etat() const{return etat ;} ; /// Returns the logical state.
	bool is_flat() const {return plat;} ; ///Is the metric a flat one?
	void set_flat(bool plate) {plat = plate;} ;///Sets the flat flag
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file

	friend ostream& operator<<(ostream& , const Qmetrique & ) ; /// Display.

    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the contravariant representation.
	 * The result is in {\tt *p\_met\_con}.
	 */
	void fait_con() const ;

	/**
	 * Calculates, if needed, the covariant representation.
	 * The result is in {\tt *p\_met\_cov}.
	 */
	void fait_cov() const ;

	/**
	 * Calculates, if needed, the Christoffel symbols.
	 * The result is in {\tt *p\_gamma}.
	 * @param dt [input] the time-step to compute time derivatives.	 
	 */
	void fait_gamma(double dt) const ;

	/**
	 * Calculates, if needed, the Ricci-curvature.
	 * The result is in {\tt *p\_ricci}.
	 * @param dt [input] the time-step to compute time derivatives.	 
	 */
	void fait_ricci(double dt) const ;

	/**
	 * Calculates, if needed, the Ricci-scalar.
	 * The result is in {\tt *p\_ricci\_scal}.
	 * @param dt [input] the time-step to compute time derivatives.	 
	 */
	void fait_ricci_scal(double dt) const ;

	/**
	 * Calculates, if needed, the determinant.
	 * The result is in {\tt *p\_determinant}.
	 */
	void fait_determinant() const ;

    // Friend classes 
    // ---------------

	friend class Qtenseur ;	/// Friend class {\tt Qtenseur}.

};
ostream& operator<<(ostream& , const Qmetrique & ) ; 

#endif
