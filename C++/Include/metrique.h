/*
 *  Definition of Lorene class Metrique
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


#ifndef __METRIQUE_H_
#define __METRIQUE_H_


/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.4  2000/02/09  19:28:54  eric
 * Ajout de l'argument triad dans le constructeur par lecture de fichier.
 *
 * Revision 2.3  2000/01/11  11:13:45  eric
 * Modif commentaires.
 *
 * Revision 2.2  2000/01/10  17:21:03  eric
 * Modif commentaires.
 *
 * Revision 2.1  1999/12/07  15:24:50  phil
 * ajout include
 *
 * Revision 2.0  1999/12/02  17:15:40  phil
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/02  17:13:23  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */
 
// Headers C++
#include <iostream.h>

// Headers Lorene 
#include "tenseur.h"

#define N_DEPEND 200

/**
 * Metric handling.
 * 
 * This class is used to described 3-dimensionnal metrics. 
 * It should be used ONLY for a metric decreasing as $\frac{1}{r^2}$ at 
 * infinity,  otherwise the calculation of the derivative terms will not be
 * correct in the external zone.
 * 
 * This class includes a dynamical calculation of the Christoffel symbols, 
 * the Ricci-curvature and the Ricci-scalar. 
 * 
 * @version #$Id$#
 */
class Metrique {

    // Data : 
    // -----
    private:
	const Map* const mp ;	/// Reference mapping.
	int etat ;  ///Logical state {\tt (ETATZERO, ETATQCQ or ETATNONDEF)}  
    	
	/**
	 * Pointer on the contravariant representation.
	 */
	mutable Tenseur_sym* p_met_con ;

	/**
	 * Pointer on the covariant representation.
	 */
	mutable Tenseur_sym* p_met_cov ;

    // Derived data : 
    // ------------
	/**
	 * Pointer on the Christoffel symbols.
	 */
	mutable Tenseur_sym* p_gamma ;

	/**
	 * Pointer on the Ricci curvature.
	 */
	mutable Tenseur_sym* p_ricci ;

	/**
	 * Pointer on the Ricci scalar.
	 */
	mutable Tenseur* p_ricci_scal ;
	
	/**
	 * Pointer on the dependancies, that means the array contains pointers
	 * on all the {\tt Tenseur} whom derivative members have been calculated
	 * using {\tt *this}.
	 */
	const Tenseur** dependances ;
	
    // Constructors - Destructor :
    // -------------------------
    public:
	/** Standard constructor.
	 *  Nothing is allocated but {\tt dependances}.
	 */
	explicit Metrique (const Map&) ;

	Metrique (const Metrique&) ;    /// Constructor by copy.

	/** Constructor from a {\tt Tenseur\_sym} of {\tt valence $= 2$}.
	 *  One representation is allocated depending on the 
	 *  type of {\tt source}.
	 */
	explicit Metrique (const Tenseur_sym& source) ;

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
	Metrique(const Map& map, const Base_vect& triad_i, FILE* fich) ;		

	~Metrique() ;		    /// Destructor
	
    // Memory management
    // -----------------
    private:
	void del_t() ;		    /// Logical destructor
	void del_gamma() ; /// Logical destructor of the derivative members.

	/**
	 * Sets all the pointer on the derivative members to zero.
	 */
	void set_der_0x0() ;

	/**
	 * Delete all the derivative members of the {\tt Tenseur} contained in
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
	
	///Assignment from another {\tt Metrique}.
	void operator= (const Metrique&) ; 

	/**
	 * Assignment from a {\tt Tenseur\_sym} of {\tt valence $=2$}.
	 * The allocated representation depends on the type of {\tt t}.
	 * All the other members are deleted.
	 */
	void operator= (const Tenseur_sym& t) ;
	
	/**
	 * Set the standard spectral basis on the allocated representations.
	 */
	void set_std_base() ;
    

    // Accessors
    // ---------
    public:
	const Tenseur_sym& con() const ; /// Returns the contravariant representation.
	const Tenseur_sym& cov() const ; /// Returns the covariant representation.
	const Tenseur_sym& gamma() const ; /// Returns the Christoffel symbols.
	const Tenseur_sym& ricci() const ; /// Returns the Ricci-curvature.
	const Tenseur& ricci_scal() const ; /// Returns the Ricci-scalar.
	
	const Map* get_mp() const{return mp ; } ; /// Returns a pointer on the mapping.
	int get_etat() const{return etat ;} ; /// Returns the logical state.
	
    // Outputs
    // -------
    public:
	void sauve(FILE *) const ;	    /// Save in a file

	friend ostream& operator<<(ostream& , const Metrique & ) ; /// Display.

    // Computation of derived members
    // ------------------------------
    private:
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
	 */
	void fait_gamma() const ;

	/**
	 * Calculates, if needed, the Ricci-curvature.
	 * The result is in {\tt *p\_ricci}.
	 */
	void fait_ricci() const ;

	/**
	 * Calculates, if needed, the Ricci-scalar.
	 * The result is in {\tt *p\_ricci\_scal}.
	 */
	void fait_ricci_scal() const ;

    // Friend classes 
    // ---------------

	friend class Tenseur ;	/// Friend class {\tt Tenseur}.

};

 /**
 * @name Utilities for {\tt Metrique}
 */
    //@{
    /**
     * Calculates the inverse of {\tt t},  being a {\tt Tenseur\_sym}
     * of {\tt valence $= 2$}.
     */
Tenseur_sym fait_inverse (const Tenseur_sym& t) ;

    //@}

#endif
