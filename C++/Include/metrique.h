/*
 *  Definition of Lorene class Metrique
 *
 */

/*
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


#ifndef __METRIQUE_H_
#define __METRIQUE_H_


/*
 * $Id$
 * $Log$
 * Revision 1.14  2014/10/13 08:52:35  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2003/10/17 13:04:18  f_limousin
 * Add new functions get_cov() and get_con().
 *
 * Revision 1.12  2003/06/20 14:22:20  f_limousin
 * The functions set_con() and set_cov() now returns a Tenseur&.
 *
 * Revision 1.11  2003/03/03 19:13:46  f_limousin
 * Add two new methods : set_con() and set_cov().
 *
 * Revision 1.10  2003/02/12 18:30:18  f_limousin
 * Added set_cov and set_con methods
 *
 * Revision 1.9  2003/02/07 15:41:46  e_gourgoulhon
 *
 * ### Added a PROVISORY set_cov method ###
 *
 * Revision 1.8  2002/10/16 14:36:29  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.7  2002/09/13 09:17:33  j_novak
 * Modif. commentaires
 *
 * Revision 1.6  2002/08/09 15:41:08  j_novak
 * New class Metconf added for conformal metric handling.
 *
 * Revision 1.5  2002/08/08 15:10:44  j_novak
 * The flag "plat" has been added to the class Metrique to show flat metrics.
 *
 * Revision 1.4  2002/08/07 16:14:10  j_novak
 * class Tenseur can now also handle tensor densities, this should be transparent to older codes
 *
 * Revision 1.3  2002/08/02 15:07:41  j_novak
 * Member function determinant has been added to the class Metrique.
 * A better handling of spectral bases is now implemented for the class Tenseur.
 *
 * Revision 1.2  2002/06/17 14:05:17  j_novak
 * friend functions are now also declared outside the class definition
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
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
 
// Headers Lorene 
#include "tenseur.h"

#define N_DEPEND 200

namespace Lorene {
/**
 * Metric handling.
 * 
 * This class is used to described 3-dimensionnal metrics. 
 * 
 * This class includes a dynamical calculation of the Christoffel symbols, 
 * the Ricci-curvature and the Ricci-scalar. 
 * 
 * @version #$Id$#
 */
class Metrique {

    // Data : 
    // -----
    protected:
	const Map* const mp ;	/// Reference mapping.
	int etat ;  ///Logical state {\tt (ETATZERO, ETATQCQ or ETATNONDEF)} 
	bool plat ; ///Flag for a flat metric
    	
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
	 * Pointer on the determinant.
	 */
	mutable Tenseur* p_determinant ;
	
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
	 *  By default, the {\tt Metrique} is not flat.
	 */
	explicit Metrique (const Map&, bool plate = false) ;

	Metrique (const Metrique&) ;    /// Constructor by copy.

	/** Constructor from a {\tt Tenseur\_sym} of {\tt valence} = 2.
	 *  One representation is allocated depending on the 
	 *  type of {\tt source}. By default, the {\tt Metrique} is not flat.
	 */
	explicit Metrique (const Tenseur_sym& source, bool plate = false) ;

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

	virtual ~Metrique() ;		    /// Destructor
	
    // Memory management
    // -----------------
    protected:
	void del_t() ;		    /// Logical destructor
	virtual void del_deriv() ; /// Logical destructor of the derivative members.

	/**
	 * Sets all the pointer on the derivative members to zero.
	 */
	virtual void set_der_0x0() ;

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
	 * Assignment from a {\tt Tenseur\_sym} of {\tt valence =2}.
	 * The allocated representation depends on the type of {\tt t}.
	 * All the other members are deleted.
	 */
	virtual void operator= (const Tenseur_sym& t) ;
	
	/**
	 * Set the standard spectral basis on the allocated representations.
	 */
	void set_std_base() ;
	
	/// Sets a component (covariant representation)
	Cmp& set_cov(int i, int j) ; 
 
	/// Sets a component (contravariant representation)
	Cmp& set_con(int i, int j) ; 

	/// Gets a component (covariant representation)
	const Cmp& get_cov(int i, int j) const ; 
 
	/// Gets a component (contravariant representation)
	const Cmp& get_con(int i, int j) const ; 
		
	/// Sets the covariant representation of the metric.
        Tenseur_sym& set_cov() ;
    
	/// Sets the contravariant representation of the metric.
	Tenseur_sym& set_con() ;

	/// Gets the covariant representation of the metric.
        const Tenseur_sym& get_cov() const ;
    
	/// Gets the contravariant representation of the metric.
	const Tenseur_sym& get_con() const ;

 
    // Accessors
    // ---------
    public:
	/// Returns the contravariant representation.
	const Tenseur_sym& con() const ; 
	/// Returns the covariant representation.
	const Tenseur_sym& cov() const ; 	
	/// Returns the Christoffel symbols.
	const Tenseur_sym& gamma() const ; 
	/// Returns the Ricci-curvature.
	const Tenseur_sym& ricci() const ; 
	/// Returns the Ricci-scalar.
	const Tenseur& ricci_scal() const ; 
	/// Returns the determinant.
	const Tenseur& determinant() const ; 
	
	/// Returns a pointer on the mapping.
	const Map* get_mp() const{return mp ; } ; 
	int get_etat() const{return etat ;} ; /// Returns the logical state.
	bool is_flat() const {return plat;} ; ///Is the metric a flat one?
	void set_flat(bool plate) {plat = plate;} ;///Sets the flat flag
    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file

	friend ostream& operator<<(ostream& , const Metrique & ) ; /// Display.

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
	 */
	virtual void fait_gamma() const ;

	/**
	 * Calculates, if needed, the Ricci-curvature.
	 * The result is in {\tt *p\_ricci}.
	 */
	virtual void fait_ricci() const ;

	/**
	 * Calculates, if needed, the Ricci-scalar.
	 * The result is in {\tt *p\_ricci\_scal}.
	 */
	virtual void fait_ricci_scal() const ;

	/**
	 * Calculates, if needed, the determinant.
	 * The result is in {\tt *p\_determinant}.
	 */
	virtual void fait_determinant() const ;

    // Friend classes 
    // ---------------

	friend class Tenseur ;	/// Friend class {\tt Tenseur}.

};
ostream& operator<<(ostream& , const Metrique & ) ; 

 /**
 * @name Utilities for {\tt Metrique}
 */
    //@{
    /**
     * Calculates the inverse of {\tt t},  being a {\tt Tenseur\_sym}
     * of {\tt valence} = 2.
     */
Tenseur_sym fait_inverse (const Tenseur_sym& t) ;

    //@}

}
#endif
