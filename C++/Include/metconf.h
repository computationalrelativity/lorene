/*
 *  Definition of Lorene class Metconf
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


#ifndef __METCONF_H_
#define __METCONF_H_

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/09/11 08:50:56  j_novak
 * Modifs. des commentaires
 *
 * Revision 1.2  2002/08/13 08:02:45  j_novak
 * Handling of spherical vector/tensor components added in the classes
 * Mg3d and Tenseur. Minor corrections for the class Metconf.
 *
 * Revision 1.1  2002/08/09 15:41:08  j_novak
 * New class Metconf added for conformal metric handling.
 *
 *
 * $Header$
 *
 */
 
// Headers C++
#include <iostream.h>

// Headers Lorene 
#include "metrique.h"

/**
 * Pseudo-metric handling.
 * 
 * This class is used to described 3-dimensionnal pseudo-metrics 
 * $\tilde{\gamma}_{ij} = \gamma^{-1/3}\gamma_{ij}$, with $\gamma_{ij}$
 * being a usual 3-metric. Thus the covariant representation is a tensor 
 * density of weight -2/3 and the contravariant a 2/3 one.
 * It is a child class of {\tt Metrique} and
 * includes a dynamical calculation of the Christoffel symbols, 
 * the Ricci-curvature and the Ricci-scalar. 
 * 
 * @version #$Id$#
 */
class Metconf: public Metrique {

    // Data : 
    // -----
    protected:
        const Metrique* gamij ;
	const Metrique* fij ;
	bool dirac ;

    // Derived data : 
    // ------------
	/**
	 * Pointer on the $\Delta^i_{jk}$'s
	 */
	mutable Tenseur_sym* p_delta ;

	/**
	 * Pointer on the divergence
	 */
	mutable Tenseur* p_Hi ;

    // Constructors - Destructor :
    // -------------------------
    public:
	/** Standard constructor.
	 *  Nothing is allocated but {\tt dependances}.
	 */
	explicit Metconf (const Map&, const Metrique& metric, 
			  const Metrique& metplat, bool jauge = false, 
			  bool plate = false) ;

	Metconf (const Metconf&) ;    /// Constructor by copy.

	/** Constructor from a {\tt Tenseur\_sym} of {\tt valence} = 2.
	 *  One representation is allocated depending on the 
	 *  type of {\tt source}.
	 */
	Metconf (const Tenseur_sym& source, const Metrique& metplat, 
	         bool jauge = false, bool plate = false) ;

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
	Metconf(const Map& map, const Base_vect& triad_i, const Metrique&
		metric, const Metrique& metplat, FILE* fich) ;

	virtual ~Metconf() ;		    /// Destructor
	
    // Memory management
    // -----------------
    protected:
	virtual void del_deriv() ; /// Logical destructor of the derivative members.

	/**
	 * Sets all the pointer on the derivative members to zero.
	 */
	virtual void set_der_0x0() ;

   // Mutators / assignment
    // ---------------------
    public:	
	///Assignment from another {\tt Metconf}.
	void operator= (const Metconf&) ; 

	/**
	 * Assignment from a {\tt Tenseur\_sym} of {\tt valence =2}.
	 * The allocated representation depends on the type of {\tt t}.
	 * All the other members are deleted.
	 */
	virtual void operator= (const Tenseur_sym& t) ;
	
    // Accessors
    // ---------
    public:
	const Tenseur_sym& delta() const ; ///Returns $\Delta^i_{jk}$
	///Returns the divergence of {\it this}, using the flat metric.
	const Tenseur& Hi() const ; 
	/// Returns a pointer on the flat metric.
	const Metrique* get_fij() const{return fij ; } ; 
	///Is the metric in the Dirac gauge?
	bool dirac_gauge() const {return dirac;} ; 
	void set_dirac_gauge(bool gauge) ; ///Sets the dirac flag
	void set_flat_metric(const Metrique& flat) ; ///Sets the flat metric

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file

    // Computation of derived members
    // ------------------------------
    protected:
	/**
	 * Calculates, if needed, the $\Delta^i_{jk}$.
	 * The result is in {\tt *p\_delta}.
	 */
	void fait_delta() const ;

	/**
	 * Calculates, if needed, the divergence of the pseudo metric.
	 * The result is in {\tt *p\_Hi}.
	 */
	void fait_Hi() const ;

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
	 * The determinant is always one for the conformal metric.
	 * This result is in {\tt *p\_determinant}.
	 */
	virtual void fait_determinant() const ;

    // Friend classes 
    // ---------------

	friend class Tenseur ;	/// Friend class {\tt Tenseur}.

};

#endif
