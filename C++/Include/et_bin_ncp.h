/*
 *  Definition of Lorene class Et_bin_ncp
 *
 */

/*
 *   Copyright (c) 2002 Francois Limousin
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

#ifndef __ET_BIN_NCP_H_ 
#define __ET_BIN_NCP_H_ 

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/12/09 10:45:06  f_limousin
 * Definition of class Et_bin_ncp
 *
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "etoile.h"
#include "metrique.h"
#include "metconf.h"

/**
 * Class for stars in binary systems with a full spatial metric. 
 * 
 * This class is a derived class from {\tt Etoile\_bin }
 *
 * @version #$Id$#
 */
class Et_bin_ncp : public Etoile_bin {

    // Data : 
    // -----
    protected:

    ///  3-metric 
    Metrique gamma ;

    ///  flat 3-metric
    Metrique flat ;

    /// Conformal 3-metric (tensor density of weight -2/3)
    Metconf gamma_tilde ;


    // Constructors - Destructor
    // -------------------------
    public:
	/** Standard constructor. 
	 * 
	 * @param mp_i Mapping on which the star will be defined
	 * @param nzet_i Number of domains occupied by the star
	 * @param relat should be {\tt true} for a relativistic
	 *			star,  {\tt false} for a Newtonian one
	 * @param eos_i Equation of state of the stellar matter
	 * @param irrot should be {\tt true} for an irrotational star, 
	 *		    {\tt false} for a corotating one
	 * @param ref_triad_i  Reference triad ("absolute frame"), with
	 *	    respect to which the components of all the member 
	 *	    {\tt Tenseur}'s are defined, except for {\tt w\_shift}
	 *	    and {\tt ssjm1\_wshift} whose components are defined
	 *	    with respect to the mapping {\tt mp} Cartesian triad. 
	 * @param flat_i flat 3-metric 
	 */
         Et_bin_ncp(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
		    bool irrot, const Base_vect& ref_triad_i, 
		    const Metrique& flat_i) ;
		  
        Et_bin_ncp(const Et_bin_ncp& ) ;		/// Copy constructor

	/** Constructor from a file (see {\tt sauve(FILE* )}). 
	 *   
	 * @param mp_i Mapping on which the star will be defined
	 * @param eos_i Equation of state of the stellar matter
	 * @param ref_triad_i  Reference triad ("absolute frame"), with
	 *	    respect to which the components of all the member 
	 *	    {\tt Tenseur}'s are defined, except for {\tt w\_shift}
	 *	    and {\tt ssjm1\_wshift} whose components are defined
	 *	    with respect to the mapping {\tt mp} Cartesian triad. 
	 * @param fich	input file (must have been created by the function
	 *	{\tt sauve})
	 */
	Et_bin_ncp(Map& mp_i, const Eos& eos_i, const Base_vect& ref_triad_i, const Metrique& flat, FILE* fich) ;    		


	virtual ~Et_bin_ncp() ;			/// Destructor
 

    // Mutators / assignment
    // ---------------------
    public:
	/// Assignment to another {\tt Et\_bin\_ncp}
	void operator=(const Et_bin_ncp& ) ;	
	
    // Accessors
    // ---------
    public:

	/// Retunrs the 3-metric
	const Metrique& get_gamma() const {return gamma;} ;

	/// Retunrs the flat 3-metric
	const Metrique& get_flat() const {return flat;} ;

	/// Returns the conformal 3-metric (tensor density of weight -2/3)
	const Metconf& get_gamma_tilde() const {return gamma_tilde;} ;

    // Outputs
    // -------
    public:
	virtual void sauve(FILE *) const ;	    /// Save in a file

   protected:
	/// Operator >> (virtual function called by the operator <<). 
	virtual ostream& operator>>(ostream& ) const ;    


};

#endif
