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
 * Revision 1.3  2003/01/20 09:31:56  f_limousin
 * Modification of the standard constructor
 *
 * Revision 1.2  2003/01/14 14:13:25  f_limousin
 * Binary NS with Nonconformally flat metric.
 *
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
		    const Metrique& flat_i,const Tenseur_sym &source) ;
		  
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


    // Computational routines
    // ----------------------
    public: 
	/** Performs the scalar product of two tensors by contracting
	 *  the last index of {\tt t1} with the first index of {\tt t2}.
	 *  Both indices are supposed to be contravariant, so that a 
	 *  multiplication by $A^2$ is performed to lower one index. 
	 *  For instance, for two vectors $V^i$ and $W^i$, this function
	 *  returns the scalar $h_{ij} V^i W^j = A^2 f_{ij} V^i W^j$.  
	 */
	virtual Tenseur sprod(const Tenseur& t1, const Tenseur& t2) const ; 


	/** Computes an equilibrium configuration.
	 * 
	 *  The values of {\tt logn\_comp}, {\tt beta\_comp}, {\tt pot\_centri}
	 *  are held fixed during the iteration. 
	 *  
	 *  @param ent_c  [input] Central enthalpy
	 *  @param mermax [input] Maximum number of steps 
	 *  @param mermax_poisson [input]   Maximum number of steps in 
	 *				    Map\_et::poisson
	 *  @param relax_poisson [input]  Relaxation factor in Map\_et::poisson
	 *  @param mermax_potvit [input]  Maximum number of steps in 
	 *				  Map\_radial::poisson\_compact
	 *  @param relax_potvit [input]   Relaxation factor in 
	 *				  Map\_radial::poisson\_compact
	 *  @param thres_adapt  [input]   Threshold on dH/dr for the adaptation 
	 *				  of the mapping
	 *  @param fact [input]    1-D {\tt Tbl} for the input
	 *                          of some factors : \\
	 *          {\tt fact(0)} : A resizing factor for the first shell
	 *  @param diff [output]   1-D {\tt Tbl} for the storage of some
	 *			    error indicators : \\
	 *	    {\tt diff(0)} : Relative change in the enthalpy field
	 *			      between two successive steps \\
	 *	    {\tt diff(1)} : Relative error returned by the routine
	 *				{\tt Etoile\_bin::velocity\_potential}  
	 *	    {\tt diff(2)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt logn\_auto} \\  
	 *	    {\tt diff(3)} : Relative error in the resolution of the
	 *			    Poisson equation for {\tt beta\_auto} \\  
	 *	    {\tt diff(4)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (x comp.) \\  
	 *	    {\tt diff(5)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (y comp.) \\  
	 *	    {\tt diff(6)} : Relative error in the resolution of the
	 *			    equation for {\tt shift\_auto} (z comp.)   
	 */
	virtual void equilibrium(double ent_c, int mermax, int mermax_poisson, 
			 double relax_poisson, int mermax_potvit, 
			 double relax_potvit, double thres_adapt, 
			 const Tbl& fact, Tbl& diff) ;	

};

#endif
