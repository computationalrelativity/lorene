/*
 *  Methods of class Isol_hor
 *
 *    (see file isol_hor.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Jose Luis Jaramillo
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

char isol_hor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/11/03 17:16:06  f_limousin
 * Change the standart constructor. Add 4 memebers : trK, trK_point,
 * gamt and gamt_point.
 * Add also a constructor from a file.
 *
 * Revision 1.3  2004/10/29 15:44:45  jl_jaramillo
 * Remove two members
 *
 * Revision 1.2  2004/09/28 16:07:16  f_limousin
 * Remove all unused functions.
 *
 * Revision 1.1  2004/09/09 14:07:26  jl_jaramillo
 * First version
 *
 * Revision 1.1  2004/03/30 14:00:31  jl_jaramillo
 * New class Isol_hor (first version).
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers

#include "time_slice.h"
#include "isol_hor.h"
#include "tensor.h"
#include "metric.h"
#include "evolution.h"



			    //--------------//
			    // Constructors //
			    //--------------//


// Constructor from conformal decomposition
// ----------------------------------------

Isol_hor::Isol_hor(const Scalar& lapse_in, const Vector& shift_in,
		   const Scalar& psi_in, const Sym_tensor& aa_in, 
		   const Metric& met_gamt, const Sym_tensor& gamt_point_in, 
		   const Scalar& trK_in, const Scalar& trK_point_in,
		   const Metric_flat& ff_in, int depth_in) 	  
    : Time_slice_conf(lapse_in, shift_in, psi_in*psi_in*psi_in*psi_in*
		      met_gamt.cov(), psi_in*psi_in*psi_in*psi_in*aa_in 
		      +1./3.*trK*met_gamt.con()/(psi_in*psi_in*psi_in*psi_in),
		      ff_in, depth_in),
      gamt(met_gamt.cov()),
      gamt_point(gamt_point_in),
      trK(trK_in),
      trK_point(trK_point_in){
}




// Copy constructor
// ----------------

Isol_hor::Isol_hor(const Isol_hor& isolhor_in) 
    : Time_slice_conf(isolhor_in),
      gamt(isolhor_in.gamt),
      gamt_point(isolhor_in.gamt_point),
      trK(isolhor_in.trK),
      trK_point(isolhor_in.trK_point){
}
/*
// Constructor from a file
// -----------------------

Isol_hor::Isol_hor(const Map& mp, const Base_vect& triad, 
		   const Metric_flat& ff_in, FILE* fich, 
		   bool partial_read, int depth_in)
    : Time_slice_conf(mp, triad, ff_in, fich, partial_read, depth_in),
      gamt(mp, triad, fich),
      gamt_point(mp, triad, fich),
      trK(mp, *(mp.get_&mg()), fich),
      trK_point(mp, *(mp.get_mg()), fich){
}

*/

			    //--------------//
			    //  Destructor  //
			    //--------------//

Isol_hor::~Isol_hor(){}


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Isol_hor::operator=(const Isol_hor& isolhor_in) {

    Time_slice_conf::operator=(isolhor_in) ; 
}


                //------------------//
                //      output      //
                //------------------//


ostream& Isol_hor::operator>>(ostream& flux) const {

    Isol_hor::operator>>(flux) ; 
    return flux ; 

}


                //--------------------------//
                //      Save in a file      //
                //--------------------------//


void Isol_hor::sauve(FILE* fich, bool partial_save) const {


    // Writing of quantities common to all derived classes of Time_slice
    // -----------------------------------------------------------------
    
    Time_slice_conf::sauve(fich, true) ; 
    
    // Writing of quantities common to all derived classes of Isol_hor
    // ---------------------------------------------------------------

    gamt.sauve(fich) ;
    gamt_point.sauve(fich) ;    
    trK.sauve(fich) ;
    trK_point.sauve(fich) ;

}
