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
		   const Sym_tensor& gamma_in, const Sym_tensor kk_in, 
		   const Metric_flat& ff_in, int depth_in) 	  
  : Time_slice_conf( lapse_in, shift_in, gamma_in, kk_in, ff_in, depth_in), 
    gam_point_evol(gamma_in, depth_in),
    gamt_point_evol(gamma_in, depth_in){

}
                 



// Copy constructor
// ----------------

Isol_hor::Isol_hor(const Isol_hor& isolhor_in) 
                    : Time_slice_conf(isolhor_in), 
                      gam_point_evol(isolhor_in.gam_point_evol), 
                      gamt_point_evol(isolhor_in.gamt_point_evol){

}
			    //--------------//
			    //  Destructor  //
			    //--------------//

Isol_hor::~Isol_hor(){

}


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Isol_hor::operator=(const Isol_hor& isolhor_in) {

    Time_slice_conf::operator=(isolhor_in) ; 
    gam_point_evol = isolhor_in.gam_point_evol ; 
    gamt_point_evol = isolhor_in.gamt_point_evol ; 
       
}


                //------------------//
                //      output      //
                //------------------//


ostream& Isol_hor::operator>>(ostream& flux) const {

    Isol_hor::operator>>(flux) ; 

    flux << "Isolated Horizon without gauge choices" << '\n' ;

    //    if (khi_evol.is_known(jtime)) {
    //        maxabs( khi_evol[jtime], "Khi", flux) ;
    //    }
    //    if (mu_evol.is_known(jtime)) {
    //        maxabs( mu_evol[jtime], "Mu", flux) ;
    //    }
    //    if (trh_evol.is_known(jtime)) {
    //        maxabs( trh_evol[jtime], "tr h", flux) ;
    //    }
    
    return flux ; 

}

