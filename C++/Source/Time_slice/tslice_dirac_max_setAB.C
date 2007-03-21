/*
 *  Methods of class Tslice_dirac_max dealing with the members potA and tildeB
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2007  Jerome Novak
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

char tslice_dirax_max_setAB_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2007/03/21 14:51:50  j_novak
 * Introduction of potentials A and tilde(B) of h^{ij} into Tslice_dirac_max.
 *
 *
 * $Header $
 *
 */

// C headers
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "param.h"

void Tslice_dirac_max::set_A_tildeB(const Scalar& potA_in, const Scalar& tildeB_in, 
				    Param* par_bc, Param* par_mat) {

    potA_evol.update(potA_in, jtime, the_time[jtime]) ; 
    tildeB_evol.update(tildeB_in, jtime, the_time[jtime]) ; 
    
    // Computation of trace h and h^{ij} to ensure det tgam_{ij} = det f_{ij} :

    hh_det_one_AB(jtime, par_bc, par_mat) ;
    
} 

void Tslice_dirac_max::hh_det_one_AB(int j0, Param* par_bc, Param* par_mat) const {

    assert (potA_evol.is_known(j0)) ;   // The starting point
    assert (tildeB_evol.is_known(j0)) ;    // of the computation 

    const Map& mp = potA_evol[j0].get_mp() ;

    // The TT part of h^{ij}, which stays unchanged during the computation :
    Sym_tensor_tt hijtt(mp, *(ff.get_triad()), ff) ;
    hijtt.set_A_tildeB(potA_evol[j0], tildeB_evol[j0], par_bc, par_mat) ;
    
    // The representation of h^{ij} as an object of class Sym_tensor_trans :
    Sym_tensor_trans hij(mp, *(ff.get_triad()), ff) ;
    hij.trace_from_det_one(hijtt) ;

    // Result set to trh_evol and hh_evol
    // ----------------------------------
    trh_evol.update(hij.the_trace(), j0, the_time[j0]) ;
    
    // The longitudinal part of h^{ij}, which is zero by virtue of Dirac gauge :
    Vector wzero(mp, CON,  *(ff.get_triad())) ; 
    wzero.set_etat_zero() ;                   

    // Temporary Sym_tensor with longitudinal part set to zero : 
    Sym_tensor hh_new(mp, CON, *(ff.get_triad())) ;
    
    hh_new.set_longit_trans(wzero, hij) ;
    
    hh_evol.update(hh_new, j0, the_time[j0]) ;
    
    if (j0 == jtime) {
        // Reset of quantities depending on h^{ij}:
        if (p_tgamma != 0x0) {
            delete p_tgamma ;
            p_tgamma = 0x0 ; 
        } 
        if (p_hdirac != 0x0) {
            delete p_hdirac ; 
            p_hdirac = 0x0 ; 
        }
        if (p_gamma != 0x0) {
            delete p_gamma ; 
            p_gamma = 0x0 ;
        }
    }
    gam_dd_evol.downdate(j0) ; 
    gam_uu_evol.downdate(j0) ;
    adm_mass_evol.downdate(j0) ;  
         
    // Test
    if (j0 == jtime) {
        maxabs(tgam().determinant() - 1, 
        "Max. of absolute value of deviation from det tgam = 1") ; 
    }
    else {
        Metric tgam_j0( ff.con() + hh_evol[j0] ) ; 
        maxabs(tgam_j0.determinant() - 1, 
        "Max. of absolute value of deviation from det tgam = 1") ; 
    }

}
