/*
 *  Methods of class Time_slice to access the various fields
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon, Jose Luis Jaramillo & Jerome Novak
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

char time_slice_access_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/03/29 12:00:16  e_gourgoulhon
 * Computation of extrinsic curvature now performed via new methods
 *  Vector::ope_killing.
 *
 * Revision 1.2  2004/03/28 21:29:45  e_gourgoulhon
 * Evolution_std's renamed with suffix "_evol"
 * Method gam() modified
 * Added special constructor for derived classes.
 *
 * Revision 1.1  2004/03/26 13:33:02  j_novak
 * New methods for accessing/updating members (nn(), beta(), gam_uu(), k_uu(), ...)
 *
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
#include "tensor.h"
#include "metric.h"
#include "evolution.h"


const Scalar& Time_slice::nn() const {

    assert( n_evol.is_known(jtime) ) ; 
    return n_evol[jtime] ; 

}


const Vector& Time_slice::beta() const {

    assert( beta_evol.is_known(jtime) ) ; 
    return beta_evol[jtime] ; 


}

const Metric& Time_slice::gam() const {

    if (p_gamma == 0x0) {
        gam_dd() ; // may force the computation of p_gamma
        if (p_gamma == 0x0) p_gamma = new Metric( gam_dd() ) ; 
    }
    
    return *p_gamma ; 

}


const Sym_tensor& Time_slice::gam_dd() const {

    if (!( gam_dd_evol.is_known(jtime)) ) {
        assert( gam_uu_evol.is_known(jtime) ) ; 
        if (p_gamma == 0x0) {
            p_gamma = new Metric( gam_uu_evol[jtime] ) ; 
        }
        
        gam_dd_evol.update(p_gamma->cov(), jtime, the_time[jtime] ) ; 
    }

    return gam_dd_evol[jtime] ;

}

const Sym_tensor& Time_slice::gam_uu() const {

    if (!( gam_uu_evol.is_known(jtime)) ) {
      assert( gam_dd_evol.is_known(jtime) ) ; 
      gam_uu_evol.update(gam().con(), jtime, the_time[jtime] ) ; 
    }

    return gam_uu_evol[jtime] ;

}



const Sym_tensor& Time_slice::k_dd() const {

    if ( ! (k_dd_evol.is_known(jtime)) ) {
       
      Vector beta_d = beta().down(0, gam()) ;

      gam_dd() ; // to make sure that gam_dd is up to date before taking its
                 // time derivative
      
      Sym_tensor resu = beta_d.ope_killing(gam()) 
                        - gam_dd_evol.time_derive(jtime, scheme_order) ; 
            
      resu = resu / (2*nn()) ;

      k_dd_evol.update(resu, jtime, the_time[jtime]) ;
        
    }

    return k_dd_evol[jtime] ;

}

const Sym_tensor& Time_slice::k_uu() const {

    if ( ! (k_uu_evol.is_known(jtime)) ) {
       
      gam_uu() ; // to make sure that gam_uu is up to date before taking its
                 // time derivative
      
      Sym_tensor resu =  beta().ope_killing(gam())
                        + gam_uu_evol.time_derive(jtime, scheme_order) ;
            
      resu = resu / (2*nn()) ;
      
      k_uu_evol.update(resu, jtime, the_time[jtime]) ;
        
    }

    return k_uu_evol[jtime] ;

}










