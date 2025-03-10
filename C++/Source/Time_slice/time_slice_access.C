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

 

/*
 * $Id$
 * $Log$
 * Revision 1.10  2016/12/05 16:18:19  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.9  2014/10/13 08:53:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:13:21  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2008/12/02 15:02:22  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.6  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.5  2004/04/05 11:52:36  j_novak
 * First operational (but not tested!) version of checks of Eintein equation.
 *
 * Revision 1.4  2004/04/01 16:09:02  j_novak
 * Trace of K_ij is now member of Time_slice (it was member of Time_slice_conf).
 * Added new methods for checking 3+1 Einstein equations (preliminary).
 *
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
#include <cassert>

// Lorene headers
#include "time_slice.h"

namespace Lorene {
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

const Scalar& Time_slice::trk() const {

  if ( ! (trk_evol.is_known(jtime)) ) {

    if ( k_uu_evol.is_known(jtime) )
      trk_evol.update( k_uu().trace(gam()), jtime, the_time[jtime] ) ;
    else 
      trk_evol.update( k_dd().trace(gam()), jtime, the_time[jtime] ) ;

  }

  return trk_evol[jtime] ; 

}










}
