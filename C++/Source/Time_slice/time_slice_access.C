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

    assert( lapse.is_known(jtime) ) ; 
    return lapse[jtime] ; 

}


const Vector& Time_slice::beta() const {

    assert( shift.is_known(jtime) ) ; 
    return shift[jtime] ; 


}


const Metric& Time_slice::gam() const {

    if (p_gamma == 0x0) {
        if ( gamma_dd.is_known(jtime) ) {
            p_gamma = new Metric( gamma_dd[jtime] ) ; 
        }
        else {
            assert( gamma_uu.is_known(jtime) ) ; 
            p_gamma = new Metric( gamma_uu[jtime] ) ; 
        }    
    }
    
    return *p_gamma ; 

}


const Sym_tensor& Time_slice::gam_dd() const {

    if (!( gamma_dd.is_known(jtime)) ) {
      assert( gamma_uu.is_known(jtime) ) ; 
      gamma_dd.update(gam().cov(), jtime, the_time[jtime] ) ; 
    }

    return gamma_dd[jtime] ;

}

const Sym_tensor& Time_slice::gam_uu() const {

    if (!( gamma_uu.is_known(jtime)) ) {
      assert( gamma_dd.is_known(jtime) ) ; 
      gamma_uu.update(gam().con(), jtime, the_time[jtime] ) ; 
    }

    return gamma_uu[jtime] ;

}



const Sym_tensor& Time_slice::k_dd() const {

    if ( ! (kk_dd.is_known(jtime)) ) {
       
      gam_dd() ;
      
      Vector beta_d = beta().down(0, gam()) ;
      const Tensor& dbeta_dd = beta_d.derive_cov(gam()) ;

      Sym_tensor resu = - gamma_dd.time_derive(jtime, scheme_order) ;
      for (int i=1; i<=3; i++) 
	for (int j=i; j<=3; j++) 
	  resu.set(i,j) += dbeta_dd(i,j) + dbeta_dd(j,i) ;
      
      resu = resu / (2*nn()) ;
      kk_dd.update(resu, jtime, the_time[jtime]) ;
        
    }

    return kk_dd[jtime] ;

}

const Sym_tensor& Time_slice::k_uu() const {

    if ( ! (kk_uu.is_known(jtime)) ) {
       
      gam_uu() ;
      
      const Tensor& dbeta_uu = beta().derive_con(gam()) ;

      Sym_tensor resu = gamma_uu.time_derive(jtime, scheme_order) ;
      for (int i=1; i<=3; i++) 
	for (int j=i; j<=3; j++) 
	  resu.set(i,j) += dbeta_uu(i,j) + dbeta_uu(j,i) ;
      
      resu = resu / (2*nn()) ;
      kk_uu.update(resu, jtime, the_time[jtime]) ;
        
    }

    return kk_uu[jtime] ;

}










