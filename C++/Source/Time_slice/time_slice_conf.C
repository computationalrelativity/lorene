/*
 *  Methods of class Time_slice_conf
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

char time_slice_conf_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/03/28 21:32:23  e_gourgoulhon
 * Corrected error in method trk().
 *
 * Revision 1.1  2004/03/28 21:30:13  e_gourgoulhon
 * First version.
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



			    //--------------//
			    // Constructors //
			    //--------------//


// Constructor from conformal decomposition
// ----------------------------------------

Time_slice_conf::Time_slice_conf(const Scalar& lapse_in, const Vector& shift_in,
            const Metric_flat& ff_in, const Scalar& psi_in, 
            const Sym_tensor& hh_in, const Sym_tensor aa_in, 
            const Scalar& trk_in, int depth_in) 
                    : Time_slice(depth_in),
                      ff(ff_in),
                      psi_evol(psi_in, depth_in), 
                      qq_evol(depth_in),
                      hh_evol(hh_in, depth_in), 
                      aa_evol(aa_in, depth_in),
                      trk_evol(trk_in, depth_in) {

    assert(hh_in.get_index_type(0) == CON) ; 
    assert(hh_in.get_index_type(1) == CON) ; 
    assert(aa_in.get_index_type(0) == CON) ; 
    assert(aa_in.get_index_type(1) == CON) ; 

    double time_init = the_time[jtime] ; 

    n_evol.update(lapse_in, jtime, time_init) ;     
    beta_evol.update(shift_in, jtime, time_init) ;    
    
    set_der_0x0() ;  
    
}
                 

// Constructor from physical metric
// --------------------------------                 

Time_slice_conf::Time_slice_conf(const Scalar& lapse_in, const Vector& shift_in,
               const Sym_tensor& gamma_in, const Sym_tensor kk_in,
               const Metric_flat& ff_in, int depth_in) 
                    : Time_slice(lapse_in, shift_in, gamma_in, kk_in, depth_in),
                      ff(ff_in),
                      psi_evol(depth_in), 
                      qq_evol(depth_in),
                      hh_evol(depth_in), 
                      aa_evol(depth_in),
                      trk_evol(depth_in) {
                    
    Scalar psi04 = pow( (Time_slice::gam()).determinant() / ff.determinant(), 
                         0.3333333333333333) ;

    psi_evol.update( pow(psi04, 0.25), jtime, the_time[jtime] ) ; 
    
    hh_evol.update( psi04 * gam().con() - ff.con(), jtime, the_time[jtime] ) ; 
    
    trk_evol.update( kk_in.trace(gam()), jtime, the_time[jtime] ) ; 
    
    aa_evol.update( psi04 *( k_uu() - 0.3333333333333333 * trk_evol[jtime] 
                                * gam().con() ), jtime, the_time[jtime] ) ; 
     
    set_der_0x0() ; 
    
    p_psi4 = new Scalar( psi04 ) ; 
                       
}


// Copy constructor
// ----------------

Time_slice_conf::Time_slice_conf(const Time_slice_conf& tin) 
                    : Time_slice(tin), 
                      ff(tin.ff),
                      psi_evol(tin.psi_evol), 
                      qq_evol(tin.qq_evol),
                      hh_evol(tin.hh_evol), 
                      aa_evol(tin.aa_evol),
                      trk_evol(tin.trk_evol) {

    set_der_0x0() ; 
                       
}
			    //--------------//
			    //  Destructor  //
			    //--------------//

Time_slice_conf::~Time_slice_conf(){

    Time_slice_conf::del_deriv() ; 

}

                    //---------------------//
                    //  Memory management  //
                    //---------------------//

void Time_slice_conf::del_deriv() const {

    if (p_tgamma != 0x0) delete p_tgamma ; 
    if (p_psi4 != 0x0) delete p_psi4 ; 
    
    set_der_0x0() ;
}


void Time_slice_conf::set_der_0x0() const {

    p_tgamma = 0x0 ; 
    p_psi4 = 0x0 ; 
    
}


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Time_slice_conf::operator=(const Time_slice_conf& tin) {

    Time_slice::operator=(tin) ; 

    psi_evol = tin.psi_evol ; 
    qq_evol = tin.qq_evol ; 
    hh_evol = tin.hh_evol ; 
    aa_evol = tin.aa_evol ; 
    trk_evol = tin.trk_evol ; 
       
    del_deriv() ; 
    
}


                //-----------------------------------------------//
                //  Update of fields from base class Time_slice  //
                //-----------------------------------------------//

const Scalar& Time_slice_conf::nn() const {

    if (!( n_evol.is_known(jtime) ) ) {

        assert( psi_evol.is_known(jtime) ) ; 
        assert( qq_evol.is_known(jtime) ) ; 
        
        n_evol.update( qq_evol[jtime] / ( psi_evol[jtime]*psi_evol[jtime] ), 
                        jtime, the_time[jtime] ) ; 
    }

    return n_evol[jtime] ;

} 



const Sym_tensor& Time_slice_conf::gam_dd() const {

    if (!( gam_dd_evol.is_known(jtime)) ) {
        gam_dd_evol.update( psi4() * tgam().cov(), jtime, the_time[jtime] ) ; 
    }

    return gam_dd_evol[jtime] ;

}


const Sym_tensor& Time_slice_conf::gam_uu() const {

    if (!( gam_uu_evol.is_known(jtime)) ) {
        gam_uu_evol.update( tgam().con() / psi4() , jtime, the_time[jtime] ) ; 
    }

    return gam_uu_evol[jtime] ;

}


const Sym_tensor& Time_slice_conf::k_dd() const {

    if ( ! (k_dd_evol.is_known(jtime)) ) {
       
        k_dd_evol.update( k_uu().up_down(gam()), jtime, the_time[jtime] ) ; 
        
    }

    return k_dd_evol[jtime] ;

}


const Sym_tensor& Time_slice_conf::k_uu() const {

    if ( ! (k_uu_evol.is_known(jtime)) ) {
       
        k_uu_evol.update( aa()/psi4() + 0.3333333333333333* trk()* gam().con(),
             jtime, the_time[jtime] ) ; 
    }

    return k_uu_evol[jtime] ;

}





                //-----------------------------------//
                //  Update of fields from this class //
                //-----------------------------------//


const Scalar& Time_slice_conf::psi() const {

    if (!( psi_evol.is_known(jtime) ) ) {

        assert( n_evol.is_known(jtime) ) ; 
        assert( qq_evol.is_known(jtime) ) ; 
        
        psi_evol.update( sqrt( qq_evol[jtime] / n_evol[jtime] ), jtime, 
                         the_time[jtime] ) ; 
    }

    return psi_evol[jtime] ;

} 

const Scalar& Time_slice_conf::psi4() const {

    if (p_psi4 == 0x0)  {

        p_psi4 = new Scalar( pow( psi(), 4.) ) ; 
    }

    return *p_psi4 ;

} 


const Scalar& Time_slice_conf::qq() const {

    if (!( qq_evol.is_known(jtime) ) ) {
        
        assert( n_evol.is_known(jtime) ) ; 
        assert( psi_evol.is_known(jtime) ) ; 

        const Scalar& psij = psi_evol[jtime] ; 
        qq_evol.update( psij*psij * n_evol[jtime], jtime, the_time[jtime] ) ; 
    }

    return qq_evol[jtime] ;

}


const Metric& Time_slice_conf::tgam() const {

    if (p_tgamma == 0x0) {
        p_tgamma = new Metric( ff.con() + hh() ) ; 
    }
    
    return *p_tgamma ; 

}


const Sym_tensor& Time_slice_conf::hh() const {

    assert( hh_evol.is_known(jtime) ) ; 
    return hh_evol[jtime] ; 

}


const Sym_tensor& Time_slice_conf::aa() const {

    assert( aa_evol.is_known(jtime) ) ; 
    return aa_evol[jtime] ; 

}


const Scalar& Time_slice_conf::trk() const {

    assert( trk_evol.is_known(jtime) ) ; 
    return trk_evol[jtime] ; 

}




                
                
                
                

