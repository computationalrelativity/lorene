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
 * Revision 1.6  2004/04/05 21:25:02  e_gourgoulhon
 * -- Added constructor as standard time slice of Minkowski spacetime.
 * -- Added some calls to Scalar::std_spectral_base() after
 *    non-arithmetical operations.
 *
 * Revision 1.5  2004/04/05 12:38:45  j_novak
 * Minor modifs to prevent some warnings.
 *
 * Revision 1.4  2004/04/01 16:09:02  j_novak
 * Trace of K_ij is now member of Time_slice (it was member of Time_slice_conf).
 * Added new methods for checking 3+1 Einstein equations (preliminary).
 *
 * Revision 1.3  2004/03/29 12:00:41  e_gourgoulhon
 * Many modifs.
 *
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
                      aa_evol(aa_in, depth_in) {

    assert(hh_in.get_index_type(0) == CON) ; 
    assert(hh_in.get_index_type(1) == CON) ; 
    assert(aa_in.get_index_type(0) == CON) ; 
    assert(aa_in.get_index_type(1) == CON) ; 

    double time_init = the_time[jtime] ; 

    n_evol.update(lapse_in, jtime, time_init) ;     
    beta_evol.update(shift_in, jtime, time_init) ; 
    trk_evol.update(trk_in, jtime, time_init) ;
    
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
                      aa_evol(depth_in) {
                    
    set_der_0x0() ; // put here in order not to erase p_psi4

    double time_init = the_time[jtime] ; 

    p_psi4 = new Scalar( pow( p_gamma->determinant() / ff.determinant(), 
                         0.3333333333333333) ) ;
    p_psi4->std_spectral_base() ;

    Scalar tmp = pow(*p_psi4, 0.25) ;
    tmp.std_spectral_base() ;
    psi_evol.update(tmp , jtime, time_init ) ; 
    
    hh_evol.update( (*p_psi4) * p_gamma->con() - ff.con(), 
                    jtime, time_init ) ; 
    
    aa_evol.update( (*p_psi4) *( Time_slice::k_uu() 
                    - 0.3333333333333333 * trk_evol[jtime] * p_gamma->con() ), 
                    jtime, time_init ) ; 
                                
}

// Constructor as standard time slice of flat spacetime (Minkowski)
// ----------------------------------------------------------------

Time_slice_conf::Time_slice_conf(const Map& mp, const Base_vect& triad, 
                                 const Metric_flat& ff_in, int depth_in) 
        : Time_slice(mp, triad, depth_in),
          ff(ff_in), 
          psi_evol(depth_in), 
          qq_evol(depth_in),
          hh_evol(depth_in), 
          aa_evol(depth_in) {
    
    double time_init = the_time[jtime] ; 
    
    // Psi identically one:
    Scalar tmp(mp) ; 
    tmp.set_etat_one() ; 
    tmp.std_spectral_base() ;
    psi_evol.update(tmp, jtime, time_init) ; 
    
    // Q identically one:
    qq_evol.update(tmp, jtime, time_init) ; 
    
    // h^{ij} identically zero:
    Sym_tensor stmp(mp, CON, triad) ; 
    stmp.set_etat_zero() ; 
    hh_evol.update(stmp, jtime, time_init) ; 
    
    // A^{ij} identically zero:
    aa_evol.update(stmp, jtime, time_init) ; 

    set_der_0x0() ; 

}


// Copy constructor
// ----------------

Time_slice_conf::Time_slice_conf(const Time_slice_conf& tin) 
                    : Time_slice(tin), 
                      ff(tin.ff),
                      psi_evol(tin.psi_evol), 
                      qq_evol(tin.qq_evol),
                      hh_evol(tin.hh_evol), 
                      aa_evol(tin.aa_evol) {

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
    if (p_ln_psi != 0x0) delete p_ln_psi ; 
    if (p_hdirac != 0x0) delete p_hdirac ; 
    
    set_der_0x0() ;

    Time_slice::del_deriv() ; 
}


void Time_slice_conf::set_der_0x0() const {

    p_tgamma = 0x0 ; 
    p_psi4 = 0x0 ; 
    p_ln_psi = 0x0 ; 
    p_hdirac = 0x0 ; 
    
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
        
        Scalar tmp = sqrt( qq_evol[jtime] / n_evol[jtime] ) ; 
        tmp.std_spectral_base() ;
        psi_evol.update(tmp, jtime, the_time[jtime] ) ; 
    }

    return psi_evol[jtime] ;

} 

const Scalar& Time_slice_conf::psi4() const {

    if (p_psi4 == 0x0)  {

        p_psi4 = new Scalar( pow( psi(), 4.) ) ; 
        p_psi4->std_spectral_base() ;
    }

    return *p_psi4 ;

} 

const Scalar& Time_slice_conf::ln_psi() const {

    if (p_ln_psi == 0x0)  {

        p_ln_psi = new Scalar( log( psi() ) ) ; 
        p_ln_psi->std_spectral_base() ;
    }

    return *p_ln_psi ;

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

    if( !(aa_evol.is_known(jtime)) ) {

        assert( hh_evol.is_known(jtime) ) ; 

        Sym_tensor resu = hh_evol.time_derive(jtime, scheme_order)
            - hh().derive_lie(beta()) 
            - 0.6666666666666666 * beta().divergence(ff)
            + beta().ope_killing_conf(ff) ; 
        resu = resu / (2*nn()) ;
        
        aa_evol.update(resu, jtime, the_time[jtime]) ;
        
    }
     
    return aa_evol[jtime] ; 

}


const Scalar& Time_slice_conf::trk() const {

    if( !(trk_evol.is_known(jtime)) ) {

        psi() ; 
        Scalar resu = beta().divergence(ff) 
                + 6.*( contract(beta(), 0, ln_psi().derive_cov(ff), 0)
                - psi_evol.time_derive(jtime, scheme_order) / psi() ) ; 
        resu = resu / nn() ; 
        
        trk_evol.update(resu, jtime, the_time[jtime]) ;

    } 
    
    return trk_evol[jtime] ; 

}


const Vector& Time_slice_conf::hdirac() const {

    if (p_hdirac == 0x0) {
        p_hdirac = new Vector( hh().divergence(ff) ) ; 
    }
    
    return *p_hdirac ; 

}


                //------------------//
                //      output      //
                //------------------//

ostream& Time_slice_conf::operator>>(ostream& flux) const {

    Time_slice::operator>>(flux) ; 

    flux << "Triad on which the components of the flat metric are defined:\n" 
        << *(ff.get_triad()) << '\n' ;  

    if (psi_evol.is_known(jtime)) {
        maxabs( psi_evol[jtime], "Psi", flux) ;
    }
    if (qq_evol.is_known(jtime)) {
        maxabs( qq_evol[jtime], "Q", flux) ;
    }
    if (hh_evol.is_known(jtime)) {
        maxabs( hh_evol[jtime], "h^{ij}", flux) ;
    }
    if (aa_evol.is_known(jtime)) {
        maxabs( aa_evol[jtime], "A^{ij}", flux) ;
    }

    if (p_tgamma != 0x0) flux << 
        "Conformal metric tilde gamma is up to date" << endl ; 
    if (p_psi4 != 0x0) maxabs( *p_psi4, "Psi^4", flux) ; 
    if (p_ln_psi != 0x0) maxabs( *p_ln_psi, "ln(Psi)", flux) ; 
    if (p_hdirac != 0x0) maxabs( *p_hdirac, "H^i", flux) ; 
    
    return flux ; 

}



                
                
                
                

