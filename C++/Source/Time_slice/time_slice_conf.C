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
 * Revision 1.15  2008/12/02 15:02:22  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.14  2004/06/24 20:36:54  e_gourgoulhon
 * Added method check_psi_dot.
 *
 * Revision 1.13  2004/05/31 09:08:18  e_gourgoulhon
 * Method sauve and constructor from binary file are now operational.
 *
 * Revision 1.12  2004/05/27 15:25:04  e_gourgoulhon
 * Added constructors from binary file, as well as corresponding
 * functions sauve and save.
 *
 * Revision 1.11  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.10  2004/05/10 09:10:05  e_gourgoulhon
 * Added "adm_mass_evol.downdate(jtime)" in methods set_*.
 *
 * Revision 1.9  2004/05/05 14:31:14  e_gourgoulhon
 * Method aa(): added *as a comment*  annulation of hh_point in the compactified
 * domain.
 *
 * Revision 1.8  2004/05/03 14:47:11  e_gourgoulhon
 * Corrected method aa().
 *
 * Revision 1.7  2004/04/08 16:43:26  e_gourgoulhon
 * Added methods set_*
 * Added test of determinant one in constructor and set_hh.
 *
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
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "utilitaires.h"



			    //--------------//
			    // Constructors //
			    //--------------//


// Constructor from conformal decomposition
// ----------------------------------------

Time_slice_conf::Time_slice_conf(const Scalar& lapse_in, const Vector& shift_in,
            const Metric_flat& ff_in, const Scalar& psi_in, 
            const Sym_tensor& hh_in, const Sym_tensor& hata_in, 
            const Scalar& trk_in, int depth_in) 
                    : Time_slice(depth_in),
                      ff(ff_in),
                      psi_evol(psi_in, depth_in), 
                      npsi_evol(depth_in),
                      hh_evol(hh_in, depth_in), 
                      hata_evol(hata_in, depth_in),
		      A_hata_evol(depth_in), B_hata_evol(depth_in){

    assert(hh_in.get_index_type(0) == CON) ; 
    assert(hh_in.get_index_type(1) == CON) ; 
    assert(hata_in.get_index_type(0) == CON) ; 
    assert(hata_in.get_index_type(1) == CON) ; 

    double time_init = the_time[jtime] ; 

    // Check whether det tgam^{ij} = det f^{ij} :
    // ----------------------------------------
    Sym_tensor tgam_in = ff_in.con() + hh_in ; 
    
    Scalar det_in = tgam_in(1, 1)*tgam_in(2, 2)*tgam_in(3, 3) 
        + tgam_in(1, 2)*tgam_in(2, 3)*tgam_in(3, 1)
        + tgam_in(1, 3)*tgam_in(2, 1)*tgam_in(3, 2) 
        - tgam_in(3, 1)*tgam_in(2, 2)*tgam_in(1, 3)
        - tgam_in(3, 2)*tgam_in(2, 3)*tgam_in(1, 1) 
        - tgam_in(3, 3)*tgam_in(2, 1)*tgam_in(1, 2) ;
    
    double diffdet = max(maxabs(det_in - 1. / ff.determinant(), 
        "Deviation of det tgam^{ij} from 1/f")) ;
    if ( diffdet > 1.e-13 ) {
        cerr << 
        "Time_slice_conf::Time_slice_conf : the input h^{ij} does not"
        << " ensure \n" << "  det tgam_{ij} = f  ! \n" 
        << "  error = " << diffdet << endl ; 
        abort() ; 
    }
          
    n_evol.update(lapse_in, jtime, time_init) ;     
    beta_evol.update(shift_in, jtime, time_init) ; 
    trk_evol.update(trk_in, jtime, time_init) ;
    A_hata() ;
    B_hata() ;
    
    set_der_0x0() ;  
    
}
                 

// Constructor from physical metric
// --------------------------------                 

Time_slice_conf::Time_slice_conf(const Scalar& lapse_in, const Vector& shift_in,
               const Sym_tensor& gamma_in, const Sym_tensor& kk_in,
               const Metric_flat& ff_in, int depth_in) 
                    : Time_slice(lapse_in, shift_in, gamma_in, kk_in, depth_in),
                      ff(ff_in),
                      psi_evol(depth_in), 
                      npsi_evol(depth_in),
                      hh_evol(depth_in), 
                      hata_evol(depth_in),
		      A_hata_evol(depth_in), B_hata_evol(depth_in) {
                    
    set_der_0x0() ; // put here in order not to erase p_psi4

    double time_init = the_time[jtime] ; 
    
    assert( p_gamma != 0x0 ) ; 
    p_psi4 = new Scalar( pow( p_gamma->determinant() / ff.determinant(), 
                         0.3333333333333333) ) ;
    p_psi4->std_spectral_base() ;

    Scalar tmp = pow(*p_psi4, 0.25) ;
    tmp.std_spectral_base() ;
    psi_evol.update(tmp , jtime, time_init ) ; 
    
    hh_evol.update( (*p_psi4) * p_gamma->con() - ff.con(), 
                    jtime, time_init ) ; 
    
    hata_evol.update( tmp*tmp*(*p_psi4)*(*p_psi4) *( Time_slice::k_uu() 
                    - 0.3333333333333333 * trk_evol[jtime] * p_gamma->con() ), 
                    jtime, time_init ) ;                   
    A_hata() ;
    B_hata() ;
}

// Constructor as standard time slice of flat spacetime (Minkowski)
// ----------------------------------------------------------------

Time_slice_conf::Time_slice_conf(const Map& mp, const Base_vect& triad, 
                                 const Metric_flat& ff_in, int depth_in) 
    : Time_slice(mp, triad, depth_in),
      ff(ff_in), 
      psi_evol(depth_in), 
      npsi_evol(depth_in),
      hh_evol(depth_in), 
      hata_evol(depth_in),
      A_hata_evol(depth_in), B_hata_evol(depth_in) {
    
    double time_init = the_time[jtime] ; 
    
    // Psi identically one:
    Scalar tmp(mp) ; 
    tmp.set_etat_one() ; 
    tmp.std_spectral_base() ;
    psi_evol.update(tmp, jtime, time_init) ; 
    
    // N Psi identically one:
    npsi_evol.update(tmp, jtime, time_init) ; 
    
    // h^{ij} identically zero:
    Sym_tensor stmp(mp, CON, triad) ; 
    stmp.set_etat_zero() ; 
    hh_evol.update(stmp, jtime, time_init) ; 
    
    // \hat{A}^{ij} identically zero:
    hata_evol.update(stmp, jtime, time_init) ; 

    tmp.set_etat_zero() ;
    A_hata_evol.update(tmp, jtime, time_init) ;
    B_hata_evol.update(tmp, jtime, time_init) ;

    set_der_0x0() ; 

}


// Constructor from binary file             
// ----------------------------

Time_slice_conf::Time_slice_conf(const Map& mp, const Base_vect& triad, 
                        const Metric_flat& ff_in, FILE* fich, 
                        bool partial_read, int depth_in) 
        : Time_slice(mp, triad, fich, true, depth_in),
          ff(ff_in), 
          psi_evol(depth_in), 
          npsi_evol(depth_in),
          hh_evol(depth_in), 
          hata_evol(depth_in),
	  A_hata_evol(depth_in), B_hata_evol(depth_in) {

    // Put here, not to destroy p_vec_X
    set_der_0x0() ; 

    // Reading of various fields
    // -------------------------
    
    int jmin = jtime - depth + 1 ; 
    int indicator ; 

    // Psi
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) {
            Scalar psi_file(mp, *(mp.get_mg()), fich) ; 
            psi_evol.update(psi_file, j, the_time[j]) ; 
        }
    } 
    // hat{A}^{ij}
    for (int j=jmin; j<=jtime; j++) {
	fread_be(&indicator, sizeof(int), 1, fich) ;	
	if (indicator == 1) {
	    Sym_tensor hat_A_file(mp, triad, fich) ; 
	    hata_evol.update( hat_A_file, j, the_time[j] ) ; 
        }
    }

    //A and B...
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) {
            Scalar A_file(mp, *(mp.get_mg()), fich) ; 
            A_hata_evol.update(A_file, j, the_time[j]) ; 
        }
    } 
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) {
            Scalar B_file(mp, *(mp.get_mg()), fich) ; 
            B_hata_evol.update(B_file, j, the_time[j]) ; 
        }
    } 

    // Case of a full reading
    // -----------------------
    if (!partial_read) {
        cout << 
        "Time_slice constructor from file: the case of full reading\n"
        << " is not ready yet !" << endl ; 
        abort() ; 
    }

} 




// Copy constructor
// ----------------

Time_slice_conf::Time_slice_conf(const Time_slice_conf& tin) 
                    : Time_slice(tin), 
                      ff(tin.ff),
                      psi_evol(tin.psi_evol), 
                      npsi_evol(tin.npsi_evol),
                      hh_evol(tin.hh_evol), 
                      hata_evol(tin.hata_evol),
		      A_hata_evol(tin.A_hata_evol),
		      B_hata_evol(tin.B_hata_evol){

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
    if (p_vec_X != 0x0) delete p_vec_X ;
    
    set_der_0x0() ;

    Time_slice::del_deriv() ; 
}


void Time_slice_conf::set_der_0x0() const {

    p_tgamma = 0x0 ; 
    p_psi4 = 0x0 ; 
    p_ln_psi = 0x0 ; 
    p_hdirac = 0x0 ;
    p_vec_X = 0x0 ;
    
}


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Time_slice_conf::operator=(const Time_slice_conf& tin) {

    Time_slice::operator=(tin) ; 

    psi_evol = tin.psi_evol ; 
    npsi_evol = tin.npsi_evol ; 
    hh_evol = tin.hh_evol ; 
    hata_evol = tin.hata_evol ; 
    A_hata_evol = tin.A_hata_evol ;
    B_hata_evol = tin.B_hata_evol ;
       
    del_deriv() ; 
    
}

void Time_slice_conf::operator=(const Time_slice& tin) {

    Time_slice::operator=(tin) ; 

    cerr << 
    "Time_slice_conf::operator=(const Time_slice& ) : not implemented yet !"
        << endl ;
    abort() ;       
    del_deriv() ; 
    
}


void Time_slice_conf::set_psi_del_npsi(const Scalar& psi_in) {

    psi_evol.update(psi_in, jtime, the_time[jtime]) ; 

    // Reset of quantities depending on Psi:
    if (p_psi4 != 0x0) {
        delete p_psi4 ; 
        p_psi4 = 0x0 ; 
    }
    if (p_ln_psi != 0x0) {
        delete p_ln_psi ; 
        p_ln_psi = 0x0 ; 
    }
    if (p_gamma != 0x0) {
        delete p_gamma ; 
        p_gamma = 0x0 ; 
    }
    npsi_evol.downdate(jtime) ; 
    gam_dd_evol.downdate(jtime) ; 
    gam_uu_evol.downdate(jtime) ; 
    adm_mass_evol.downdate(jtime) ;  
     
}

void Time_slice_conf::set_psi_del_n(const Scalar& psi_in) {

    psi_evol.update(psi_in, jtime, the_time[jtime]) ; 

    // Reset of quantities depending on Psi:
    if (p_psi4 != 0x0) {
        delete p_psi4 ; 
        p_psi4 = 0x0 ; 
    }
    if (p_ln_psi != 0x0) {
        delete p_ln_psi ; 
        p_ln_psi = 0x0 ; 
    }
    if (p_gamma != 0x0) {
        delete p_gamma ; 
        p_gamma = 0x0 ; 
    }
    n_evol.downdate(jtime) ; 
    gam_dd_evol.downdate(jtime) ; 
    gam_uu_evol.downdate(jtime) ; 
    adm_mass_evol.downdate(jtime) ;  
     
}


void Time_slice_conf::set_npsi_del_psi(const Scalar& npsi_in) {

    npsi_evol.update(npsi_in, jtime, the_time[jtime]) ; 

    // Reset of quantities depending on N Psi:
    psi_evol.downdate(jtime) ; 
    if (p_psi4 != 0x0) {
        delete p_psi4 ; 
        p_psi4 = 0x0 ; 
    }
    if (p_ln_psi != 0x0) {
        delete p_ln_psi ; 
        p_ln_psi = 0x0 ; 
    }
    if (p_gamma != 0x0) {
        delete p_gamma ; 
        p_gamma = 0x0 ; 
    }
    gam_dd_evol.downdate(jtime) ; 
    gam_uu_evol.downdate(jtime) ; 
    adm_mass_evol.downdate(jtime) ;  
     
}


void Time_slice_conf::set_npsi_del_n(const Scalar& npsi_in) {

    npsi_evol.update(npsi_in, jtime, the_time[jtime]) ; 

    // Reset of quantities depending on N Psi:
    n_evol.downdate(jtime) ; 
    adm_mass_evol.downdate(jtime) ;  
     
}


void Time_slice_conf::set_hh(const Sym_tensor& hh_in) {

    // Check whether det tgam^{ij} = det f^{ij} :
    // ----------------------------------------
    Sym_tensor tgam_in = ff.con() + hh_in ; 
    
    Scalar det_in = tgam_in(1, 1)*tgam_in(2, 2)*tgam_in(3, 3) 
        + tgam_in(1, 2)*tgam_in(2, 3)*tgam_in(3, 1)
        + tgam_in(1, 3)*tgam_in(2, 1)*tgam_in(3, 2) 
        - tgam_in(3, 1)*tgam_in(2, 2)*tgam_in(1, 3)
        - tgam_in(3, 2)*tgam_in(2, 3)*tgam_in(1, 1) 
        - tgam_in(3, 3)*tgam_in(2, 1)*tgam_in(1, 2) ;
    
    double diffdet = max(maxabs(det_in - 1. / ff.determinant(), 
        "Deviation of det tgam^{ij} from 1/f")) ;
    if ( diffdet > 1.e-13 ) {
        cerr << 
        "Time_slice_conf::set_hh : the input h^{ij} does not"
        << " ensure \n" << "  det tgam_{ij} = f  ! \n" 
        << "  error = " << diffdet << endl ; 
        abort() ; 
    }
          
    hh_evol.update(hh_in, jtime, the_time[jtime]) ; 
    
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
    gam_dd_evol.downdate(jtime) ; 
    gam_uu_evol.downdate(jtime) ; 
    adm_mass_evol.downdate(jtime) ;  
     
}


void Time_slice_conf::set_hata(const Sym_tensor& hata_in) {

    hata_evol.update(hata_in, jtime, the_time[jtime]) ; 

    // Reset of quantities depending on A^{ij}:
    A_hata_evol.downdate(jtime) ;
    B_hata_evol.downdate(jtime) ;
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 

}

void Time_slice_conf::set_hata_TT(const Sym_tensor_tt& hata_tt) {
    
    A_hata_evol.update(hata_tt.compute_A(true), jtime, the_time[jtime]) ;
    B_hata_evol.update(hata_tt.compute_tilde_B_tt(true), jtime, the_time[jtime]) ;

    hata_evol.downdate(jtime) ;
    // Reset of quantities depending on A^{ij}:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 
}

void Time_slice_conf::set_hata_from_XAB(Param* par_bc, Param* par_mat) {

    assert (p_vec_X != 0x0) ;
    assert (A_hata_evol.is_known(jtime)) ;
    assert (B_hata_evol.is_known(jtime)) ;

    const Map& mp = p_vec_X->get_mp() ;

    Sym_tensor_tt hata_tt(mp, mp.get_bvect_spher(), ff) ;
    hata_tt.set_A_tildeB(A_hata_evol[jtime], B_hata_evol[jtime], par_bc, par_mat) ;

    hata_evol.update( hata_tt + p_vec_X->ope_killing_conf(ff), jtime, the_time[jtime]) ;

    // Reset of quantities depending on A^{ij}:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ; 
    
}


                //-----------------------------------------------//
                //  Update of fields from base class Time_slice  //
                //-----------------------------------------------//

const Scalar& Time_slice_conf::nn() const {

    if (!( n_evol.is_known(jtime) ) ) {

        assert( psi_evol.is_known(jtime) ) ; 
        assert( npsi_evol.is_known(jtime) ) ; 
        
        n_evol.update( npsi_evol[jtime] / psi_evol[jtime] , 
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
       
        k_uu_evol.update( hata()/(psi4()*psi4()*psi()*psi()) 
			  + 0.3333333333333333* trk()* gam().con(),
			  jtime, the_time[jtime] ) ; 
    }

    return k_uu_evol[jtime] ;

}





                //-----------------------------------//
                //  Update of fields from this class //
                //-----------------------------------//

const Scalar& Time_slice_conf::A_hata() const {

    if ( !( A_hata_evol.is_known(jtime) ) ) {

	assert( hata_evol.is_known(jtime) ) ;

	A_hata_evol.update( hata_evol[jtime].compute_A(true), jtime, the_time[jtime] ) ;
    }
    return A_hata_evol[jtime] ;
}

const Scalar& Time_slice_conf::B_hata() const {

    if ( !( B_hata_evol.is_known(jtime) ) ) {

	assert( hata_evol.is_known(jtime) ) ;

	B_hata_evol.update( hata_evol[jtime].compute_tilde_B_tt(true), jtime, 
			    the_time[jtime] ) ;
    }
    return A_hata_evol[jtime] ;
}


const Scalar& Time_slice_conf::psi() const {

    if (!( psi_evol.is_known(jtime) ) ) {

        assert( n_evol.is_known(jtime) ) ; 
        assert( npsi_evol.is_known(jtime) ) ; 
        
        psi_evol.update( npsi_evol[jtime] / n_evol[jtime] , jtime, the_time[jtime] ) ; 
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


const Scalar& Time_slice_conf::npsi() const {

    if (!( npsi_evol.is_known(jtime) ) ) {
        
        assert( n_evol.is_known(jtime) ) ; 
        assert( psi_evol.is_known(jtime) ) ; 

        npsi_evol.update( psi_evol[jtime] * n_evol[jtime], jtime, the_time[jtime] ) ; 
    }

    return npsi_evol[jtime] ;

}


const Metric& Time_slice_conf::tgam() const {

    if (p_tgamma == 0x0) {
        p_tgamma = new Metric( ff.con() + hh() ) ; 
    }
    
    return *p_tgamma ; 

}


const Sym_tensor& Time_slice_conf::hh(Param*, Param*) const {

    assert( hh_evol.is_known(jtime) ) ; 
    return hh_evol[jtime] ; 

}

Sym_tensor Time_slice_conf::aa() const {

    return hata()/(psi4()*psi()*psi()) ; 

}


const Sym_tensor& Time_slice_conf::hata() const {

    if( !(hata_evol.is_known(jtime)) ) {

        assert( hh_evol.is_known(jtime) ) ; 

        Sym_tensor hh_point = hh_evol.time_derive(jtime, scheme_order) ; 
        hh_point.inc_dzpuis(2) ; // dzpuis : 0 -> 2
    
        Sym_tensor resu = hh_point - hh().derive_lie(beta()) 
            - 0.6666666666666666 * beta().divergence(ff) * hh()
            + beta().ope_killing_conf(ff) ; 

        resu = psi4()*psi()*psi() * resu / (2*nn()) ;
        
        hata_evol.update(resu, jtime, the_time[jtime]) ;
        
    }
     
    return hata_evol[jtime] ; 

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

const Vector& Time_slice_conf::vec_X(int methode_poisson) const {

    if (p_vec_X == 0x0) {
	assert ( hata_evol.is_known(jtime) ) ;
	Vector source = hata_evol[jtime].divergence(ff) ;
        p_vec_X = new Vector( source.poisson( 1./3., methode_poisson ) ) ;
    }
    
    return *p_vec_X ;     
}

void Time_slice_conf::compute_X_from_momentum_constraint
(const Vector& hat_S, Param* par_bc, Param* par_mat, int iter_max, double precis,
 double relax, int meth_poisson) {

    // Some checks
    assert( hat_S.get_index_type(0) == CON ) ;
    assert( A_hata_evol.is_known(jtime) ) ;
    assert( B_hata_evol.is_known(jtime) ) ;
    assert( hh_evol.is_known(jtime) ) ;

    // Initializations
    //----------------
    const Map& mp = hat_S.get_mp() ;
    const Tensor_sym& delta = tgam().connect().get_delta() ;
    Sym_tensor_tt hata_tt(mp, mp.get_bvect_cart(), ff) ;
    hata_tt.set_A_tildeB( A_hata_evol[jtime], B_hata_evol[jtime], par_bc, par_mat ) ;
    if (p_vec_X != 0x0) {
	delete p_vec_X ;
	p_vec_X = 0x0 ;
    }

    // Constant part of the source
    Vector source = hat_S - contract( delta, 1, 2, hata_tt, 0, 1 ) 
	- 2./3.*psi4()*psi()*psi()*contract( tgam().con(), 0, trk().derive_cov(ff), 0 );

    p_vec_X = new Vector( source.poisson( 1./3., meth_poisson ) ) ;

    // Iteration on the vector X
    //--------------------------
    for (int it=0; it<iter_max; it++) {

	Vector fuente = source 
	    - contract( delta, 1, 2, p_vec_X->ope_killing_conf(ff), 0, 1 ) ;

	Vector X_new = fuente.poisson( 1./3., meth_poisson) ;//### check the dzpuis!

	// Control of the convergence
	double diff = max( max( abs( X_new - *p_vec_X ) ) ) ;

	// Relaxation
	(*p_vec_X) = relax*X_new + ( 1. - relax )*(*p_vec_X) ;

	// If converged, gets out of the loop
	if (diff < precis) break ;
    }

    // Update of \hat{A}^{ij} and reset of related quantities
    hata_evol.update( hata_tt + p_vec_X->ope_killing_conf(ff), jtime, the_time[jtime] ) ;

    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ;     
}

void Time_slice_conf::set_AB_hata(const Scalar& A_in, const Scalar& B_in) {

    A_hata_evol.update(A_in, jtime, the_time[jtime]) ;
    B_hata_evol.update(B_in, jtime, the_time[jtime]) ;

    hata_evol.downdate(jtime) ;
    // Reset of quantities depending on A^{ij}:
    k_dd_evol.downdate(jtime) ; 
    k_uu_evol.downdate(jtime) ;     

}


void Time_slice_conf::check_psi_dot(Tbl& tlnpsi_dot, Tbl& tdiff, Tbl& tdiff_rel) const {

    // Computation of d/dt ln Psi :
     
    Scalar lnpsi_dot(psi().get_mp()) ; 
    if ( psi_evol.is_known(jtime-1) ) {       
        lnpsi_dot = psi_evol.time_derive(jtime, scheme_order) / psi() ; 
    }
    else {
        lnpsi_dot = npsi_evol.time_derive(jtime, scheme_order) / npsi()
                           - n_evol.time_derive(jtime, scheme_order) / nn()  ; 
    }
    
    tlnpsi_dot = max(abs(lnpsi_dot)) ; 
        
    // Error on the d/dt ln Psi relation : 
    
    Scalar diff = contract(beta(),0, ln_psi().derive_cov(ff),0) 
        + 0.1666666666666666 * ( beta().divergence(ff) - nn() * trk() ) ;

    diff.dec_dzpuis(2) ; 
    
    Tbl tref = max(abs(diff)) + tlnpsi_dot ; 
    
    diff -= lnpsi_dot ; 
    
    tdiff = max(abs(diff)) ;
            
    tdiff_rel = tdiff / tref  ; 
    
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
    if (npsi_evol.is_known(jtime)) {
        maxabs( npsi_evol[jtime], "N Psi", flux) ;
    }
    if (hh_evol.is_known(jtime)) {
        maxabs( hh_evol[jtime], "h^{ij}", flux) ;
    }
    if (hata_evol.is_known(jtime)) {
        maxabs( hata_evol[jtime], "hat{A}^{ij}", flux) ;
    }

    if (p_tgamma != 0x0) flux << 
        "Conformal metric tilde gamma is up to date" << endl ; 
    if (p_psi4 != 0x0) maxabs( *p_psi4, "Psi^4", flux) ; 
    if (p_ln_psi != 0x0) maxabs( *p_ln_psi, "ln(Psi)", flux) ; 
    if (p_hdirac != 0x0) maxabs( *p_hdirac, "H^i", flux) ; 
    if (p_vec_X != 0x0) maxabs( *p_vec_X, "X^i", flux) ; 
    
    return flux ; 

}



void Time_slice_conf::sauve(FILE* fich, bool partial_save) const {
    
    // Writing of quantities common to all derived classes of Time_slice
    // -----------------------------------------------------------------
    
    Time_slice::sauve(fich, true) ; 
    
    // Writing of quantities common to all derived classes of Time_slice_conf
    // ----------------------------------------------------------------------
    
    int jmin = jtime - depth + 1 ; 

    // Psi
    psi() ;     // forces the update at the current time step
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (psi_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) psi_evol[j].sauve(fich) ;
    }

    // hat{A}
    hata() ;
    for (int j=jmin; j<=jtime; j++) {
	int indicator = ( hata_evol.is_known(j) ? 1 : 0 ) ;
	fwrite_be(&indicator, sizeof(int), 1, fich) ;
	if (indicator == 1) hata_evol[j].sauve(fich) ;
    }

    //A and B...
    A_hata() ;
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (A_hata_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) A_hata_evol[j].sauve(fich) ;
    }

    B_hata() ;
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (B_hata_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) B_hata_evol[j].sauve(fich) ;
    }	

    // Case of a complete save
    // -----------------------
    if (!partial_save) {
        cout << "Time_slice_conf::sauve: the full writing is not ready yet !" 
             << endl ; 
        abort() ; 
    }
    
}
