/*
 *  Methods of class Tslice_dirac_max
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

char tslice_dirac_max_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.6  2004/05/10 09:16:32  e_gourgoulhon
 * -- Method initial_data_cts: added a call to del_deriv() at the end.
 * -- Methods set_trh and hh_det_one: added "adm_mass_evol.downdate(jtime)".
 * -- Method trh() : the update is now performed via a call to hh_det_one().
 *
 * Revision 1.5  2004/05/06 15:23:55  e_gourgoulhon
 * Added method initial_data_cts.
 *
 * Revision 1.4  2004/05/03 08:15:48  e_gourgoulhon
 * Method hh_det_one(): added check at the end (deviation from det = 1).
 *
 * Revision 1.3  2004/04/08 16:44:19  e_gourgoulhon
 * Added methods set_* and hh_det_one().
 *
 * Revision 1.2  2004/04/05 21:22:49  e_gourgoulhon
 * Added constructor as standard time slice of Minkowski spacetime.
 *
 * Revision 1.1  2004/03/30 14:00:31  j_novak
 * New class Tslide_dirac_max (first version).
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
#include "utilitaires.h"



			    //--------------//
			    // Constructors //
			    //--------------//


// Constructor from conformal decomposition
// ----------------------------------------

Tslice_dirac_max::Tslice_dirac_max(const Scalar& lapse_in, const Vector& shift_in,
            const Metric_flat& ff_in, const Scalar& psi_in, 
            const Sym_tensor_trans& hh_in, const Sym_tensor aa_in, 
            int depth_in) 
  : Time_slice_conf( lapse_in, shift_in, ff_in, psi_in, hh_in, aa_in, 
		     0*lapse_in, depth_in), 
    khi_evol(hh_in.tt_part().khi(), depth_in), 
    mu_evol(hh_in.tt_part().mu(), depth_in),
    trh_evol(hh_in.the_trace(), depth_in) { }
                 

// Constructor as standard time slice of flat spacetime (Minkowski) 
// ----------------------------------------------------------------

Tslice_dirac_max::Tslice_dirac_max(const Map& mp, const Base_vect& triad, 
                                   const Metric_flat& ff_in, int depth_in) 
        : Time_slice_conf(mp, triad, ff_in, depth_in),
          khi_evol(depth_in),   
          mu_evol(depth_in),   
          trh_evol(depth_in) {

    double time_init = the_time[jtime] ; 
    
    // khi identically zero:
    Scalar tmp(mp) ; 
    tmp.set_etat_zero() ; 
    khi_evol.update(tmp, jtime, time_init) ; 
    
    // mu identically zero:
    mu_evol.update(tmp, jtime, time_init) ; 
    
    // tr h identically zero:
    trh_evol.update(tmp, jtime, time_init) ; 
    
}   



// Copy constructor
// ----------------

Tslice_dirac_max::Tslice_dirac_max(const Tslice_dirac_max& tin) 
                    : Time_slice_conf(tin), 
                      khi_evol(tin.khi_evol), 
                      mu_evol(tin.mu_evol),
                      trh_evol(tin.trh_evol) { }
                      
                      
			    //--------------//
			    //  Destructor  //
			    //--------------//

Tslice_dirac_max::~Tslice_dirac_max(){ }


                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Tslice_dirac_max::operator=(const Tslice_dirac_max& tin) {

    Time_slice_conf::operator=(tin) ; 

    khi_evol = tin.khi_evol ; 
    mu_evol = tin.mu_evol ; 
    trh_evol = tin.trh_evol ; 
       
}


void Tslice_dirac_max::set_hh(const Sym_tensor& hh_in) {

    Time_slice_conf::set_hh(hh_in) ; 

    // Reset of quantities depending on h^{ij}:
    khi_evol.downdate(jtime) ; 
    mu_evol.downdate(jtime) ; 
    trh_evol.downdate(jtime) ; 
         
}


void Tslice_dirac_max::initial_data_cts(const Sym_tensor& uu, 
                const Scalar& trk_in, const Scalar& trk_point, 
                double pdt, double precis,
                const Scalar* p_ener_dens, const Vector* p_mom_dens, 
                const Scalar* p_trace_stress) {

    
    Time_slice_conf::initial_data_cts(uu, trk_in, trk_point, pdt, precis,
                p_ener_dens, p_mom_dens, p_trace_stress) ;

    for (int j = jtime-depth+1 ; j < jtime; j++) {
            
        Sym_tensor hhtmp = hh_evol[j] ;
        int nz = uu.get_mp().get_mg()->get_nzone() ; 
        hhtmp.annule_domain(nz-1) ; 
        for (int i=1; i<=3; i++) {
            for (int k=i; k<=3; k++) {
                hhtmp.set(i,k).set_dzpuis(4) ; 
            }
        }
        
        Scalar tmp = hhtmp.transverse(ff).tt_part().khi() ;
        tmp.annule_domain(nz-1) ;
//##        khi_evol.update(tmp, j, the_time[j]) ;

// Formula valid for vanishing time derivative only : 
        khi_evol.update(hh_evol[jtime].transverse(ff).tt_part().khi(),
                           j, the_time[j]) ;

        tmp = hhtmp.transverse(ff).tt_part().mu() ;
        tmp.annule_domain(nz-1) ;
//##        mu_evol.update(tmp, j, the_time[j]) ;

// Formula valid for vanishing time derivative only : 

        mu_evol.update(hh_evol[jtime].transverse(ff).tt_part().mu(),
                           j, the_time[j]) ;
    }
  
    khi_evol.update(hh_evol[jtime].transverse(ff).tt_part().khi(),
                           jtime, the_time[jtime]) ;
    mu_evol.update(hh_evol[jtime].transverse(ff).tt_part().mu(),
                           jtime, the_time[jtime]) ;
    
    maxabs(khi_evol[jtime] - khi_evol[jtime-1], "khi^J - khi^{J-1}") ; 
    
    maxabs(mu_evol[jtime] - mu_evol[jtime-1], "mu^J - mu^{J-1}") ; 
    
    // Reset of derived quantities (at the new time step jtime)
    // ---------------------------
    del_deriv() ; 
    
}


void Tslice_dirac_max::set_khi_mu(const Scalar& khi_in, const Scalar& mu_in) {

    khi_evol.update(khi_in, jtime, the_time[jtime]) ; 
    mu_evol.update(mu_in, jtime, the_time[jtime]) ; 
    
    // Computation of trace h and h^{ij} to ensure det tgam_{ij} = det f_{ij} :

    hh_det_one() ;
    
    
} 


void Tslice_dirac_max::set_trh(const Scalar& trh_in) {

    trh_evol.update(trh_in, jtime, the_time[jtime]) ; 
    cout << "Tslice_dirac_max::set_trh : #### WARNING : \n"
        << "   this method does not check whether det(tilde gamma) = 1"
        << endl ; 
        
    // Reset of quantities depending on mu:
    hh_evol.downdate(jtime) ; 
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


void Tslice_dirac_max::hh_det_one() const {

    assert (khi_evol.is_known(jtime)) ;   // The starting point
    assert (mu_evol.is_known(jtime)) ;    // of the computation 

    int it_max = 100 ;
    double precis = 1.e-14 ;
  
    const Map& mp = khi().get_mp() ;
    
    // The TT part of h^{ij}, which stays unchanged during the computation :
    Sym_tensor_tt hijtt(mp, *(ff.get_triad()), ff) ;
    hijtt.set_khi_mu(khi(), mu(), 2) ;
    
    // The representation of h^{ij} as an object of class Sym_tensor_trans :
    Sym_tensor_trans hij = hijtt ; 

    // The trace h = f_{ij} h^{ij} :
    Scalar htrace(mp) ;
        
    // Value of h at previous step of the iterative procedure below :
    Scalar htrace_prev(mp) ;
    htrace_prev.set_etat_zero() ;   // initialisation to zero
    
    int it ;
    for (it=0; it<it_max; it++) {
      
        // Trace h from the condition det(f^{ij} + h^{ij}) = det f^{ij} :
      
        htrace = hij(1,1) * hij(2,3) * hij(2,3) 
	    + hij(2,2) * hij(1,3) * hij(1,3) + hij(3,3) * hij(1,2) * hij(1,2)
	    - 2.* hij(1,2) * hij(1,3) * hij(2,3) 
            - hij(1,1) * hij(2,2) * hij(3,3) ;
        
        htrace.dec_dzpuis(2) ; // dzpuis: 6 --> 4
        
	htrace += hij(1,2) * hij(1,2) + hij(1,3) * hij(1,3) 
                    + hij(2,3) * hij(2,3) - hij(1,1) * hij(2,2) 
                    - hij(1,1) * hij(3,3) - hij(2,2) * hij(3,3) ;

        // New value of hij from htrace and hijtt (obtained by solving 
        //    the Poisson equation for Phi) : 

        hij.set_tt_trace(hijtt, htrace) ; 

        double diff = max(max(abs(htrace - htrace_prev))) ;
        cout << "Tslide_dirac_max::hh_det_one() : " 
	     << "iteration : " << it << " difference : " << diff << endl ;
        if (diff < precis) break ;
        else htrace_prev = htrace ;

    }
    
    if (it == it_max) {
        cerr << "Tslide_dirac_max::hh_det_one() : convergence not reached \n" ;
        cerr << "  for the required accuracy (" << precis << ") ! " << endl ;
        abort() ;
    }

    // Updates
    // -------
    trh_evol.update(htrace, jtime, the_time[jtime]) ;
    
    // The longitudinal part of h^{ij}, which is zero by virtue of Dirac gauge :
    Vector wzero(mp, CON,  *(ff.get_triad())) ; 
    wzero.set_etat_zero() ;                   

    // Temporary Sym_tensor with longitudinal part set to zero : 
    Sym_tensor hh_new(mp, CON, *(ff.get_triad())) ;
    
    hh_new.set_longit_trans(wzero, hij) ;
    
    hh_evol.update(hh_new, jtime, the_time[jtime]) ;
    
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
         
    maxabs(tgam().determinant() - 1, 
    "Max. of absolute value of deviation from det tgam = 1") ; 

}

                //----------------------------------------------------//
                //  Update of fields from base class Time_slice_conf  //
                //----------------------------------------------------//


const Sym_tensor& Tslice_dirac_max::hh() const {

    if (!( hh_evol.is_known(jtime) ) ) {

        assert (khi_evol.is_known(jtime)) ;
        assert (mu_evol.is_known(jtime)) ;
    
        // Computation of h^{ij} to ensure det tgam_{ij} = det f_{ij} :

        hh_det_one() ;
    }
  
    return hh_evol[jtime] ; 

}


const Scalar& Tslice_dirac_max::trk() const {

    if( !(trk_evol.is_known(jtime)) ) {

      Scalar resu(ff.get_mp()) ;
      resu.set_etat_zero() ;
      
      trk_evol.update(resu, jtime, the_time[jtime]) ;

    } 
    
    return trk_evol[jtime] ; 

}


const Vector& Tslice_dirac_max::hdirac() const {

    if (p_hdirac == 0x0) {
        p_hdirac = new Vector(ff.get_mp(), CON, ff.get_triad() ) ;
	p_hdirac->set_etat_zero() ;
    }
    
    return *p_hdirac ; 

}




                //-----------------------------------//
                //  Update of fields from this class //
                //-----------------------------------//


const Scalar& Tslice_dirac_max::khi() const {

    if (!( khi_evol.is_known(jtime) ) ) {

        assert( hh_evol.is_known(jtime) ) ; 
        
        khi_evol.update( hh().transverse(ff).tt_part().khi(), 
			 jtime, the_time[jtime] ) ; 
    }

    return khi_evol[jtime] ;

} 

const Scalar& Tslice_dirac_max::mu() const {

    if (!( mu_evol.is_known(jtime) ) ) {

        assert( hh_evol.is_known(jtime) ) ; 
      
        mu_evol.update( hh().transverse(ff).tt_part().mu(), 
			jtime, the_time[jtime] ) ; 
    }

    return mu_evol[jtime] ;

}


const Scalar& Tslice_dirac_max::trh() const {

    if( !(trh_evol.is_known(jtime)) ) {
    
        // Computation of tr(h) to ensure det tgam_{ij} = det f_{ij} :
        hh_det_one() ;
        
    }
    
    return trh_evol[jtime] ; 

}



                //------------------//
                //      output      //
                //------------------//

ostream& Tslice_dirac_max::operator>>(ostream& flux) const {

    Time_slice_conf::operator>>(flux) ; 

    flux << "Dirac gauge and maximal slicing" << '\n' ;

    if (khi_evol.is_known(jtime)) {
        maxabs( khi_evol[jtime], "Khi", flux) ;
    }
    if (mu_evol.is_known(jtime)) {
        maxabs( mu_evol[jtime], "Mu", flux) ;
    }
    if (trh_evol.is_known(jtime)) {
        maxabs( trh_evol[jtime], "tr h", flux) ;
    }
    
    return flux ; 

}



                
                
                
                

