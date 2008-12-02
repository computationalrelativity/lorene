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
 * Revision 1.23  2008/12/02 15:02:22  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
 * Revision 1.22  2007/11/06 14:47:07  j_novak
 * New constructor from a rotating star in Dirac gauge (class Star_rot_Dirac).
 * Evolution can take into account matter terms.
 *
 * Revision 1.21  2007/09/25 16:54:11  j_novak
 * *** empty log message ***
 *
 * Revision 1.20  2007/09/25 16:52:15  j_novak
 * *** empty log message ***
 *
 * Revision 1.19  2007/06/05 07:38:37  j_novak
 * Better treatment of dzpuis for A and tilde(B) potentials. Some errors in the bases manipulation have been also corrected.
 *
 * Revision 1.18  2007/04/25 15:21:01  j_novak
 * Corrected an error in the initialization of tildeB in
 * Tslice_dirac_max::initial_dat_cts. + New method for solve_hij_AB.
 *
 * Revision 1.17  2007/03/21 14:51:50  j_novak
 * Introduction of potentials A and tilde(B) of h^{ij} into Tslice_dirac_max.
 *
 * Revision 1.16  2004/12/28 14:21:48  j_novak
 * Added the method Sym_tensor_trans::trace_from_det_one
 *
 * Revision 1.15  2004/07/08 12:29:01  j_novak
 * use of new method Tensor::annule_extern_cn
 *
 * Revision 1.14  2004/06/30 08:02:40  j_novak
 * Added filtering in l of khi_new and mu_new. ki_source is forced to go to
 * zero at least as r^2.
 *
 * Revision 1.13  2004/06/17 06:59:41  e_gourgoulhon
 * -- Method initial_data_cts: re-organized treatment of vanishing uu.
 * -- Method hh_det_one: replaced the attenuation with tempo by a call
 *    to the new method Tensor::annule_extern_c2.
 *
 * Revision 1.12  2004/06/08 14:05:06  j_novak
 * Added the attenuation of khi and mu in the last domain in ::det_one(). They are set to zero in the CED.
 *
 * Revision 1.11  2004/05/31 20:31:31  e_gourgoulhon
 * -- Method hh_det_one takes now a time step as argument, to compute
 *    h^{ij} from khi and mu at some arbitrary time step and not only at
 *    the latest one.
 * -- h^{ij} is no longer saved in binary files (method sauve);
 *    accordingly, the constructor from file calls the new version of
 *    hh_det_one to restore h^{ij}.
 *
 * Revision 1.10  2004/05/31 09:08:18  e_gourgoulhon
 * Method sauve and constructor from binary file are now operational.
 *
 * Revision 1.9  2004/05/27 15:25:04  e_gourgoulhon
 * Added constructors from binary file, as well as corresponding
 * functions sauve and save.
 *
 * Revision 1.8  2004/05/17 19:54:10  e_gourgoulhon
 * Method initial_data_cts: added arguments graph_device and method_poisson_vect.
 *
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
            const Sym_tensor_trans& hh_in, const Sym_tensor& hata_in, 
            int depth_in) 
  : Time_slice_conf( lapse_in, shift_in, ff_in, psi_in, hh_in, hata_in, 
		     0*lapse_in, depth_in), 
    A_hh_evol(depth_in), B_hh_evol(depth_in), source_A_hh_evol(depth_in), 
    source_B_hh_evol(depth_in), source_A_hata_evol(depth_in) ,
    source_B_hata_evol(depth_in), trh_evol(hh_in.the_trace(), depth_in) 
{ }              

// Constructor as standard time slice of flat spacetime (Minkowski) 
// ----------------------------------------------------------------

Tslice_dirac_max::Tslice_dirac_max(const Map& mp, const Base_vect& triad, 
                                   const Metric_flat& ff_in, int depth_in) 
        : Time_slice_conf(mp, triad, ff_in, depth_in),
          A_hh_evol(depth_in), B_hh_evol(depth_in),   
          source_A_hh_evol(depth_in), source_B_hh_evol(depth_in),   
          source_A_hata_evol(depth_in), source_B_hata_evol(depth_in),   
          trh_evol(depth_in) {

    double time_init = the_time[jtime] ; 
    
    // All potentials identically zero:
    Scalar tmp(mp) ; 
    tmp.set_etat_zero() ; 
    A_hh_evol.update(tmp, jtime, time_init) ; 
    B_hh_evol.update(tmp, jtime, time_init) ; 

    source_A_hh_evol.update(tmp, jtime, time_init) ; 
    source_B_hh_evol.update(tmp, jtime, time_init) ; 

    source_A_hata_evol.update(tmp, jtime, time_init) ; 
    source_B_hata_evol.update(tmp, jtime, time_init) ; 
    
    // tr h identically zero:
    trh_evol.update(tmp, jtime, time_init) ;     
}   


// Constructor from binary file             
// ----------------------------

Tslice_dirac_max::Tslice_dirac_max(const Map& mp, const Base_vect& triad, 
                        const Metric_flat& ff_in, FILE* fich, 
                        bool partial_read, int depth_in) 
        : Time_slice_conf(mp, triad, ff_in, fich, true, depth_in),
          A_hh_evol(depth_in), B_hh_evol(depth_in),   
          source_A_hh_evol(depth_in), source_B_hh_evol(depth_in),   
          source_A_hata_evol(depth_in), source_B_hata_evol(depth_in),   
          trh_evol(depth_in) {

    if (partial_read) {
        cout << 
        "Constructor of Tslice_dirac_max from file: the case of partial reading\n"
        << "  is not ready yet !"
            << endl ; 
        abort() ; 
    }
    
    // Reading of various fields
    // -------------------------
    
    int jmin = jtime - depth + 1 ; 
    int indicator ; 

    // h^{ij}
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;	
        if (indicator == 1) {
            Sym_tensor hh_file(mp, triad, fich) ; 
            hh_evol.update(hh_file, j, the_time[j]) ; 
        }
    }

    // A - hh
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;	
        if (indicator == 1) {
            Scalar A_hh_file(mp, *(mp.get_mg()), fich) ; 
            A_hh_evol.update(A_hh_file, j, the_time[j]) ; 
        }
    }

    // B - hh
    for (int j=jmin; j<=jtime; j++) {
        fread_be(&indicator, sizeof(int), 1, fich) ;	
        if (indicator == 1) {
            Scalar B_hh_file(mp, *(mp.get_mg()), fich) ; 
            B_hh_evol.update(B_hh_file, j, the_time[j]) ; 
        }
    }    
}

// Constructor from a rotating star             
// --------------------------------

Tslice_dirac_max::Tslice_dirac_max(const Star_rot_Dirac& star, double pdt, int depth_in) 
    : Time_slice_conf(star.get_nn(), star.get_beta(), star.get_mp().flat_met_spher(),
		      0.*star.get_nn(), star.get_hh(), 0.*star.get_aa(),
		      0.*star.get_nn(), depth_in),
          A_hh_evol(depth_in), B_hh_evol(depth_in),   
          source_A_hh_evol(depth_in), source_B_hh_evol(depth_in),   
          source_A_hata_evol(depth_in), source_B_hata_evol(depth_in),   
          trh_evol(depth_in) {

    Scalar tmp = exp(star.get_ln_psi()) ;
    tmp.std_spectral_base() ;
    psi_evol.downdate(jtime) ;
    psi_evol.update(tmp, jtime, the_time[jtime]) ;

    Sym_tensor tmp2 = psi4()*psi()*psi()*star.get_aa() ;
    hata_evol.downdate(jtime) ;
    hata_evol.update(tmp2, jtime, the_time[jtime]) ;

    A_hh() ;
    B_hh() ;
    A_hata() ;
    B_hata() ;
    compute_sources() ;

    // Update of various fields
    // -------------------------
    double ttime1 = the_time[jtime] ; 
    int jtime1 = jtime ; 
    for (int j=1; j < depth; j++) {
        jtime1++ ; 
        ttime1 += pdt ; 
        psi_evol.update(psi_evol[jtime], jtime1, ttime1) ;  
        n_evol.update(n_evol[jtime], jtime1, ttime1) ;  
        beta_evol.update(beta_evol[jtime], jtime1, ttime1) ;  
        hh_evol.update(hh_evol[jtime], jtime1, ttime1) ;
        trk_evol.update(trk_evol[jtime], jtime1, ttime1) ;
	A_hh_evol.update(A_hh_evol[jtime], jtime1, ttime1) ;
	B_hh_evol.update(B_hh_evol[jtime], jtime1, ttime1) ;
	A_hata_evol.update(A_hata_evol[jtime], jtime1, ttime1) ;
	B_hata_evol.update(B_hata_evol[jtime], jtime1, ttime1) ;
	trh_evol.update(trh_evol[jtime], jtime1, ttime1) ;
	k_dd_evol.update(k_dd_evol[jtime], jtime1, ttime1) ;
        the_time.update(ttime1, jtime1, ttime1) ;         
    } 
    jtime += depth - 1 ;

    initialize_sources_copy() ;
}

// Copy constructor
// ----------------

Tslice_dirac_max::Tslice_dirac_max(const Tslice_dirac_max& tin) 
                    : Time_slice_conf(tin), 
                      A_hh_evol(tin.A_hh_evol), 
                      B_hh_evol(tin.B_hh_evol),
                      source_A_hh_evol(tin.source_A_hh_evol), 
                      source_B_hh_evol(tin.source_B_hh_evol),
                      source_A_hata_evol(tin.source_A_hata_evol), 
                      source_B_hata_evol(tin.source_B_hata_evol),
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

    A_hh_evol = tin.A_hh_evol ; 
    B_hh_evol = tin.B_hh_evol ;
    source_A_hh_evol = tin.source_A_hh_evol ; 
    source_B_hh_evol = tin.source_B_hh_evol ;
    source_A_hata_evol = tin.source_A_hata_evol ; 
    source_B_hata_evol = tin.source_B_hata_evol ;
    trh_evol = tin.trh_evol ; 
       
}


void Tslice_dirac_max::set_hh(const Sym_tensor& hh_in) {

    Time_slice_conf::set_hh(hh_in) ; 

    // Reset of quantities depending on h^{ij}:
    A_hh_evol.downdate(jtime) ; 
    B_hh_evol.downdate(jtime) ; 
    source_A_hh_evol.downdate(jtime) ; 
    source_B_hh_evol.downdate(jtime) ; 
    source_A_hata_evol.downdate(jtime) ; 
    source_B_hata_evol.downdate(jtime) ; 
    trh_evol.downdate(jtime) ; 
         
}


void Tslice_dirac_max::initial_data_cts(const Sym_tensor& uu, 
                const Scalar& trk_in, const Scalar& trk_point, 
                double pdt, double precis, int method_poisson_vect,
                const char* graph_device, const Scalar* p_ener_dens, 
                const Vector* p_mom_dens, const Scalar* p_trace_stress) {

    
    Time_slice_conf::initial_data_cts(uu, trk_in, trk_point, pdt, precis,
                                      method_poisson_vect, graph_device, 
                                      p_ener_dens, p_mom_dens, p_trace_stress) ;

    int nz = trk_in.get_mp().get_mg()->get_nzone() ;
    // Setting khi and mu for j <= jtime, taking into account u^{ij} = dh^{ij}/dt
    //--------------------------------------------------------------------------
    for (int j = jtime-depth+1 ; j <= jtime; j++) {
            
	// A and tildeB are computed from the value of hh
	Scalar tmp = hh_evol[j].compute_A(true) ;
	assert (tmp.get_etat() != ETATNONDEF) ;
	if (tmp.get_etat() != ETATZERO) {
	    assert(tmp.get_mp().get_mg()->get_type_r(nz-1) == UNSURR) ;
	    tmp.annule_domain(nz-1) ;
	}
	tmp.set_dzpuis(0) ;
	A_hh_evol.update(tmp, j, the_time[j]) ;

	tmp = hh_evol[jtime].compute_tilde_B_tt(true) ;
	assert (tmp.get_etat() != ETATNONDEF) ;
	if (tmp.get_etat() != ETATZERO) {
	    assert(tmp.get_mp().get_mg()->get_type_r(nz-1) == UNSURR) ;
	    tmp.annule_domain(nz-1) ;
	}
	tmp.set_dzpuis(0) ;
	B_hh_evol.update(tmp, j, the_time[j]) ;
    }
  
    cout << endl << 
    "Tslice_dirac_max::initial_data_cts : variation of A and tilde(B) for J = " 
    << jtime << " :\n" ;  
    maxabs(A_hh_evol[jtime] - A_hh_evol[jtime-1], "A(h)^J - A(h)^{J-1}") ; 
    
    maxabs(B_hh_evol[jtime] - B_hh_evol[jtime-1], "B(h)^J - B(h)^{J-1}") ; 
    
    maxabs(A_hata_evol[jtime] - A_hata_evol[jtime-1], "A(hat{A})^J - A(hat{A})^{J-1}") ; 
    
    maxabs(B_hata_evol[jtime] - B_hata_evol[jtime-1], "B(hat{A})^J - B(hat{A})^{J-1}") ; 
    
    // Reset of derived quantities (at the new time step jtime)
    // ---------------------------
    del_deriv() ; 
    
    compute_sources() ;
    initialize_sources_copy() ; //## should be a call to the Runge-Kutta integrator
}


void Tslice_dirac_max::set_khi_mu(const Scalar& khi_in, const Scalar& mu_in) {

    const Map& mp = khi_in.get_mp() ;
    Sym_tensor_tt hh_tt(mp, mp.get_bvect_spher(), mp.flat_met_spher());
    hh_tt.set_khi_mu(khi_in, mu_in) ;
    Sym_tensor_trans hh_tmp(mp, mp.get_bvect_spher(), mp.flat_met_spher());
    hh_tmp.trace_from_det_one(hh_tt) ;

    set_hh( hh_tmp ) ;
    
} 

void Tslice_dirac_max::set_trh(const Scalar& trh_in) {

    trh_evol.update(trh_in, jtime, the_time[jtime]) ; 
    cout << "Tslice_dirac_max::set_trh : #### WARNING : \n"
        << "   this method does not check whether det(tilde gamma) = 1"
        << endl ; 
        
    // Reset of quantities depending on the trace:
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
    source_A_hh_evol.downdate(jtime) ; 
    source_B_hh_evol.downdate(jtime) ; 
    source_A_hata_evol.downdate(jtime) ; 
    source_B_hata_evol.downdate(jtime) ; 
    gam_dd_evol.downdate(jtime) ; 
    gam_uu_evol.downdate(jtime) ; 
    adm_mass_evol.downdate(jtime) ;  
     
} 


                //----------------------------------------------------//
                //  Update of fields from base class Time_slice_conf  //
                //----------------------------------------------------//


const Sym_tensor& Tslice_dirac_max::hh(Param* par_bc, Param* par_mat) const {

    if (!( hh_evol.is_known(jtime) ) ) {

        assert (A_hh_evol.is_known(jtime)) ;
        assert (B_hh_evol.is_known(jtime)) ;
    
        // Computation of h^{ij} to ensure det tgam_{ij} = det f_{ij} :

        hh_det_one(jtime, par_bc, par_mat) ;
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

const Scalar& Tslice_dirac_max::A_hh() const {

    if (!( A_hh_evol.is_known(jtime) ) ) {
	assert( hh_evol.is_known(jtime) ) ;

	A_hh_evol.update( hh_evol[jtime].compute_A(true), jtime, the_time[jtime] ) ;

    }
    return A_hh_evol[jtime] ;
} 

const Scalar& Tslice_dirac_max::B_hh() const {

    if (!( B_hh_evol.is_known(jtime) ) ) {
	assert( hh_evol.is_known(jtime) ) ;

	B_hh_evol.update( hh_evol[jtime].compute_tilde_B_tt(true), jtime, the_time[jtime] ) ;

    }
    return B_hh_evol[jtime] ;
} 

const Scalar& Tslice_dirac_max::trh() const {

    if( !(trh_evol.is_known(jtime)) ) {
    
        // Computation of tr(h) to ensure det tgam_{ij} = det f_{ij} :
        hh_det_one(jtime) ;
        
    }
    
    return trh_evol[jtime] ; 

}



                //------------------//
                //      output      //
                //------------------//

ostream& Tslice_dirac_max::operator>>(ostream& flux) const {

    Time_slice_conf::operator>>(flux) ; 

    flux << "Dirac gauge and maximal slicing" << '\n' ;

    if (A_hh_evol.is_known(jtime)) {
        maxabs( A_hh_evol[jtime], "A_hh", flux) ;
    }
    if (B_hh_evol.is_known(jtime)) {
        maxabs( B_hh_evol[jtime], "B_hh", flux) ;
    }
    if (trh_evol.is_known(jtime)) {
        maxabs( trh_evol[jtime], "tr h", flux) ;
    }
    
    return flux ; 

}


void Tslice_dirac_max::sauve(FILE* fich, bool partial_save) const {
    
    if (partial_save) {
        cout << 
        "Tslice_dirac_max::sauve : the partial_save case is not ready yet !"
            << endl ; 
        abort() ; 
    }
    
    // Writing of quantities common to all derived classes of Time_slice_conf
    // ----------------------------------------------------------------------
    
    Time_slice_conf::sauve(fich, true) ; 
    
    // Writing of the other fields
    // ---------------------------
        
    int jmin = jtime - depth + 1 ; 

    // h^{ij}
    assert( hh_evol.is_known(jtime) ) ;
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (hh_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) hh_evol[j].sauve(fich) ; 
    }	

    // A_hh
    A_hh() ;     // forces the update at the current time step
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (A_hh_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) A_hh_evol[j].sauve(fich) ; 
    }

    // B_hh
    B_hh() ;     // forces the update at the current time step
    for (int j=jmin; j<=jtime; j++) {
        int indicator = (B_hh_evol.is_known(j)) ? 1 : 0 ; 
        fwrite_be(&indicator, sizeof(int), 1, fich) ;
        if (indicator == 1) B_hh_evol[j].sauve(fich) ; 
    }
}


                
                
                
                

