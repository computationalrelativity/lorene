 /*
 *  Methods of class Iso_hor to compute sources for Psi, N y beta
 *
 *    (see file hor_isol.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
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

char source_hor_C[] = "$Header$" ;

/*
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <assert.h>

// Lorene headers
#include "time_slice.h"
#include "isol_hor.h"
#include "metric.h"
#include "evolution.h"
#include "unites.h"
#include "graphique.h"
#include "utilitaires.h"

Scalar Isol_hor::source_psi_hor( const Scalar* p_ener_dens, const Vector* p_mom_dens, 
                const Scalar* p_trace_stress) {

    using namespace Unites ;
   

    // Initialisations
    // ---------------
    //    double ttime = the_time[jtime] ; 
         
    //    trk_evol.update(trk_in, jtime, ttime) ; 

    const Map& map = ff.get_mp() ; 
    const Base_vect& triad = *(ff.get_triad()) ;
    
    Scalar tmp(map) ;
    Scalar tmp_scal(map) ; 
    Sym_tensor tmp_sym(map, CON, triad) ;


    // Reset of quantities depending on K:
    //    k_dd_evol.downdate(jtime) ; 
    //    k_uu_evol.downdate(jtime) ; 

    //    tmp_sym = ( beta().ope_killing_conf(tgam()) + uu ) 
    //                                / (2.* nn()) ; 
        
    //    set_aa(uu / (2.* nn()) ) ; // porque beta es cero, supongo

    Scalar ener_dens(map) ; 
    if (p_ener_dens != 0x0) ener_dens = *(p_ener_dens) ; 
    else ener_dens.set_etat_zero() ; 
    
    Vector mom_dens(map, CON, triad) ; 
    if (p_mom_dens != 0x0) mom_dens = *(p_mom_dens) ; 
    else mom_dens.set_etat_zero() ; 
    
    Scalar trace_stress(map) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
       

    Scalar source_psi(map) ; 
       
    //===============================================
    //  Computations of the source for Psi 
    //===============================================
    
    const Vector& dpsi = psi().derive_cov(ff) ;       // D_i Psi
    //    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    //    const Vector& dnn = nn().derive_cov(ff) ;         // D_i N
    
    Sym_tensor taa = aa().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 

    //    cout<<"Dzpuis de aa_quad:   "<<aa_quad.get_dzpuis() <<endl ;
        
    //    arrete() ;

    // Source for Psi 
    // --------------
    tmp = 0.125* psi() * tgam().ricci_scal() 
      - contract(hh(), 0, 1, dpsi.derive_cov(ff), 0, 1 ) ;
    tmp.inc_dzpuis() ; // dzpuis : 3 -> 4
    
    tmp -= contract(hdirac(), 0, dpsi, 0) ;  
                
    source_psi = tmp - psi()*psi4()* ( 0.5*qpig* ener_dens 
				       + 0.125* aa_quad 
				       - 8.33333333333333e-2* trk()*trk() ) ;

    return source_psi ;

}


Scalar Isol_hor::source_nn_hor( const Scalar& trk_point, const Scalar* p_ener_dens, const Vector* p_mom_dens, 
                const Scalar* p_trace_stress) {

    using namespace Unites ;
   

    // Initialisations
    // ---------------
    //    double ttime = the_time[jtime] ; 
         
    //    trk_evol.update(trk_in, jtime, ttime) ; 

    const Map& map = ff.get_mp() ; 
    const Base_vect& triad = *(ff.get_triad()) ;
    
    Scalar tmp(map) ;
    Scalar tmp_scal(map) ; 
    Sym_tensor tmp_sym(map, CON, triad) ;



    // Reset of quantities depending on K:
    //    k_dd_evol.downdate(jtime) ; 
    //    k_uu_evol.downdate(jtime) ; 

    //    tmp_sym = ( beta().ope_killing_conf(tgam()) + uu ) 
    //                                / (2.* nn()) ; 
        
    //    set_aa(uu / (2.* nn()) ) ; // porque beta es cero, supongo

    Scalar ener_dens(map) ; 
    if (p_ener_dens != 0x0) ener_dens = *(p_ener_dens) ; 
    else ener_dens.set_etat_zero() ; 
    
    Vector mom_dens(map, CON, triad) ; 
    if (p_mom_dens != 0x0) mom_dens = *(p_mom_dens) ; 
    else mom_dens.set_etat_zero() ; 
    
    Scalar trace_stress(map) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
       

    Scalar source_nn(map) ; 
       
    //===============================================
    //  Computations of the source for NN 
    //===============================================
    
    //    const Vector& dpsi = psi().derive_cov(ff) ;       // D_i Psi
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& dnn = nn().derive_cov(ff) ;         // D_i N
    
    Sym_tensor taa = aa().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 

    // Source for N 
    // ------------
        
    //    tmp_scal = 0.3333333333333333* trk()*trk() ;
    //    tmp_scal.inc_dzpuis(4) ;
    

    source_nn = psi4()*( nn()*( qpig* (ener_dens + trace_stress) + aa_quad
				- 0.3333333333333333* trk()*trk() )
			 - trk_point ) 
      - 2.* contract(dln_psi, 0, nn().derive_con(tgam()), 0)  
      - contract(hdirac(), 0, dnn, 0) ; 
        
    tmp = psi4()* contract(beta(), 0, trk().derive_cov(ff), 0) 
      - contract( hh(), 0, 1, dnn.derive_cov(ff), 0, 1 ) ;
        
    tmp.inc_dzpuis() ; // dzpuis: 3 -> 4
        
    source_nn += tmp ;


    return source_nn ;

}



Vector Isol_hor::source_beta_hor( const Scalar* p_ener_dens, const Vector* p_mom_dens, 
                const Scalar* p_trace_stress) {

    using namespace Unites ;
   

    // Initialisations
    // ---------------
    //    double ttime = the_time[jtime] ; 
         
    //    trk_evol.update(trk_in, jtime, ttime) ; 

    const Map& map = ff.get_mp() ; 
    const Base_vect& triad = *(ff.get_triad()) ;
    
    Scalar tmp(map) ;
    Scalar tmp_scal(map) ; 
    Sym_tensor tmp_sym(map, CON, triad) ;
    Vector tmp_vect(map, CON, triad) ;



    // Reset of quantities depending on K:
    //    k_dd_evol.downdate(jtime) ; 
    //    k_uu_evol.downdate(jtime) ; 

    //    tmp_sym = ( beta().ope_killing_conf(tgam()) + uu ) 
    //                                / (2.* nn()) ; 
        
    //    set_aa(uu / (2.* nn()) ) ; // porque beta es cero, supongo

    Scalar ener_dens(map) ; 
    if (p_ener_dens != 0x0) ener_dens = *(p_ener_dens) ; 
    else ener_dens.set_etat_zero() ; 
    
    Vector mom_dens(map, CON, triad) ; 
    if (p_mom_dens != 0x0) mom_dens = *(p_mom_dens) ; 
    else mom_dens.set_etat_zero() ; 
    
    Scalar trace_stress(map) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
       

    Vector source_beta(map, CON, triad) ; 

    //===============================================
    //  Computations of the source for beta 
    //===============================================
    
    //    const Vector& dpsi = psi().derive_cov(ff) ;       // D_i Psi
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& dnn = nn().derive_cov(ff) ;         // D_i N
    
    Sym_tensor taa = aa().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 

    Sym_tensor uu =  2.* nn() * aa() -  beta().ope_killing_conf(tgam())  ;


    // Source for beta 
    // ---------------
    
    source_beta = 2.* contract(aa(), 1, 
			       dnn - 6.*nn() * dln_psi, 0) ;
                
    tmp_vect = 0.66666666666666666* trk().derive_con(tgam()) ;
    tmp_vect.inc_dzpuis() ;

    source_beta += 2.* nn() * ( 2.*qpig* psi4() * mom_dens 
				+ tmp_vect
				- contract(tgam().connect().get_delta(), 1, 2, 
					   aa(), 0, 1) ) ;
            
    Vector vtmp = contract(hh(), 0, 1, 
                           beta().derive_cov(ff).derive_cov(ff), 1, 2)
      + 0.3333333333333333*
      contract(hh(), 1, beta().divergence(ff).derive_cov(ff), 0) 
      - hdirac().derive_lie(beta()) 
      + uu.divergence(ff) ;                       // zero in the Driac gauge
    vtmp.inc_dzpuis() ; // dzpuis: 3 -> 4
    
    source_beta -= vtmp ; 
        
    source_beta += 0.66666666666666666* beta().divergence(ff) * hdirac() ;



    return source_beta ;

}




  
     
        

