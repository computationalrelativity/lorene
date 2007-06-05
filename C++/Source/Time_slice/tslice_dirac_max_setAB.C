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
 * Revision 1.4  2007/06/05 07:38:37  j_novak
 * Better treatment of dzpuis for A and tilde(B) potentials. Some errors in the bases manipulation have been also corrected.
 *
 * Revision 1.3  2007/05/24 12:10:41  j_novak
 * Update of khi_evol and mu_evol.
 *
 * Revision 1.2  2007/04/25 15:21:01  j_novak
 * Corrected an error in the initialization of tildeB in
 * Tslice_dirac_max::initial_dat_cts. + New method for solve_hij_AB.
 *
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
#include "unites.h"
#include "graphique.h"

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

    // The representation of h^{ij} as an object of class Sym_tensor_trans :
    Sym_tensor_trans hij(mp, *(ff.get_triad()), ff) ;
    const Scalar* ptrace = 0x0 ;
    if (trh_evol.is_known(j0-1)) ptrace = &trh_evol[j0-1] ;
    hij.set_AtBtt_det_one(potA_evol[j0], tildeB_evol[j0], ptrace, par_bc, par_mat) ;

    Scalar khi_new = hij.tt_part()(1,1) ;
    khi_new.mult_r() ;
    khi_new.mult_r() ;
    khi_evol.update(khi_new, j0, the_time[j0]) ;
    mu_evol.update(hij.tt_part().mu(), j0, the_time[j0]) ;

    // Result set to trh_evol and hh_evol
    // ----------------------------------
    trh_evol.update(hij.the_trace(), j0, the_time[j0]) ;
    
    // The longitudinal part of h^{ij}, which is zero by virtue of Dirac gauge :
    Vector wzero(mp, CON,  *(ff.get_triad())) ; 
    wzero.set_etat_zero() ;                   

    // Temporary Sym_tensor with longitudinal part set to zero : 
    Sym_tensor hh_new(mp, CON, *(ff.get_triad())) ;
        hh_new = hij ;

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
                    //-------------------------------//
                    //      Equation for h^{ij}      //
                    //-------------------------------//

void Tslice_dirac_max::solve_hij_AB(Param& par_A, Param& par_B,
                                 Scalar& A_new, Scalar& B_new,
                                 const char* graph_device, 
                                 const Sym_tensor* p_strain_tens) const 
{
  using namespace Unites ;
  
  const Map& map = hh().get_mp() ; 
  const Base_vect& otriad = *hh().get_triad() ;
    
  // For graphical outputs:
    int ngraph0 = 80 ;  // index of the first graphic device to be used
    int nz = map.get_mg()->get_nzone() ; 
    double ray_des = 1.25*map.val_r(nz-2, 1., 0., 0.) ; // outermost radius
                                                          // for plots

  Sym_tensor strain_tens(map, CON, otriad) ; 
  if (p_strain_tens != 0x0) strain_tens = *(p_strain_tens) ; 
  else strain_tens.set_etat_zero() ; 

  // Time derivatives:
  Scalar nn_point = n_evol.time_derive(jtime, scheme_order) ; 
  nn_point.inc_dzpuis(2) ; // dzpuis : 0 -> 2
  
  Vector beta_point = beta_evol.time_derive(jtime, scheme_order) ; 
  beta_point.inc_dzpuis(2) ; // dzpuis : 0 -> 2
  
  Sym_tensor hh_point = hh_evol.time_derive(jtime, scheme_order) ; 
  
  hh_point.inc_dzpuis(2) ; // dzpuis : 0 -> 2
        
  des_meridian(hh_point(1,1), 0., ray_des, "dot h\\urr\\d", ngraph0, 
		graph_device) ; 
  // des_meridian(hh_point(2,3), 0., ray_des, "dot h\\u\\gh\\gf\\d", ngraph0+1,
  //            graph_device) ; 
  // des_meridian(hh_point(3,3), 0., ray_des, "dot h\\u\\gf\\gf\\d", ngraph0+2,
  //            graph_device) ; 
  
  //==================================
  // Source for hij
  //==================================
        
  const Sym_tensor& tgam_dd = tgam().cov() ;    // {\tilde \gamma}_{ij}
  const Sym_tensor& tgam_uu = tgam().con() ;    // {\tilde \gamma}^{ij}
  const Tensor_sym& dtgam = tgam_dd.derive_cov(ff) ;// D_k {\tilde \gamma}_{ij}
  const Tensor_sym& dhh = hh().derive_cov(ff) ; // D_k h^{ij}
  const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
  const Vector& tdln_psi_u = ln_psi().derive_con(tgam()) ; // tD^i ln(Psi)
  const Vector& tdnn_u = nn().derive_con(tgam()) ;       // tD^i N
  const Vector& dqq = qq().derive_cov(ff) ;         // D_i Q
  const Scalar& div_beta = beta().divergence(ff) ;  // D_k beta^k

  Scalar tmp(map) ;
  Sym_tensor sym_tmp(map, CON, otriad) ; 

  // Quadratic part of the Ricci tensor of gam_tilde 
  // ------------------------------------------------
        
  Sym_tensor ricci_star(map, CON, otriad) ; 
        
  ricci_star = contract(hh(), 0, 1, dhh.derive_cov(ff), 2, 3) ; 

  ricci_star.inc_dzpuis() ;   // dzpuis : 3 --> 4

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  tmp += dhh(i,k,l) * dhh(j,l,k) ; 
	}
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
  ricci_star -= sym_tmp ;

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  for (int m=1; m<=3; m++) {
	    for (int n=1; n<=3; n++) {
                            
     tmp += 0.5 * tgam_uu(i,k)* tgam_uu(j,l) 
       * dhh(m,n,k) * dtgam(m,n,l)
       + tgam_dd(n,l) * dhh(m,n,k) 
       * (tgam_uu(i,k) * dhh(j,l,m) + tgam_uu(j,k) *  dhh(i,l,m) )
       - tgam_dd(k,l) *tgam_uu(m,n) * dhh(i,k,m) * dhh(j,l,n) ;
	    }
	  } 
	}
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
  ricci_star += sym_tmp ;

  ricci_star = 0.5 * ricci_star ; 
        
  // Curvature scalar of conformal metric :
  // -------------------------------------
        
  Scalar tricci_scal = 
    0.25 * contract(tgam_uu, 0, 1,
		    contract(dhh, 0, 1, dtgam, 0, 1), 0, 1 ) 
    - 0.5  * contract(tgam_uu, 0, 1,
		      contract(dhh, 0, 1, dtgam, 0, 2), 0, 1 ) ;  
                                                       
  // Full quadratic part of source for h : S^{ij}
  // --------------------------------------------
        
  Sym_tensor ss(map, CON, otriad) ; 
        
  sym_tmp = nn() * (ricci_star + 8.* tdln_psi_u * tdln_psi_u)
    + 4.*( tdln_psi_u * tdnn_u + tdnn_u * tdln_psi_u ) 
    - 0.3333333333333333 * 
    ( nn() * (tricci_scal  + 8.* contract(dln_psi, 0, tdln_psi_u, 0) )
      + 8.* contract(dln_psi, 0, tdnn_u, 0) ) *tgam_uu ;

  ss = sym_tmp / psi4()  ;
        
  sym_tmp = contract(tgam_uu, 1, 
		     contract(tgam_uu, 1, dqq.derive_cov(ff), 0), 1) ;
                            
  sym_tmp.inc_dzpuis() ; // dzpuis : 3 --> 4
        
  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  tmp += ( hh()(i,k)*dhh(l,j,k) + hh()(k,j)*dhh(i,l,k)
		   - hh()(k,l)*dhh(i,j,k) ) * dqq(l) ; 
	}
      }
      sym_tmp.set(i,j) += 0.5 * tmp ; 
    }
  }
        
  tmp = qq().derive_con(tgam()).divergence(tgam()) ; 
  tmp.inc_dzpuis() ; // dzpuis : 3 --> 4
        
  sym_tmp -= 0.3333333333333333 * tmp *tgam_uu ; 
                    
  ss -= sym_tmp / (psi4()*psi()*psi()) ; 

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  tmp += tgam_dd(k,l) * aa()(i,k) * aa()(j,l) ; 
	}
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
        
  tmp = psi4() * strain_tens.trace(tgam()) ; // S = S_i^i 

  ss += (2.*nn()) * ( sym_tmp - qpig*( psi4()* strain_tens 
                                       - 0.3333333333333333 * tmp * tgam_uu ) 
                    )   ; 

  maxabs(ss, "ss tot") ; 
  
  // Source for h^{ij} 
  // -----------------
                 
  Sym_tensor lbh = hh().derive_lie(beta()) ; 

  Sym_tensor source_hh = (nn()*nn()/psi4() - 1.) 
    * hh().derive_con(ff).divergence(ff) 
    + 2.* hh_point.derive_lie(beta()) - lbh.derive_lie(beta()) ;
  source_hh.inc_dzpuis() ; 
        
  source_hh += 2.* nn() * ss ;
              
  //## Provisory: waiting for the Lie derivative to allow
  //  derivation with respect to a vector with dzpuis != 0
  Vector vtmp = beta_point ; 
  vtmp.dec_dzpuis(2) ; 
  sym_tmp = hh().derive_lie(vtmp) ; 
  sym_tmp.inc_dzpuis(2) ;             

  source_hh += sym_tmp 
    + 1.3333333333333333 * div_beta* (hh_point - lbh)
    + 2. * (nn_point - nn().derive_lie(beta())) * aa()  ;
              

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	tmp += ( hh().derive_con(ff)(k,j,i) 
		 + hh().derive_con(ff)(i,k,j) 
		 - hh().derive_con(ff)(i,j,k) ) * dqq(k) ;
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
            
  source_hh -= nn() / (psi4()*psi()*psi()) * sym_tmp ; 
         
  tmp =  beta_point.divergence(ff) - div_beta.derive_lie(beta()) ; 
  tmp.inc_dzpuis() ; 
  source_hh += 0.6666666666666666* 
    ( tmp - 0.6666666666666666* div_beta * div_beta ) * hh() ; 
               
        
  // Term (d/dt - Lie_beta) (L beta)^{ij}--> sym_tmp        
  // ---------------------------------------
  Sym_tensor l_beta = beta().ope_killing_conf(ff) ; 

  sym_tmp = beta_point.ope_killing_conf(ff) - l_beta.derive_lie(beta()) ;
  
  sym_tmp.inc_dzpuis() ; 
  
  // Final source:
  // ------------
  source_hh += 0.6666666666666666* div_beta * l_beta - sym_tmp ; 
           
  maxabs(hh(), "h^{ij}") ;
  maxabs(source_hh, "Maxabs source_hh") ; 

  maxabs( source_hh.divergence(ff), "Divergence of source_hh") ; 
                
    //=============================================
    // Resolution of wave equation for h
    //=============================================
    
  Scalar A_source = source_hh.compute_A(true) ; 
//  A_source.annule_extern_cn(nz-2, 1) ;  
//    filtre_l(A_source, 0, 1, true) ;
//    A_source.filtre_r(nfiltre) ;
  A_new = potA_evol[jtime].avance_dalembert(par_A, potA_evol[jtime-1], A_source) ;
//    filtre_l(A_new, 0, 1) ;
//  A_new.set_spectral_va().ylm_i() ;
  maxabs(A_new - potA_evol[jtime], "Variation of A") ;  
  
  Scalar B_source = source_hh.compute_tilde_B_tt(true) ;      
//  filtre_l(B_source, 0, 1) ;
  B_new = tildeB_evol[jtime].avance_dalembert(par_B, tildeB_evol[jtime-1], B_source) ;
//    filtre_l(B_new, 0, 1) ;
//  B_new.set_spectral_va().ylm_i() ;
  maxabs(B_new - tildeB_evol[jtime], "Variation of tilde(B)") ;  
                                        
}


