/*
 *  Methods of class Tslice_dirac_max for solving Einstein equations
 *
 *    (see file time_slice.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon & Jerome Novak
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

char tslice_dirac_max_solve_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.10  2004/06/14 20:47:31  e_gourgoulhon
 * Added argument method_poisson to method solve_hij.
 *
 * Revision 1.9  2004/06/03 10:02:45  j_novak
 * Some filtering is done on source_khi and khi_new.
 *
 * Revision 1.8  2004/05/24 21:00:44  e_gourgoulhon
 * Method solve_hij: khi and mu.smooth_decay(2,2) --> smooth_decay(2,1) ;
 *   added exponential_decay(khi) and exponential_decay(mu) after the
 *   call to smooth_decay. Method exponential_decay is provisory defined
 *   in this file.
 *
 * Revision 1.7  2004/05/17 19:56:25  e_gourgoulhon
 * -- Method solve_beta: added argument method
 * -- Method solve_hij: added argument graph_device
 *
 * Revision 1.6  2004/05/12 15:24:20  e_gourgoulhon
 * Reorganized the #include 's, taking into account that
 * time_slice.h contains now an #include "metric.h".
 *
 * Revision 1.5  2004/05/05 14:47:05  e_gourgoulhon
 * Modified text and graphical outputs.
 *
 * Revision 1.4  2004/05/03 15:06:27  e_gourgoulhon
 * Added matter source in solve_hij.
 *
 * Revision 1.3  2004/05/03 14:50:00  e_gourgoulhon
 * Finished the implementation of method solve_hij.
 *
 * Revision 1.2  2004/04/30 14:36:15  j_novak
 * Added the method Tslice_dirac_max::solve_hij(...)
 * NOT READY YET!!!
 *
 * Revision 1.1  2004/04/30 10:52:14  e_gourgoulhon
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
#include "unites.h"
#include "graphique.h"
#include "proto.h"

void exponential_decay(Scalar& ) ;
void filtre_l(Scalar& , int, int) ;

                    //--------------------------//
                    //      Equation for N      //
                    //--------------------------//

Scalar Tslice_dirac_max::solve_n(const Scalar* p_ener_dens,
                                 const Scalar* p_trace_stress) const {

    using namespace Unites ;

    const Map& map = nn().get_mp() ; 
    Scalar ener_dens(map) ; 
    if (p_ener_dens != 0x0) ener_dens = *(p_ener_dens) ; 
    else ener_dens.set_etat_zero() ; 
    
    Scalar trace_stress(map) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
    
    // Source for N 
    // ------------
        
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& dnn = nn().derive_cov(ff) ;         // D_i N

    Sym_tensor taa = aa().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 
        
    Scalar source_nn = psi4()* nn()*( qpig* (ener_dens + trace_stress) 
                                + aa_quad )
                      - 2.* contract(dln_psi, 0, nn().derive_con(tgam()), 0) ; 
         
    Scalar tmp = - contract( hh(), 0, 1, dnn.derive_cov(ff), 0, 1 ) ;
        
    tmp.inc_dzpuis() ; // dzpuis: 3 -> 4
        
    source_nn += tmp ;
        
    // Resolution of the Poisson equation for the lapse
    // ------------------------------------------------
        
    Scalar nn_new = source_nn.poisson() + 1. ; 

    if (nn_new.get_etat() == ETATUN) nn_new.std_spectral_base() ; 

    // Test:
    maxabs(nn_new.laplacian() - source_nn,
                "Absolute error in the resolution of the equation for N") ;  

    return nn_new ; 

}

                
                    //--------------------------//
                    //      Equation for Q      //
                    //--------------------------//


Scalar Tslice_dirac_max::solve_q(const Scalar* p_trace_stress) const {

    using namespace Unites ;

    Scalar trace_stress(qq().get_mp()) ; 
    if (p_trace_stress != 0x0) trace_stress = *(p_trace_stress) ; 
    else trace_stress.set_etat_zero() ; 
    
    // Source for Q
    // ------------
        
    const Vector& dqq = qq().derive_cov(ff) ;           // D_i Q
    const Vector& dln_psi = ln_psi().derive_cov(ff) ;   // D_i ln(Psi)
    const Tensor_sym& dhh = hh().derive_cov(ff) ;       // D_k h^{ij}
    const Tensor_sym& dtgam = tgam().cov().derive_cov(ff) ;    
                                                    // D_k {\tilde \gamma}_{ij}

    Sym_tensor taa = aa().up_down(tgam()) ; 
        
    Scalar aa_quad = contract(taa, 0, 1, aa(), 0, 1) ; 
        
    Scalar source_qq = psi4() * qq() * ( qpig* trace_stress + 0.75* aa_quad ) ;
        
    Scalar tmp = contract( hh(), 0, 1, dqq.derive_cov(ff), 0, 1 ) ;             
    tmp.inc_dzpuis() ; 
        
    source_qq -= tmp ;  
                        
    tmp = 0.0625 * contract( dhh, 0, 1, dtgam, 0, 1 ).trace(tgam()) 
          - 0.125 * contract( dhh, 0, 1, dtgam, 0, 2 ).trace(tgam()) 
          + 2.* contract( dln_psi, 0, ln_psi().derive_con(tgam()), 0) ; 
     
    source_qq += psi()*psi() * ( nn() * tmp 
                + 2*contract(dln_psi, 0, nn().derive_con(tgam()), 0) ) ; 
                             
               
    // Resolution of the Poisson equation for Q
    // -----------------------------------------
        
    Scalar qq_new = source_qq.poisson() + 1. ; 

    if (qq_new.get_etat() == ETATUN) qq_new.std_spectral_base() ; 

    // Test:
    maxabs(qq_new.laplacian() - source_qq,
                "Absolute error in the resolution of the equation for Q") ;  

    return qq_new ; 

}
                


                    //--------------------------//
                    //      Equation for beta   //
                    //--------------------------//


Vector Tslice_dirac_max::solve_beta(const Vector* p_mom_dens, int method) 
    const {

    using namespace Unites ;

    Vector mom_dens(beta().get_mp(), CON, beta().get_triad()) ; 
    if (p_mom_dens != 0x0) mom_dens = *(p_mom_dens) ; 
    else mom_dens.set_etat_zero() ; 

    // Source for beta
    // ---------------
        
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& dnn = nn().derive_cov(ff) ;         // D_i N
 
    Vector source_beta = 2.* contract(aa(), 1, 
                                   dnn - 6.*nn() * dln_psi, 0) ;
                
    source_beta += 2.* nn() * ( 2.*qpig* psi4() * mom_dens 
                        - contract(tgam().connect().get_delta(), 1, 2, 
                                   aa(), 0, 1) ) ;
            
    Vector vtmp = contract(hh(), 0, 1, 
                           beta().derive_cov(ff).derive_cov(ff), 1, 2)
                + 0.3333333333333333*
                  contract(hh(), 1, beta().divergence(ff).derive_cov(ff), 0) ; 
    vtmp.inc_dzpuis() ; // dzpuis: 3 -> 4
                    
    source_beta -= vtmp ; 
    
    //##
    //for (int i=1; i<=3; i++) {
    //    if (source_beta(i).get_etat() != ETATZERO) {
    //        const Mg3d& mg = *(beta().get_mp().get_mg()) ;
    //        int nz = mg.get_nzone() ;
    //        int nr = mg.get_nr(nz-1) ; 
    //        int nt = mg.get_nt(nz-1) ; 
    //        int np = mg.get_np(nz-1) ; 
    //        Tbl tb(np, nt, nr) ;
    //        tb.annule_hard() ; 
    //        source_beta.set(i).set_domain(nz-1) = tb ; 
    //    }
    // }

    // Resolution of the vector Poisson equation 
    //------------------------------------------
    
    Vector beta_new = source_beta.poisson(0.3333333333333333, ff, method) ; 
        
    // Test:
    Vector test_beta = (beta_new.derive_con(ff)).divergence(ff)
            +  0.3333333333333333 * (beta_new.divergence(ff)).derive_con(ff) ;
    test_beta.inc_dzpuis() ;  
    maxabs(test_beta - source_beta,
                "Absolute error in the resolution of the equation for beta") ; 

    return beta_new ; 

}
                

                    //-------------------------------//
                    //      Equation for h^{ij}      //
                    //-------------------------------//

void Tslice_dirac_max::solve_hij(Param& par_khi, Param& par_mu,
                                 Scalar& khi_new, Scalar& mu_new,
                                 int method_poisson, 
                                 const char* graph_device, 
                                 const Sym_tensor* p_strain_tens) const 
{
  using namespace Unites ;
  
  const Map& map = hh().get_mp() ; 
  const Base_vect& otriad = *hh().get_triad() ;
    
  // For graphical outputs:
  //  int ngraph0 = 40 ;  // index of the first graphic device to be used
    int nz = map.get_mg()->get_nzone() ; 
    //   double ray_des = map.val_r(nz-2, 1., 0., 0.) ; // outermost radius
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
        
  // des_meridian(hh_point(1,1), 0., ray_des, "dot h\\urr\\d", ngraph0, 
  //            graph_device) ; 
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
    maxabs( source_hh.transverse(ff, 0x0, method_poisson).divergence(ff), 
                "Divergence of source_hh_transverse") ;
    maxabs( source_hh.transverse(ff).trace(ff), 
                "Trace of source_hh_transverse") ; 


    //=============================================
    // Resolution of wave equation for h
    //=============================================
    
    const Sym_tensor_tt& source_htt = source_hh.transverse(ff, 0x0, 
                                                method_poisson).tt_part() ;
         
    maxabs( source_htt.divergence(ff), "Divergence of source_htt") ; 
    maxabs( source_htt.trace(ff), "Trace of source_hhtt") ; 

    int *nfiltre = new int[nz] ;
    nfiltre[0] = 0 ;
    nfiltre[nz-1] = 0 ;
    for (int lz=1; lz<nz-1; lz++)
      nfiltre[lz] = map.get_mg()->get_nr(lz) / 3 + 1 ;
    
    Scalar khi_source = source_htt.khi() ; 
    filtre_l(khi_source, 0, 1) ;
    khi_source.filtre_r(nfiltre) ;

    const Scalar& mu_source = source_htt.mu() ; 
                            
    khi_new = khi_evol[jtime].avance_dalembert(par_khi,
                                         khi_evol[jtime-1], khi_source) ;
    khi_new.smooth_decay(2,1) ; 
    exponential_decay(khi_new) ; 
    khi_new.filtre_r(nfiltre) ;
        
    maxabs(khi_new - khi_evol[jtime], "Variation of khi") ;  
        
    mu_new = mu_evol[jtime].avance_dalembert(par_mu,  
                                         mu_evol[jtime-1], mu_source) ;
    mu_new.smooth_decay(2,1) ; 
    exponential_decay(mu_new) ; 
                                        

}


void exponential_decay(Scalar& uu) {

    const Map& mp = uu.get_mp() ; 
    const Mg3d& mg = *(mp.get_mg()) ; 
    int nz = mg.get_nzone() ;
    int nzm1 = nz - 1 ;  

    const Map_af* mapaf  = dynamic_cast<const Map_af*>(&mp) ;
    if (mapaf == 0x0) {
	    cout << "exponential_decay: present version supports only \n" 
         << "  affine mappings !" << endl ;
	    abort() ;
    }
    
    double rbound = mapaf->get_alpha()[nzm1-1] + mapaf->get_beta()[nzm1-1] ; 

    Mtbl xx = mp.r - rbound ; 
    Scalar tmp(mp) ; 
    tmp = exp( - xx*xx ) ; 
    tmp.annule(0, nz-2) ; 
    tmp *= uu ; 
    uu.set_domain(nzm1) = tmp.domain(nzm1) ; 
    
}

void filtre_l(Scalar& uu, int l_min, int l_max) {

    const Map& mp = uu.get_mp() ; 
    const Mg3d& mg = *(mp.get_mg()) ; 
    int nz = mg.get_nzone() ;

    Valeur& uuva = uu.set_spectral_va() ;
    if (uuva.get_etat() != ETATZERO) {
      assert (uuva.get_etat() == ETATQCQ) ;
      uuva.ylm() ;
      for (int lz=0; lz<nz; lz++) {
	int np = mg.get_np(lz) ;
	int nt = mg.get_nt(lz) ;
	int nr = mg.get_nr(lz) ;
	if (uuva.c_cf->operator()(lz).get_etat() != ETATZERO)
	  for (int k=0; k<np+1; k++) 
	    for (int j=0; j<nt; j++) {
	      int l_q, m_q, base_r ;
	      if (nullite_plm(j, nt, k, np, uuva.base) == 1) {
		donne_lm(nz, lz, j, k, uuva.base, m_q, l_q, base_r) ;
		if ( (l_q>=l_min) && (l_q <=l_max) ) 
		  for (int i=0; i<nr; i++) 
		    uuva.c_cf->set(lz, k, j, i) = 0. ;
	      }
	    } 
      }
      if (uuva.c != 0x0) {
	delete uuva.c ;
	uuva.c = 0x0 ;
      }
    }
}
