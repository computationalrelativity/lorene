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

 

/*
 * $Id$
 * $Log$
 * Revision 1.12  2016/12/05 16:18:19  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.11  2014/10/13 08:53:48  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.10  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.9  2012/02/06 12:59:07  j_novak
 * Correction of some errors.
 *
 * Revision 1.8  2011/07/22 13:21:02  j_novak
 * Corrected an error on BC treatment.
 *
 * Revision 1.7  2010/10/20 07:58:10  j_novak
 * Better implementation of the explicit time-integration. Not fully-tested yet.
 *
 * Revision 1.6  2008/12/04 18:22:49  j_novak
 * Enhancement of the dzpuis treatment + various bug fixes.
 *
 * Revision 1.5  2008/12/02 15:02:22  j_novak
 * Implementation of the new constrained formalism, following Cordero et al. 2009
 * paper. The evolution eqs. are solved as a first-order system. Not tested yet!
 *
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
#include <cassert>

// Lorene headers
#include "time_slice.h"
#include "param.h"
#include "unites.h"
#include "proto.h"
#include "graphique.h"

namespace Lorene {
void Tslice_dirac_max::set_AB_hh(const Scalar& A_in, const Scalar& B_in) {

    A_hh_evol.update(A_in, jtime, the_time[jtime]) ; 
    B_hh_evol.update(B_in, jtime, the_time[jtime]) ; 

    // Reset of quantities depending on h^{ij}:
    hh_evol.downdate(jtime) ;
    trh_evol.downdate(jtime) ;
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

void Tslice_dirac_max::hh_det_one(int j0, Param* par_bc, Param* par_mat) const {

    assert (A_hh_evol.is_known(j0)) ;   // The starting point
    assert (B_hh_evol.is_known(j0)) ;    // of the computation 

    const Map& mp = A_hh_evol[j0].get_mp() ;

    // The representation of h^{ij} as an object of class Sym_tensor_trans :
    Sym_tensor_trans hij(mp, *(ff.get_triad()), ff) ;
    const Scalar* ptrace = 0x0 ;
    if (trh_evol.is_known(j0-1)) ptrace = &trh_evol[j0-1] ;
    hij.set_AtBtt_det_one(A_hh_evol[j0], B_hh_evol[j0], ptrace, par_bc, par_mat) ;

    // Result set to trh_evol and hh_evol
    // ----------------------------------
    trh_evol.update(hij.the_trace(), j0, the_time[j0]) ;
    
    // The longitudinal part of h^{ij}, which is zero by virtue of Dirac gauge :
    Vector wzero(mp, CON,  *(ff.get_triad())) ; 
    wzero.set_etat_zero() ;                   

    // Temporary Sym_tensor with longitudinal part set to zero : 
    Sym_tensor hh_new(mp, CON, *(ff.get_triad())) ;
    hh_new.set_longit_trans(wzero, hij) ;
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

#ifndef NDEBUG         
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
#endif
}

void Tslice_dirac_max::hh_det_one(const Sym_tensor_tt& hijtt, Param* par_mat) 
    const {

    const Map& mp = hijtt.get_mp() ;

    // The representation of h^{ij} as an object of class Sym_tensor_trans :
    Sym_tensor_trans hij(mp, *(ff.get_triad()), ff) ;
    const Scalar* ptrace = 0x0 ;
    if ( trh_evol.is_known( jtime - 1 ) ) ptrace = &trh_evol[jtime-1] ;
    hij.set_tt_part_det_one(hijtt, ptrace, par_mat) ;

    // Result set to trh_evol and hh_evol
    // ----------------------------------
    trh_evol.update(hij.the_trace(), jtime, the_time[jtime]) ;
    
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

#ifndef NDEBUG         
    // Test
    maxabs(tgam().determinant() - 1, 
	   "Max. of absolute value of deviation from det tgam = 1") ; 
#endif
}
                 //----------------------------------------------//
                 //   Equations for h^{ij} and \hat{A}^{ij}      //
                 //----------------------------------------------//

void Tslice_dirac_max::compute_sources( const Sym_tensor* p_strain_tens) const {

    using namespace Unites ;
    
    const Map& map = hh().get_mp() ; 
    const Base_vect& otriad = *hh().get_triad() ;
    int nz = map.get_mg()->get_nzone() ;
    
    Sym_tensor strain_tens(map, CON, otriad) ; 
    if (p_strain_tens != 0x0) 
	strain_tens = *(p_strain_tens) ; 
    else 
	strain_tens.set_etat_zero() ; 

    Sym_tensor aij = aa() ;
    aij.annule_domain(nz-1) ;
    
    const Sym_tensor& tgam_dd = tgam().cov() ;    // {\tilde \gamma}_{ij}
    const Sym_tensor& tgam_uu = tgam().con() ;    // {\tilde \gamma}^{ij}
    const Tensor_sym& dtgam = tgam_dd.derive_cov(ff) ;// D_k {\tilde \gamma}_{ij}
    const Tensor_sym& dhh = hh().derive_cov(ff) ; // D_k h^{ij}
    const Vector& dln_psi = ln_psi().derive_cov(ff) ; // D_i ln(Psi)
    const Vector& tdln_psi_u = ln_psi().derive_con(tgam()) ; // tD^i ln(Psi)
    Scalar log_N = log(nn()) ;
    log_N.std_spectral_base() ;
    const Vector& tdlnn_u = log_N.derive_con(tgam()) ;       // tD^i ln(N)
    const Scalar& div_beta = beta().divergence(ff) ;  // D_k beta^k
    Scalar qq = nn()*psi()*psi() ;
    qq.annule_domain(nz-1) ;
    const Vector& dqq = qq.derive_cov(ff) ;         // D_i Q
    Sym_tensor a_hat = hata() ;
    a_hat.annule_domain(nz-1) ;
    Scalar psi6 = psi4()*psi()*psi() ;
    Sym_tensor sym_tmp(map, CON, otriad) ; 
    Scalar tmp(map) ;

    //==================================
    // Source for hij
    //==================================
    
    Sym_tensor source_hij = hh().derive_lie(beta()) + 2*(nn()/psi6 - 1.)*a_hat 
      - beta().ope_killing_conf(ff) + 0.6666666666666667*div_beta*hh() ;
    source_hij.annule_domain(nz-1) ;
    for (int i=1; i<=3; i++)
 	for (int j=i; j<=3; j++)
 	    source_hij.set( i, j ).set_dzpuis(0) ;

    tmp = 2.*A_hata_evol[jtime] + source_hij.compute_A(true) ;
    tmp.set_spectral_va().ylm() ;
    tmp.annule_domain(nz-1) ;
    tmp.set_dzpuis(0) ;
    source_A_hh_evol.update( tmp, jtime, the_time[jtime] ) ;

    tmp = 2.*B_hata_evol[jtime] + source_hij.compute_tilde_B_tt(true) ;
    tmp.set_spectral_va().ylm() ;
    tmp.annule_domain(nz-1) ;
    tmp.set_dzpuis(0) ;
    source_B_hh_evol.update(tmp, jtime, the_time[jtime] ) ;
    
    //==================================
    // Source for \hat{A}^{ij}
    //==================================
    
    Sym_tensor source_aij = a_hat.derive_lie(beta())  
    	+ 1.666666666666667*a_hat*div_beta ;
    
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
    	    ricci_star.set(i,j) -= tmp ; 
    	}
    }

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
    ricci_star += sym_tmp ; // a factor 1/2 is still missing, shall be put later
    
    // Curvature scalar of conformal metric :
    // -------------------------------------
        
    Scalar tricci_scal = 
    	0.25 * contract(tgam_uu, 0, 1,
    			contract(dhh, 0, 1, dtgam, 0, 1), 0, 1 ) 
    	- 0.5  * contract(tgam_uu, 0, 1,
    			  contract(dhh, 0, 1, dtgam, 0, 2), 0, 1 ) ;  

    Scalar lap_A = A_hh_evol[jtime].laplacian(2) ;
    Scalar tilde_lap_B(map) ;
    tilde_laplacian( B_hh_evol[jtime], tilde_lap_B) ;
    Sym_tensor_tt laplace_h(map, otriad, ff) ;
    laplace_h.set_A_tildeB(lap_A, tilde_lap_B) ;
    laplace_h.annule_domain(nz-1) ;

    //   sym_tmp.inc_dzpuis() ; // dzpuis : 3 --> 4

    source_aij += (0.5*(qq - 1.))*laplace_h + qq*(0.5*ricci_star + 8.*tdln_psi_u * tdln_psi_u 
  	+ 4.*( tdln_psi_u * tdlnn_u + tdlnn_u * tdln_psi_u )
  	-  0.3333333333333333 * (tricci_scal + 8.*(contract(dln_psi, 0, tdln_psi_u, 0) 
  						   + contract(dln_psi, 0, tdlnn_u, 0) ) 
  	    )*tgam_uu
  	) ;
			   
    sym_tmp = contract(tgam_uu, 1, contract(tgam_uu, 1, dqq.derive_cov(ff), 0), 1) ;
    
    for (int i=1; i<=3; i++) {
  	for (int j=1; j<=i; j++) {
  	    tmp = 0 ; 
  	    for (int k=1; k<=3; k++) {
  		for (int l=1; l<=3; l++) {
  		    tmp += ( tgam_uu(i,k)*dhh(l,j,k) + tgam_uu(k,j)*dhh(i,l,k)
  			     - tgam_uu(k,l)*dhh(i,j,k) ) * dqq(l) ; 
  		}
  	    }
  	    sym_tmp.set(i,j) += 0.5 * tmp ; 
  	}
    }
        
  source_aij -= sym_tmp 
    - ( 0.3333333333333333*qq.derive_con(tgam()).divergence(tgam()) ) *tgam_uu ; 
                    
  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
  	for (int l=1; l<=3; l++) {
  	  tmp += tgam_dd(k,l) * a_hat(i,k) * aij(j,l) ; 
  	}
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
        
  tmp = psi4() * strain_tens.trace(tgam()) ; // S = S_i^i 

  source_aij += (2.*nn()) 
    * ( 
       sym_tmp - qpig*psi6*( psi4()* strain_tens - (0.3333333333333333 * tmp) * tgam_uu ) 
	)   ; 

  source_aij.annule_domain(nz-1) ;
  for (int i=1; i<=3; i++)
      for (int j=i; j<=3; j++)
	  source_aij.set(i,j).set_dzpuis(0) ;
#ifndef NDEBUG
  maxabs(source_aij, "source_aij tot") ; 
#endif

  tmp = 0.5*lap_A + source_aij.compute_A(true) ;
  tmp.annule_domain(nz-1) ;
  tmp.set_dzpuis(0) ;
  source_A_hata_evol.update( tmp, jtime, the_time[jtime] ) ;
  tmp = 0.5*tilde_lap_B + source_aij.compute_tilde_B_tt(true) ;
  tmp.annule_domain(nz-1) ;
  tmp.set_dzpuis(0) ;
  // Scalar dess = tmp - tilde_lap_B ;
  // dess.set_spectral_va().ylm_i() ;
  // des_profile(dess, 0, 8., 1, 1) ;
  source_B_hata_evol.update( 0.5*tilde_lap_B, jtime, the_time[jtime] ) ;
}

void Tslice_dirac_max::initialize_sources_copy() const {

    assert( source_A_hh_evol.is_known(jtime) ) ;
    assert( source_B_hh_evol.is_known(jtime) ) ;
    assert( source_A_hata_evol.is_known(jtime) ) ;
    assert( source_B_hata_evol.is_known(jtime) ) ;

    Scalar tmp_Ah = source_A_hh_evol[jtime] ;
    Scalar tmp_Bh = source_B_hh_evol[jtime] ;
    Scalar tmp_Aa = source_A_hata_evol[jtime] ;
    Scalar tmp_Ba = source_B_hata_evol[jtime] ;

    source_A_hh_evol.downdate(jtime) ;
    source_B_hh_evol.downdate(jtime) ;
    source_A_hata_evol.downdate(jtime) ;
    source_B_hata_evol.downdate(jtime) ;

    int jtime1 = jtime - depth + 1; 
    for (int j=jtime1; j <= jtime; j++) {
        source_A_hh_evol.update( tmp_Ah, j, the_time[j] ) ;  
        source_B_hh_evol.update( tmp_Bh, j, the_time[j] ) ;  
        source_A_hata_evol.update( tmp_Aa, j, the_time[j] ) ;
        source_B_hata_evol.update( tmp_Ba, j, the_time[j] ) ;
    } 
}

}
