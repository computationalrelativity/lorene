/*
 *  Solution of the tensor Poisson equation for rotating stars in Dirac gauge.
 *
 *    (see file star_rot_dirac.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005 Lap-Ming Lin & Jerome Novak
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

char strot_dirac_solvehij_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2005/04/20 14:26:29  j_novak
 * Removed some outputs.
 *
 * Revision 1.5  2005/04/05 09:24:05  j_novak
 * minor modifs
 *
 * Revision 1.4  2005/02/17 17:32:16  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.3  2005/02/16 12:51:49  j_novak
 * Corrected an error on the matter source terms.
 *
 * Revision 1.2  2005/02/09 13:34:56  lm_lin
 *
 * Remove the Laplacian of hij from the term source_hh and fix an overall
 * minus error.
 *
 * Revision 1.1  2005/01/31 08:51:48  j_novak
 * New files for rotating stars in Dirac gauge (still under developement).
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "star_rot_dirac.h"
#include "unites.h"

void Star_rot_Dirac::solve_hij(Sym_tensor_trans& hij_new) const {

    using namespace Unites ;

    const Metric_flat& mets = mp.flat_met_spher() ;
    const Base_vect_spher& bspher = mp.get_bvect_spher() ;
    
  //==================================
  // Source for hij
  //==================================
    
    const Vector& dln_psi = ln_psi.derive_cov(mets) ; // D_i ln(Psi)
    const Vector& dqq = qqq.derive_cov(mets) ;           // D_i Q
    const Tensor_sym& dhh = hh.derive_cov(mets) ;       // D_k h^{ij}
    const Tensor_sym& dtgam = tgamma.cov().derive_cov(mets) ;    
    const Sym_tensor& tgam_dd = tgamma.cov() ;    // {\tilde \gamma}_{ij}
    const Sym_tensor& tgam_uu = tgamma.con() ;    // {\tilde \gamma}^{ij}
    const Vector& tdln_psi_u = ln_psi.derive_con(tgamma) ; // tD^i ln(Psi)
    const Vector& tdnn_u = nn.derive_con(tgamma) ;       // tD^i N
//    const Scalar& div_beta = beta.divergence(mets) ;  // D_k beta^k
    Scalar div_beta(mp) ; //##
    div_beta.set_etat_zero() ;

    Scalar tmp(mp) ;
    Sym_tensor sym_tmp(mp, CON, bspher) ; 

  // Quadratic part of the Ricci tensor of gam_tilde 
  // ------------------------------------------------
        
  Sym_tensor ricci_star(mp, CON, bspher) ; 
        
  ricci_star = contract(hh, 0, 1, dhh.derive_cov(mets), 2, 3) ; 

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
        
  Sym_tensor ss(mp, CON, bspher) ; 
        
  sym_tmp = nn * (ricci_star + 8.* tdln_psi_u * tdln_psi_u)
    + 4.*( tdln_psi_u * tdnn_u + tdnn_u * tdln_psi_u ) 
    - 0.3333333333333333 * 
    ( nn * (tricci_scal  + 8.* contract(dln_psi, 0, tdln_psi_u, 0) )
      + 8.* contract(dln_psi, 0, tdnn_u, 0) ) *tgam_uu ;

  ss = sym_tmp / psi4  ;
        
  sym_tmp = contract(tgam_uu, 1, 
		     contract(tgam_uu, 1, dqq.derive_cov(mets), 0), 1) ;
                            
  sym_tmp.inc_dzpuis() ; // dzpuis : 3 --> 4
        
  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  tmp += ( hh(i,k)*dhh(l,j,k) + hh(k,j)*dhh(i,l,k)
		   - hh(k,l)*dhh(i,j,k) ) * dqq(l) ; 
	}
      }
      sym_tmp.set(i,j) += 0.5 * tmp ; 
    }
  }
        
  tmp = qqq.derive_con(tgamma).divergence(tgamma) ; 
  tmp.inc_dzpuis() ; // dzpuis : 3 --> 4
        
  sym_tmp -= 0.3333333333333333 * tmp *tgam_uu ; 
                    
  ss -= sym_tmp / (psi4*psi2) ; 

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	for (int l=1; l<=3; l++) {
	  tmp += tgam_dd(k,l) * aa(i,k) * aa(j,l) ; 
	}
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
        
  ss += (2.*nn) * ( sym_tmp - qpig*( psi4* stress_euler
                                       - 0.3333333333333333 * s_euler * tgam_uu ) 
                    )   ; 

//  maxabs(ss, "ss tot") ; 
  
  // Source for h^{ij} 
  // -----------------
                 
  Sym_tensor lbh = hh.derive_lie(beta) ; 

  Sym_tensor source_hh = - lbh.derive_lie(beta) ;
  source_hh.inc_dzpuis() ; 
        
  source_hh += 2.* nn * ss ;
              
  source_hh += - 1.3333333333333333 * div_beta* lbh
    - 2. * nn.derive_lie(beta) * aa  ;
              

  for (int i=1; i<=3; i++) {
    for (int j=1; j<=i; j++) {
      tmp = 0 ; 
      for (int k=1; k<=3; k++) {
	tmp += ( hh.derive_con(mets)(k,j,i) 
		 + hh.derive_con(mets)(i,k,j) 
		 - hh.derive_con(mets)(i,j,k) ) * dqq(k) ;
      }
      sym_tmp.set(i,j) = tmp ; 
    }
  }
            
  source_hh -= nn / (psi4*psi2) * sym_tmp ; 
         
  tmp = - div_beta.derive_lie(beta) ; 
  tmp.inc_dzpuis() ; 
  source_hh += 0.6666666666666666* 
    ( tmp - 0.6666666666666666* div_beta * div_beta ) * hh ; 
               
        
  // Term (d/dt - Lie_beta) (L beta)^{ij}--> sym_tmp        
  // ---------------------------------------
  Sym_tensor l_beta = beta.ope_killing_conf(mets) ; 

  sym_tmp = - l_beta.derive_lie(beta) ;

  sym_tmp.inc_dzpuis() ; 
  
  // Final source:
  // ------------
  source_hh += 0.6666666666666666* div_beta * l_beta - sym_tmp ; 

  source_hh = - ( psi4/nn/nn )*source_hh ;

  cout << " Max( trace of source_hh ) " << endl ;
  cout << max(abs(source_hh.trace(mets))) ;

  for (int i=1; i<=3; i++)
      for (int j=i; j<=3; j++) {
	  source_hh.set(i,j).set_dzpuis(4) ;
      }
  
    Sym_tensor_trans source_hht(mp, bspher, mets) ;
    source_hht = source_hh ;
 //  const Sym_tensor_tt& source_htt = source_hh.transverse(mets, 0x0, 4).tt_part() ;
    const Sym_tensor_tt& source_htt = source_hht.tt_part() ;

  maxabs(source_htt, "source_htt tot") ; 
  
  cout << " Max( divergence of source_hh ) " << endl ;
  for (int i=1; i<=3; i++) 
      cout << max(abs(source_htt.divergence(mets)(i))) ;


  hij_new.trace_from_det_one(source_htt.poisson()) ;
  hij_new.dec_dzpuis(2) ;

}
