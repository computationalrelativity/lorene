/*
 * Methods Star_bin::update_metric_der_comp
 *
 * (see file star.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004 Francois Limousin
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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
 * Revision 1.18  2016/12/05 16:18:15  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.17  2014/10/13 08:53:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.16  2014/10/06 15:13:17  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.15  2006/05/24 16:52:50  f_limousin
 * New computation of tkij_comp
 *
 * Revision 1.14  2005/09/13 19:38:31  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.13  2005/02/24 16:26:46  f_limousin
 * Add the name of the variable for the companion which is now used
 * to compute dlogn_comp, dlnq_comp...
 *
 * Revision 1.12  2005/02/24 16:07:23  f_limousin
 * Improve the computation of dlogn, dlnq and dlnpsi. We import the part
 * coming from the companion and add the auto part.
 *
 * Revision 1.11  2005/02/18 13:14:18  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
 * Revision 1.10  2005/02/17 17:34:28  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.9  2004/06/22 12:52:47  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
 * Revision 1.8  2004/06/07 16:25:14  f_limousin
 * Minor modif.
 *
 * Revision 1.7  2004/04/08 16:33:32  f_limousin
 * The new variable is ln(Q) instead of Q=psi^2*N. It improves the
 * convergence of the code.
 *
 * Revision 1.6  2004/03/23 10:00:09  f_limousin
 * We now make the derivation with respect to the metric tilde
 * instead of the flat metric for the computation of dshift_comp.
 *
 * Revision 1.5  2004/02/27 09:53:14  f_limousin
 * Correction of an error on the computation of kcar_comp.
 *
 * Revision 1.4  2004/02/18 18:47:01  e_gourgoulhon
 * divshift_comp now computed via Tensor::divergence, the
 * method Tensor::scontract having disappeared.
 *
 * Revision 1.3  2004/01/20 15:20:23  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
#include <cmath>

// Headers Lorene
#include "star.h"
#include "utilitaires.h"
#include "graphique.h"

namespace Lorene {
void Star_bin::update_metric_der_comp(const Star_bin& comp, double om) {
  
  // Derivatives of metric coefficients
  // ----------------------------------
  
    // dcov_logn 
    Vector temp ((comp.logn_auto).derive_cov(comp.get_flat())) ;
    temp.dec_dzpuis(2) ;
    temp.change_triad(temp.get_mp().get_bvect_cart()) ;
    temp.change_triad(mp.get_bvect_cart()) ;
    Base_val sauve_base1 (temp(1).get_spectral_va().get_base()) ;
    Base_val sauve_base2 (temp(2).get_spectral_va().get_base()) ;
    Base_val sauve_base3 (temp(3).get_spectral_va().get_base()) ;
    assert ( *(temp.get_triad()) == *(dcov_logn.get_triad())) ;
	
    dcov_logn.set(1).import(temp(1)) ;
    dcov_logn.set(2).import(temp(2)) ;
    dcov_logn.set(3).import(temp(3)) ;

    dcov_logn.set(1).set_spectral_va().set_base(sauve_base1) ;
    dcov_logn.set(2).set_spectral_va().set_base(sauve_base2) ;
    dcov_logn.set(3).set_spectral_va().set_base(sauve_base3) ;
    dcov_logn.inc_dzpuis(2) ;

    dcov_logn = dcov_logn + logn_auto.derive_cov(flat) ;

    // dcon_logn 
    Vector temp_con((comp.logn_auto).derive_con(comp.get_flat())) ;
    temp_con.dec_dzpuis(2) ;
    temp_con.change_triad(temp_con.get_mp().get_bvect_cart()) ;
    temp_con.change_triad(mp.get_bvect_cart()) ;
    sauve_base1 = temp_con(1).get_spectral_va().get_base() ;
    sauve_base2 = temp_con(2).get_spectral_va().get_base() ;
    sauve_base3 = temp_con(3).get_spectral_va().get_base() ;
    assert ( *(temp_con.get_triad()) == *(dcon_logn.get_triad())) ;
	
    dcon_logn.set(1).import(temp_con(1)) ;
    dcon_logn.set(2).import(temp_con(2)) ;
    dcon_logn.set(3).import(temp_con(3)) ;

    dcon_logn.set(1).set_spectral_va().set_base(sauve_base1) ;
    dcon_logn.set(2).set_spectral_va().set_base(sauve_base2) ;
    dcon_logn.set(3).set_spectral_va().set_base(sauve_base3) ;
    dcon_logn.inc_dzpuis(2) ;

    dcon_logn = dcon_logn + logn_auto.derive_con(flat) ;

    // dcov_phi 
    temp = (0.5*(comp.lnq_auto-comp.logn_auto)).derive_cov(comp.get_flat()) ;
    temp.dec_dzpuis(2) ;
    temp.change_triad(temp.get_mp().get_bvect_cart()) ;
    temp.change_triad(mp.get_bvect_cart()) ;
    sauve_base1 = temp(1).get_spectral_va().get_base() ;
    sauve_base2 = temp(2).get_spectral_va().get_base() ;
    sauve_base3 = temp(3).get_spectral_va().get_base() ;
    assert ( *(temp.get_triad()) == *(dcov_phi.get_triad())) ;
	
    dcov_phi.set(1).import(temp(1)) ;
    dcov_phi.set(2).import(temp(2)) ;
    dcov_phi.set(3).import(temp(3)) ;

    dcov_phi.set(1).set_spectral_va().set_base(sauve_base1) ;
    dcov_phi.set(2).set_spectral_va().set_base(sauve_base2) ;
    dcov_phi.set(3).set_spectral_va().set_base(sauve_base3) ;
    dcov_phi.inc_dzpuis(2) ;

    dcov_phi = dcov_phi + 0.5*(lnq_auto - logn_auto).derive_cov(flat) ;

    // dcon_phi 
    temp_con = (0.5*(comp.lnq_auto-comp.logn_auto))
	.derive_con(comp.get_flat()) ;
    temp_con.dec_dzpuis(2) ;
    temp_con.change_triad(temp_con.get_mp().get_bvect_cart()) ;
    temp_con.change_triad(mp.get_bvect_cart()) ;
    sauve_base1 = temp_con(1).get_spectral_va().get_base() ;
    sauve_base2 = temp_con(2).get_spectral_va().get_base() ;
    sauve_base3 = temp_con(3).get_spectral_va().get_base() ;
    assert ( *(temp_con.get_triad()) == *(dcon_phi.get_triad())) ;
	
    dcon_phi.set(1).import(temp_con(1)) ;
    dcon_phi.set(2).import(temp_con(2)) ;
    dcon_phi.set(3).import(temp_con(3)) ;

    dcon_phi.set(1).set_spectral_va().set_base(sauve_base1) ;
    dcon_phi.set(2).set_spectral_va().set_base(sauve_base2) ;
    dcon_phi.set(3).set_spectral_va().set_base(sauve_base3) ;
    dcon_phi.inc_dzpuis(2) ;

    dcon_phi = dcon_phi + 0.5*(lnq_auto - logn_auto).derive_con(flat) ;

    
  // Construction of Omega d/dphi
  // ----------------------------
  
  const Mg3d* mg = mp.get_mg() ; 
  int nz = mg->get_nzone() ;	    // total number of domains
  Vector omdsdp (mp, CON, mp.get_bvect_cart()) ;
  Scalar yya (mp) ;
  yya = mp.ya ;
  Scalar xxa (mp) ;
  xxa = mp.xa ;
  
  if (fabs(mp.get_rot_phi()) < 1e-10){ 
    omdsdp.set(1) = - om * yya ;
    omdsdp.set(2) = om * xxa ;
    omdsdp.set(3).annule_hard() ;
  }
  else{
    omdsdp.set(1) = om * yya ;
    omdsdp.set(2) = - om * xxa ;
    omdsdp.set(3).annule_hard() ;
  }
  
  omdsdp.set(1).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[0])) ;
  omdsdp.set(2).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[1])) ;
  omdsdp.set(3).set_spectral_va()
    .set_base(*(mp.get_mg()->std_base_vect_cart()[2])) ;
  
  omdsdp.annule_domain(nz-1) ;
  
  // Computation of tkij_comp
  // ------------------------
  
  // Gradient tilde (partial derivatives with respect to
  //           the cartesian coordinates of the mapping)
  // D~^j beta^i
    
  const Tensor& dbeta_comp = beta_comp.derive_con(gtilde) ;
                        
  // Trace of D~_j beta^i  :
  Scalar divbeta_comp = beta_comp.divergence(gtilde) ;
		  
  // Computation of K^{ij}
  // -------------------------
  
  for (int i=1; i<=3; i++) 
    for (int j=i; j<=3; j++) {
      
      tkij_comp.set(i, j) = dbeta_comp(i, j) + dbeta_comp(j, i) - 
	double(2) /double(3) * divbeta_comp * (gtilde.con())(i,j) ; 
    }
  
  // Addition (or not !) of u^{ij}
  tkij_comp = tkij_comp - 0*hij_comp.derive_lie(omdsdp) ;

  tkij_comp = 0.5 * tkij_comp / nn ;
  
  /*
  Sym_tensor aa_comp (mp, CON, mp.get_bvect_cart()) ;
  aa_comp.set_etat_qcq() ;
  Sym_tensor comp_aa (comp.get_tkij_auto()) ; 
  comp_aa.dec_dzpuis(2) ;
  comp_aa.change_triad(mp.get_bvect_cart()) ;
 
  
  assert(*(aa_comp.get_triad()) == *(comp_aa.get_triad())) ;
  // importations :
  for (int i=1 ; i<=3 ; i++)
    for (int j=i ; j<=3 ; j++) {
      aa_comp.set(i, j).import(comp_aa(i, j)) ;
      aa_comp.set(i, j).set_spectral_va().set_base(comp_aa(i, j).
					  get_spectral_va().get_base()) ;
    }
  aa_comp.inc_dzpuis(2) ;

  for (int i=1 ; i<=3 ; i++)
    for (int j=i ; j<=3 ; j++)
      for (int l=0 ; l<=nz-2 ; l++) 
	tkij_comp.set(i,j).set_domain(l) = aa_comp(i,j).domain(l) ;
  */


  // Computation of kcar_comp
  // ------------------------
  
  Sym_tensor tkij_auto_cov = tkij_auto.up_down(gtilde) ;
  
  kcar_comp = contract(tkij_auto_cov, 0, 1, tkij_comp, 0, 1,true) ; 
  
  // The derived quantities are obsolete
  // -----------------------------------
  
  del_deriv() ;
  
}      

}
