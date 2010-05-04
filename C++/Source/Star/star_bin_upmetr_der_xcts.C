/*
 * Methods Star_bin_xcts::update_metric_der_comp
 * (see file star.h for documentation)
 */

/*
 *   Copyright (c) 2010 Michal Bejger
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

char star_bin_upmetr_der_xcts_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2010/05/04 07:51:05  m_bejger
 * Initial version
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Headers Lorene
#include "star.h"
#include "utilitaires.h"
#include "graphique.h"

void Star_bin_xcts::update_metric_der_comp(const Star_bin_xcts& comp) {
  
  // Derivatives of metric coefficients
  // ----------------------------------
      
    // dcov_Psi
    
    Vector temp = (comp.Psi_auto).derive_cov(comp.flat) ;
    temp.dec_dzpuis(2) ;
    
    temp.change_triad(temp.get_mp().get_bvect_cart()) ;
    temp.change_triad(mp.get_bvect_cart()) ;
    Base_val sauve_base1 (temp(1).get_spectral_va().get_base()) ;
    Base_val sauve_base2 (temp(2).get_spectral_va().get_base()) ;
    Base_val sauve_base3 (temp(3).get_spectral_va().get_base()) ;
    assert ( *(temp.get_triad()) == *(dcov_Psi.get_triad())) ;
	
	for(int i=1; i<=3; i++) dcov_Psi.set(i).import(temp(i)) ;

    dcov_Psi.set(1).set_spectral_va().set_base(sauve_base1) ;
    dcov_Psi.set(2).set_spectral_va().set_base(sauve_base2) ;
    dcov_Psi.set(3).set_spectral_va().set_base(sauve_base3) ;
    dcov_Psi.inc_dzpuis(2) ;

    dcov_Psi += Psi_auto.derive_cov(flat) ;
    
    // dcov_chi
    temp = ((comp.chi_auto).derive_cov(comp.flat)) ;
    temp.dec_dzpuis(2) ;
    
    temp.change_triad(temp.get_mp().get_bvect_cart()) ;
    temp.change_triad(mp.get_bvect_cart()) ;
    sauve_base1 = (temp(1).get_spectral_va().get_base()) ;
    sauve_base2 = (temp(2).get_spectral_va().get_base()) ;
    sauve_base3 = (temp(3).get_spectral_va().get_base()) ;
    assert ( *(temp.get_triad()) == *(dcov_chi.get_triad())) ;

    for(int i=1; i<=3; i++) dcov_chi.set(i).import(temp(i)) ;

    dcov_chi.set(1).set_spectral_va().set_base(sauve_base1) ;
    dcov_chi.set(2).set_spectral_va().set_base(sauve_base2) ;
    dcov_chi.set(3).set_spectral_va().set_base(sauve_base3) ;
    dcov_chi.inc_dzpuis(2) ;

    dcov_chi += chi_auto.derive_cov(flat) ;
  
  // Computation of \hat{A}^{ij}_{comp}
  // ----------------------------------

  // D^j beta^i    
  const Tensor& dbeta_comp = beta_comp.derive_con(flat) ;
  
  /*
  Tensor temp_beta = (comp.beta_auto).derive_cov(comp.flat) ;
  temp_beta.dec_dzpuis(2) ;
  Tensor dbeta_comp(mp, CON, CON, mp.get_bvect_cart()) ;  
   
    temp_beta.change_triad(mp.get_bvect_cart()) ;
    Base_val sauve_base11 = (temp_beta(1,1).get_spectral_va().get_base()) ;
    Base_val sauve_base12 = (temp_beta(1,2).get_spectral_va().get_base()) ;
    Base_val sauve_base13 = (temp_beta(1,3).get_spectral_va().get_base()) ;
    Base_val sauve_base21 = (temp_beta(2,1).get_spectral_va().get_base()) ;
    Base_val sauve_base22 = (temp_beta(2,2).get_spectral_va().get_base()) ;
    Base_val sauve_base23 = (temp_beta(2,3).get_spectral_va().get_base()) ;
    Base_val sauve_base31 = (temp_beta(3,1).get_spectral_va().get_base()) ;
    Base_val sauve_base32 = (temp_beta(3,2).get_spectral_va().get_base()) ;
    Base_val sauve_base33 = (temp_beta(3,3).get_spectral_va().get_base()) ;    
    
    assert ( *(temp_beta.get_triad()) == *(dbeta_comp.get_triad())) ;
	
	dbeta_comp = temp_beta ; 
	
	for(int i=1; i<=3; i++) 
	for(int j=1; j<=3; j++) dbeta_comp.set(i,j).import(temp_beta(i,j)) ;

    dbeta_comp.set(1,1).set_spectral_va().set_base(sauve_base11) ;
    dbeta_comp.set(1,2).set_spectral_va().set_base(sauve_base12) ;
    dbeta_comp.set(1,3).set_spectral_va().set_base(sauve_base13) ;
    dbeta_comp.set(2,1).set_spectral_va().set_base(sauve_base21) ;
    dbeta_comp.set(2,2).set_spectral_va().set_base(sauve_base22) ;
    dbeta_comp.set(2,3).set_spectral_va().set_base(sauve_base23) ;
    dbeta_comp.set(3,1).set_spectral_va().set_base(sauve_base31) ;
    dbeta_comp.set(3,2).set_spectral_va().set_base(sauve_base32) ;
    dbeta_comp.set(3,3).set_spectral_va().set_base(sauve_base33) ;    
    
    dbeta_comp.inc_dzpuis(2) ;
    */
                          
  // Trace of D_j beta^i  :
  Scalar divbeta_comp = beta_comp.divergence(flat) ;
		  
  for (int i=1; i<=3; i++) 
    for (int j=i; j<=3; j++) {
      
      haij_comp.set(i, j) = dbeta_comp(i, j) + dbeta_comp(j, i) - 
	double(2) /double(3) * divbeta_comp * (flat.con())(i,j) ; 
    }

  haij_comp = 0.5 * pow(Psi, 7.) * haij_comp / chi ;   
  
  // Computation of (\hat{A}_{ij}\hat{A}^{ij})_{comp}
  // ------------------------------------------------
  
  Sym_tensor haij_auto_cov = haij_auto.up_down(flat) ;
  
  hacar_comp = contract(haij_auto_cov, 0, 1, haij_comp, 0, 1, true) ; 
  
  // The derived quantities are obsolete
  // -----------------------------------
  
  del_deriv() ;
  
}      

