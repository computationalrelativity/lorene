/*
 * Methods Et_bin_ncp::update_metric_der_comp
 *
 * (see file et_bin_ncp.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Francois Limousin
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


char et_bin_ncp_upmetrder_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// Headers Lorene
#include "et_bin_ncp.h"
#include "utilitaires.h"
#include "graphique.h"

void Et_bin_ncp::update_metric_der_comp(const Et_bin_ncp& comp) {
  
  
  // Derivatives of metric coefficients
  // ----------------------------------
  
  int nz = mp.get_mg()->get_nzone() ;

  dcov_logn_auto = logn_auto.derive_cov(gtilde) ;
  dcov_logn = (logn_auto + logn_comp).derive_cov(gtilde) ;
  dcon_logn = (logn_auto + logn_comp).derive_con(gtilde) ;
  
  // dcondcov_nnn
  
  Tenseur dcov_nnn (nnn.derive_cov(gtilde)) ;
  dcov_nnn.dec2_dzpuis() ;
  dcondcov_nnn = contract(dcov_nnn.derive_con(gtilde), 0, 1) ;
  dcondcov_nnn.inc2_dzpuis() ;
  
  dcondcov_nnn.annule(nz-1) ;
  
  // dcovdcov_logn_auto
  
  Tenseur dcov_logn_auto2 = dcov_logn_auto ;
  dcov_logn_auto2.dec2_dzpuis() ;
  dcovdcov_logn_auto = dcov_logn_auto2.derive_cov(gtilde) ;
  dcovdcov_logn_auto.inc2_dzpuis() ;
 
  dcovdcov_logn_auto.annule(nz-1) ;
 
  // dcondcov_logn_auto and laplacien
  
  dcondcov_logn_auto = contract(dcov_logn_auto2.derive_con(gtilde), 0, 1) ;
  dcondcov_logn_auto.inc2_dzpuis() ;
   
  lap_logn_auto = logn_auto().laplacien() ;
 
  dcondcov_logn_auto.annule(nz-1) ;
  lap_logn_auto.annule(nz-1) ;

  cout << "dcondcov_logn_auto" << endl << norme(dcondcov_logn_auto()) << endl ;
  cout << "lap_logn_auto" << endl << norme(lap_logn_auto()) << endl ;
  

  dcov_acar_auto = a_car_auto.derive_cov(gtilde) ;
  dcon_acar_auto = a_car_auto.derive_con(gtilde) ;
  
  dcov_acar = (a_car_auto + a_car_comp).derive_cov(gtilde) ;
  dcon_acar = (a_car_auto + a_car_comp).derive_con(gtilde) ;
  
  dcov_acar.annule(nz-1) ;

/*
  cout << "dcov_acar avant" << endl << norme(dcov_acar(0)) << endl ;
  cout << "dcon_acar avant" << endl << norme(dcon_acar(0)) << endl ;

  dcov_acar.dec2_dzpuis() ;
  dcov_acar.inc2_dzpuis() ;
  dcon_acar.dec2_dzpuis() ;
  dcon_acar.inc2_dzpuis() ;
  cout << "dcov_acar apres" << endl << norme(dcov_acar(0)) << endl ;
  cout << "dcon_acar apres" << endl << norme(dcon_acar(0)) << endl ;
  */

  // dcovdcov_acar_auto, dcondcov_acar_auto and laplacien
 
  Tenseur dcov_acar_auto2 = dcov_acar_auto ;
  dcov_acar_auto2.dec2_dzpuis() ;
 
  dcovdcov_acar_auto = dcov_acar_auto2.derive_cov(gtilde) ;
  dcondcov_acar_auto = contract(dcov_acar_auto2.derive_con(gtilde), 0, 1) ;
 
  dcovdcov_acar_auto.inc2_dzpuis() ;
  dcondcov_acar_auto.inc2_dzpuis() ;
  
  lap_acar_auto = a_car_auto().laplacien() ;
  
  dcovdcov_acar_auto.annule(nz-1) ; 
  dcondcov_acar_auto.annule(nz-1) ;
  lap_acar_auto.annule(nz-1) ;

  //  des_coef_xi(lap_acar_auto().va, 2, 0, 0); 
  //  des_coef_xi(dcondcov_acar_auto().va, 2, 0, 0); 

 
  cout << "dcondcov_acar_auto" << endl << norme(dcondcov_acar_auto()) << endl ;
  cout << "lap_acar_auto" << endl << norme(lap_acar_auto()) << endl ;
 
  /*
  // Derivatives of logn
  // -------------------

  dcov_logn_auto = logn_auto.derive_cov(gtilde) ;
 
  // dcov_logn

  Tenseur temp0 = (comp.logn_auto).derive_cov(comp.gtilde) ;
  temp0.dec2_dzpuis() ;
  dcov_logn.set_etat_qcq() ;
  (dcov_logn.set(0)).import(temp0(0)) ;
  (dcov_logn.set(1)).import(temp0(1)) ;
  (dcov_logn.set(2)).import(temp0(2)) ;
  dcov_logn.set_std_base() ;
  dcov_logn.inc2_dzpuis() ;
  // dcov_logn = dcov_logn + dcov_logn_auto ;

  cout << "dcov_logn_comp" << endl << norme(dcov_logn(0)) ;
  cout << "dcov_logn_comp" << endl << norme(logn_comp.derive_cov(gtilde)(0)) ;

  // dcon_logn

  Tenseur temp1 = (comp.logn_auto).derive_con(comp.gtilde) ;
  temp1.dec2_dzpuis() ;
  dcon_logn.set_etat_qcq() ;
  (dcon_logn.set(0)).import(temp1(0)) ;
  (dcon_logn.set(1)).import(temp1(1)) ;
  (dcon_logn.set(2)).import(temp1(2)) ;
  dcon_logn.set_std_base() ;
  dcon_logn.inc2_dzpuis() ;
  dcon_logn = dcon_logn + logn_auto.derive_con(gtilde) ;

  /* 
  cout << "dcon_logn avant" << endl << norme(dcon_logn(0)) << endl ;
  dcon_logn.dec2_dzpuis() ;
  dcon_logn.inc2_dzpuis() ;
  cout << "dcon_logn apres" << endl << norme(dcon_logn(0)) << endl ;
  
  Tenseur dcon_logn2 = (logn_auto + logn_comp).derive_con(gtilde) ;
  cout << "dcon_logn2 avant" << endl << norme(dcon_logn2(0)) << endl ;
  dcon_logn2.dec2_dzpuis() ;
  dcon_logn2.inc2_dzpuis() ;
  cout << "dcon_logn2 apres" << endl << norme(dcon_logn2(0)) << endl ;
  */
  /*
  // dcondcov_nnn

  Tenseur dcov_nnn (nnn.derive_cov(gtilde)) ;
  dcov_nnn.dec2_dzpuis() ;
  dcondcov_nnn = contract(dcov_nnn.derive_con(gtilde), 0, 1) ;
  dcondcov_nnn.inc2_dzpuis() ;
 
  // dcovdcov_logn_auto

  dcov_logn_auto.dec2_dzpuis() ;
  dcovdcov_logn_auto = dcov_logn_auto.derive_cov(gtilde) ;
  dcov_logn_auto.inc2_dzpuis() ;
  dcovdcov_logn_auto.inc2_dzpuis() ;

  // dcondcov_logn_auto and laplacien

  dcov_logn_auto.dec2_dzpuis() ;
  dcondcov_logn_auto = contract(dcov_logn_auto.derive_con(gtilde), 0, 1) ;
  dcov_logn_auto.inc2_dzpuis() ;
  dcondcov_logn_auto.inc2_dzpuis() ;

  lap_logn_auto = logn_auto().laplacien() ;
  
  //cout << "dcondcov_logn_auto" << endl << norme(dcondcov_logn_auto()) << endl ;
  //cout << "lap_logn_auto" << endl << norme(lap_logn_auto()) << endl ;


  // Derivatives of a_car
  // --------------------

  dcov_acar_auto = a_car_auto.derive_cov(gtilde) ;
  dcon_acar_auto = a_car_auto.derive_con(gtilde) ;

  // dcov_acar

  Tenseur temp2 = (comp.a_car_auto).derive_cov(comp.gtilde) ;
  temp2.dec2_dzpuis() ;
  dcov_acar.set_etat_qcq() ;
  (dcov_acar.set(0)).import(temp2(0)) ;
  (dcov_acar.set(1)).import(temp2(1)) ;
  (dcov_acar.set(2)).import(temp2(2)) ;
  dcov_acar.set_std_base() ;
  dcov_acar.inc2_dzpuis() ;
  dcov_acar = dcov_acar + dcov_acar_auto ;

  // dcon_acar

  Tenseur temp3 = (comp.a_car_auto).derive_con(comp.gtilde) ;
  temp3.dec2_dzpuis() ;
  dcon_acar.set_etat_qcq() ;
  (dcon_acar.set(0)).import(temp3(0)) ;
  (dcon_acar.set(1)).import(temp3(1)) ;
  (dcon_acar.set(2)).import(temp3(2)) ;
  dcon_acar.set_std_base() ;
  dcon_acar.inc2_dzpuis() ;
  dcon_acar = dcon_acar + dcon_acar_auto ;
 
  // dcondcov_acar

  Tenseur temp4 = contract(temp2.derive_con(comp.gtilde), 0, 1) ;
  temp4.dec2_dzpuis() ;
  dcondcov_acar.set_etat_qcq() ;
  (dcondcov_acar.set()).import(temp4()) ;
  dcondcov_acar.set_std_base() ;
  dcondcov_acar.inc2_dzpuis() ;
  dcondcov_acar.inc2_dzpuis() ;
  Tenseur temp4b = dcov_acar_auto ;
  temp4b.dec2_dzpuis() ;
  Tenseur dcondcov_acar2 = contract(temp4b.derive_con(gtilde), 0, 1) ;
  dcondcov_acar2.inc2_dzpuis() ;
  dcondcov_acar = dcondcov_acar + dcondcov_acar2 ;
  
  // dcovdcov_acar_auto, dcondcov_acar_auto and laplacien

  dcov_acar_auto.dec2_dzpuis() ;
  dcovdcov_acar_auto = dcov_acar_auto.derive_cov(gtilde) ;
  dcondcov_acar_auto = contract(dcov_acar_auto.derive_con(gtilde), 0, 1) ;
  dcov_acar_auto.inc2_dzpuis() ;
  dcovdcov_acar_auto.inc2_dzpuis() ;
  dcondcov_acar_auto.inc2_dzpuis() ;
  
  lap_acar_auto = a_car_auto().laplacien() ;
  
*/




  // Derivatives of shift
  // --------------------

  // dcondcovj_shiftj
  
  Tenseur dshift_auto(contract(shift_auto.derive_cov(gtilde), 0, 1)) ;
  dshift_auto.dec2_dzpuis() ;

  Tenseur dcondcovj_shiftj_auto(dshift_auto.derive_con(gtilde)) ;
  dcondcovj_shiftj_auto.inc2_dzpuis() ;

  // diffdidj_shift_autoj
  
  Tenseur dflat_shiftauto(contract(shift_auto.derive_cov(flat), 0, 1)) ;
  dflat_shiftauto.dec2_dzpuis() ;

  Tenseur ddflat_shiftauto(dflat_shiftauto.derive_con(flat)) ;
  ddflat_shiftauto.inc2_dzpuis() ;
  
  diffdidj_shift_autoj = ddflat_shiftauto - dcondcovj_shiftj_auto ;

  cout << "diffshift" << endl << norme(diffdidj_shift_autoj(0)) << endl ;
  
  // dcovdcon_shift_auto and laplacien
  
  Tenseur dcon_shift_auto (shift_auto.derive_con(gtilde)) ;
  dcon_shift_auto.dec2_dzpuis() ;
  dcovdcon_shift_auto = contract(dcon_shift_auto.derive_cov(gtilde), 0, 1) ;
  dcovdcon_shift_auto.inc2_dzpuis() ;


  lap_shift_auto.set_etat_qcq() ;
  for(int i=0; i<=2; i++) {
    lap_shift_auto.set(i) = shift_auto(i).laplacien() ;
  }

  for(int i=0; i<=2; i++) {
    lap_shift_auto.set(i).annule(nz-1) ;
    dcovdcon_shift_auto.set(i).annule(nz-1) ;
    diffdidj_shift_autoj.set(i).annule(nz-1) ;
  }


  lap_gtilde_auto.set_etat_qcq() ;
  for(int i=0; i<=2; i++) {
    for(int j=i; j<=2; j++) {
      lap_gtilde_auto.set(i,j) = ((gtilde_auto.cov())(i,j)).laplacien() ;
      lap_gtilde_auto.set(i,j).annule(nz-1) ;
    }
  }
  /*
  des_profile(a_car_auto(), 0, 20, 0, 0) ;
  des_coef_xi(a_car_auto().va, 0, 0, 0) ;
  des_profile(dcov_acar_auto(0), 0, 20, 0, 0) ;
  des_coef_xi(dcov_acar_auto(0).va, 0, 0, 0) ;
  des_profile(dcovdcov_acar_auto(0,0), 0, 20, 0, 0) ;
  des_coef_xi(dcovdcov_acar_auto(0,0).va, 0, 0, 0) ;
  des_profile(logn_auto(), 0, 20, 0, 0) ;
  des_coef_xi(logn_auto().va, 0, 0, 0) ;
  des_profile(dcovdcov_logn_auto(0,0), 0, 20, 0, 0) ;
  des_coef_xi(dcovdcov_logn_auto(0,0).va, 0, 0, 0) ;
  */


  /* 
  Cmp source_tot = + pow(gamma(),-1./3.)*(dcov_acar(2)*dcov_logn_auto(2) 
					  + dcov_acar(2)*dcov_logn_auto(2)) ;
  source_tot.std_base_scal() ;
  des_profile(source_tot, 0, 20, 0, 0) ;
  des_coef_xi(source_tot.va, 0, 0, 0) ;
  des_coef_xi(source_tot.va, 1, 0, 0) ;
  des_coef_xi(source_tot.va, 2, 0, 0) ;
  */
  


  // Computation of tkij_comp
  // ------------------------
  
    if ( (comp.tkij_auto).get_etat() == ETATZERO ) {
	tkij_comp.set_etat_zero() ;
    }
    else{

      // Components of shift_comp with respect to the Cartesian triad
      //  (d/dx, d/dy, d/dz) of the mapping :
      Tenseur shift_comp_local = shift_comp ;
 
      // Gradient tilde (partial derivatives with respect to
      //           the Cartesian coordinates of the mapping)
      // D~_j beta^i

      Tenseur dshift_comp = shift_comp_local.derive_con(gtilde) ;

      // Trace of D~_j beta^i  :
      Tenseur divshift_comp = contract(shift_comp_local.derive_cov(gtilde), 0, 1) ;

      // Computation of K^{ij}
      // -------------------------
      tkij_comp.set_etat_qcq() ;

      for (int i=0; i<3; i++) {
	for (int j=i; j<3; j++) {
	  tkij_comp.set(i, j) = dshift_comp(i, j) + dshift_comp(j, i) - 
	    double(2) /double(3) * divshift_comp() * (gtilde.con())(i,j) ; 
	  
	}
      }
      tkij_comp = - 0.5 * tkij_comp / nnn ;
      tkij_comp.set_std_base() ;
      
      // Computation of kcar_comp
      // -------------------------
      
      kcar_comp.set_etat_qcq() ;
      
      kcar_comp.set() = 0 ;
      
      for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {

	  Tenseur tkij_auto_cov = tkij_auto ;
	  tkij_auto_cov.set(i,j) = contract(contract(gtilde.cov() * 				     gtilde.cov() * tkij_auto, 1, 4), 2, 3)(i,j) ;
	  
	  kcar_comp.set() += tkij_auto_cov(i,j) % tkij_comp(i,j) ; 
	}
      }
      
      kcar_comp.set() =  a_car() % kcar_comp() ; 
      kcar_comp.set_std_base() ;

      
      // Computation of tkij_auto_cov, tkij_comp_cov and kcar_cov
      // --------------------------------------------------------
      
      
      Tenseur tkij_auto_cov (contract(gtilde.cov() * contract(gtilde.cov() 
      				      * tkij_auto, 1, 3), 1, 3)) ;
      
      Tenseur tkij_comp_cov (contract(gtilde.cov() * contract(gtilde.cov() 
	       	        * tkij_comp, 1, 3), 1, 3)) ;
      
      tkij_auto_cov.set_std_base() ;
      tkij_comp_cov.set_std_base() ;
      
      kcar_cov.set_etat_qcq() ;
      kcar_cov = pow(gamma, -2./3.) * contract( contract(gtilde.con() 
	  * tkij_auto_cov, 0, 3) * (tkij_auto_cov + tkij_comp_cov), 0, 3) ;
      
      kcar_cov.set_std_base() ;
      
      
      // The derived quantities are obsolete
      // -----------------------------------
      
      del_deriv() ;

      
    }
}
