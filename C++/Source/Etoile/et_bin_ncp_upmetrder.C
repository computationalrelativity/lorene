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
  
    cout << "update metric derivatives 1" << endl ;

  int nr = mp.get_mg()->get_nr(0) ;
  int nt = mp.get_mg()->get_nt(0) ;
  int np = mp.get_mg()->get_np(0) ;
  
  int nz = mp.get_mg()->get_nzone() ;

  dcov_logn_auto = logn_auto.derive_cov(flat) ;
  dcov_logn = (logn_auto + logn_comp).derive_cov(flat) ;
  dcon_logn = (logn_auto + logn_comp).derive_con(flat) ;
    
  // dcovdcov_logn_auto
  
  Tenseur dcov_logn_auto2 = logn_auto.derive_cov(flat) ;
  dcov_logn_auto2.dec2_dzpuis() ;
  dcovdcov_logn_auto = dcov_logn_auto2.derive_cov(flat) ;
  dcovdcov_logn_auto.inc2_dzpuis() ;
 
    cout << "update metric derivatives 2" << endl ;


  // dcov_acar...
 
  dcov_acar_auto = a_car_auto.derive_cov(flat) ;
  dcov_acar = (a_car_auto + a_car_comp).derive_cov(flat) ;
 

/*
  cout << "dcov_acar avant" << endl << norme(dcov_acar(0)) << endl ;

  dcov_acar.dec2_dzpuis() ;
  dcov_acar.inc2_dzpuis() ;
  cout << "dcov_acar apres" << endl << norme(dcov_acar(0)) << endl ;
*/

  // dcovdcov_acar_auto

  Tenseur dcov_acar_auto2 = a_car_auto.derive_cov(flat) ;
  dcov_acar_auto2.dec2_dzpuis() ;
  dcovdcov_acar_auto = dcov_acar_auto2.derive_cov(flat) ;
  dcovdcov_acar_auto.inc2_dzpuis() ;


    cout << "update metric derivatives 3" << endl ;


  // Derivatives of gtilde
  // ---------------------

  const Tenseur& dcov_gtilde_auto = (gtilde_auto.set_cov()).derive_cov(flat) ;
  const Tenseur& dcov_gtilde = (gtilde.set_cov()).derive_cov(flat) ;
  const Tenseur& dcov_gtilde_con = (gtilde.set_con()).derive_cov(flat) ;
  
  Tenseur dcov_gtilde_auto2 = (gtilde_auto.set_cov()).derive_cov(flat) ;
  dcov_gtilde_auto2.dec2_dzpuis() ; 
  Tenseur dcovdcov_gtilde_auto = dcov_gtilde_auto2.derive_cov(flat) ;
  dcovdcov_gtilde_auto.inc2_dzpuis() ;

  // Derivatives of hij
  // ------------------

  const Tenseur& dcov_hij = hij.derive_cov(flat) ;

    cout << "update metric derivatives 4" << endl ;


  // deltakij and deltakij_auto

  deltakij_auto.set_etat_qcq() ;
  deltakij.set_etat_qcq() ;

  const Tenseur& gtilde_con = gtilde.con() ; 
  Tenseur temp1 = contract(gtilde_con,1, dcov_gtilde_auto,1) ;
  Tenseur temp2 = contract(gtilde_con,1, dcov_gtilde_auto,2) ;
  Tenseur temp3 = contract(gtilde_con,1, dcov_gtilde_auto,0) ;
  
  Tenseur temp4 = contract(gtilde_con,1, dcov_gtilde,1) ;
  Tenseur temp5 = contract(gtilde_con,1, dcov_gtilde,2) ;
  Tenseur temp6 = contract(gtilde_con,1, dcov_gtilde,0) ;
 
    
  for (int k=0; k<3; k++)
    for (int i=0; i<3; i++)
      for (int j=i; j<3; j++){
   
	deltakij_auto.set(k,i,j) = ( temp1(k,i,j) + temp2(k,j,i) 
				     - temp3(k,i,j) ) * 0.5 ;

	deltakij.set(k,i,j) = ( temp4(k,i,j) + temp5(k,j,i) 
				     - temp6(k,i,j) ) * 0.5 ;

      }

  for (int k=0; k<3; k++)
    for (int i=1; i<3; i++)
      for (int j=0; j<i; j++){
  
	  deltakij_auto.set(k,i,j) = deltakij_auto(k,j,i) ;
	  deltakij.set(k,i,j) = deltakij(k,j,i) ;

      }

  temp1 = 0 ;
  temp2 = 0 ;
  temp3 = 0 ;
  temp4 = 0 ;
  temp5 = 0 ;
  temp6 = 0 ;


    cout << "update metric derivatives 5" << endl ;


  // Ricci
  
    ricci_auto.set_etat_qcq() ;

    Tenseur temp7 = contract(contract(gtilde.con(),0,
				    dcovdcov_gtilde_auto,0), 0, 1) ;
    Tenseur temp8 = contract(contract(dcov_gtilde_con,1,
				    dcov_gtilde_auto,0), 1, 2) ;
    Tenseur temp9 = contract(contract(dcov_gtilde_con,1,
				    dcov_gtilde_auto,0), 1, 3) ;
    Tenseur temp10 = contract(contract(deltakij,0, deltakij,2), 1, 2) ;
  

  for (int i=0; i<3; i++)
    for (int j=i; j<3; j++){  
      ricci_auto.set(i,j) = - ( temp7(i,j) + temp8(i,j) + temp9(j,i) ) * 0.5 
    - temp10(i,j) ;
    }

  temp7 = 0 ;
  temp8 = 0 ;
  temp9 = 0 ;
  temp10 = 0 ;

    cout << "update metric derivatives 6" << endl ;

  // Ricci_scal

  ricci_scal = contract(contract(contract(contract(gtilde.con(),0, dcov_hij
       ,0),0, dcov_gtilde,0), 0, 2), 0, 1)*0.25 
       - contract(contract(contract(contract(gtilde.con(),0, dcov_hij
       ,0),0, dcov_gtilde,2), 0, 3), 0, 1)*0.5 ;

    cout << "update metric derivatives 7" << endl ;


  /*
  des_profile(a_car_auto(), 0, 20, 0, 0) ;
  des_coef_xi(a_car_auto().va, 0, 0, 0) ;
  des_profile(dcov_acar_auto(0), 0, 20, 0, 0) ;
  des_coef_xi(dcov_acar_auto(0).va, 0, 0, 0) ;
  des_profile(dcovdcov_acar_auto(0,0), 0, 20, 0, 0) ;
  des_coef_xi(dcovdcov_acar_auto(0,0).va, 0, 0, 0) ;
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

      const Tenseur& dshift_comp = shift_comp_local.derive_con(gtilde) ;

      // Trace of D~_j beta^i  :
      Tenseur divshift_comp = contract(shift_comp_local.derive_cov(gtilde), 0, 1) ;

      // Computation of K^{ij}
      // -------------------------

      tkij_comp.set_etat_qcq() ;
  
      for (int i=0; i<3; i++) {
	for (int j=i; j<3; j++) {

	  tkij_comp.set(i, j) = dshift_comp(i, j) + dshift_comp(j, i) - 
	    double(2) /double(3) * divshift_comp() % (gtilde.con())(i,j) ; 	  
	}
      }

      tkij_comp = - 0.5 * tkij_comp / nnn ;
      tkij_comp.set_std_base() ;
      
      // Computation of kcar_comp
      // -------------------------
      
      kcar_comp.set_etat_qcq() ;     
      kcar_comp.set() = 0 ;
      
      Tenseur temp11 = contract(gtilde.cov(),1, contract( gtilde.cov(),1, 
				             tkij_auto,1),1) ;

      Tenseur tkij_auto_cov = tkij_auto ;

      for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	  
	  tkij_auto_cov.set(i,j) = temp11(i,j) ;
	  
	  kcar_comp.set() += tkij_auto_cov(i,j) % tkij_comp(i,j) ; 
	}
      }
      
      kcar_comp.set() =  a_car() % kcar_comp() ; 
      kcar_comp.set_std_base() ;

      
      // Computation of kcar_cov
      // -----------------------
            
      kcar_con.set_etat_qcq() ;
      kcar_con = contract( contract(gtilde.cov(),0,
	  tkij_auto,1),0, tkij_auto + tkij_comp,1) ;
      
      kcar_con.set_std_base() ;
      
      
      // The derived quantities are obsolete
      // -----------------------------------
      
      del_deriv() ;

          cout << "update metric derivatives 8" << endl ;

    }
}
