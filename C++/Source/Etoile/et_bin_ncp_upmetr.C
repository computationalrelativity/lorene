/*
 * Methods of Et_bin_ncp::update_metric
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


char et_bin_ncp_upmetr_C[] = "$Header$" ;

/*
 * $Header$ *
 */

// Headers Lorene
#include "et_bin_ncp.h"
#include "graphique.h"

//----------------------------------//
//	 Version without relaxation //
//----------------------------------//

void Et_bin_ncp::update_metric(const Et_bin_ncp& comp) {

  // Computation of quantities coming from the companion
  // ---------------------------------------------------

  if ( (comp.logn_auto).get_etat() == ETATZERO ) {
    logn_comp.set_etat_zero() ;
  }
  else{
    logn_comp.set_etat_qcq() ;
    (logn_comp.set()).import_symy( comp.logn_auto() ) ;
    logn_comp.set_std_base() ;   // set the bases for spectral expansions
  }


  if ( (comp.shift_auto).get_etat() == ETATZERO ) {
    shift_comp.set_etat_zero() ; 
  }
  else{  
    shift_comp.set_etat_qcq() ; 

    shift_comp.set_triad( *((comp.shift_auto).get_triad()) ) ;  
    assert ( *(shift_comp.get_triad()) == *((comp.shift_auto).get_triad())) ;
	
    (shift_comp.set(0)).import_asymy( comp.shift_auto(0) ) ;  // N^x antisym
    (shift_comp.set(1)).import_symy( comp.shift_auto(1) ) ;   // N^y sym.
    (shift_comp.set(2)).import_asymy( comp.shift_auto(2) ) ;  // N^z anisym


    shift_comp.change_triad( mp.get_bvect_cart() ) ;    	
    shift_comp.set_std_base() ;   // set the bases for spectral expansions
  }
    


  if ( (comp.a_car_auto).get_etat()  == ETATZERO ) {
    a_car_comp.set_etat_zero() ;
  }
  else{
    a_car_comp.set_etat_qcq() ;
  
    (a_car_comp.set()).import_symy( comp.a_car_auto() ) ;
    a_car_comp.set_std_base() ;   // set the bases for spectral expansions
  }	


  if (!conf_flat){
      
    (gtilde_comp.set_cov()).set_triad( *(((comp.gtilde_auto).cov()).get_triad()) ) ;  

    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
	
	if ( ((comp.gtilde_auto).cov())(i,j).get_etat()  == ETATZERO ) {
	  (gtilde_comp.set_cov(i,j)).set_etat_zero() ;
	}
	else{
	  (gtilde_comp.set_cov(i,j)).set_etat_qcq() ;
	  (gtilde_comp.set_cov(i,j)).import( ((comp.gtilde_auto).cov())(i,j) ) ;
	}
      }
    }
    
    
    cout << "gtilde_auto" << endl << norme(gtilde_auto.cov()(0,0)) << endl ;
    cout << "gtilde_comp" << endl << norme(gtilde_comp.cov()(0,0)) << endl ;
    
    (gtilde_comp.set_cov()).change_triad( mp.get_bvect_cart() ) ;
    gtilde_comp.set_std_base() ;   // set the bases for spectral expansions
    
  }
  
  // Lapse function N
  // ----------------

  Tenseur logn_total = logn_auto + logn_comp ; 

  nnn = exp( unsurc2 * logn_total ) ;

  nnn.set_std_base() ;   // set the bases for spectral expansions



  // Shift vector beta^i
  // -------------------

  shift = shift_auto + shift_comp ;

 
  // Conformal factor A^2
  // ----------------------------------

  a_car = a_car_auto + a_car_comp ;

  gamma = pow(a_car, 3.) ;
  gamma.set_poids(2.) ;
  gamma.set_std_base() ;


  if (!conf_flat){
      
    // Coefficients of the 3-metric tilde
    // ----------------------------------
    cout << "gtilde" << endl << norme(gtilde.cov()(0,0)) << endl << norme(gtilde.cov()(0,1)) << endl << norme(gtilde.cov()(1,1)) << endl ;
     
    assert ( *(gtilde_auto.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
      
    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
   
	gtilde.set_cov(i,j) = (gtilde_auto.cov())(i,j) + (gtilde_comp.cov())(i,j) ;
      }
    }
    cout << "gtilde" << endl << norme(gtilde.cov()(0,0)) << endl << norme(gtilde.cov()(0,1)) << endl << norme(gtilde.cov()(1,1)) << endl ;


    gtilde.set_std_base() ;
      
    // Determinant of gtilde

    Tenseur det_gtilde(mp) ;
    det_gtilde = gtilde.cov()(0, 0)*gtilde.cov()(1, 1)*gtilde.cov()(2, 2) 
      + gtilde.cov()(0, 1)*gtilde.cov()(1, 2)*gtilde.cov()(2, 0)
      + gtilde.cov()(0, 2)*gtilde.cov()(1, 0)*gtilde.cov()(2, 1) 
      - gtilde.cov()(2, 0)*gtilde.cov()(1, 1)*gtilde.cov()(0, 2)
      - gtilde.cov()(2, 1)*gtilde.cov()(1, 2)*gtilde.cov()(0, 0) 
      - gtilde.cov()(2, 2)*gtilde.cov()(1, 0)*gtilde.cov()(0, 1) ;

    int nzone = mp.get_mg()->get_nzone() ;
    int nr = mp.get_mg()->get_nr(0);
    int nt = mp.get_mg()->get_nt(0);
    int np = mp.get_mg()->get_np(0);
      
    double max_det[nzone] ;
    double min_det[nzone] ;
    double moy_det[nzone] ;
      
    for (int i=0; i<nzone; i++){
      min_det[i] = 2 ;
      moy_det[i] = 0 ;
      max_det[i] = 0 ;
    }

    for (int l=0; l<nzone; l++)
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
	  for (int i=0; i<nr; i++){
	      
	    moy_det[l] = moy_det[l] + det_gtilde()(l,k,j,i) ;
	    if (det_gtilde()(l,k,j,i) > max_det[l]){
	      max_det[l] = det_gtilde()(l,k,j,i) ;
	    }
	    if (det_gtilde()(l,k,j,i) < min_det[l]){
	      min_det[l] = det_gtilde()(l,k,j,i) ;
	    }
	  }
     
    cout << "average determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nzone; l++){
      cout << moy_det[l]/(nr*nt*np) << "  " ;
    }
    cout << endl ;

      
    cout << "maximum of the determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nzone; l++){
      cout << max_det[l] << "  " ;
    }
    cout << endl ;

    cout << "minimum of the determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nzone; l++){
      cout << min_det[l] << "  " ;
    }
    cout << endl ;
      

       
    // Coefficients of the 3-metric
    // ----------------------------
      
    assert ( *(met_gamma.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    assert ( *(metgamma_auto.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    assert ( *(metgamma_comp.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
      
      
    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
	  
	met_gamma.set_cov(i,j) = a_car() * gtilde.cov()(i,j) ;
	metgamma_auto.set_cov(i,j) = a_car() * gtilde_auto.cov()(i,j) ;
	metgamma_comp.set_cov(i,j) = a_car() * gtilde_comp.cov()(i,j) ;
      }
    }  
    met_gamma.set_std_base() ;
    metgamma_auto.set_std_base() ;
    metgamma_comp.set_std_base() ;
      
  }
  else {
      
    // for conformally flat metric
      
    for (int i=0; i<3; i++) {
      gtilde_auto.set_cov(i,i) = a_car_auto() / a_car()  ; 
      gtilde_comp.set_cov(i,i) = a_car_comp() / a_car()  ; 
      met_gamma.set_cov(i,i) = a_car() ;
      metgamma_auto.set_cov(i,i) = a_car_auto() ;
      metgamma_comp.set_cov(i,i) = a_car_comp() ;
    }
      
    gtilde_auto.set_std_base() ;
    gtilde_comp.set_std_base() ;
    met_gamma.set_std_base() ;
    metgamma_auto.set_std_base() ;
    metgamma_comp.set_std_base() ;
      
  }
    
    
  if (relativistic) {
    // ... extrinsic curvature (tkij_auto and kcar_auto)
    extrinsic_curvature() ;
  }
   
  // The derived quantities are obsolete
  // -----------------------------------

  del_deriv() ;


}



//----------------------------------//
//	  Version with relaxation   //
//----------------------------------//

void Et_bin_ncp::update_metric(const Et_bin_ncp& comp,
			       const Et_bin_ncp& star_jm1, double relax) {


  // Computation of quantities coming from the companion
  // ---------------------------------------------------

  if ( (comp.logn_auto).get_etat() == ETATZERO ) {
    logn_comp.set_etat_zero() ;
  }
  else{
    logn_comp.set_etat_qcq() ;
    (logn_comp.set()).import_symy( comp.logn_auto() ) ;
    logn_comp.set_std_base() ;   // set the bases for spectral expansions
  }


  if ( (comp.shift_auto).get_etat() == ETATZERO ) {
    shift_comp.set_etat_zero() ; 
  }
  else{  
    shift_comp.set_etat_qcq() ; 

    shift_comp.set_triad( *((comp.shift_auto).get_triad()) ) ;  
    assert ( *(shift_comp.get_triad()) == *((comp.shift_auto).get_triad())) ;

    (shift_comp.set(0)).import_asymy( comp.shift_auto(0) ) ;  // N^x antisym
    (shift_comp.set(1)).import_symy( comp.shift_auto(1) ) ;   // N^y sym.
    (shift_comp.set(2)).import_asymy( comp.shift_auto(2) ) ;  // N^z anisym

    shift_comp.change_triad( mp.get_bvect_cart() ) ;    	
    shift_comp.set_std_base() ;   // set the bases for spectral expansions
  }


  if ( (comp.a_car_auto).get_etat()  == ETATZERO ) {
    a_car_comp.set_etat_zero() ;
  }
  else{
    a_car_comp.set_etat_qcq() ;

    (a_car_comp.set()).import_symy( comp.a_car_auto() ) ;
    a_car_comp.set_std_base() ;   // set the bases for spectral expansions
  }	

    cout << "gtilde_auto" << endl << norme(gtilde_auto.cov()(0,0)) << endl ;
    cout << "gtilde_comp" << endl << norme((comp.gtilde_auto).cov()(0,0)) << endl ;
	
     
  if (!conf_flat){
    
    (gtilde_comp.set_cov()).set_triad( *(((comp.gtilde_auto).cov()).get_triad()) ) ;  
      
    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
	
	if ( ((comp.gtilde_auto).cov())(i,j).get_etat()  == ETATZERO ) {
	  (gtilde_comp.set_cov(i,j)).set_etat_zero() ;
	}
	else{
	  (gtilde_comp.set_cov(i,j)).set_etat_qcq() ;
	  (gtilde_comp.set_cov(i,j)).import( ((comp.gtilde_auto).cov())(i,j) ) ; 	  
	  gtilde_comp.set_std_base() ; // set the bases for spectral expansions
	}   
      }	
    }
    
    (gtilde_comp.set_cov()).change_triad( mp.get_bvect_cart() ) ;
    
  }
    
  // Relaxation on logn_comp, shift_comp, loggamma_comp, gtilde_comp
  // ---------------------------------------------------------------
  double relaxjm1 = 1. - relax ; 
    
  logn_comp = relax * logn_comp + relaxjm1 * (star_jm1.get_logn_comp()) ; 
    
  shift_comp = relax * shift_comp + relaxjm1 * (star_jm1.get_shift_comp()) ; 

  a_car_comp = relax * a_car_comp + relaxjm1 * (star_jm1.get_acar_comp()) ;

  if (!conf_flat){
      
    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
	gtilde_comp.set_cov(i,j) = relax * (gtilde_comp.cov())(i,j) + relaxjm1 * ((star_jm1.get_gtilde_comp()).cov())(i,j) ;
      }
    }
  }  
  
  // Lapse function N
  // ----------------
    
  Tenseur logn_total = logn_auto + logn_comp ; 
    
  nnn = exp( unsurc2 * logn_total ) ; 
    
  nnn.set_std_base() ;   // set the bases for spectral expansions
    
    
  // Shift vector N^i
  // ----------------
    
  shift = shift_auto + shift_comp ; 


  // Conformal factor A^2
  // ------------------

  a_car = a_car_auto + a_car_comp ;

  gamma = pow(a_car, 3.) ;
  gamma.set_poids(2.) ;
  gamma.set_std_base() ;
 
  if (!conf_flat){
      
    // Coefficients of the 3-metric tilde
    // ----------------------------------
      
    assert ( *(gtilde_auto.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
      
    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
	  
	gtilde.set_cov(i,j) = (gtilde_auto.cov())(i,j) 
	  + (gtilde_comp.cov())(i,j) ;
      }
    }
    gtilde.set_std_base() ;

    // Determinant of gtilde

    Tenseur det_gtilde(mp) ;
    det_gtilde = gtilde.cov()(0, 0)*gtilde.cov()(1, 1)*gtilde.cov()(2, 2) 
      + gtilde.cov()(0, 1)*gtilde.cov()(1, 2)*gtilde.cov()(2, 0)
      + gtilde.cov()(0, 2)*gtilde.cov()(1, 0)*gtilde.cov()(2, 1) 
      - gtilde.cov()(2, 0)*gtilde.cov()(1, 1)*gtilde.cov()(0, 2)
      - gtilde.cov()(2, 1)*gtilde.cov()(1, 2)*gtilde.cov()(0, 0) 
      - gtilde.cov()(2, 2)*gtilde.cov()(1, 0)*gtilde.cov()(0, 1) ;


    int nzone = mp.get_mg()->get_nzone() ;
    int nr = mp.get_mg()->get_nr(0);
    int nt = mp.get_mg()->get_nt(0);
    int np = mp.get_mg()->get_np(0);
      
    double max_det[nzone] ;
    double min_det[nzone] ;
    double moy_det[nzone] ;
      
    for (int i=0; i<nzone; i++){
      min_det[i] = 2 ;
      moy_det[i] = 0 ;
      max_det[i] = 0 ;
    }

    for (int l=0; l<nzone; l++)
      for (int k=0; k<np; k++)
	for (int j=0; j<nt; j++)
	  for (int i=0; i<nr; i++){
	      
	    moy_det[l] = moy_det[l] + det_gtilde()(l,k,j,i) ;
	    if (det_gtilde()(l,k,j,i) > max_det[l]){
	      max_det[l] = det_gtilde()(l,k,j,i) ;
	    }
	    if (det_gtilde()(l,k,j,i) < min_det[l]){
	      min_det[l] = det_gtilde()(l,k,j,i) ;
	    }
	  }
     
    cout << "average determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nzone; l++){
      cout << moy_det[l]/(nr*nt*np) << "  " ;
    }
    cout << endl ;

      
    cout << "maximum of the determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nzone; l++){
      cout << max_det[l] << "  " ;
    }
    cout << endl ;

    cout << "minimum of the determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nzone; l++){
      cout << min_det[l] << "  " ;
    }
    cout << endl << endl ;
      
      
    // Coefficients of the 3-metric
    // ----------------------------
      
    assert ( *(met_gamma.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    assert ( *(metgamma_auto.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    assert ( *(metgamma_comp.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
	  
	met_gamma.set_cov(i,j) = a_car() * gtilde.cov()(i,j) ;
	metgamma_auto.set_cov(i,j) = a_car() * gtilde_auto.cov()(i,j) ;
	metgamma_comp.set_cov(i,j) = a_car() * gtilde_comp.cov()(i,j) ;
      }
    }    
    metgamma_auto.set_std_base() ;
    metgamma_comp.set_std_base() ;
    met_gamma.set_std_base() ;
    
  }
  else {
      
    // for conformally flat metric
      
    for (int i=0; i<3; i++) {
      gtilde_auto.set_cov(i,i) = a_car_auto() / a_car()  ; 
      gtilde_comp.set_cov(i,i) = a_car_comp() / a_car()  ; 
      met_gamma.set_cov(i,i) = a_car() ;
      metgamma_auto.set_cov(i,i) = a_car_auto() ;
      metgamma_comp.set_cov(i,i) = a_car_comp() ;
    }
      
    gtilde_auto.set_std_base() ;
    gtilde_comp.set_std_base() ;
    met_gamma.set_std_base() ;
    metgamma_auto.set_std_base() ;
    metgamma_comp.set_std_base() ;
      
  }
      
  // ... extrinsic curvature (tkij_auto and kcar_auto)
  extrinsic_curvature() ; 
    
  // The derived quantities are obsolete
  // -----------------------------------
    
  del_deriv() ;                
    
			      
} 


void Et_bin_ncp::update_metric_init(const Et_bin_ncp& comp) {

  // Computation of quantities coming from the companion
  // ---------------------------------------------------

  if ( (comp.logn_auto).get_etat() == ETATZERO ) {
    logn_comp.set_etat_zero() ;
  }
  else{
    logn_comp.set_etat_qcq() ;
    (logn_comp.set()).import_symy( comp.logn_auto() ) ;
    logn_comp.set_std_base() ;   // set the bases for spectral expansions
  }

  a_car_auto = metgamma_auto.cov()(0,0) ;
  a_car_auto = a_car_auto() - 1.+ decouple ;
     
  if ( ((comp.metgamma_auto).cov()(0,0)).get_etat() == ETATZERO ) {
    a_car_comp.set_etat_zero() ;
  }
  else{
    a_car_comp.set_etat_qcq() ;
    (a_car_comp.set()).import_symy( (comp.metgamma_auto).cov()(0,0) - 1.
				    + comp.decouple) ;
    a_car_comp.set_std_base() ;   // set the bases for spectral expansions
  }


  // Lapse function N
  // ----------------

  Tenseur logn_total = logn_auto + logn_comp ; 

  nnn = exp( unsurc2 * logn_total ) ;

  nnn.set_std_base() ;   // set the bases for spectral expansions

  // Determination of a_car 
  // ----------------------

  a_car = a_car_auto + a_car_comp ;

  a_car.set_std_base() ;   // set the bases for spectral expansions

  // Determinant of the metric gamma.
  //---------------------------------

  gamma = pow(a_car, 3.) ;
  gamma.set_poids(2.) ;
  gamma.set_std_base() ;
   

  for (int i=0; i<3; i++) {
    gtilde_auto.set_cov(i,i) = decouple   ; 
    //	met_gamma.set_cov(i,i) = a_car() ;
    //	metgamma_auto.set_cov(i,i) = a_car_auto() ;
    //	metgamma_comp.set_cov(i,i) = a_car_comp() ;
  }
  gtilde_auto.set_std_base() ;
  gtilde_comp.set_std_base() ;
  met_gamma.set_std_base() ;
  metgamma_auto.set_std_base() ;
  metgamma_comp.set_std_base() ;

  cout << "gtilde_auto" << endl << norme(gtilde_auto.cov()(0,0)) << endl ;
 
}
