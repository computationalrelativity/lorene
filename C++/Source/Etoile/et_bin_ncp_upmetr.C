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

		    //----------------------------------//
		    //	 Version without relaxation	//
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

	shift_comp.set_triad( mp.get_bvect_cart() ) ;    	
	shift_comp.set_std_base() ;   // set the bases for spectral expansions
    }
    


    if ( (comp.loggamma_auto).get_etat()  == ETATZERO ) {
      loggamma_comp.set_etat_zero() ;
    }
    else{
      loggamma_comp.set_etat_qcq() ;
      (loggamma_comp.set()).import_symy( comp.loggamma_auto() ) ;
      loggamma_comp.set_std_base() ;   // set the bases for spectral expansions
    }	


    (gtilde_comp.set_cov()).set_triad( *(((comp.gtilde_auto).cov()).get_triad()) ) ;  
     assert ( *((gtilde_comp.cov()).get_triad()) == *(((comp.gtilde_auto).cov()).get_triad())) ; 
    
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

    (gtilde_comp.set_cov()).set_triad( mp.get_bvect_cart() ) ;
    gtilde_comp.set_std_base() ;   // set the bases for spectral expansions



    // Lapse function N
    // ----------------

    Tenseur logn_total = logn_auto + logn_comp ; 

    nnn = exp( unsurc2 * logn_total ) ;

    nnn.set_std_base() ;   // set the bases for spectral expansions



    // Shift vector beta^i
    // -------------------

    shift = shift_auto + shift_comp ;

    // log (determinant ) and determinant
    // ----------------------------------

    loggamma.set() = loggamma_auto() + loggamma_comp() ;

    gamma = exp (loggamma) ;
    gamma.set_poids(2.) ;
    gamma.set_std_base() ;


    // Conformal factor A^2
    // ---------------------
 
    a_car = pow(gamma(), 1./6.) ;
    a_car.set_std_base() ;   // set the bases for spectral expansions    

   // Coefficients of the 3-metric tilde
    // ----------------------------------

    assert ( *(gtilde_auto.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
 
    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
   
	gtilde.set_cov(i,j) = (gtilde_auto.cov())(i,j) + (gtilde_comp.cov())(i,j) ;
      }
    }
 

   gtilde.set_std_base() ;

    // Coefficients of the 3-metric
    // ----------------------------

    assert ( *(met_gamma.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    assert ( *(metgamma_auto.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    assert ( *(metgamma_comp.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;


    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
     
	met_gamma.set_cov(i,j) = pow(gamma(), 1./3.) * gtilde.cov()(i,j) ;
      }
    }  
    met_gamma.set_std_base() ;
    

    // Coefficients of the 3-metric auto and comp
    // ------------------------------------------

    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
     
	metgamma_auto.set_cov(i,j) = pow(gamma(), 1./3.) * gtilde_auto.cov()(i,j) ;
      }
    }    
    metgamma_auto.set_std_base() ;

    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
     
	metgamma_comp.set_cov(i,j) = pow(gamma(), 1./3.) * gtilde_comp.cov()(i,j) ;
      }
    }    
    metgamma_comp.set_std_base() ;


    // Derivatives of metric coefficients
    // ----------------------------------

    // ... (d/dX,d/dY,d/dZ)(logn_auto) :
    d_logn_auto_regu = logn_auto_regu.gradient() ;    // (d/dx, d/dy, d/dz)
    d_logn_auto_regu.change_triad(mp.get_bvect_cart()) ;   // -->  (d/dX, d/dY, d/dZ)

 	// Change the basis from spherical coordinate to Cartesian one
	d_logn_auto_div.change_triad( mp.get_bvect_cart() ) ;

	
    d_logn_auto = d_logn_auto_regu + d_logn_auto_div ;

    // ... (d/dX,d/dY,d/dZ)(beta_auto) :
    d_beta_auto = beta_auto.gradient() ;    // (d/dx, d/dy, d/dz)
    d_beta_auto.change_triad(mp.get_bvect_cart()) ;   // -->  (d/dX, d/dY, d/dZ)

    if (relativistic) {
	// ... extrinsic curvature (tkij_auto and kcar_auto)
      extrinsic_curvature() ;
    }
    
    cout << "update metric" << endl ;

    // The derived quantities are obsolete
    // -----------------------------------

    del_deriv() ;


}



		    //----------------------------------//
		    //	  Version with relaxation       //
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

	shift_comp.set_triad( mp.get_bvect_cart() ) ;    	
	shift_comp.set_std_base() ;   // set the bases for spectral expansions
    }


    if ( (comp.loggamma_auto).get_etat()  == ETATZERO ) {
      loggamma_comp.set_etat_zero() ;
    }
    else{
      loggamma_comp.set_etat_qcq() ;
      (loggamma_comp.set()).import_symy( comp.loggamma_auto() ) ;
      loggamma_comp.set_std_base() ;   // set the bases for spectral expansions
    }	

	 

    (gtilde_comp.set_cov()).set_triad( *(((comp.gtilde_auto).cov()).get_triad()) ) ;  
    assert ( *(gtilde_comp.cov().get_triad()) == *(((comp.gtilde_auto).cov()).get_triad())) ; 

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
    (gtilde_comp.set_cov()).set_triad( mp.get_bvect_cart() ) ;


    // Relaxation on logn_comp, shift_comp, loggamma_comp, gtilde_comp
    // ---------------------------------------------------------------
    double relaxjm1 = 1. - relax ; 
    
    logn_comp = relax * logn_comp + relaxjm1 * (star_jm1.get_logn_comp()) ; 
    
    shift_comp = relax * shift_comp + relaxjm1 * (star_jm1.get_shift_comp()) ; 

    loggamma_comp = relax * loggamma_comp + relaxjm1 * (star_jm1.get_loggamma_comp()) ;

    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
	gtilde_comp.set_cov(i,j) = relax * (gtilde_comp.cov())(i,j) + relaxjm1 * ((star_jm1.get_gtilde_comp()).cov())(i,j) ;
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


    // log (determinant )
    // ------------------

    loggamma = loggamma_auto + loggamma_comp ;

    gamma = exp (loggamma) ;
    gamma.set_poids(2.) ;
    gamma.set_std_base() ;
 
  
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
    
    // Coefficients of the 3-metric
    // ----------------------------

    assert ( *(met_gamma.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    assert ( *(metgamma_auto.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    assert ( *(metgamma_comp.cov().get_triad())== *(gtilde.cov().get_triad()) ) ;
    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {

	met_gamma.set_cov(i,j) = pow(gamma(), 1./3.) * gtilde.cov()(i,j) ;
      }
    }    
    met_gamma.set_std_base() ;

    // Coefficients of the 3-metric auto and comp
    // ------------------------------------------

    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
     
	metgamma_auto.set_cov(i,j) = pow(gamma(), 1./3.) * gtilde_auto.cov()(i,j) ;

      }
    }    
    metgamma_auto.set_std_base() ;

    for(int i=0; i<=2; i++) {
      for(int j=i; j<=2; j++) {
     
	metgamma_comp.set_cov(i,j) = pow(gamma(), 1./3.) * gtilde_comp.cov()(i,j) ;

      }
    }    
    metgamma_comp.set_std_base() ;


    // Derivatives of metric coefficients
    // ----------------------------------
    
    // ... (d/dX,d/dY,d/dZ)(logn_auto) : 
    d_logn_auto_regu = logn_auto_regu.gradient() ;    // (d/dx, d/dy, d/dz)
    d_logn_auto_regu.change_triad(mp.get_bvect_cart()) ;   // -->  (d/dX, d/dY, d/dZ)  
    

	// Change the basis from spherical coordinate to Cartesian one
	d_logn_auto_div.change_triad( mp.get_bvect_cart() ) ;


    d_logn_auto = d_logn_auto_regu + d_logn_auto_div ;      
    
    // ... (d/dX,d/dY,d/dZ)(beta_auto) : 
    d_beta_auto = beta_auto.gradient() ;    // (d/dx, d/dy, d/dz)
    d_beta_auto.change_triad(mp.get_bvect_cart()) ;   // -->  (d/dX, d/dY, d/dZ)        

    // ... extrinsic curvature (tkij_auto and kcar_auto)
    extrinsic_curvature() ; 
    
    // The derived quantities are obsolete
    // -----------------------------------
    
    del_deriv() ;                
    
			      
} 





