/*
 * Methods of Star_bin::update_metric
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


char star_binupmetr_C[] = "$Header$" ;

/*
 * $Header$ *
 */

// Headers Lorene
#include "star.h"
#include "graphique.h"

//----------------------------------//
//	 Version without relaxation //
//----------------------------------//

void Star_bin::update_metric(const Star_bin& comp) {

    // Computation of quantities coming from the companion
    // ---------------------------------------------------

    int nzone = mp.get_mg()->get_nzone() ;
    int nr = mp.get_mg()->get_nr(0);
    int nt = mp.get_mg()->get_nt(0);
    int np = mp.get_mg()->get_np(0);
  
    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	logn_comp.import_symy( comp.logn_auto ) ;
	logn_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }


    shift_comp.set_etat_qcq() ; 
  
    shift_comp.set_triad( *((comp.shift_auto).get_triad()) ) ;  
    assert ( *(shift_comp.get_triad()) == *((comp.shift_auto).get_triad())) ;
  
    (shift_comp.set(0)).import_asymy( comp.shift_auto(0) ) ;  // N^x antisym
    (shift_comp.set(1)).import_symy( comp.shift_auto(1) ) ;   // N^y sym.
    (shift_comp.set(2)).import_asymy( comp.shift_auto(2) ) ;  // N^z anisym
  
    shift_comp.change_triad( mp.get_bvect_spher() ) ;    	
    shift_comp.std_spectral_base() ;   // set the bases for spectral expansions




    if ( (comp.qq_auto).get_etat()  == ETATZERO ) {
	qq_comp.set_etat_zero() ;
    }
    else{
	qq_comp.set_etat_qcq() ;
  
	qq_comp.import_symy( comp.qq_auto ) ;
	qq_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }	


    hij_comp.set_triad( *(comp.hij_auto.get_triad()) ) ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    
	    hij_comp.set(i,j).set_etat_qcq() ;
	    hij_comp.set(i,j).import( (comp.hij_auto)(i,j) ) ;
	}
 
    hij_comp.change_triad( mp.get_bvect_cart() ) ;
    hij_comp.std_spectral_base() ;   // set the bases for spectral expansions
  
// Lapse function N
// ----------------

    Scalar logn_total = logn_auto + logn_comp ; 

    nnn = exp( logn_total ) ;

    nnn.std_spectral_base() ;   // set the bases for spectral expansions

// Shift vector beta^i
// -------------------

    shift = shift_auto + shift_comp ;

 
// Quantity qq = psi^2*N
// ----------------------

    qq = qq_auto + qq_comp ;

    psi4 = qq * qq / (nnn*nnn) ;
    psi4.std_spectral_base() ;

// Coefficients of the 3-metric tilde
// ----------------------------------
 
    Sym_tensor gtilde_con (mp, CON, mp.get_bvect_spher()) ; 
    
     for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
   
	    hij.set(i,j) = hij_auto(i,j) + hij_comp(i,j) ;
	    gtilde_con.set(i,j) = hij(i,j) + flat.con()(i,j) ;
	}
    
    gtilde_con.std_spectral_base() ;
    gtilde = gtilde_con ;
    hij.std_spectral_base() ;
      
// Determinant of gtilde

    Scalar det_gtilde(mp) ;
    
    Sym_tensor gtilde_cov (gtilde.cov()) ;
    det_gtilde = gtilde_cov(0, 0)*gtilde_cov(1, 1)*gtilde_cov(2, 2) 
	+ gtilde_cov(0, 1)*gtilde_cov(1, 2)*gtilde_cov(2, 0)
	+ gtilde_cov(0, 2)*gtilde_cov(1, 0)*gtilde_cov(2, 1) 
	- gtilde_cov(2, 0)*gtilde_cov(1, 1)*gtilde_cov(0, 2)
	- gtilde_cov(2, 1)*gtilde_cov(1, 2)*gtilde_cov(0, 0) 
	- gtilde_cov(2, 2)*gtilde_cov(1, 0)*gtilde_cov(0, 1) ;

       
    double* max_det = new double[nzone] ;
    double* min_det = new double[nzone] ;
    double* moy_det = new double[nzone] ;
      
    for (int i=0; i<nzone; i++){
	min_det[i] = 2 ;
	moy_det[i] = 0 ;
	max_det[i] = 0 ;
    }

    for (int l=0; l<nzone; l++)
	for (int k=0; k<np; k++)
	    for (int j=0; j<nt; j++)
		for (int i=0; i<nr; i++){
	      
		    moy_det[l] = moy_det[l] + det_gtilde.point(l,k,j,i) ;
		    if (det_gtilde.point(l,k,j,i) > max_det[l]){
			max_det[l] = det_gtilde.point(l,k,j,i) ;
		    }
		    if (det_gtilde.point(l,k,j,i) < min_det[l]){
			min_det[l] = det_gtilde.point(l,k,j,i) ;
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
          
    // ... extrinsic curvature (tkij_auto and kcar_auto)
    extrinsic_curvature() ;

   
// The derived quantities are obsolete
// -----------------------------------

    del_deriv() ;


}



//----------------------------------//
//	  Version with relaxation   //
//----------------------------------//

void Star_bin::update_metric(const Star_bin& comp,
			       const Star_bin& star_jm1, double relax) {


     // Computation of quantities coming from the companion
    // ---------------------------------------------------

    int nzone = mp.get_mg()->get_nzone() ;
    int nr = mp.get_mg()->get_nr(0);
    int nt = mp.get_mg()->get_nt(0);
    int np = mp.get_mg()->get_np(0);
  
    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	logn_comp.import_symy( comp.logn_auto ) ;
	logn_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }


    shift_comp.set_etat_qcq() ; 
  
    shift_comp.set_triad( *((comp.shift_auto).get_triad()) ) ;  
    assert ( *(shift_comp.get_triad()) == *((comp.shift_auto).get_triad())) ;
  
    (shift_comp.set(0)).import_asymy( comp.shift_auto(0) ) ;  // N^x antisym
    (shift_comp.set(1)).import_symy( comp.shift_auto(1) ) ;   // N^y sym.
    (shift_comp.set(2)).import_asymy( comp.shift_auto(2) ) ;  // N^z anisym
  
    shift_comp.change_triad( mp.get_bvect_spher() ) ;    	
    shift_comp.std_spectral_base() ;   // set the bases for spectral expansions




    if ( (comp.qq_auto).get_etat()  == ETATZERO ) {
	qq_comp.set_etat_zero() ;
    }
    else{
	qq_comp.set_etat_qcq() ;
  
	qq_comp.import_symy( comp.qq_auto ) ;
	qq_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }	


    hij_comp.set_triad( *(comp.hij_auto.get_triad()) ) ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {

	    hij_comp.set(i,j).set_etat_qcq() ;
	    hij_comp.set(i,j).import( (comp.hij_auto)(i,j) ) ;
	}
 
    hij_comp.change_triad( mp.get_bvect_cart() ) ;
    hij_comp.std_spectral_base() ;   // set the bases for spectral expansions
  
// Relaxation on logn_comp, shift_comp, qq_comp, hij_comp
// ---------------------------------------------------------------
    double relaxjm1 = 1. - relax ; 
    
    logn_comp = relax * logn_comp + relaxjm1 * (star_jm1.logn_comp) ; 
    
    shift_comp = relax * shift_comp + relaxjm1 * (star_jm1.shift_comp) ; 

    qq_comp = relax * qq_comp + relaxjm1 * (star_jm1.qq_comp) ;

       
    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {

	    hij_comp.set(i,j) = relax * hij_comp(i,j) 
		+ relaxjm1 * (star_jm1.hij_comp)(i,j) ; 
	
	}

// Lapse function N
// ----------------

    Scalar logn_total = logn_auto + logn_comp ; 

    nnn = exp( logn_total ) ;

    nnn.std_spectral_base() ;   // set the bases for spectral expansions

// Shift vector beta^i
// -------------------

    shift = shift_auto + shift_comp ;

 
// Quantity qq = psi^2*N
// ----------------------

    qq = qq_auto + qq_comp ;

    psi4 = qq * qq / (nnn*nnn) ;
    psi4.std_spectral_base() ;

// Coefficients of the 3-metric tilde
// ----------------------------------
     
    Sym_tensor gtilde_con(mp, CON, mp.get_bvect_spher()) ;
  
    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
   
	    hij.set(i,j) = hij_auto(i,j) + hij_comp(i,j) ;
	    gtilde_con.set(i,j) = hij(i,j) + flat.con()(i,j) ;
	}
    
    gtilde_con.std_spectral_base() ;
    gtilde = gtilde_con ;
    hij.std_spectral_base() ;
      
// Determinant of gtilde

    Scalar det_gtilde(mp) ;
    Sym_tensor gtilde_cov(gtilde.cov()) ;
   
    det_gtilde = gtilde_cov(0, 0)*gtilde_cov(1, 1)*gtilde_cov(2, 2) 
	+ gtilde_cov(0, 1)*gtilde_cov(1, 2)*gtilde_cov(2, 0)
	+ gtilde_cov(0, 2)*gtilde_cov(1, 0)*gtilde_cov(2, 1) 
	- gtilde_cov(2, 0)*gtilde_cov(1, 1)*gtilde_cov(0, 2)
	- gtilde_cov(2, 1)*gtilde_cov(1, 2)*gtilde_cov(0, 0) 
	- gtilde_cov(2, 2)*gtilde_cov(1, 0)*gtilde_cov(0, 1) ;

       
    double* max_det = new double[nzone] ;
    double* min_det = new double[nzone] ;
    double* moy_det = new double[nzone] ;
      
    for (int i=0; i<nzone; i++){
	min_det[i] = 2 ;
	moy_det[i] = 0 ;
	max_det[i] = 0 ;
    }

    for (int l=0; l<nzone; l++)
	for (int k=0; k<np; k++)
	    for (int j=0; j<nt; j++)
		for (int i=0; i<nr; i++){
	      
		    moy_det[l] = moy_det[l] + det_gtilde.point(l,k,j,i) ;
		    if (det_gtilde.point(l,k,j,i) > max_det[l]){
			max_det[l] = det_gtilde.point(l,k,j,i) ;
		    }
		    if (det_gtilde.point(l,k,j,i) < min_det[l]){
			min_det[l] = det_gtilde.point(l,k,j,i) ;
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
          
    // ... extrinsic curvature (tkij_auto and kcar_auto)
    extrinsic_curvature() ;

   
// The derived quantities are obsolete
// -----------------------------------

    del_deriv() ;


}

void Star_bin::update_metric_init(const Star_bin& comp) {

// Computation of quantities coming from the companion
// ---------------------------------------------------

    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	logn_comp.import_symy( comp.logn_auto ) ;
	logn_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }

    qq_auto = qq_auto - 1 + decouple ;

// Lapse function N
// ----------------

    logn = logn_auto + logn_comp ; 

    nnn = exp( logn ) ;

    nnn.std_spectral_base() ;   // set the bases for spectral expansions

}


