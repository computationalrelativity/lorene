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
 * $Id$
 * $Log$
 * Revision 1.8  2004/06/07 16:23:52  f_limousin
 * New treatment for conformally flat metrics.
 *
 * Revision 1.7  2004/04/08 16:33:16  f_limousin
 * The new variable is ln(Q) instead of Q=psi^2*N. It improves the
 * convergence of the code.
 *
 * Revision 1.6  2004/03/23 09:58:55  f_limousin
 * Add function Star::update_decouple()
 *
 * Revision 1.5  2004/02/27 10:54:27  f_limousin
 * To avoid errors when merging versions of Lorene.
 *
 * Revision 1.4  2004/02/27 09:56:42  f_limousin
 * Many minor changes.
 *
 * Revision 1.3  2004/02/21 17:05:13  e_gourgoulhon
 * Method Scalar::point renamed Scalar::val_grid_point.
 * Method Scalar::set_point renamed Scalar::set_grid_point.
 *
 * Revision 1.2  2004/01/20 15:20:08  f_limousin
 * First version
 *
 *
 * $Header$ *
 */

// Headers Lorene
#include "cmp.h"
#include "star.h"
#include "graphique.h"
#include "utilitaires.h"

//----------------------------------//
//	 Version without relaxation //
//----------------------------------//

void Star_bin::update_metric(const Star_bin& comp) {

    // Computation of quantities coming from the companion
    // ---------------------------------------------------

    int nz = mp.get_mg()->get_nzone() ;
    int nr = mp.get_mg()->get_nr(0);
    int nt = mp.get_mg()->get_nt(0);
    int np = mp.get_mg()->get_np(0);

    const Map& mp_comp (comp.get_mp()) ;
  
    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	logn_comp.import_symy( comp.logn_auto ) ;
	logn_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }


    shift_comp.set_etat_qcq() ; 
    shift_comp.set_triad(mp.get_bvect_cart()) ;

    Vector comp_shift(comp.shift_auto) ;
    comp_shift.change_triad(mp_comp.get_bvect_cart()) ;
    comp_shift.change_triad(mp.get_bvect_cart()) ;

    assert ( *(shift_comp.get_triad()) == *(comp_shift.get_triad())) ;

    (shift_comp.set(1)).import( comp_shift(1) ) ;  
    (shift_comp.set(2)).import( comp_shift(2) ) ;  
    (shift_comp.set(3)).import( comp_shift(3) ) ;  

    shift_comp.std_spectral_base() ;   
    shift_comp.change_triad(mp.get_bvect_spher()) ;
 

    if ( (comp.qq_auto).get_etat()  == ETATZERO ) {
	qq_comp.set_etat_zero() ;
    }
    else{
	qq_comp.set_etat_qcq() ;  
	qq_comp.import_symy( comp.qq_auto ) ;
	qq_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }	


    hij_comp.set_triad(mp.get_bvect_cart()) ;
    Tensor comp_hij(comp.hij_auto) ;
    comp_hij.change_triad(mp_comp.get_bvect_cart()) ;
    comp_hij.change_triad(mp.get_bvect_cart()) ;

    assert ( *(hij_comp.get_triad()) == *(comp_hij.get_triad())) ;


    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    
	    hij_comp.set(i,j).set_etat_qcq() ;
	    hij_comp.set(i,j).import( (comp_hij)(i,j) ) ;
	}
 
    hij_comp.std_spectral_base() ;   // set the bases for spectral expansions
    hij_comp.change_triad( mp.get_bvect_spher() ) ;
   
// Lapse function N
// ----------------

    logn = logn_auto + logn_comp ; 

    nnn = exp( logn ) ;

    nnn.std_spectral_base() ;   // set the bases for spectral expansions

// Quantity qq = log(psi^2*N)
// ----------------------

    qq = qq_auto + qq_comp ;
    
    psi4 = exp(2*qq) / (nnn * nnn) ;
    psi4.std_spectral_base() ;

// Shift vector 
// -------------

    shift = shift_auto + shift_comp ;

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

    Sym_tensor tens_gamma = gtilde_con / psi4 ;
    gamma = tens_gamma ;


    // For conformally flat metrics
    // ----------------------------

    if (conf_flat) {
	hij_auto.set_etat_zero() ; 
	hij_comp.set_etat_zero() ; 
	hij.set_etat_zero() ; 
	gtilde = flat ;
	tens_gamma = flat.con() / psi4 ;
	gamma = tens_gamma ;
    }


    // Determinant of gtilde

    Scalar det_gtilde (gtilde.determinant()) ;
       
    double* max_det = new double[nz] ;
    double* min_det = new double[nz] ;
    double* moy_det = new double[nz] ;
      
    for (int i=0; i<nz; i++){
	min_det[i] = 2 ;
	moy_det[i] = 0 ;
	max_det[i] = 0 ;
    }

    for (int l=0; l<nz; l++)
	for (int k=0; k<np; k++)
	    for (int j=0; j<nt; j++)
		for (int i=0; i<nr; i++){
	      
		    moy_det[l] = moy_det[l] + det_gtilde.val_grid_point(l,k,j,i) ;
		    if (det_gtilde.val_grid_point(l,k,j,i) > max_det[l]){
			max_det[l] = det_gtilde.val_grid_point(l,k,j,i) ;
		    }
		    if (det_gtilde.val_grid_point(l,k,j,i) < min_det[l]){
			min_det[l] = det_gtilde.val_grid_point(l,k,j,i) ;
		    }
		}
     
    cout << "average determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nz; l++){
	cout << moy_det[l]/(nr*nt*np) << "  " ;
    }
    cout << endl ;

      
    cout << "maximum of the determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nz; l++){
	cout << max_det[l] << "  " ;
    }
    cout << endl ;

    cout << "minimum of the determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nz; l++){
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

    int nz = mp.get_mg()->get_nzone() ;
    int nr = mp.get_mg()->get_nr(0);
    int nt = mp.get_mg()->get_nt(0);
    int np = mp.get_mg()->get_np(0);
  
    const Map& mp_comp (comp.get_mp()) ;

    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	logn_comp.import_symy( comp.logn_auto ) ;
	logn_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }


    shift_comp.set_etat_qcq() ; 
    shift_comp.set_triad(mp.get_bvect_cart()) ;

    Vector comp_shift(comp.shift_auto) ;
    comp_shift.change_triad(mp_comp.get_bvect_cart()) ;
    comp_shift.change_triad(mp.get_bvect_cart()) ;

    assert ( *(shift_comp.get_triad()) == *(comp_shift.get_triad())) ;

    (shift_comp.set(1)).import( comp_shift(1) ) ;  
    (shift_comp.set(2)).import( comp_shift(2) ) ;  
    (shift_comp.set(3)).import( comp_shift(3) ) ;  

    shift_comp.std_spectral_base() ;   
    shift_comp.change_triad(mp.get_bvect_spher()) ;
 
   if ( (comp.qq_auto).get_etat()  == ETATZERO ) {
	qq_comp.set_etat_zero() ;
    }
    else{
	qq_comp.set_etat_qcq() ;
	qq_comp.import_symy( comp.qq_auto ) ;
	qq_comp.std_spectral_base() ;   // set the bases for spectral expansions
    }	

    hij_comp.set_triad(mp.get_bvect_cart()) ;

    Tensor comp_hij(comp.hij_auto) ;
    comp_hij.change_triad(mp_comp.get_bvect_cart()) ;
    comp_hij.change_triad(mp.get_bvect_cart()) ;

    assert ( *(hij_comp.get_triad()) == *(comp_hij.get_triad())) ;


    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    
	    hij_comp.set(i,j).set_etat_qcq() ;
	    hij_comp.set(i,j).import( (comp_hij)(i,j) ) ;
	}
 
    hij_comp.std_spectral_base() ;
    hij_comp.change_triad( mp.get_bvect_spher() ) ;


  
// Relaxation on logn_comp, shift_comp, qq_comp, hij_comp
// ---------------------------------------------------------------
    double relaxjm1 = 1. - relax ; 
    
    logn_comp = relax * logn_comp + relaxjm1 * (star_jm1.logn_comp) ; 
    
    shift_comp = relax * shift_comp + relaxjm1 
	                               * (star_jm1.shift_comp) ; 

    qq_comp = relax * qq_comp + relaxjm1 * (star_jm1.qq_comp) ;

       
    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {

	    hij_comp.set(i,j) = relax * hij_comp(i,j) 
		+ relaxjm1 * (star_jm1.hij_comp)(i,j) ; 
	
	}

// Lapse function N
// ----------------

    logn = logn_auto + logn_comp ; 

    nnn = exp( logn ) ;

    nnn.std_spectral_base() ;   // set the bases for spectral expansions


// Quantity qq = log(psi^2*N)
// --------------------------

    qq = qq_auto + qq_comp ;
    
    psi4 = exp(2*qq) / (nnn * nnn) ;
    psi4.std_spectral_base() ;

// Shift vector
// ------------

    shift = shift_auto + shift_comp ;
    
// Coefficients of the 3-metric tilde
// ----------------------------------
     
    Sym_tensor gtilde_con(mp, CON, mp.get_bvect_spher()) ;
    hij_auto.std_spectral_base() ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
   
	    hij.set(i,j) = hij_auto(i,j) + hij_comp(i,j) ;
	    gtilde_con.set(i,j) = hij(i,j) + flat.con()(i,j) ;
	}
    
    gtilde_con.std_spectral_base() ;
    gtilde = gtilde_con ;
    Tensor tens_gamma(gtilde_con / psi4) ;
    gamma = tens_gamma ;
      

    // For conformally flat metrics
    // ----------------------------

    if (conf_flat) {
	hij_auto.set_etat_zero() ; 
	hij_comp.set_etat_zero() ; 
	hij.set_etat_zero() ; 
	gtilde = flat ;
	tens_gamma = flat.con() / psi4 ;
	gamma = tens_gamma ;
    }


/*
    // Determination of h33 in order to have det(gtilde) = 1
    // -----------------------------------------------------

    Scalar gtilde11 (gtilde.cov()(1,1)) ;
    Scalar gtilde21 (gtilde.cov()(2,1)) ;
    Scalar gtilde31 (gtilde.cov()(3,1)) ;
    Scalar gtilde22 (gtilde.cov()(2,2)) ;
    Scalar gtilde32 (gtilde.cov()(3,2)) ;
    Scalar gtilde33 (gtilde.cov()(3,3)) ;

    gtilde33 = ( 1 - (gtilde21*gtilde32*gtilde31 + gtilde31*gtilde21*gtilde32
	       - gtilde31*gtilde22*gtilde31 - gtilde11*gtilde32*gtilde32) ) /
	  ( gtilde11*gtilde22 - gtilde21*gtilde21 ) ;

    Sym_tensor gtilde_cov(mp, COV, mp.get_bvect_spher()) ;
    gtilde_cov.set(1,1) = gtilde11 ;
    gtilde_cov.set(2,1) = gtilde21 ;
    gtilde_cov.set(3,1) = gtilde31 ;
    gtilde_cov.set(2,2) = gtilde22 ;
    gtilde_cov.set(3,2) = gtilde32 ;
    gtilde_cov.set(3,3) = gtilde33 ;

    gtilde_cov.std_spectral_base() ;
    gtilde = gtilde_cov ;
    tens_gamma = gtilde_con / psi4 ;
    gamma = tens_gamma ;
 
    hij.set(3,3) = gtilde.con()(3,3) - 1 ;
    hij.std_spectral_base() ;
*/
// Determinant of gtilde

    Scalar det_gtilde (gtilde.determinant()) ;
       
    double* max_det = new double[nz] ;
    double* min_det = new double[nz] ;
    double* moy_det = new double[nz] ;
      
    for (int i=0; i<nz; i++){
	min_det[i] = 2 ;
	moy_det[i] = 0 ;
	max_det[i] = 0 ;
    }

    for (int l=0; l<nz; l++)
	for (int k=0; k<np; k++)
	    for (int j=0; j<nt; j++)
		for (int i=0; i<nr; i++){
	      
		    moy_det[l] = moy_det[l] + det_gtilde.val_grid_point(l,k,j,i) ;
		    if (det_gtilde.val_grid_point(l,k,j,i) > max_det[l]){
			max_det[l] = det_gtilde.val_grid_point(l,k,j,i) ;
		    }
		    if (det_gtilde.val_grid_point(l,k,j,i) < min_det[l]){
			min_det[l] = det_gtilde.val_grid_point(l,k,j,i) ;
		    }
		}
     
    cout << "average determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nz; l++){
	cout << moy_det[l]/(nr*nt*np) << "  " ;
    }
    cout << endl ;

      
    cout << "maximum of the determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nz; l++){
	cout << max_det[l] << "  " ;
    }
    cout << endl ;

    cout << "minimum of the determinant of gtilde in each zone : " << endl ; 
    for (int l=0; l<nz; l++){
	cout << min_det[l] << "  " ;
    }
    cout << endl ;
            
    //## Juste pour des test de convergence

//    tens_gamma = flat.con() / psi4 ;
//   gamma = tens_gamma ;



    // ... extrinsic curvature (tkij_auto and kcar_auto)
    extrinsic_curvature() ;

   
// The derived quantities are obsolete
// -----------------------------------

    del_deriv() ;


}

void Star_bin::update_metric_init1() {

    logn_auto = logn ;
    qq_auto = qq ;

}

void Star_bin::update_metric_init2(const Star_bin& comp) {

    logn_auto = logn ;
    
    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	logn_comp.import_symy( comp.logn_auto ) ;
	logn_comp.std_spectral_base() ;  // set the bases for spectral expansions
    }
    
    logn = logn_auto + logn_comp ; 
}

void Star_bin::update_decouple(const Star_bin& comp) {

    int nr = mp.get_mg()->get_nr(0);
    int nt = mp.get_mg()->get_nt(0);
    int np = mp.get_mg()->get_np(0);

//    decouple = (logn_auto - 1.e-12) / (logn - 2.e-12) ;

    decouple = 0.5 ;

    Cmp decouple_cmp (decouple) ;

    cout << "decouple" << endl << norme(decouple/(nr*nt*np)) << endl ;

//    des_profile(decouple, 0., 40., 1., 1.5, 0.) ;

}
