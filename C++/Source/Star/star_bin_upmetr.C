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
 * Revision 1.10  2005/02/17 17:34:10  f_limousin
 * Change the name of some quantities to be consistent with other classes
 * (for instance nnn is changed to nn, shift to beta, beta to lnq...)
 *
 * Revision 1.9  2004/06/22 12:52:26  f_limousin
 * Change qq, qq_auto and qq_comp to beta, beta_auto and beta_comp.
 *
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
#include "param.h"

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
    }


    beta_comp.set_etat_qcq() ; 
    beta_comp.set_triad(mp.get_bvect_cart()) ;

    Vector comp_beta(comp.beta_auto) ;
    comp_beta.change_triad(mp_comp.get_bvect_cart()) ;
    comp_beta.change_triad(mp.get_bvect_cart()) ;

    assert ( *(beta_comp.get_triad()) == *(comp_beta.get_triad())) ;

    (beta_comp.set(1)).import( comp_beta(1) ) ;  
    (beta_comp.set(2)).import( comp_beta(2) ) ;  
    (beta_comp.set(3)).import( comp_beta(3) ) ;  

    beta_comp.std_spectral_base() ;   
    beta_comp.change_triad(mp.get_bvect_spher()) ;
 

    if ( (comp.lnq_auto).get_etat()  == ETATZERO ) {
	lnq_comp.set_etat_zero() ;
    }
    else{
	lnq_comp.set_etat_qcq() ;  
	lnq_comp.import_symy( comp.lnq_auto ) ;
    }	


    hh_comp.set_triad(mp.get_bvect_cart()) ;
    Tensor comp_hh(comp.hh_auto) ;
    comp_hh.change_triad(mp_comp.get_bvect_cart()) ;
    comp_hh.change_triad(mp.get_bvect_cart()) ;

    assert ( *(hh_comp.get_triad()) == *(comp_hh.get_triad())) ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    
	    hh_comp.set(i,j).set_etat_qcq() ;
	    hh_comp.set(i,j).import( (comp_hh)(i,j) ) ;
	}
 
    hh_comp.std_spectral_base() ;   // set the bases for spectral expansions
    hh_comp.change_triad( mp.get_bvect_spher() ) ;
   
// Lapse function N
// ----------------

    logn = logn_auto + logn_comp ; 

    nn = exp( logn ) ;

    nn.std_spectral_base() ;   // set the bases for spectral expansions

// Quantity lnq = log(psi^2*N)
// ----------------------

    lnq = lnq_auto + lnq_comp ;
    
    psi4 = exp(2*lnq) / (nn * nn) ;
    psi4.std_spectral_base() ;

// Beta vector 
// -------------

    beta = beta_auto + beta_comp ;

// Coefficients of the 3-metric tilde
// ----------------------------------
 
    Sym_tensor gtilde_con (mp, CON, mp.get_bvect_spher()) ; 
    
     for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
   
	    hh.set(i,j) = hh_auto(i,j) + hh_comp(i,j) ;
	    gtilde_con.set(i,j) = hh(i,j) + flat.con()(i,j) ;
	}

     gtilde = gtilde_con ;

     Sym_tensor tens_gamma = gtilde_con / psi4 ;
     gamma = tens_gamma ;


    // For conformally flat metrics
    // ----------------------------

    if (conf_flat) {
	hh_auto.set_etat_zero() ; 
	hh_comp.set_etat_zero() ; 
	hh.set_etat_zero() ; 
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
         
    // ... extrinsic curvature (aa_auto and aa_quad_auto)
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
    }


    beta_comp.set_etat_qcq() ; 
    beta_comp.set_triad(mp.get_bvect_cart()) ;

    Vector comp_beta(comp.beta_auto) ;
    comp_beta.change_triad(mp_comp.get_bvect_cart()) ;
    comp_beta.change_triad(mp.get_bvect_cart()) ;

    assert ( *(beta_comp.get_triad()) == *(comp_beta.get_triad())) ;

    (beta_comp.set(1)).import( comp_beta(1) ) ;  
    (beta_comp.set(2)).import( comp_beta(2) ) ;  
    (beta_comp.set(3)).import( comp_beta(3) ) ;  

    beta_comp.std_spectral_base() ;   
    beta_comp.change_triad(mp.get_bvect_spher()) ;

 
   if ( (comp.lnq_auto).get_etat()  == ETATZERO ) {
	lnq_comp.set_etat_zero() ;
    }
    else{
	lnq_comp.set_etat_qcq() ;
	lnq_comp.import_symy( comp.lnq_auto ) ;
    }	

    hh_comp.set_triad(mp.get_bvect_cart()) ;

    Tensor comp_hh(comp.hh_auto) ;
    comp_hh.change_triad(mp_comp.get_bvect_cart()) ;
    comp_hh.change_triad(mp.get_bvect_cart()) ;

    assert ( *(hh_comp.get_triad()) == *(comp_hh.get_triad())) ;


    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    
	    hh_comp.set(i,j).set_etat_qcq() ;
	    hh_comp.set(i,j).import( (comp_hh)(i,j) ) ;
	}
 
    hh_comp.std_spectral_base() ;
    hh_comp.change_triad( mp.get_bvect_spher() ) ;


  
// Relaxation on logn_comp, beta_comp, lnq_comp, hh_comp
// ---------------------------------------------------------------
    double relaxjm1 = 1. - relax ; 
    
    logn_comp = relax * logn_comp + relaxjm1 * (star_jm1.logn_comp) ; 
    
    beta_comp = relax * beta_comp + relaxjm1 
	                               * (star_jm1.beta_comp) ; 

    lnq_comp = relax * lnq_comp + relaxjm1 * (star_jm1.lnq_comp) ;

       
    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {

	    hh_comp.set(i,j) = relax * hh_comp(i,j) 
		+ relaxjm1 * (star_jm1.hh_comp)(i,j) ; 
	
	}

// Lapse function N
// ----------------

    logn = logn_auto + logn_comp ; 

    nn = exp( logn ) ;

    nn.std_spectral_base() ;   // set the bases for spectral expansions


// Quantity lnq = log(psi^2 * N)
// --------------------------

    lnq = lnq_auto + lnq_comp ;
    
    psi4 = exp(2*lnq) / (nn * nn) ;
    psi4.std_spectral_base() ;

// Beta vector
// ------------

    beta = beta_auto + beta_comp ;
    

// Computation of hh (at this point we only have h^{ij} TT))
// ---------------------------------------------------
     
    hh_auto.std_spectral_base() ;
    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) 
	    hh.set(i,j) = hh_auto(i,j) + hh_comp(i,j) ;
	
     
    // Computation of h ( See eq. 116 of BGGN )
    // ----------------------------------------
    
    int it_max = 200 ;
    double precis = 1.e-4 ;
	    
    // The trace h = f_{ij} h^{ij} :
    Scalar htrace(mp) ;
    
    hh.inc_dzpuis(2) ;

    Sym_tensor_tt hh_tt (mp, mp.get_bvect_spher(), flat) ;
    hh_tt = hh ;
    
    // Value of h at previous step of the iterative procedure below :
    Scalar htrace_prev(mp) ;
    htrace_prev.set_etat_zero() ;   // initialisation to zero
     // Parameters for the poisson equation in set_tt_trace
    Param par ;
    int mermax_poisson = 4 ;
    double relax_poisson = 0.5 ;
    double precis_poisson = 1.e-15 ;
    int niter ;
    Cmp ssjm1_htrace (mp) ;
    ssjm1_htrace = 0 ;

    par.add_int(mermax_poisson,  0) ;  // maximum number of iterations
    par.add_double(relax_poisson,  0) ; // relaxation parameter
    par.add_double(precis_poisson, 1) ; // required precision
    par.add_int_mod(niter, 0) ; // number of iterations actually used 
    par.add_cmp_mod( ssjm1_htrace ) ; 


    int it ;
    for (it=0; it<it_max; it++) {
      
	// Trace h from the condition det(f^{ij} + h^{ij}) = det f^{ij} :
      
	htrace = hh(1,1) * hh(2,3) * hh(2,3) 
	    + hh(2,2) * hh(1,3) * hh(1,3) 
	    + hh(3,3) * hh(1,2) * hh(1,2)
	    - 2.* hh(1,2) * hh(1,3) * hh(2,3) 
	    - hh(1,1) * hh(2,2) * hh(3,3) ;
     
        
	htrace.dec_dzpuis(2) ; // dzpuis: 6 --> 4
        
	htrace += hh(1,2) * hh(1,2) 
	    + hh(1,3) * hh(1,3) 
	    + hh(2,3) * hh(2,3) 
	    - hh(1,1) * hh(2,2) 
	    - hh(1,1) * hh(3,3) 
	    - hh(2,2) * hh(3,3) ;

	// New value of hh from htrace and hh_tt 
	// (obtained by solving 
	// the Poisson equation for Phi) : 
		
	htrace.std_spectral_base() ;
	Sym_tensor_trans hh_trans (mp, mp.get_bvect_spher(), flat) ;
	if (htrace.get_etat() == ETATQCQ){
	    hh_trans.set_tt_trace(hh_tt, htrace, &par) ; 
	    hh = hh_trans ;
	}
	else
	    hh.set_etat_zero() ;

	    cout << "hh_trans :" << endl ;
	    for (int i=1; i<=3; i++)
		for (int j=1; j<=i; j++) {
		    cout << "  Comp. " << i << " " << j << " :  " ;
		    for (int l=0; l<nz; l++){
			cout << norme(hh(i,j)/(nr*nt*np))(l) << " " ;
		    }
		    cout << endl ;
		}
	    cout << endl ;
	    cout << "h trace" << endl << norme(htrace) << endl ;

	double diff = max(max(abs(htrace - htrace_prev))) ;
	cout << " det gtilde = 1: "  
	     << "iteration : " << it << " difference : " << diff << endl ;
	if (diff < precis) break ;
	else htrace_prev = htrace ;

    }

    if (it == it_max) {
	cerr << "hh_auto det = 1 : convergence not reached \n" ;
	cerr << "  for the required accuracy (" << precis << ") ! " << endl ;
	abort() ;
    }

    hh.dec_dzpuis(2) ;
    

// Coefficients of the 3-metric tilde
// ----------------------------------
     
    Sym_tensor gtilde_con(mp, CON, mp.get_bvect_spher()) ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
   
	    hh.set(i,j) = hh_auto(i,j) + hh_comp(i,j) ;
	    gtilde_con.set(i,j) = hh(i,j) + flat.con()(i,j) ;
	}
    
    gtilde = gtilde_con ;
    Tensor tens_gamma(gtilde_con / psi4) ;
    gamma = tens_gamma ;
      

    // For conformally flat metrics
    // ----------------------------

    if (conf_flat) {
	hh_auto.set_etat_zero() ; 
	hh_comp.set_etat_zero() ; 
	hh.set_etat_zero() ; 
	gtilde = flat ;
	tens_gamma = flat.con() / psi4 ;
	gamma = tens_gamma ;
    }





// Computation of det(gtilde)

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
            

    // ... extrinsic curvature (aa_auto and aa_quad_auto)
    extrinsic_curvature() ;

   
// The derived quantities are obsolete
// -----------------------------------

    del_deriv() ;


}

void Star_bin::update_metric_init1() {

    logn_auto = logn ;
    lnq_auto = lnq ;

}

void Star_bin::update_metric_init2(const Star_bin& comp) {

    logn_auto = logn ;
    
    if ( (comp.logn_auto).get_etat() == ETATZERO ) {
	logn_comp.set_etat_zero() ;
    }
    else{
	logn_comp.set_etat_qcq() ;
	logn_comp.import_symy( comp.logn_auto ) ;
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
