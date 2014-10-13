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
 * Revision 1.14  2014/10/13 08:53:38  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2005/09/13 19:38:31  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.12  2005/02/24 16:05:14  f_limousin
 * Change the name of some variables (for instance dcov_logn --> dlogn).
 *
 * Revision 1.11  2005/02/18 13:14:18  j_novak
 * Changing of malloc/free to new/delete + suppression of some unused variables
 * (trying to avoid compilation warnings).
 *
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

//----------------------------------//
//	 Version without relaxation //
//----------------------------------//

namespace Lorene {
void Star_bin::update_metric(const Star_bin& comp, double om) {

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
	logn_comp.import( comp.logn_auto ) ;
	logn_comp.std_spectral_base() ;
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

    if ( (comp.lnq_auto).get_etat()  == ETATZERO ) {
	lnq_comp.set_etat_zero() ;
    }
    else{
	lnq_comp.set_etat_qcq() ;  
	lnq_comp.import( comp.lnq_auto ) ;
	lnq_comp.std_spectral_base() ;
    }	


    hij_comp.set_triad(mp.get_bvect_cart()) ;
    Sym_tensor comp_hij(comp.hij_auto) ;
    comp_hij.change_triad(mp_comp.get_bvect_cart()) ;
    comp_hij.change_triad(mp.get_bvect_cart()) ;

    assert ( *(hij_comp.get_triad()) == *(comp_hij.get_triad())) ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    
	    hij_comp.set(i,j).set_etat_qcq() ;
	    hij_comp.set(i,j).import( (comp_hij)(i,j) ) ;
	}
 
    hij_comp.std_spectral_base() ;   // set the bases for spectral expansions
   
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

// Shift vector 
// -------------

    beta = beta_auto + beta_comp ;

// Coefficients of the 3-metric tilde
// ----------------------------------
 
    Sym_tensor gtilde_con (mp, CON, mp.get_bvect_cart()) ; 
    
     for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
   
	    hij.set(i,j) = hij_auto(i,j) + hij_comp(i,j) ;
	    gtilde_con.set(i,j) = hij(i,j) + flat.con()(i,j) ;
	}

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
    extrinsic_curvature(om) ;

   
// The derived quantities are obsolete
// -----------------------------------

    del_deriv() ;

    delete max_det ;
    delete moy_det ;
    delete min_det ;
}



//----------------------------------//
//	  Version with relaxation   //
//----------------------------------//

void Star_bin::update_metric(const Star_bin& comp,
			     const Star_bin& star_jm1, 
			     double relax, double om) {


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
	logn_comp.import( comp.logn_auto ) ;
	logn_comp.std_spectral_base() ;
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

 
    if ( (comp.lnq_auto).get_etat()  == ETATZERO ) {
      lnq_comp.set_etat_zero() ;
    }
    else{
	lnq_comp.set_etat_qcq() ;
	lnq_comp.import( comp.lnq_auto ) ;
 	lnq_comp.std_spectral_base() ;
   }	

    hij_comp.set_triad(mp.get_bvect_cart()) ;

    Sym_tensor comp_hij(comp.hij_auto) ;
    comp_hij.change_triad(mp_comp.get_bvect_cart()) ;
    comp_hij.change_triad(mp.get_bvect_cart()) ;

    assert ( *(hij_comp.get_triad()) == *(comp_hij.get_triad())) ;


    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
    
	    hij_comp.set(i,j).set_etat_qcq() ;
	    hij_comp.set(i,j).import( (comp_hij)(i,j) ) ;
	}
 
    hij_comp.std_spectral_base() ;
  
// Relaxation on logn_comp, beta_comp, lnq_comp, hij_comp
// ---------------------------------------------------------------
    double relaxjm1 = 1. - relax ; 
    
    logn_comp = relax * logn_comp + relaxjm1 * (star_jm1.logn_comp) ; 
    
    beta_comp = relax * beta_comp + relaxjm1 
	                               * (star_jm1.beta_comp) ; 

    lnq_comp = relax * lnq_comp + relaxjm1 * (star_jm1.lnq_comp) ;

       
    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {

	    hij_comp.set(i,j) = relax * hij_comp(i,j) 
		+ relaxjm1 * (star_jm1.hij_comp)(i,j) ; 
	
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

// Shift vector
// ------------

    beta = beta_auto + beta_comp ;

// Coefficients of the 3-metric tilde
// ----------------------------------
     
    Sym_tensor gtilde_con(mp, CON, mp.get_bvect_cart()) ;

    for(int i=1; i<=3; i++) 
	for(int j=i; j<=3; j++) {
   
	    hij.set(i,j) = hij_auto(i,j) + hij_comp(i,j) ;
	    gtilde_con.set(i,j) = hij(i,j) + flat.con()(i,j) ;
	}

    
    gtilde = gtilde_con ;
    Sym_tensor tens_gamma(gtilde_con / psi4) ;
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
            

    // ... extrinsic curvature (tkij_auto and kcar_auto)
    extrinsic_curvature(om) ;

   
// The derived quantities are obsolete
// -----------------------------------

    del_deriv() ;

    delete max_det ;
    delete moy_det ;
    delete min_det ;

}

}
