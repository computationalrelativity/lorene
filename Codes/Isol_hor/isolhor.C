/*
 *  Main code for Isolated Horizon in arbitrary gauge
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
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

char isolhor_C[] = "$Header$" ;

/* 
 * $Id$
 * $Log$
 * Revision 1.4  2004/09/28 15:57:45  f_limousin
 * Add the 2 lines  $Id and $Log to see the comments
 *
 *
 * Revision 1.1  2004/02/18 19:16:28  jl_jaramillo
 * First version: c'est loin d'etre pret tout ca !!!
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Lorene headers
#include "tenseur.h"
#include "metric.h"
#include "evolution.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "time_slice.h"
#include "isol_hor.h"


int main() {

    //======================================================================
    //      Construction and initialization of the various objects
    //======================================================================

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------

    int nz, nt, np, nr ;

    ifstream fpar("par_hor.d") ;
    fpar.ignore(1000, '\n') ;
    fpar.ignore(1000, '\n') ;
    fpar >> nz; fpar.ignore(1000, '\n');
    fpar >> nt; fpar.ignore(1000, '\n');
    fpar >> np; fpar.ignore(1000, '\n');
    fpar >> nr; fpar.ignore(1000, '\n');

    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified
  
    double radius, relax, seuil, niter ;
    fpar >> radius; fpar.ignore(1000, '\n');
    fpar >> relax; fpar.ignore(1000, '\n');
    fpar >> seuil; fpar.ignore(1000, '\n');
    fpar >> niter; fpar.ignore(1000, '\n');
    

    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    //    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid -->physical space(Lorene class Map_af)
    // -----------------------------------------------------------------------

    // radial boundaries of each domain:
    double* bornes = new double[nz+1];
    bornes[0] = 0. ;
    for (int l=1; l<nz; l++) 
	bornes[l] = pow(2, l-1) * radius ;
	
    bornes[nz] = __infinity ; 
  
    Map_af map(mgrid, bornes) ;   // Mapping construction
  	
    //    cout << map << endl ;  

   // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    const Coord& cost = map.cost ;  // cos(theta) field
    const Coord& sint = map.sint ;  // sin(theta) field
    const Coord& cosp = map.cosp ;  // cos(phi) field
    const Coord& sinp = map.sinp ;  // sin(phi) field

 
    // Flat metric f
    // -------------

    const Metric_flat& ff = map.flat_met_spher() ; 
    
    // Triad orthonormal with respect to the flat metric f
    // ----------------------------------------------------

    const Base_vect_spher& otriad = map.get_bvect_spher() ;
    

    // Working stuff
    // -------------
    
    Scalar tmp_scal(map) ;
    Vector tmp_vect(map, CON, otriad) ;
    Sym_tensor tmp_sym(map, CON, otriad) ; 

   
    // Physical Parameters
    //--------------------

     // Key function
    Mtbl usr = 1 / r ;

    Scalar unsr(map) ;
    unsr = usr ;
    
    unsr.set_domain(0) = 1 ; // scalar set to 1 in the nucleus  ("inside" the horizon)
    unsr.std_spectral_base() ;

     Mtbl expr = exp(-r) ;

    Scalar expmr(map) ;
    expmr = expr ;

    
    // Set up of field Psi 
    // -------------------

    Scalar psi_init(map) ; 
    Scalar tmp(map) ; 
    
    psi_init =  pow(1+unsr, 1./4.) ;
     
    psi_init.std_spectral_base() ;    // sets standard spectral bases

    
    // Set up of the 3-metric tilde
    // ----------------------------
    
    Sym_tensor gammat_dd_init(map, COV, map.get_bvect_spher()) ;
    gammat_dd_init.set(1,1) = 1. ;
    gammat_dd_init.set(2,2) = 1/(1+unsr) ;
    gammat_dd_init.set(3,3) = 1/(1+unsr) ;
    gammat_dd_init.set(2,1) = 0. ;
    gammat_dd_init.set(3,1) = 0. ;
    gammat_dd_init.set(3,2) = 0. ;


    // Set up of the 3-metric 

    Metric gamma_init(gammat_dd_init * psi_init* psi_init* psi_init* psi_init);
    Sym_tensor gamma_dd_init = gamma_init.cov() ;

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(map, CON, otriad ) ; 
    beta_init.set_etat_zero() ; 

    tmp_scal = unsr / (1+unsr) ;
    tmp_scal.std_spectral_base() ;
    
    beta_init.set(1) = tmp_scal ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;
    beta_init.annule_domain(0) ;
    
    beta_init.std_spectral_base() ;
     
    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(map) ; 

    nn_init = pow(1+unsr, -0.5) ;

    nn_init.std_spectral_base() ;    // sets standard spectral bases

  
    // Set up of the conformal extrinsic curvature KK
    //-----------------------------------------------
   	
    Sym_tensor kk_init(map, CON, otriad) ;   // A^{ij}
    
    for (int i=1; i<=3; i++) {
      for (int j=1; j<=i; j++) {
	kk_init.set(i,j) = 0.5 * (beta_init.derive_con(gamma_init)(i,j)  
	  + beta_init.derive_con(gamma_init)(j,i)) ;
      }
    }


    kk_init = kk_init/nn_init ;

    //-------------------------------------
    //     Construction of the space-time
    //-------------------------------------
     
    Isol_hor isolhor_tmp(nn_init, beta_init, gamma_dd_init, kk_init, ff, 3) ;

    Scalar trk = isolhor_tmp.trk() ;


    beta_init.set(1) = 1 / (1 + 4*unsr) ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;
    beta_init.annule_domain(0) ;
    
    beta_init.std_spectral_base() ;

/*
    nn_init = pow(1+unsr, -0.5) ;

    nn_init.std_spectral_base() ;    // sets standard spectral bases
*/


    Isol_hor isolhor(nn_init, beta_init, gamma_dd_init, kk_init, ff, 3) ;



/*
    Vector source (map, CON, map.get_bvect_spher()) ;

    Vector beta_j (map, CON, map.get_bvect_spher()) ;
    beta_j = beta_init ;

    Vector beta_old (map, CON, map.get_bvect_spher()) ;
    beta_old = beta_init ;
    
    Vector bound (map, CON, map.get_bvect_spher()) ;
    bound.set(1) = unsr ;
    bound.set(2) = 0. ;
    bound.set(3) = 0. ;
    bound.std_spectral_base() ;
    bound.change_triad(map.get_bvect_cart()) ;
    
    
    Valeur boundary_x (map.get_mg()->get_angu() )  ;
    Valeur boundary_y (map.get_mg()->get_angu() )  ;
    Valeur boundary_z (map.get_mg()->get_angu() )  ;
    
    boundary_x = 1. ;
    boundary_y = 1. ;
    boundary_z = 1. ;

    int nnp = map.get_mg()->get_np(1) ;
    int nnt = map.get_mg()->get_nt(1) ;

    for (int k=0 ; k<nnp ; k++)
	for (int j=0 ; j<nnt ; j++){
	    boundary_x.set(0, k, j, 0) = bound(1).val_grid_point(1, k, j, 0) ;
	    boundary_y.set(0, k, j, 0) = bound(2).val_grid_point(1, k, j, 0) ;
	    boundary_z.set(0, k, j, 0) = bound(3).val_grid_point(1, k, j, 0) ;
	}
    boundary_x.set_base(bound(1).get_spectral_va().get_base()) ; 
    boundary_y.set_base(bound(2).get_spectral_va().get_base()) ; 
    boundary_z.set_base(bound(3).get_spectral_va().get_base()) ; 

  
    for (int i=0; i<niter; i++) {
   
	beta_old = beta_j ;
	des_profile(beta_j(1), 1.00001, 10, 0., 0., "beta (1)") ;
	
	double lambda = 0. ;
	source = 0.66666666666666666 * (isolhor.trk()).derive_con(ff) ;
	source += -(1./3.- lambda) * beta_j.divergence(ff).derive_con(ff) ;
	Vector source2 = beta_init.derive_con(ff).divergence(ff) ;
	cout << "source" << endl << norme(source(1)) << endl ;
	cout << "source2" << endl << norme(source2(1)) << endl ;
	cout << "source2 - source" << endl << norme(source2(1)-source(1)) << endl ;	

	source.inc_dzpuis() ;
	
	double precision = 1e-8 ;
	poisson_vect_boundary(lambda, source, beta_j, boundary_x, boundary_y, 
			     boundary_z, 0, precision, 20) ;


//	beta_j = relax * beta_j + (1 - relax) * beta_old ;

	source.dec_dzpuis() ;
   	maxabs(beta_j.derive_con(ff).divergence(ff) 
	       + lambda * beta_j.divergence(ff).derive_con(ff) - source,
	     "Absolute error in the resolution of the equation for beta") ;  
	cout << endl ;
	source.inc_dzpuis() ;

        double diff_beta = max( maxabs(beta_j - beta_old) ) ; 
	cout << "  diff_beta = " << diff_beta << endl ;
	if ( (diff_beta < seuil) )
	break ; 
   }
 
    Scalar beta_ana = unsr ;
    beta_ana = pow(beta_ana, 0.5) ;
    beta_ana.set_spectral_va().set_base(isolhor.beta()(1).get_spectral_va().get_base()) ;

    Scalar diff_beta = beta_j(1) - beta_ana ;

     des_profile(beta_j(1), 1., 10, 0., 0., "beta") ;
     des_profile(beta_ana, 1., 10, 0., 0., "beta_ana") ;
     des_profile(diff_beta, 1., 10, 0., 0., "diff_beta") ;

    abort() ;
*/

    //-----------------------------------------
    //          "Physical Parameters"
    //-----------------------------------------

    
    tmp_scal.set_etat_zero() ;
    tmp_sym.set_etat_zero() ;
    tmp_scal.set_dzpuis(4) ;

    isolhor.init_data_schwar(tmp_sym, trk, tmp_scal, seuil, relax, 
			     niter) ;

    isolhor.check_hamiltonian_constraint() ;
    isolhor.check_momentum_constraint() ;
 

    des_profile(isolhor.nn(), 1., 10, 0., 0., "nn") ;
    des_profile(pow(1+unsr, -0.5), 1., 10, 0., 0., "nn ana") ;
    des_profile(pow(1+unsr, -0.5)-isolhor.nn(), 1., 10, 0., 0., "diff nn") ;
    des_profile(isolhor.psi(), 1., 10, 0., 0., "psi") ;
    des_profile(isolhor.beta()(1), 1., 10, 0., 0., "beta") ;
    des_profile(unsr/(1+unsr), 1., 10, 0., 0., "beta ana") ;
    des_profile(unsr/(1+unsr) - isolhor.beta()(1), 1., 10, 0., 0., "diff beta") ;



/*
    des_profile(isolhor.nn(), 1., 10, 0., 0., "nn") ;
    des_profile(isolhor.psi(), 1., 10, 0., 0., "psi") ;


    Mtbl nn_anaa = (1 - 1 / r) / (1 + 1 / r) ;
    Scalar nn_ana (map) ;
    nn_ana = nn_anaa ;
    nn_ana.std_spectral_base() ;
    nn_ana.set_domain(0) = 0 ;

    Mtbl psi_anaa (1 + 1 / r) ;
    Scalar psi_ana (map) ;
    psi_ana = psi_anaa ;
    psi_ana.std_spectral_base() ;
    psi_ana.set_domain(0) = 0 ;

    Scalar diff_nn = isolhor.nn() - nn_ana ;
    Scalar diff_psi = isolhor.psi() - psi_ana ; 

    des_profile(diff_nn, 1.00001, 10, 0., 0., "diff nn") ;
    des_profile(diff_psi, 1.00001, 10, 0., 0., "diff psi") ;
*/

/*
    Scalar beta_ana = unsr ;
    beta_ana = pow(beta_ana, 0.5) ;
    beta_ana.set_spectral_va().set_base(isolhor.beta()(1).get_spectral_va().get_base()) ;

    Scalar diff_beta = isolhor.beta()(1) - beta_ana ;

     des_profile(isolhor.beta()(1), 1., 10, 0., 0., "beta") ;
     des_profile(beta_ana, 1., 10, 0., 0., "beta_ana") ;
     des_profile(diff_beta, 1., 10, 0., 0., "diff_beta") ;
*/

    
    double rr_hor =  isolhor.radius_hor() ;
    cout<<"Radius of the horizon = "<< rr_hor <<endl ;
    
    Vector rad = isolhor.radial_vect_hor() ;
 
    double jj_hor =  isolhor.ang_mom_hor() ;
    cout<<"Angular momentum of the horizon = "<< jj_hor <<endl ; 

    double mm_hor = isolhor.mass_hor() ;
    cout<<"Mass of the horizon = "<< mm_hor <<endl ;  

    double kappa_hor = isolhor.kappa_hor() ;
    cout<<"Surface gravity of the horizon = "<< kappa_hor <<endl ; 

    double omega_hor = isolhor.omega_hor() ;
    cout<<"Orbital velocity of the horizon = "<< omega_hor <<endl ; 


    //--------------------------------------
    //--------------------------------------

    cout<<"Tout va bien !!!"<<endl ;



    return EXIT_SUCCESS ; 
}


    
    
