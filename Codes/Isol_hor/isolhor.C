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
    //      Parameters of the computation
    //======================================================================
  /*
    double pdt = 0.01 ; 
    int jmax = 1000 ; 
    int jstop = jmax ; 
    bool compute_source = true ; 

    double relativistic_init = 0.5;     // 0 = flat space
    double ampli_h_init = 0.001 ;     // 0 = flat space
  */    

    //======================================================================
    //      Construction and initialization of the various objects
    //======================================================================

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int nz = 4 ; 	// Number of domains
    int nr = 17 ; 	// Number of collocation points in r in each domain
    int nt = 9 ; 	// Number of collocation points in theta in each domain
    int np = 8 ; 	// Number of collocation points in phi in each domain
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified
  
    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    //    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 1., 2., 4., __infinity} ; 
    assert( nz == 4 ) ;  // since the above array described only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    //    cout << map << endl ;  

   // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    const Coord& cost = map.cost ;  // cos(theta) field
    const Coord& sint = map.sint ;  // sin(theta) field

 
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

    double mm = 1. ;
     
    // Key function
    Mtbl mmr = 2. * mm / r ;

    Scalar Mr(map) ;
    Mr = mmr ;
    
    Mr.set_domain(0) = 1. ; // scalar set to 1 in the nucleus  ("inside" the horizon)
    Mr.std_spectral_base() ;

    
    // Set up of the metric
    // --------------------
    
    Sym_tensor gamma_dd_init = ff.cov() ;

    gamma_dd_init.set(1,1) = 1. + Mr ;

    gamma_dd_init.std_spectral_base() ;

    Metric gamma_init(gamma_dd_init) ; 

    /* 
    // Set up of field Psi 
    // -------------------
    
    Scalar psi_init(map) ; 
    Scalar tmp(map) ; 
    
    psi_init = 1. ;
     
    psi_init.std_spectral_base() ;    // sets standard spectral bases

    */

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(map, CON, otriad ) ; 
    //    beta_init.set_etat_zero() ; 

    tmp_scal = Mr/(1 + Mr) ;
    tmp_scal.std_spectral_base() ;
    //    tmp_scal.div_r() ;
    //    tmp_scal = sqrt(tmp_scal) ;
    //    tmp_scal.std_spectral_base() ;
    
    
    beta_init.set(1) = tmp_scal ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;

    beta_init.std_spectral_base() ;
        
    /*  
    for (int i=1 ; i<4 ; i++) {
      cout<<"component: ("<<i<<") of beta"<<endl ;      
	des_profile(beta_init(i), 1.000001, 10., M_PI/2., 0., "vector comp beta_init ", "radial distance") ;
    }

    */
    
    

    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(map) ; 

    nn_init = 1./sqrt(1. + Mr) ;

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
     
    Isol_hor isolhor(nn_init, beta_init, gamma_dd_init, kk_init, ff, 3) ;

  
    //---------------------------------------------
    //         Conformal/ Physical param
    //---------------------------------------------

    Scalar trk = isolhor.trk() ;
    trk.dec_dzpuis(2) ;

    //    des_profile (trk,  1., 10., 1.57, 0., "K", "radial distance") ;

    Scalar alpha = isolhor.nn() ;

    Vector beta = isolhor.beta() ;
    
    Scalar psi = isolhor.psi() ;

    
    //    des_profile (psi,  1., 10., 1.57, 0., "Psi", "radial distance") ;
 
    //---------------------------------------------
    //                 Sources
    //---------------------------------------------
 
    tmp_scal.set_etat_zero() ;

    Scalar source_psi = isolhor.source_psi_hor_ih() ;

    Scalar source_nn = isolhor.source_nn_hor_ih(tmp_scal) ;
    
    Vector source_beta = isolhor.source_beta_hor_ih() ;
    
  



    

    //---------------------------------------------
    //       Test a las ecuaciones de Einstein
    //---------------------------------------------

    Tbl ham = isolhor.check_hamiltonian_constraint () ;


    Tbl mom = isolhor.check_momentum_constraint () ;
 	
    
    /*

    //---------------------------------------------
    //         Boundary condition for Psi
    //---------------------------------------------

 
    Scalar div_beta = contract(beta.derive_cov(ff), 0, 1) ;
    div_beta.dec_dzpuis(2) ;
    
    des_profile (div_beta,  1., 10., 1.57, 0., "div_beta", "radial distance") ;
 

    Scalar diff = div_beta - alpha * trk ;

    des_profile (diff,  1., 10., 1.57, 0., "diff", "radial distance") ;
 

    
    Scalar bc_psi = contract(beta, 0, psi.derive_cov(ff), 0) ;
    bc_psi.dec_dzpuis(2) ;
      
    bc_psi = bc_psi + psi * diff/ 6. ;

    des_profile (bc_psi,  1., 10., 1.57, 0., "bc_psi", "radial distance") ;
 

    cout<<"bc_psi:   "<<bc_psi<<endl ;
    
    */

    /*    
    //---------------------------------------------
    //         Boundary condition for lapse
    //---------------------------------------------

    

    //  Radial vector
    //---------------
    
    Vector radial(map, CON, otriad) ;

    for (int i=1 ; i<4 ; i++){
      radial.set(i) = kerr_sch.gam_uu()(1,i)/sqrt(kerr_sch.gam_uu()(1,1)) ;
    }
    radial.std_spectral_base() ;

    */

    //  Extrinsic curvature
    //---------------------

    //    Sym_tensor k_dd = kerr_sch.k_dd()    ;

    
    // Surface gravity
    //----------------

    /*

    Scalar kappa(map) ;
    kappa = 0.5 / r ;
    
    kappa.set_domain(0) = 1. ;

    // Test to WIH equation
    //---------------------

    Scalar wih(map) ;

    wih = contract(radial, 0, alpha.derive_cov(gamma_init), 0) 
		   - contract(radial, 0, contract(radial, 0, k_dd,0), 0) * alpha ;

    wih.dec_dzpuis(2) ;
    
    wih = wih - kappa ;

    
       
    des_profile (wih,  1.0000000000001, 10., 1.57, 0., "WIH", "radial distance") ;
 

    */
    







    //-----------------------------------------
    //          "Physical Parameters"
    //-----------------------------------------

    
    double rr_hor =  isolhor.radius_hor_ih() ;
    cout<<"Radius of the horizon = "<< rr_hor <<endl ;
    
    Vector rad = isolhor.radial_vect_hor_ih() ;
 
    double jj_hor =  isolhor.ang_mom_hor_ih() ;
    cout<<"Angular momentum of the horizon = "<< jj_hor <<endl ; 

    double mm_hor = isolhor.mass_hor_ih() ;
    cout<<"Mass of the horizon = "<< mm_hor <<endl ;  

    double kappa_hor = isolhor.kappa_hor_ih() ;
    cout<<"Surface gravity of the horizon = "<< kappa_hor <<endl ; 

    double omega_hor = isolhor.omega_hor_ih() ;
    cout<<"Orbital velocity of the horizon = "<< omega_hor <<endl ; 
    

 
    /* 

    //-----------------------------------------
    //          "Boundary conditions"
    //-----------------------------------------

    //    Valeur bound_psi_Dir = horizon.boundary_psi_Dir() ;
    //    Valeur bound_psi_Neu = horizon.boundary_psi_Neu() ;
    Vector beta_cart = horizon.beta_bound_cart() ;
    Valeur bound_beta_r = horizon.boundary_beta_r() ;
    Valeur bound_beta_theta = horizon.boundary_beta_theta() ;
    Valeur bound_beta_phi = horizon.boundary_beta_phi() ;
    Valeur bound_beta_x = horizon.boundary_beta_x() ;
    Valeur bound_beta_y = horizon.boundary_beta_y() ;
    Valeur bound_beta_z = horizon.boundary_beta_z() ;
    //    Valeur bound_nn_Dir = horizon.boundary_nn_Dir_kk() ;
    Valeur bound_nn_Neu = horizon.boundary_nn_Neu_kk() ;
    Valeur bound_psi_Neu = horizon.boundary_psi_Neu() ;

    */
    

 
    //   horizon.init_data_rot(tmp_sym, tmp_scal,tmp_scal) ;



    //--------------------------------------
    //--------------------------------------

    cout<<"Todo bien"<<endl ;



    return EXIT_SUCCESS ; 
}


    
    
