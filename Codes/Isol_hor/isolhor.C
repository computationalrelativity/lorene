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
  
    double radius, relax, seuil ;
    fpar >> radius; fpar.ignore(1000, '\n');
    fpar >> relax; fpar.ignore(1000, '\n');
    fpar >> seuil; fpar.ignore(1000, '\n');
    

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
    Mtbl mmr = mm / r ;

    Scalar msr(map) ;
    msr = mmr ;
    
    msr.set_domain(0) = 1. ; // scalar set to 1 in the nucleus  ("inside" the horizon)
    msr.std_spectral_base() ;

    
    // Set up of field Psi 
    // -------------------
    
    Scalar psi_init(map) ; 
    Scalar tmp(map) ; 
    
    psi_init = 0.5*msr + 1. ;
     
    psi_init.std_spectral_base() ;    // sets standard spectral bases

    
    // Set up of the metric tilde
    // --------------------------
    
    Sym_tensor gammat_dd_init = ff.cov() ;

    // Setup of the physical metric

    Metric gamma_init(gammat_dd_init * psi_init* psi_init* psi_init* psi_init);
    Sym_tensor gamma_dd_init (gamma_init.cov()) ;
    

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(map, CON, otriad ) ; 
    beta_init.set_etat_zero() ; 
/*
    tmp_scal = dmsr/(1 + dmsr) ;
    tmp_scal.std_spectral_base() ;
    //    tmp_scal.div_r() ;
    //    tmp_scal = sqrt(tmp_scal) ;
    //    tmp_scal.std_spectral_base() ;
    
    
    beta_init.set(1) = tmp_scal ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;

    beta_init.std_spectral_base() ;
*/      
     
    

    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(map) ; 

    nn_init = - 0.5 * msr + 1. ;

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

    Scalar alpha = isolhor.nn() ;

    Vector beta = isolhor.beta() ;
    
    Scalar psi = isolhor.psi() ;

    
    //    des_profile (psi,  1., 10., 1.57, 0., "Psi", "radial distance") ;
 
    //---------------------------------------------
    //                 Sources
    //---------------------------------------------
 
    tmp_scal.set_etat_zero() ;

    Scalar source_psi = isolhor.source_psi_hor() ;

    Scalar source_nn = isolhor.source_nn_hor(tmp_scal) ;
    
    Vector source_beta = isolhor.source_beta_hor() ;
    
     
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

    
    tmp_scal.set_etat_zero() ;
    tmp_sym.set_etat_zero() ;

    isolhor.init_data_schwar(tmp_sym, tmp_scal,tmp_scal, seuil, relax) ;

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

    cout<<"Todo bien"<<endl ;



    return EXIT_SUCCESS ; 
}


    
    
