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
 * Revision 1.7  2004/11/02 16:27:07  f_limousin
 * Add new parameter ang_vel in function init_data(...).
 *
 * Revision 1.6  2004/10/29 15:41:02  jl_jaramillo
 * ADM angular momentum added
 *
 * Revision 1.5  2004/10/01 16:48:47  f_limousin
 * *** empty log message ***
 *
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
	bornes[l] = pow(2., l-1) * radius ;
	
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

    
    
    // Initialization of the TrK 
    //--------------------------

    // Set up of field Psi 
    // -------------------

    Scalar psi_init(map) ; 
    Scalar tmp(map) ; 
    
    psi_init =  1 + unsr ;
     
    psi_init.std_spectral_base() ;    // sets standard spectral bases

    
    // Set up of the 3-metric tilde
    // ----------------------------
    
    Sym_tensor gammat_dd_init(map, COV, map.get_bvect_spher()) ;
    gammat_dd_init.set(1,1) = 1. ;
    gammat_dd_init.set(2,2) = 1. ;
    gammat_dd_init.set(3,3) = 1. ;
    gammat_dd_init.set(2,1) = 0. ;
    gammat_dd_init.set(3,1) = 0.001 * unsr ;
    gammat_dd_init.set(3,2) = 0. ;
    gammat_dd_init.std_spectral_base() ;

    // Set up of the 3-metric 

    Metric gamma_init(gammat_dd_init * psi_init* psi_init* psi_init* psi_init);
    Sym_tensor gamma_dd_init = gamma_init.cov() ;

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(map, CON, otriad ) ; 
    beta_init.set_etat_zero() ; 

    beta_init.set(1) = 0. ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;
    beta_init.annule_domain(0) ;
    
    beta_init.std_spectral_base() ;
     
    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(map) ; 

    nn_init = (1 - unsr)/(1+unsr) ;

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
    cout << "trk" << endl << norme(trk) << endl ;

    
    // Actual initialization of the 3+1 decomposition
    //-----------------------------------------------

    psi_init = 1. + unsr*unsr ;

    beta_init.set(1) = 0. ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;
    beta_init.annule_domain(0) ;  
    beta_init.std_spectral_base() ;


    nn_init = 1. + unsr*unsr ;
    
    nn_init.std_spectral_base() ;    // sets standard spectral bases
        
    for (int i=1; i<=3; i++) {
      for (int j=1; j<=i; j++) {
	kk_init.set(i,j) = 0.5 * (beta_init.derive_con(gamma_init)(i,j)  
	  + beta_init.derive_con(gamma_init)(j,i)) ;
      }
    }

    kk_init = kk_init/nn_init ;

    gammat_dd_init.set(1,1) = 1. ;
    gammat_dd_init.set(2,2) = 1. ;
    gammat_dd_init.set(3,3) = 1. ;
    gammat_dd_init.set(2,1) = 0. ;
    gammat_dd_init.set(3,1) = 0.001*unsr ;
    gammat_dd_init.set(3,2) = 0. ;
    gammat_dd_init.std_spectral_base() ;



    Isol_hor isolhor(nn_init, beta_init, gamma_dd_init, kk_init, ff, 3) ;
 
 
    //-----------------------------------------
    //          "Call to init_data.C"
    //-----------------------------------------
    
    tmp_scal.set_etat_zero() ;
    tmp_sym.set_etat_zero() ;
    tmp_scal.set_dzpuis(4) ;

    double ang_vel = 0.01 ;
    isolhor.init_data(tmp_sym, trk, tmp_scal, seuil, relax, 
			     niter, ang_vel) ;

    //    FILE* fresu = fopen("resu.d", "w") ; 
    //    isolhor.sauve(fresu) ;
    //    fclose(fresu) ;     
    

    // Test of the constraints
    //------------------------------------

    cout<<"  ----------------------------------------" <<endl ;
    
    isolhor.check_hamiltonian_constraint() ;
    isolhor.check_momentum_constraint() ;

    cout<<"  ----------------------------------------" <<endl ;
 
    // Graphic output of the different fields
    //---------------------------------------

    des_profile(isolhor.nn(), 1.00001, 10, 1., 1., "nn") ;
    
    des_profile(isolhor.psi(), 1.00001, 10, 1., 1., "psi") ;
     
    des_profile(isolhor.beta()(1), 1.00001, 10, 1., 1., "beta_r") ;
    des_profile(isolhor.beta()(3), 1.00001, 10, 1., 1., "beta_phi en 1,1") ;
    des_profile(isolhor.beta()(3), 1.00001, 10, M_PI/2., 0., "beta_phi en pi/2") ;


    
    // Physical parameters of the Black Hole
    //--------------------------------------
    
    cout<<"------------------------------------------------"<<endl;
    cout<<"      Physical parameters of the Black Hole     "<<endl;
    cout<<"------------------------------------------------"<<endl;
    
    double rr_hor =  isolhor.radius_hor() ;
    cout<<"Radius of the horizon = "<< rr_hor <<endl ;
    
    double jj_hor =  isolhor.ang_mom_hor() ;
    cout<<"Angular momentum of the horizon = "<< jj_hor <<endl ; 

    double mm_hor = isolhor.mass_hor() ;
    cout<<"Mass of the horizon = "<< mm_hor <<endl ;  

    double kappa_hor = isolhor.kappa_hor() ;
    cout<<"Surface gravity of the horizon = "<< kappa_hor <<endl ; 

    double omega_hor = isolhor.omega_hor() ;
    cout<<"Orbital velocity of the horizon = "<< omega_hor <<endl ; 


    // Physical parameters of the Bulk
    //--------------------------------
    cout<<""<<endl;
    cout<<"------------------------------------------------"<<endl;
    cout<<"      Physical parameters of the Bulk     "<<endl;
    cout<<"------------------------------------------------"<<endl;
    
    double mm_adm = isolhor.adm_mass() ;
    cout<<"ADM mass= "<< mm_adm <<endl ;  

    double jj_adm = isolhor.ang_mom_adm() ;
    cout<<"ADM angular momentum= "<< jj_adm <<endl ;  

    








    //--------------------------------------
    //--------------------------------------

    cout<<"Tout va bien boudiou / Todo bien!!! (Viva Cai!)"<<endl ;

    return EXIT_SUCCESS ; 
}


    
    
