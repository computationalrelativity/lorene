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
 * Revision 1.10  2004/11/05 10:15:22  f_limousin
 * The angular velocity is now a parameter (in par_hor.d).
 *
 * Revision 1.8  2004/11/02 17:42:00  f_limousin
 * New method sauve(...) to save in a binary file.
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
  
    double radius, relax, seuil, niter, ang_vel ;
    fpar >> radius; fpar.ignore(1000, '\n');
    fpar >> relax; fpar.ignore(1000, '\n');
    fpar >> seuil; fpar.ignore(1000, '\n');
    fpar >> niter; fpar.ignore(1000, '\n');
    fpar >> ang_vel; fpar.ignore(1000, '\n');
    

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

    // Key function
    Mtbl usr = 1 / r ;
    Scalar unsr(map) ;
    unsr = usr ;
    
    unsr.set_domain(0) = 1 ; // scalar set to 1 in the nucleus 
    unsr.std_spectral_base() ;

    Mtbl expr = exp(-r) ;
    Scalar expmr(map) ;
    expmr = expr ;

    
    // Physical Parameters
    //--------------------
    
    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(map) ; 
    nn_init = 1 + unsr ;
    nn_init.std_spectral_base() ;    // sets standard spectral bases

    // Set up of field Psi 
    // -------------------

    Scalar psi_init(map) ; 
    psi_init =  1 + unsr*unsr ;
    psi_init.std_spectral_base() ;    // sets standard spectral bases

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(map, CON, otriad ) ; 
    beta_init.set_etat_zero() ; 

    beta_init.set(1) = 0. ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;
    beta_init.annule_domain(0) ;
    
    beta_init.std_spectral_base() ;
        

    // TrK, TrK_point
    // --------------

    Scalar trK (map) ;
    trK = 0. ;
    trK.std_spectral_base() ;

    Scalar trK_point (map) ;
    trK_point = 0. ;
    trK_point.std_spectral_base() ;
	
    // gamt, gamt_point
    // ----------------

    Sym_tensor gamt(map, COV, map.get_bvect_spher()) ;
    gamt.set(1,1) = 1. ;
    gamt.set(2,2) = 1. ;
    gamt.set(3,3) = 1. ;
    gamt.set(2,1) = 0. ;
    gamt.set(3,1) = 0. ;
    gamt.set(3,2) = 0. ;
    gamt.std_spectral_base() ;
    Metric met_gamt (gamt) ;     

    Sym_tensor gamt_point(map, CON, map.get_bvect_spher()) ;
    gamt_point.set_etat_zero() ;

    // Set up of extrinsic curvature
    // -----------------------------
    
    Metric met_gam(psi_init*psi_init*psi_init*psi_init*gamt) ;
    Sym_tensor kk_init (map,  CON, map.get_bvect_spher()) ;
    for (int i=1; i<=3; i++) {
      for (int j=1; j<=i; j++) {
	kk_init.set(i,j) = 0.5 * (beta_init.derive_con(met_gam)(i,j)  
	  + beta_init.derive_con(met_gam)(j,i)) ;
      }
    }
    kk_init = kk_init/nn_init ;

    Sym_tensor aa_init (map, CON, map.get_bvect_spher()) ;
    aa_init = psi_init*psi_init*psi_init*psi_init*kk_init
	- 1./3. * trK * met_gamt.con() ;
   

    //-------------------------------------
    //     Construction of the space-time
    //-------------------------------------

    Isol_hor isolhor(nn_init, beta_init, psi_init, aa_init, met_gamt,
		     gamt_point, trK, trK_point, ff, 3) ;
 
 
    //-----------------------------------------
    //          "Call to init_data.C"
    //-----------------------------------------
    
    tmp_scal.set_etat_zero() ;
    tmp_sym.set_etat_zero() ;
    tmp_scal.set_dzpuis(4) ;

    isolhor.init_data(seuil, relax, niter, ang_vel) ;

    // Save in a file
    // --------------
    
    FILE* fresu = fopen("resu.d", "w") ;
    isolhor.sauve(fresu, true) ;
    fclose(fresu) ;     
    
    // Test of the constraints
    //------------------------

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


    
    
