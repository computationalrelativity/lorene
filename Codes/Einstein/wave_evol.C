/*
 *  Main code for time evolution of a wave packet within maximal slicing
 *  + Dirac gauge.
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon & Jerome Novak
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

char wave_evol_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/04/30 10:53:32  e_gourgoulhon
 * Added resolution of elliptic Einstein equations (new methods
 * Tslice_dirac_max::solve_*) for tests at the end.
 *
 * Revision 1.4  2004/04/29 17:13:08  e_gourgoulhon
 * New argument pdt to Time_slice_conf::initial_data_cts.
 *
 * Revision 1.3  2004/04/08 16:47:09  e_gourgoulhon
 * Many changes.
 *
 * Revision 1.2  2004/04/07 07:59:22  e_gourgoulhon
 * Added check of constraints at the end.
 *
 * Revision 1.1  2004/04/05 21:26:25  e_gourgoulhon
 * First version (not ready yet !).
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
#include "time_slice.h"
#include "metric.h"
#include "evolution.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

int main() {

    //======================================================================
    //      Parameters of the computation
    //======================================================================

    double pdt = 0.01 ; 
    int jmax = 1000 ; 
    int jstop = jmax ; 

    double ampli_h_init = 1. ;     // 0 = flat space
        

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
	
    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------

    // radial boundaries of each domain:
    double r_limits[] = {0., 1., 2., 4., __infinity} ; 
    assert( nz == 4 ) ;  // since the above array described only 3 domains
  
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << map << endl ;  
    
    // Flat metric f
    // -------------

    const Metric_flat& ff = map.flat_met_spher() ; 
    
    // Triad orthonormal with respect to the flat metric f
    // ----------------------------------------------------

    const Base_vect_spher& otriad = map.get_bvect_spher() ;
    
    
    // Construction of a time slice with maximal slicing and Dirac gauge
    // -----------------------------------------------------------------

    Tslice_dirac_max sigmat(map, otriad, ff) ;  

    // Set up of potentials khi and mu
    // -------------------------------
    
    const Coord& x = map.x ; 
    const Coord& y = map.y ; 
    // const Coord& z = map.z ; 
    const Coord& r = map.r ; 
    //const Coord& cost = map.cost ; 
    //const Coord& sint = map.sint ; 
    //const Coord& cosp = map.cosp ; 
    // const Coord& sinp = map.sinp ; 
    
    Scalar khi_init(map) ; 

    khi_init = ampli_h_init * exp( - r*r ) * x*y ;
    khi_init.set_outer_boundary(nz-1, 0.) ; 
    
    //khi_init = ampli_h_init * (3*cost*cost-1) / 
    //   ( (r*r + 1./(r*r)) * sqrt(1.+r*r) ) ; 
    khi_init.std_spectral_base() ; 
    
    khi_init.spectral_display("khi_init") ;   
    //if (khi_init.get_etat() == ETATQCQ) 
    //    des_meridian(khi_init, 0., 5., "khi_init", 1) ; 

    //## khi_init.smooth_decay(2, 1) ; 
    
    Scalar mu_init(map) ; 
    mu_init = 0. * ampli_h_init / (1+r*r*r*r*r*r) ; 
    mu_init.std_spectral_base() ; 
    mu_init.mult_r() ; 
    mu_init.mult_r() ; 
    mu_init.mult_r() ; 
    mu_init.mult_cost() ; 
    
    mu_init.spectral_display("mu_init") ;   
    if (mu_init.get_etat() == ETATQCQ) 
        des_meridian(mu_init, 0., 5., "mu_init", 1) ; 
    
    
                  
    // The potentials khi and mu are used to construct h^{ij}:
    // ------------------------------------------------------
        
    sigmat.set_khi_mu(khi_init, mu_init) ; // the trace h = f_{ij} h^{ij]
                                           // is computed to ensure
                                           // det tgam_{ij} = f

    // Resolution of the initial data equations within 
    // the conformal thin sandwich framework
    // ------------------------------------------------

    // u^{ij} = d/dt h^{ij}
    Sym_tensor_trans uu_init(map, otriad, ff) ;  
    uu_init.set_etat_zero() ; 
    uu_init = - 0.5 * ( sigmat.hh() 
        - 0.33333333333333333 * sigmat.hh().trace(sigmat.tgam())
            * sigmat.tgam().con() ) ;
    uu_init.inc_dzpuis(2) ;  
    
    // tr K = K
    Scalar tmp(map) ; 
    tmp.set_etat_zero() ; 
    
    sigmat.initial_data_cts(uu_init, tmp, tmp, pdt, 1.e-10) ;
        
    cout << "sigmat : " << sigmat << endl ;  
    
    cout << "Test upon khi : difference between khi and khi_init : " << endl ; 
    Scalar diff_khi = sigmat.khi() - khi_init ;
    maxabs(diff_khi, "diff_khi") ; 
    arrete() ; 
    
    diff_khi.spectral_display("diff_khi", 1.e-14) ; 
    
    // des_meridian(sigmat.hh(), 0., 5., "h") ; 
    
    // sigmat.trh().visu_section ('x', 0., -4., 4., -4., 4., "h in x=0 plane", "h_x") ;
    
    // sigmat.trh().visu_section ('z', 0., -4., 4., -4., 4., "h in z=0 plane", "h_z") ;
    
    // sigmat.psi().visu_section ('z', 0., -4., 4., -4., 4., "Psi in z=0 plane",
    //                          "psi_z") ;
    
    // sigmat.trh().visu_box(-4., 4.,-4., 4.,-4., 4., "h") ; 
    
    // Check of constraints:
    sigmat.check_hamiltonian_constraint() ;    
    sigmat.check_momentum_constraint() ; 
    
    // Extra check
    Scalar diffr = sigmat.gam().ricci_scal() 
            - ( sigmat.tgam().ricci_scal()
     - 8.* sigmat.psi().derive_con(sigmat.tgam()).divergence(sigmat.tgam())
        / sigmat.psi()
            ) / sigmat.psi4() ;  
    maxabs(diffr, 
    "Error in the relation (involving psi) between the Ricci scalars of gam and tgam") ; 
    
    maxabs(sigmat.gam().cov()(1,1) / sigmat.tgam().cov()(1,1) -
            sigmat.psi4(), "Difference between the conformal factor and psi4") ; 
    
    arrete() ; 

    Scalar nn_new = sigmat.solve_n() ; 
    maxabs(nn_new - sigmat.nn(), "Difference between N and N_new") ;     
    
    Scalar qq_new = sigmat.solve_q() ; 
    maxabs(qq_new - sigmat.qq(), "Difference between Q and Q_new") ;     
    
    Vector beta_new = sigmat.solve_beta() ; 
    maxabs(beta_new - sigmat.beta(), "Difference between beta and beta_new") ;     
    
    
    return EXIT_SUCCESS ; 
}

