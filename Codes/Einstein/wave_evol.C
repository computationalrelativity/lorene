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
 * Revision 1.8  2004/05/17 12:59:55  e_gourgoulhon
 * Parameters of the computation are now read in file par_wave_evol.d.
 *
 * Revision 1.7  2004/05/05 14:51:48  e_gourgoulhon
 * Introduced parameters ampli_init_khi and ampli_init_mu.
 * Added some checks regarding khi and mu.
 *
 * Revision 1.6  2004/05/03 14:50:38  e_gourgoulhon
 * First full version (time evolution).
 *
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
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

int main() {

    //======================================================================
    //     Reading the parameters of the computation
    //======================================================================

    ifstream fpar("par_wave_evol.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem with opening the file par_wave_evol.d ! " << endl ;
        abort() ;
    }
    
    // Reading of physical paramaters
    double ampli_init_khi, ampli_init_mu ;
    fpar.ignore(1000,'\n') ;    // skip title
    fpar >> ampli_init_khi ; fpar.ignore(1000,'\n') ;
    fpar >> ampli_init_mu ; fpar.ignore(1000,'\n') ;
    
    // Reading of computational paramaters
    int nb_time_steps, niter_elliptic, graph, graph_init ;
    double pdt, relax_elliptic ; 
    fpar.ignore(1000,'\n') ;    // skip title
    fpar >> pdt ; fpar.ignore(1000,'\n') ;
    fpar >> nb_time_steps ; fpar.ignore(1000,'\n') ;
    fpar >> niter_elliptic ; fpar.ignore(1000,'\n') ;
    fpar >> relax_elliptic ; fpar.ignore(1000,'\n') ;
    fpar >> graph ; fpar.ignore(1000,'\n') ;
    fpar >> graph_init ; fpar.ignore(1000,'\n') ;

    char graph_device[40] ;
    if (graph == 0) strcpy(graph_device, "/n") ;
    else if (graph == 1) strcpy(graph_device, "/xwin") ;
        else if (graph == 2) strcpy(graph_device, "?") ;
            else {
                cerr << 
                "Unexpected value of input parameter graph: graph = " << graph
                << " !" << endl ; 
                abort() ; 
            }   
    
    char graph_device_init[40] ;
    if (graph_init == 0) strcpy(graph_device_init, "/n") ;
    else if (graph_init == 1) strcpy(graph_device_init, "/xwin") ;
        else if (graph_init == 2) strcpy(graph_device_init, "?") ;
            else {
                cerr << 
                "Unexpected value of input parameter graph: graph_init = " 
                << graph_init << " !" << endl ; 
                abort() ; 
            }   
    

    // Reading of multi-domain grid parameters
    int symmetry_phi0, nz, nr, nt, np ;
    fpar.ignore(1000,'\n') ;    // skip title
    fpar >> symmetry_phi0 ; fpar.ignore(1000,'\n') ;
    fpar >> nz ; fpar.ignore(1000,'\n') ;
    fpar >> nr ; fpar.ignore(1000,'\n') ;
    fpar >> nt ; fpar.ignore(1000,'\n') ;
    fpar >> np ; fpar.ignore(1000,'\n') ;
    fpar.ignore(1000,'\n') ;    // skip title
    double* r_limits = new double[nz+1] ; 
    for (int l=0; l<nz; l++) {
        fpar >> r_limits[l] ; fpar.ignore(1000,'\n') ;
    }
    r_limits[nz] = __infinity ; 
    
    // fpar >> ; fpar.ignore(1000,'\n') ;
    fpar.close() ; 

    cout << "Physical parameters: \n"
         << "-------------------  \n" ; 
    cout << "   ampli_init_khi = " <<  ampli_init_khi << endl ;        
    cout << "   ampli_init_mu = " <<  ampli_init_mu << endl ;        
    cout << "Computational parameters: \n"
         << "------------------------ \n" ; 
    cout << "   pdt = " << pdt << endl ; 
    cout << "   nb_time_steps = " << nb_time_steps << endl ; 
    cout << "   niter_elliptic = " << niter_elliptic << endl ; 
    cout << "   relax_elliptic = " << relax_elliptic << endl ; 
    cout << "   graph_device = " << graph_device << endl ;
    cout << "   graph_device_init = " << graph_device_init << endl ;

    //======================================================================
    //      Construction and initialization of the various objects
    //======================================================================

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = (symmetry_phi0 == 1) ? SYM : NONSYM ; //  symmetry in phi
    bool compact = true ; // external domain is compactified
  
    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    cout << "Computational grid :\n" 
         << "------------------ \n" 
         << "  " << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------
    
    Map_af map(mgrid, r_limits) ;   // Mapping construction
  	
    cout << "Mapping computational grid --> physical space :\n" 
         << "---------------------------------------------\n" 
         << "  " << map << endl ;  
    
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
    const Coord& r = map.r ; 
    
    Scalar khi_init(map) ; 

    khi_init = ampli_init_khi * exp( - r*r ) * x*y ;
    khi_init.set_outer_boundary(nz-1, 0.) ;     // zero at spatial infinity
    
    khi_init.std_spectral_base() ; 
    
    //## khi_init.smooth_decay(2, 1) ; 

    khi_init.spectral_display("khi_init") ;   
    if (khi_init.get_etat() == ETATQCQ) 
        des_meridian(khi_init, 0., 5., "khi_init", 1) ; 
    
    Scalar mu_init(map) ; 
    mu_init = ampli_init_mu * x*y* exp( - r*r ) ;
    mu_init.set_outer_boundary(nz-1, 0.) ; 
    mu_init.std_spectral_base() ; 
    mu_init.mult_r() ; 
    mu_init.set_outer_boundary(nz-1, 0.) ; 
    mu_init.mult_cost() ; 
    
    mu_init.spectral_display("mu_init") ;   
    // if (mu_init.get_etat() == ETATQCQ) 
    //    des_meridian(mu_init, 0., 5., "mu_init", 1) ; 
    
    
                  
    // The potentials khi and mu are used to construct h^{ij}:
    // ------------------------------------------------------
        
    sigmat.set_khi_mu(khi_init, mu_init) ; // the trace h = f_{ij} h^{ij]
                                           // is computed to ensure
                                           // det tgam_{ij} = f
    
    // sigmat.hh().transverse(ff).tt_part().khi().spectral_display("khi") ; 
    // sigmat.hh().transverse(ff).tt_part().eta().spectral_display("eta") ; 
    // sigmat.hh().transverse(ff).tt_part().mu().spectral_display("mu") ; 
                                           
    //======================================================================
    //      Resolution of the initial data equations within 
    //      the conformal thin sandwich framework
    //======================================================================

    // u^{ij} = d/dt h^{ij}
    Sym_tensor_trans uu_init(map, otriad, ff) ;  
    uu_init.set_etat_zero() ; 
    
    // uu_init = 0.5 * ( sigmat.hh() 
    //    - 0.33333333333333333 * sigmat.hh().trace(sigmat.tgam())
    //        * sigmat.tgam().con() ) ;
    //  uu_init.inc_dzpuis(2) ;  
    
    // tr K = K
    Scalar tmp(map) ; 
    tmp.set_etat_zero() ; 
    
    sigmat.initial_data_cts(uu_init, tmp, tmp, pdt, 1.e-10) ;
        
    sigmat.khi() ;  // forces updates 
    sigmat.mu() ;   //
    sigmat.trh() ;  //   
        
    cout << "Initial data : " << sigmat << endl ;  
    cout << "ADM mass : " << sigmat.adm_mass() << endl ; 
    
    cout << "Test upon khi : difference between khi and khi_init : " << endl ; 
    Scalar diff_khi = sigmat.khi() - khi_init ;
    maxabs(diff_khi, "diff_khi") ; 
    
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
        
    // Check of khi and mu

    Sym_tensor_tt htest(map, otriad, ff) ;
    htest.set_khi_mu(khi_init, mu_init) ; 

    Sym_tensor_tt htest2 = sigmat.hh().transverse(ff).tt_part() ; 
    htest2.dec_dzpuis(2) ; 
    
    maxabs(htest - htest2, "difference htest - htest2") ; 
    
    maxabs(htest.khi() - khi_init, "htest.khi() - khi_init") ; 
    maxabs(htest2.khi() - khi_init, "htest2.khi() - khi_init") ; 
    
    maxabs(htest.mu() - mu_init, "htest.mu() - mu_init") ; 
    maxabs(htest2.mu() - mu_init, "htest2.mu() - mu_init") ; 
    
    // sigmat.set_hh(htest) ; 
    
    maxabs(sigmat.psi() - 1., "Psi - 1") ;   
    
    // Scalar psitest(map) ; 
    // psitest = 1. ; 
    // sigmat.set_psi_del_q(psitest) ; 

    arrete() ; 

    //======================================================================
    //          Time evolution 
    //======================================================================    
    
    sigmat.evolve(pdt, nb_time_steps, niter_elliptic, relax_elliptic) ; 
    
    // Freeing dynamically allocated memory
    // ------------------------------------
    delete [] r_limits ; 
    
    return EXIT_SUCCESS ; 
}

