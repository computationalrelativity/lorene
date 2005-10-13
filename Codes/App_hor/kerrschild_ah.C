/*
 *  Test of the AH finder on Kerr-Schild metric
 *
 *    (see file ah_finder.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005  Lap-Ming Lin & Jerome Novak
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

char kerrschild_ah_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2005/10/13 08:53:22  j_novak
 * Main program to test the apparent horizon finder.
 *
 *
 * $Header$
 *
 */

// C headers
#include <math.h>

// Lorene headers
#include "app_hor.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"


int main() {

  const int nz = 3 ; // Number of domains
  int nr = 33 ; // Number of collocation points in r in each domain
  int nt = 33 ; // Number of collocation points in theta in each domain
  int np = 4 ; // 4 points in phi for rotation to/from cartesian components
      
  int symmetry_theta = SYM ; // symmetry w.r.t. the equatoial plane
                             // ##Bug report: symmetry_theta=NONSYM does not work correctly!
                             //     (related to change_triad....)
  int symmetry_phi = SYM ; // symmetry in phi

  int nbr[] = {nr, nr, nr} ;
  int nbt[] = {nt, nt, nt} ;
  int nbp[] = {np, np, np} ;

  int type_r_std[] = {FIN, FIN, UNSURR} ; // sampling for the standard mapping
  int type_r_2[] = {FIN, FIN, FIN} ; // sampling for the 2nd mapping

  // Construct 3D grid for the standard mapping (with nucleus)
  Mg3d mgrid(nz, nbr, type_r_std, nbt, symmetry_theta, nbp, symmetry_phi) ;
  

  // Construct 3D grid for the 2nd mapping (without nucleus)
  // The apparent horizon is defined on the 2nd mapping
  Mg3d mgrid_2(nz, nbr, type_r_2, nbt, symmetry_theta, nbp, symmetry_phi) ;


  // Construct angular grid for h(theta,phi) 
  const Mg3d& g_angu = *mgrid_2.get_angu() ;

  // orig_lim_1(2) is used to cut out the signular region....
  // note: For the purpose of interpolating the metric defined on the 1st mapping to 
  //       that on the 2nd one, you must choose orig_lim_2 >= orig_lim_1
  double orig_lim_1 = 0.55 ;
  double orig_lim_2 = 0.55 ; 

  double r_limits_std[] = {orig_lim_1, 1.5, 3.5, __infinity} ; // for the 1st mapping
  double r_limits_2[] = {orig_lim_2, 1., 2., 3.} ; // for the 2nd mapping
  //note: the outer boundary of the 2nd mapping should not be in the compactified domain 
  //      of the 1st mapping (for the purpose of interpolation)
  

  // Standard mapping
  Map_af map(mgrid, r_limits_std) ;

  // 2nd mapping 
  Map_af map_2(mgrid_2, r_limits_2) ; 

  const Base_vect_spher& bspher = map.get_bvect_spher() ;
  const Base_vect_spher& bspher2 = map_2.get_bvect_spher() ;

  // Some helpful stuff...
//  const Coord& rr = map_2.r ;
  const Coord& theta = map_2.tet ;
  const Coord& phi = map_2.phi ;
  const Coord& costh = map_2.cost ;
  const Coord& cosph = map_2.cosp ;
  const Coord& sinth = map_2.sint ;
  const Coord& sinph = map_2.sinp ;


  // User's input of \gamma_{ij} & K_{ij}
  Sym_tensor gam_dd_in(map, COV, bspher) ;
  Sym_tensor k_dd_in(map, COV, bspher) ;

  // The \gamma_{ij} & K_{ij} defined on the 2nd mapping 
  Sym_tensor gam_dd(map_2, COV, bspher2) ;
  Sym_tensor k_dd(map_2, COV, bspher2) ;


  
  //Parameters to control the iteration

  double tol = 1.e-12 ;
  int it_max  = 100 ;  
  int it_relax = 100 ;
  double relax_fac = 1. ;

  // Parameters for the initial guess to the 2-surface h
  // "aa=bb=cc=R" corresponds to a sphere of radius R (usually good enough)
  // If you want aa != bb, remember to set "symmetry_phi = NONSYM"
  double aa = 1.5 ; 
  double bb = 1.5 ;   
  double cc = 1.5 ;


  // Parameters for the analytic model input:

  double mass = 0.5 ;
  double kerr_a = 0.3 ; 

  
  Valeur h_correct(g_angu) ;  // Analytic result h(\theta,phi) for comparison 
  h_correct.annule_hard() ;



      //***** Kerr-Schild data (original Cartesian form) ********//


    const Coord& xx = map.x ;
    const Coord& yy = map.y ;
    const Coord& zz = map.z ;
    const Coord& rr1 = map.r ;

    Scalar x(map) ;
    x = xx ;

    Scalar y(map) ;
    y = yy ;

    Scalar z(map) ;
    z = zz ;

    Scalar r_rho(map) ;
    r_rho = rr1 ;


    // Define Boyer-Lindquist radial coordinate
    Scalar rbl(map) ;
    rbl = 0.5*(r_rho*r_rho - kerr_a*kerr_a) + 
      sqrt( 0.25*pow(r_rho*r_rho - kerr_a*kerr_a , 2) + kerr_a*kerr_a*z*z ) ;

    rbl = sqrt(rbl) ;
    rbl.std_spectral_base() ;


    Scalar hhh(map) ;
    hhh = mass*pow(rbl,3) / ( pow(rbl,4) + kerr_a*kerr_a*z*z ) ;    
    hhh.std_spectral_base() ;


    // Define Ingoing null vector 
    const Base_vect_cart& bcart = map.get_bvect_cart() ; //Cartesian base vector

    Vector lmu(map, COV, bcart) ; 

    lmu.set(1) = (rbl*x + kerr_a*y) / (rbl*rbl + kerr_a*kerr_a) ;
    lmu.set(2) = (rbl*y - kerr_a*x) / (rbl*rbl + kerr_a*kerr_a) ;
    lmu.set(3) = z / rbl ;
    lmu.std_spectral_base() ;


    // Lapse 
    Scalar lapse(map) ;

    lapse = pow( 1. + 2.*hhh , -0.5 ) ;
    lapse.std_spectral_base() ;


    //Shift
    Vector beta(map, COV, bcart) ;
    
    beta.set(1) = 2.*hhh*lmu.set(1) ;
    beta.set(2) = 2.*hhh*lmu.set(2) ;
    beta.set(3) = 2.*hhh*lmu.set(3) ;
    beta.std_spectral_base() ;
    

    // Define 3-metric (in cartesian basis)

    Sym_tensor g_tmp(map, COV, bcart) ; 

    g_tmp.set(1,1) = 1. + 2.*hhh*lmu.set(1)*lmu.set(1) ; 
    g_tmp.set(1,2) = 2.*hhh*lmu.set(1)*lmu.set(2) ; 
    g_tmp.set(1,3) = 2.*hhh*lmu.set(1)*lmu.set(3) ;
    g_tmp.set(2,2) = 1. + 2.*hhh*lmu.set(2)*lmu.set(2) ;
    g_tmp.set(2,3) = 2.*hhh*lmu.set(2)*lmu.set(3) ;
    g_tmp.set(3,3) = 1. + 2.*hhh*lmu.set(3)*lmu.set(3) ;

    g_tmp.std_spectral_base() ;

    
    //Change to spherical orthonormal basis

    g_tmp.change_triad(bspher) ; 
    beta.change_triad(bspher) ;

    gam_dd_in = g_tmp ;      // gam_2 is used in the iteration loop to define the 3-metric


    const Metric gam_met(g_tmp) ;

    k_dd_in = 0.5*beta.ope_killing(gam_met) / lapse ;    

    gam_dd_in.std_spectral_base() ;
    k_dd_in.std_spectral_base() ;

    // Done with input data
    

    //-------------------------------------------------------------
    // Calculate the analytic AH surface for Kerr-Schild data 

      Scalar ct(map_2) ;  
      ct = costh ;
      ct.std_spectral_base() ;

      Scalar st(map_2) ;
      st = sinth ;
      st.std_spectral_base() ;

      Scalar cp(map_2) ;
      cp = cosph ;
      cp.std_spectral_base() ;

      Scalar sp(map_2) ;
      sp = sinph ;
      sp.std_spectral_base() ;


      Scalar h_tmp(map_2) ;

      double r_analy = mass + sqrt( mass*mass - kerr_a*kerr_a ) ;


      h_tmp = st*st / ( r_analy*r_analy + kerr_a*kerr_a ) + ct*ct/(r_analy*r_analy) ;
      h_tmp = 1./ sqrt( h_tmp ) ;

      h_tmp.std_spectral_base() ;


      for (int l=0; l<nz; l++) {
	  for (int k=0; k<np; k++) {
	      for (int j=0; j<nt; j++) {
		  for (int i=0; i<nr; i++) {
		      
		      h_correct.set(l,k,j,0)  = h_tmp.val_grid_point(l,k,j,i) ; 
		      
		  }
	      }
	  }
      }

      h_correct.std_base_scal() ;


      // End calculate the correct AH surface
      //-------------------------------------------------------------

    //******** End Kerr-Schild data (original Cartesian form) *********//


      // Interpolate the input data from the mapping "map" to "map_2"
      // on which the apparent horizon is to be found. 
      // Interpolate gam_dd_in (k_dd_in) to gam_dd (k_dd)
      // (This part is to be removed later.....)

      k_dd_in.dec_dzpuis(2) ; //The dzpuis of the Scalar to be imported 
                          //must be zero!

      gam_dd.set(1,1).import( gam_dd_in.set(1,1) ) ;
      gam_dd.set(1,2).import( gam_dd_in.set(1,2) ) ;
      gam_dd.set(1,3).import( gam_dd_in.set(1,3) ) ;
      gam_dd.set(2,2).import( gam_dd_in.set(2,2) ) ;
      gam_dd.set(2,3).import( gam_dd_in.set(2,3) ) ;
      gam_dd.set(3,3).import( gam_dd_in.set(3,3) ) ;
      gam_dd.std_spectral_base() ;

      Metric gamma(gam_dd) ;   // Set up the 3-metric 

      k_dd.set(1,1).import( k_dd_in.set(1,1) ) ;
      k_dd.set(1,2).import( k_dd_in.set(1,2) ) ;
      k_dd.set(1,3).import( k_dd_in.set(1,3) ) ;
      k_dd.set(2,2).import( k_dd_in.set(2,2) ) ;
      k_dd.set(2,3).import( k_dd_in.set(2,3) ) ;
      k_dd.set(3,3).import( k_dd_in.set(3,3) ) ;
      k_dd.std_spectral_base() ;

      k_dd_in.inc_dzpuis(2) ; //Increase the dzpuis of k_dd_in back to two

      // End importing data


      // Call the AH finder (return the 2-surface h & expansion function)

      Valeur h(g_angu) ;  // 2-surface of the apparent horizon
      h.annule_hard() ;
      Scalar exp_fcn(map_2) ; // expansion function
      
//      ah_finder(gamma, k_dd, h, exp_fcn, flag, aa, bb, cc) ;

      cout << " AH flag = " 
	   << ah_finder(gamma, k_dd, h, exp_fcn, aa, bb, cc, tol, it_max, 
			it_relax, relax_fac) << endl ;


      // Test and analysis.....

      cout << "########################" << endl ;
      cout << "Done with the AH finder" << endl ;
      cout << "########################" << endl ;
      cout << " " << endl ;
      cout << "Print the 2-surface h(theta,phi) of the apparent horizon" << endl ;
      cout << h(0) << endl ;

      cout << "########################" << endl ;
      cout << " Print (h_analytic - h) / h_analytic " << endl ;
      cout <<  ( h_correct(0) - h(0) ) / h_correct(0) << endl ;


      exp_fcn.std_spectral_base() ;

      Valeur exp_AH(g_angu) ;   // expansion function evaulated on the apparent 
                                 // horizon (which should be zero)
      exp_AH.annule_hard() ;

    for (int l=0; l<nz; l++) {
      for (int k=0; k<np; k++) {
	for (int j=0; j<nt; j++) {

	  exp_AH.set(l,k,j,0) = exp_fcn.val_point(h(l,k,j,0),(+theta)(l,k,j,0)
						     ,(+phi)(l,k,j,0))  ;  
	}
      }
    }

    cout << "##############################################" << endl ;
    cout << "Test: evaluate the expansion function on the apparent horizon " << endl ;
    cout << "   " << endl ; 

    cout << exp_AH(0) << endl ;

    des_profile(exp_fcn, orig_lim_2, 2., 0., 0.) ;


} // End main


