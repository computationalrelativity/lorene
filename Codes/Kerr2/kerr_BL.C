/*
 *  Kerr metric in Boyer-Lindquist coordinates
 */

/*
 *   Copyright (c) 2011 Eric Gourgoulhon
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

char kerr_BL_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2011/11/25 16:44:42  e_gourgoulhon
 * First version; not fully checked yet
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cmath>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "cmp.h"
#include "proto.h"
#include "graphique.h"

int main() {

	// Parameters of the computation
	// -----------------------------
	
    ifstream fpar("par_kerr_BL.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_kerr_BL.d ! " << endl ;
        abort() ;
    }
    
	double aa ; // Kerr parameter a/M
    fpar >> aa ; fpar.ignore(1000,'\n') ;
	double aa2 = aa*aa ; 

	int nr ; // Number of collocation points in r in each domain
    fpar >> nr; fpar.ignore(1000,'\n') ;

	int nt ; // Number of collocation points in theta in each domain
    fpar >> nt; fpar.ignore(1000,'\n') ;

    int np ; // Number of collocation points in phi in each domain
	fpar >> np; fpar.ignore(1000,'\n') ;
	
	int nz ; // Number of domains
    fpar >> nz ; fpar.ignore(1000,'\n') ;
    int nzm1 = nz - 1 ; // Index of outermost domain

    fpar.ignore(1000,'\n') ; // skip title
    double* r_limits = new double[nz+1];  // radial boundaries of each domain in units of M      
    for (int l=0; l<nz; l++) {
		fpar >> r_limits[l]; 
    }
    r_limits[nz] = __infinity ;

    fpar.close();

    //## check
    cout << "r_limits : " << endl ; 
    for (int l=0; l<nz+1; l++) {
      cout << r_limits[l] << "  " ; 
    }	
    cout << endl ; 
    
    // value of coordinate r at the event horizon: 
    double r_hor = 1 + sqrt(1 - aa2) ; 
    cout << "Value of coordinate r at the event horizon : " << r_hor << " M" << endl ; 
    if (r_limits[1] <= r_hor) {
      cerr << "Inner boundary of domain no. 1 below the horizon : " << endl ; 
      cerr << "  r_limits[1] : " << r_limits[1] << endl ; 
      cerr << "  r_hor :       " << r_hor << endl ; 
      abort() ; 
    }

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------
  
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = NONSYM ; // no symmetry in phi
    bool compact = true ; // external domain is compactified

    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
    cout << mgrid << endl ; 

  
    // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
    // --------------------------------------------------------------------------
  
    Map_af map(mgrid, r_limits) ;
  	
    // cout << map << endl ;  
    
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    const Coord& cost = map.cost ;  // cos(theta) field
    const Coord& sint = map.sint ;  // sin(theta) field
    Mtbl r2 = r*r ; 
    Mtbl cost2 = cost*cost ; 
    Mtbl sint2 = sint*sint ; 

    // rho^2
    Scalar rho2(map) ; 
    rho2 = r2 + aa2 * cost2 ; 
    rho2.std_spectral_base() ; 
    
    // rho^2 / r^2
    Scalar rho2_ovr2(map) ; 
    rho2_ovr2 = 1 + aa2 * cost2 / r2 ; 
    rho2_ovr2.std_spectral_base() ; 
    
    // Delta / r^2
    Scalar delta_ovr2(map) ; 
    delta_ovr2 = 1 - 2/r + aa2/r2 ; 
    delta_ovr2.std_spectral_base() ; 
    
    // B^2 = Sigma^2 / (r*rho)^2 
    Scalar bb2(map) ; 
    bb2 = 1 +  aa2/r2 + 2*aa2*sint2 / (r*rho2) ; 
    bb2.std_spectral_base() ;     
    
    // Lapse
    Scalar nn = sqrt( delta_ovr2 / bb2 ) ; 
    nn.set_domain(0) = 1 ; 
    nn.std_spectral_base() ; 
     
    // Shift vector
    Vector beta(map, CON, map.get_bvect_spher()) ;
    beta.set(1) = 0 ;
    beta.set(2) = 0 ;
    beta.set(3) = - 2 *aa * sint / (rho2*(1 +  aa2/r2) + 2*aa2*sint2/r) ; 
    beta.set(3).annule_domain(0) ; 
    beta.set(3).std_spectral_base() ; 
    
    // 3-metric
    Sym_tensor gamma(map, COV, map.get_bvect_spher()) ;
    gamma.set(1,1) = rho2_ovr2 / delta_ovr2 ; 
    gamma.set(1,1).set_domain(0) = 1 ; 
    gamma.set(1,2) = 0 ; 
    gamma.set(1,3) = 0 ; 
    gamma.set(2,2) = rho2_ovr2 ; 
    gamma.set(2,2).set_domain(0) = 1 ; 
    gamma.set(2,3) = 0 ; 
    gamma.set(3,3) = bb2 ; 
    gamma.set(3,3).set_domain(0) = 1 ; 
    
    Sym_tensor inv_gamma(map, CON, map.get_bvect_spher()) ;
    inv_gamma.set(1,1) = 1 / gamma(1,1) ; 
    inv_gamma.set(1,2) = 0 ; 
    inv_gamma.set(1,3) = 0 ; 
    inv_gamma.set(2,2) = 1 / gamma(2,2) ; 
    inv_gamma.set(2,3) = 0 ; 
    inv_gamma.set(3,3) = 1 / gamma(3,3) ; 
    
    // Extrinsic curvature
    Scalar beta_phi(map) ; 
    beta_phi = - 2*aa / (r*rho2*(1 +  aa2/r2) + 2*aa2*sint2) ; 
    beta_phi.annule_domain(0) ;
    beta_phi.std_spectral_base() ; 
    Sym_tensor kk(map, COV, map.get_bvect_spher()) ;
    kk.set(1,1) = 0 ; 
    kk.set(1,2) = 0 ; 
    Scalar tmp = 0.5 * bb2 * beta_phi.dsdr() / nn ;
    tmp.mult_rsint() ;
    kk.set(1,3) = tmp ;
    kk.set(2,2) = 0 ; 
    tmp = 0.5 * bb2 * beta_phi.dsdt() / nn ;
    tmp.mult_sint() ;
    tmp.inc_dzpuis(2) ;
    kk.set(2,3) = tmp ;
    kk.set(3,3) = 0 ;

    
    // Drawings
    
    des_meridian(nn, 0, 1.5*r_limits[nzm1], "lapse", 1) ; 

    des_meridian(beta(3), 0, 1.5*r_limits[nzm1], "shift (phi)", 2) ; 

    des_meridian(gamma(1,1), 0, 1.5*r_limits[nzm1], "gamma_11", 3) ; 
    des_meridian(gamma(2,2), 0, 1.5*r_limits[nzm1], "gamma_22", 4) ; 
    des_meridian(gamma(3,3), 0, 1.5*r_limits[nzm1], "gamma_33", 5) ; 
    
    des_meridian(kk(1,3), 0, r_limits[nzm1], "K_13", 6) ; 
    des_meridian(kk(2,3), 0, r_limits[nzm1], "K_23", 7) ; 
    
    arrete() ; 

    // Output file
    //------------
    FILE* file_out = fopen("gyoto_kerr_BL.d", "w") ;
    double total_time = 0. ; // for compatibility

    fwrite_be(&total_time, sizeof(double), 1, file_out) ;
    mgrid.sauve(file_out) ;
    map.sauve(file_out) ;
    nn.sauve(file_out) ;
    beta.sauve(file_out) ;
    gamma.sauve(file_out) ;
    inv_gamma.sauve(file_out) ;
    kk.sauve(file_out) ;

    fclose(file_out) ;    
    
    return EXIT_SUCCESS ; 

}
