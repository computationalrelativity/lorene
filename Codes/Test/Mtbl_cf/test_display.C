/*
 *  Code for testing the new Mtbl_cf display
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon 
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

char test_display_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/10/19 20:03:49  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"


int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	int nz = 3 ; 	// Number of domains
	int nr = 9 ; 	// Number of collocation points in r in each domain
	int nt = 5 ; 	// Number of collocation points in theta in each domain
	int np = 8 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi
	bool compact = true ; // external domain is compactified
  
	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double r_limits[] = {0., 1., 2., __infinity} ; 
  	assert( nz == 3 ) ;  // since the above array describes only 3 domains
  
	Map_af map(mgrid, r_limits) ; 
  

	// Construction of a scalar field (Scalar)
	// ---------------------------------------

	const Coord& x = map.x ; 
	const Coord& y = map.y ; 
	const Coord& z = map.z ; 
	const Coord& r = map.r ; 
	const Coord& sint = map.sint ; 
	const Coord& sinp = map.sinp ; 

	Scalar uu(map) ; 
	Scalar tmp(map) ; 

	uu = 1 + x + 2*y - z*z + x*y + x*x*y - z*z*y  ; 	
	tmp = sint * sinp / r ;
	uu.set_domain(nz-1) = tmp.domain(nz-1) ; // y/r^2 in the external domain
	uu.set_domain(1) = 0 ; 
	
	uu.std_spectral_base() ;   // sets the standard spectral basis for 
									// expansion of a scalar field
									
	Valeur va = uu.get_spectral_va() ;
	va.coef() ; 
	va.ylm() ; 
	va.c_cf->affiche_seuil(cout) ; 
	va.c_cf->display() ; 
	
	arrete() ; 
	cout << "uu : " << endl ; 
	uu.spectral_display() ; 
	arrete() ; 
	
	// Construction of a flat metric
	// -----------------------------

	Metric_flat mets(map, map.get_bvect_spher()) ; // spherical representation

	// Gradient of the scalar field
	// ----------------------------
	
	Vector vv = uu.derive_cov(mets)  ;  
	vv.set(1) = 0 ; 

	cout << "vv : " << endl ; 
	vv.spectral_display() ; 	
	
	
	return EXIT_SUCCESS ; 
}
