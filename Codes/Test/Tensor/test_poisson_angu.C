/*
 *  Code for testing the angular Poisson equation.
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

char test_poisson_angu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/10/15 21:15:25  e_gourgoulhon
 * Test for Scalar::poisson_angu().
 *
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
	int nr = 5 ; 	// Number of collocation points in r in each domain
	int nt = 5 ; 	// Number of collocation points in theta in each domain
	int np = 12 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi
	bool compact = true ; // external domain is compactified
  
	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double r_limits[] = {0., 1., 2., __infinity} ; 
  	assert( nz == 3 ) ;  // since the above array described only 3 domains
  
	Map_af map(mgrid, r_limits) ; 
  
	

	cout << endl << 
	"================ TEST FOR A SCALAR FIELD =================\n" ; 
	

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

	uu = x + y + x*y + z*z*x - pow(z,4)*y*x ; 	
	tmp = sint * sinp / r ;
	uu.set_domain(nz-1) = tmp.domain(nz-1) ; // y/r^2 in the external domain
	
	uu.std_spectral_base() ;   // sets the standard spectral basis for 
									// expansion of a scalar field
									
	cout << "uu : " << uu << endl ; 
	uu.spectral_display(cout) ; 
	arrete() ; 
	
	Scalar lap = uu.lapang() ; 

	cout << "lap : " << endl ; 
	lap.spectral_display(cout) ; 
	arrete() ; 

	Scalar uu1 = lap.poisson_angu() ; 
	
	cout << "uu1 : " << endl ; 
	uu1.spectral_display(cout) ; 
	arrete() ; 

	Scalar diff = uu - uu1 ; 
	
	cout << "Norm of diff: " << norme(diff) << endl ; 	
	
	return EXIT_SUCCESS ; 
}
