/*
 *  Code for testing the class Vector_divfree.
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
 * Revision 1.2  2003/10/17 16:33:36  e_gourgoulhon
 * Added more tests.
 *
 * Revision 1.1  2003/10/16 21:39:43  e_gourgoulhon
 * Test code for class Vector_divfree.
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
  	

	// Construction of a flat metric
	// -----------------------------

	Metric_flat mets(map, map.get_bvect_spher()) ; // spherical representation
	Metric_flat metc(map, map.get_bvect_cart()) ;  // Cartesian representation


	// Construction of a divergence free vector field 
	// ----------------------------------------------

	/* const Coord& x = map.x ; 
	const Coord& y = map.y ; 
	const Coord& z = map.z ; 
	const Coord& r = map.r ; 
	const Coord& sint = map.sint ; 
	const Coord& sinp = map.sinp ; 
	*/
	
	Vector_divfree vvc(map, map.get_bvect_cart(), metc ) ; 
	
	vvc.set(1) = 1 ; 
	vvc.set(2) = 1 ; 
	vvc.set(3) = 0 ; 
	vvc.std_spectral_base() ; 

	Vector_divfree vvs = vvc ; 
	vvs.change_triad( map.get_bvect_spher() ) ; 
	
	cout << "vvs : " << endl ;
	for (int i=1; i<=3; i++) {
		cout << "Component " << i << " : " << endl ; 
		vvs(i).spectral_display(cout) ; 
		arrete() ; 
	}
	

	cout << "eta : " << endl ; 
	vvs.eta().spectral_display(cout) ; 
	arrete() ; 

	cout << "mu : " << endl ; 
	vvs.mu().spectral_display(cout) ; 
	
	arrete() ; 
	cout << "divergence vvc : " << vvc.divergence(metc) << endl ;
	cout << "  norme: " << norme( vvc.divergence(metc) ) << endl ; 
	arrete() ; 
	cout << "divergence vvs : " << vvs.divergence(mets) << endl ;
	cout << "  norme: " << norme( vvs.divergence(mets) ) << endl ; 
	
	Vector_divfree vvs2(map, map.get_bvect_spher(), mets ) ;
	vvs2.set_vr_eta_mu(vvs(1), vvs.eta(), vvs.mu()) ; 
	Vector diff = vvs - vvs2 ; 
	cout << "diff : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << "Component " << i << " : " << endl ; 
		diff(i).spectral_display(cout) ; 
		cout << "  norme: " << norme(diff(i)) << endl ; 
		arrete() ; 
	}
	

	
	return EXIT_SUCCESS ; 
}
