/*
 *  Code for testing the divergence-free vector Poisson equation.
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

char test_vdf_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/10/20 14:46:22  e_gourgoulhon
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
	int nzm1 = nz - 1 ;  
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

	const Coord& x = map.x ; 
	const Coord& y = map.y ; 
	const Coord& z = map.z ; 
	const Coord& r = map.r ; 
	const Coord& cost = map.cost ; 
	const Coord& sint = map.sint ; 
	const Coord& cosp = map.cosp ; 
	const Coord& sinp = map.sinp ; 
	

	cout << "========================================================" << endl ;
	cout << "                Test with a constant vector" << endl ;
	cout << "                Cartesian comp.: V^i = (1, 1, 0)" << endl ; 
	cout << "========================================================" << endl ;

	Vector_divfree vvc(map, map.get_bvect_cart(), metc ) ; 
	
	vvc.set(1) = 1 ; 
	vvc.set(1).set_dzpuis(4) ; 
	vvc.set(2) = 1 ; 
	vvc.set(2).set_dzpuis(4) ; 
	vvc.set(3) = 0 ; 
	vvc.std_spectral_base() ; 

	Vector_divfree vvs = vvc ; 
	vvs.change_triad( map.get_bvect_spher() ) ; 
	
	cout << "Spherical components : vvs : " << endl ;
	vvs.spectral_display() ; 
	arrete() ; 

	cout << "mu : " << endl ; 
	cout << "----" << endl ; 
	vvs.mu().spectral_display() ; 
	
	cout << "Norme divergence vvc : " << norme( vvc.divergence(metc) ) << endl ; 
	cout << "Norme divergence vvs : " << norme( vvs.divergence(mets) ) << endl ; 
	cout << "Max divergence vvc : " << max( abs(vvc.divergence(metc)) ) << endl ; 
	cout << "Max divergence vvs : " << max( abs(vvs.divergence(mets)) ) << endl ; 
	arrete() ; 
		
	Vector_divfree wws = vvs.poisson() ; 
	
	cout << "Solution wws : " << endl ;
	wws.spectral_display() ; 
	arrete() ; 

	cout << "========================================================" << endl ;
	cout << "                Test with a rotation vector" << endl ;
	cout << "                Cartesian comp.: V^i = (-y, x, 0)" << endl ; 
	cout << "========================================================" << endl ;

	vvc.set(1) = -y ; 
	vvc.set(2) = x ; 
	vvc.set(3) = 0 ; 
	vvc.annule_domain(nzm1) ; 
	vvc.std_spectral_base() ; 

	vvs = vvc ; 
	vvs.change_triad( map.get_bvect_spher() ) ; 
	
	cout << "Cartesian components : vvc : " << endl ;
	vvc.spectral_display() ; 
	arrete() ; 

	cout << "Spherical components : vvs : " << endl ;
	vvs.spectral_display() ; 
	arrete() ; 

	cout << "eta : " << endl ; 
	cout << "----" << endl ; 
	vvs.eta().spectral_display() ; 
	arrete() ; 

	cout << "mu : " << endl ; 
	cout << "----" << endl ; 
	vvs.mu().spectral_display() ; 
	
	wws = vvs.poisson() ; 
	
	cout << "Solution wws : " << endl ;
	wws.spectral_display() ; 
	arrete() ; 
	
	return EXIT_SUCCESS ; 
}
