/*
 *  Code for testing the covariant derivatives through the Connection class.
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

char test_connect_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2003/10/06 20:53:16  e_gourgoulhon
 * New version: constructs flat_metric and calls Tensor::derive_cov.
 *
 * Revision 1.3  2003/10/05 21:18:08  e_gourgoulhon
 * Added test onto a vector field.
 *
 * Revision 1.2  2003/10/03 14:27:51  e_gourgoulhon
 * First non trivial test (successfull !).
 *
 * Revision 1.1  2003/10/02 21:33:02  e_gourgoulhon
 * Test code for Connection.
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

	int arret = 1 ; 

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


	// Construction of a flat connection
	// ---------------------------------
	
	// Representation on a spherical orthonormal basis
	Connection_fspher ders(map, map.get_bvect_spher()) ; 
	
	// Representation on a Cartesian orthonormal basis
	 Connection_fcart derc(map, map.get_bvect_cart()) ; 


	cout << endl << 
	"================ TEST FOR A SCALAR FIELD =================\n" ; 
	

	// Construction of a scalar field (Scalar)
	// ---------------------------------------

	const Coord& x = map.x ; 
	const Coord& z = map.z ; 
	const Coord& r = map.r ; 
	const Coord& sint = map.sint ; 
	const Coord& sinp = map.sinp ; 

	Scalar uu(map) ; 
	Scalar tmp(map) ; 

	uu = x ; 	
	tmp = sint * sinp / r ;
	uu.set_domain(nz-1) = tmp.domain(nz-1) ; // y/r^2 in the external domain
	
	uu.std_spectral_base() ;   // sets the standard spectral basis for 
									// expansion of a scalar field
									
	cout << "uu : " << uu << endl ; 
	uu.spectral_display(cout) ; 
	arrete(arret) ; 
	
	// Gradient of the scalar field
	// ----------------------------
	
	Vector duuc = uu.derive_cov(metc)  ;  

	cout << "duuc : " << duuc << endl ; 

	arrete(arret) ; 

	Vector duus = uu.derive_cov(mets) ; 
	
	cout << "duus : " << duus << endl ; 

	arrete(arret) ; 
	
	// Test
	// ----
	
	Vector duus_c = duus ; 
	duus_c.change_triad( map.get_bvect_cart() ) ; 
	
	Vector diffc = duus_c - duuc ; 
	
	cout << "diffc : " << diffc << endl ; 
	cout << "Norm of diffc: " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << norme(diffc(i)) << endl ; 
	}
	arrete(arret) ; 
	
	
	Vector duuc_s = duuc ; 
	duuc_s.change_triad( map.get_bvect_spher() ) ; 
	
	Vector diffs = duuc_s - duus ; 

	cout << "diffs : " << diffs << endl ; 
	cout << "Norm of diffs: " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << norme(diffs(i)) << endl ; 
	}
	arrete(arret) ; 
	
	cout << endl << 
	"================ TEST FOR A VECTOR FIELD =================\n" ; 
	
	Vector vvc(map, COV, map.get_bvect_cart()) ; 
	
	vvc.set(1) = uu ; 
	vvc.set(2) = uu*uu ; 
	vvc.set(3) = z + x * z ; 
	vvc.set(3).set_domain(nz-1) = duuc(3).domain(nz-1) ; 
	
	vvc.std_spectral_base() ; 
	Tensor dvvc = vvc.derive_cov(metc) ;
	
	// cout << "dvvc : " << dvvc << endl ; 
	
	Vector vvs = vvc ; 
	vvs.change_triad( map.get_bvect_spher() ) ; 
	
	Tensor dvvs = ders.derive_cov(vvs) ; 
	
	Tensor dvvs_c = dvvs ; 
	dvvs_c.change_triad( map.get_bvect_cart() ) ; 
	
	Tensor diffvvc = dvvs_c - dvvc ; 
	
	cout << "dvvc(1,1) : "  << endl ;
	dvvc(1,1).spectral_display(cout) ;  
	arrete() ; 
	
	cout << "dvvs_c(1,1) : " << endl ; 
	dvvs_c(1,1).spectral_display(cout) ;  
	arrete() ; 
	
	cout << "Norm of diffvvc: " << endl ; 
	for (int i=1; i<=3; i++) {
		for (int j=1; j<=3; j++) {
			cout << i << " " << j << " : " << norme(diffvvc(i,j)) << endl ; 
		}
	}
	arrete(arret) ; 



	Tensor dvvc_s = dvvc ; 
	dvvc_s.change_triad( map.get_bvect_spher() ) ; 

	Tensor diffvvs = dvvc_s - dvvs ; 
	
	cout << "Norm of diffvvs: " << endl ; 
	for (int i=1; i<=3; i++) {
		for (int j=1; j<=3; j++) {
			cout << i << " " << j << " : " << norme(diffvvs(i,j)) << endl ; 
		}
	}
	arrete(arret) ; 
	
	return EXIT_SUCCESS ; 
}
