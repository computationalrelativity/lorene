/*
 *  Code for testing the classes Sym_tensor_trans and Sym_tensor_tt.
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

char test_sym_tensor_tt_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2003/10/30 17:29:07  e_gourgoulhon
 * new version
 *
 * Revision 1.4  2003/10/29 13:15:26  e_gourgoulhon
 * Change of method name: Scalar::laplacien --> Scalar::laplacian.
 *
 * Revision 1.3  2003/10/28 21:37:22  e_gourgoulhon
 * New version.
 *
 * Revision 1.2  2003/10/28 12:36:52  e_gourgoulhon
 * improved version
 *
 * Revision 1.1  2003/10/27 10:56:57  e_gourgoulhon
 * Test code for class Sym_tensor_*.
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
	int nr = 9 ; 	// Number of collocation points in r in each domain
	int nt = 9 ; 	// Number of collocation points in theta in each domain
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


	// Construction of a divergence free tensor field 
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
	cout << "                Test with a constant tensor" << endl ;
	cout << "   Cart. comp. h^{ij} = (0, 1, 2), (2,-1,3), (2,3,1) " << endl ; 
	cout << "========================================================" << endl ;

	Sym_tensor_trans hhc(map, map.get_bvect_cart(), metc ) ; 
	
	hhc.set(1,1) = 1 ; 
	hhc.set(1,2) = 2 ; 
	hhc.set(1,3) = 0 ; 
	hhc.set(2,2) = -4 ; 
	hhc.set(2,3) = 0 ; 
	hhc.set(3,3) = 3 ; 
	hhc.std_spectral_base() ; 

	cout << "Cartesian components : hhc : " << hhc << endl ;
	arrete() ; 

	Tensor tmp = hhc ; 
	tmp.change_triad( map.get_bvect_spher() ) ; 
	Sym_tensor_trans hhs(map, map.get_bvect_spher(), mets ) ; 
	hhs = tmp ; 
	
	cout << "Spherical components : hhs : " << endl ;
	hhs.spectral_display() ; 
	arrete() ; 
	
	hhc.divergence(metc) ; 
	cout << "Norme divergence hhc : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << norme( hhc.divergence(metc)(i) ) << endl ; 
	}
	cout << "Norme divergence hhs : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << norme( hhs.divergence(mets)(i) ) << endl ; 
	}

	cout << "Max divergence hhc : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << max( hhc.divergence(metc)(i) ) << endl ; 
	}
	
	cout << "Max divergence hhs : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << max( hhs.divergence(mets)(i) ) << endl ; 
	}
		
	arrete() ;
	 
	cout << "Trace of hhc : " << hhc.trace() << endl ; 	
	arrete() ; 
		
	cout << "Trace of hhs : " << hhs.trace() << endl ; 	
	

	cout << "========================================================" << endl ;
	cout << "                Test with the tensor" << endl ;
	cout << "   Cart. comp. h^{ij} = d_i d_j Phi  with Lap(Phi) = 0 " << endl ; 
	cout << "========================================================" << endl ;
	
	arrete() ; 
	
	Scalar pot(map) ; 
	pot = x*y ; 
	pot.annule_domain(nzm1) ; 
	Mtbl potced = 1 / r ; 
	pot.set_domain(nzm1) = potced(nzm1) ; 
	pot.std_spectral_base() ; 

	cout << "Potential : " << endl ; 
	pot.spectral_display() ; 
	arrete() ; 
	
	cout << "Laplacian of potential : " << endl ; 
	pot.laplacian().spectral_display() ; 
	arrete() ; 
	
	hhc.set(1,1) = pot.dsdx(0).dsdx(0) ; 
	hhc.set(1,2) = pot.dsdx(0).dsdy(0) ; 
	hhc.set(1,3) = pot.dsdx(0).dsdz(0) ; 
	hhc.set(2,2) = pot.dsdy(0).dsdy(0) ; 
	hhc.set(2,3) = pot.dsdy(0).dsdz(0) ; 
	hhc.set(3,3) = pot.dsdz(0).dsdz(0) ; 

	cout << "Cartesian components : hhc : " <<  endl ;
	hhc.spectral_display() ; 
	arrete() ; 

	tmp = hhc ; 
	tmp.change_triad( map.get_bvect_spher() ) ; 
	hhs = tmp ; 
	
	cout << "Spherical components : hhs : " << endl ;
	hhs.spectral_display() ; 
	arrete() ; 
	
	Vector divc = hhc.divergence(metc) ; 
	cout << "Divergence of hhc : " << endl ; 
	divc.spectral_display() ; 
	
	cout << "Max divergence hhc : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << max( abs(divc(i)) ) << endl ; 
	}
	
	arrete() ;
	
	Vector divs = hhs.divergence(mets) ; 
	cout << "Divergence of hhs : " << endl ; 
	divs.spectral_display() ; 
	
	cout << "Max divergence hhs : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << max( abs(divs(i)) ) << endl ; 
	}
		
	cout << "Max of trace of hhc : " << endl ; 
	cout << max( abs(hhc.trace()) ) << endl ; 	
		
	cout << "Max of trace of hhs : " << endl ; 
	cout << max( abs(hhs.trace()) ) << endl ; 	
		
		
	return EXIT_SUCCESS ; 
}
