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

 

/*
 * $Id$
 * $Log$
 * Revision 1.8  2016/12/05 16:18:29  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:54:02  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:12:55  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2003/10/29 13:16:11  e_gourgoulhon
 * Change of method name: Scalar::laplacien --> Scalar::laplacian.
 *
 * Revision 1.4  2003/10/29 11:06:11  e_gourgoulhon
 * inc2_dzpuis() replaced by inc_dzpuis(2).
 *
 * Revision 1.3  2003/10/21 13:59:36  e_gourgoulhon
 * new version
 *
 * Revision 1.2  2003/10/20 19:46:41  e_gourgoulhon
 * First successful version.
 *
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
#include <cstdlib>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"


using namespace Lorene ;

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	int nz = 3 ; 	// Number of domains
	int nzm1 = nz - 1 ;  
	int nr = 17 ; 	// Number of collocation points in r in each domain
	int nt = 17; 	// Number of collocation points in theta in each domain
	int np = 12 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi
	bool compact = true ; // external domain is compactified
  
	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double r_limits[] = {0., 0.5, 1., __infinity} ; 
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
	cout << "                Test with a pretty general vector" << endl ;
	cout << "                           V = curl(A) " << endl ; 
	cout << "========================================================" << endl ;

	Vector aa(map, CON, map.get_bvect_cart()) ; 
	aa.set(1) = z * (x + x*y - 3*z*z) ; 
	aa.set(2) = z * ( z*z*x - 2 *y + 1 ); 
	aa.set(3) =  1 + x*x - y + x + x*y + z*z  ; 
	aa.annule_domain(nzm1) ; 

	Mtbl tced = cost * ( sint*cosp/r + sint*sint*cosp*sinp - 3*cost*cost ) 
				/ (r*r*r) ; 
	tced = cost / (r*r*r) ; 
	aa.set(1).set_domain(nzm1) = tced(nzm1) ; 

	tced =  cost * ( cost*cost*sint*cosp - 2*sint*sinp / r + 1) / (r*r*r) ; 
//	aa.set(2).set_domain(nzm1) = tced(nzm1) ; 

//	tced = (1 + sint*sint*cosp*cosp / r - sint*sinp + sint*cosp*sint*sinp
//			+ cost*cost) / (r*r*r)  ; 
	tced = cost*cost*sint*sinp / (r*r*r)  ; 
	aa.set(3).set_domain(nzm1) = tced(nzm1) ; 


	aa.std_spectral_base() ; 
	aa.set(1).set_spectral_va().set_base_r(0, R_CHEBPIM_I) ; 
	aa.set(1).set_spectral_va().set_base_t(T_COSSIN_CI) ; 
	aa.set(2).set_spectral_va().set_base_r(0, R_CHEBPIM_I) ; 
	aa.set(2).set_spectral_va().set_base_t(T_COSSIN_CI) ; 
	aa.set(3).set_spectral_va().set_base_r(0, R_CHEBPIM_P) ; 
	aa.set(3).set_spectral_va().set_base_t(T_COSSIN_CP) ; 
	
	cout << "aa : " << endl ; 
	aa.spectral_display() ; 
	arrete() ; 


	// Curl of aa:
	Vector_divfree vvc(map, map.get_bvect_cart(), metc ) ; 	
	vvc.set(1) = aa(3).dsdy() - aa(2).dsdz() ; 
	vvc.set(2) = aa(1).dsdz() - aa(3).dsdx() ; 
	vvc.set(3) = aa(2).dsdx() - aa(1).dsdy() ; 
	
	cout << "Cartesian components : vvc : " << endl ;
	vvc.spectral_display() ; 
	arrete() ; 

	Vector_divfree vvs = vvc ; 
	vvs.change_triad( map.get_bvect_spher() ) ; 
	
	cout << "Spherical components : vvs : " << endl ;
	vvs.spectral_display() ; 
	arrete() ; 

	vvs.inc_dzpuis(2) ; 
	
	cout << "vvs after inc_dzpuis(2) : " << endl ;
	vvs.spectral_display() ; 

	cout << "mu : " << endl ; 
	cout << "----" << endl ; 
	vvs.mu().spectral_display() ; 

	Vector_divfree wws = vvs.poisson() ; 
	
	cout << "Solution wws : " << endl ;
	wws.spectral_display() ; 
	arrete() ; 

	Vector wwc = wws ;
	wwc.change_triad( map.get_bvect_cart() ) ;

	cout << "Max divergence wwc : " << max( abs(wwc.divergence(metc)) ) << endl ; 
	cout << "Max divergence wws : " << max( abs(wws.divergence(mets)) ) << endl ; 
	arrete() ; 
	
	Vector vdiff(map, CON, map.get_bvect_cart() ) ;
	vdiff.set(1) = wwc(1).laplacian(2) - vvc(1) ; 		// dzpuis = 2
	vdiff.set(2) = wwc(2).laplacian(2) - vvc(2) ; 
	vdiff.set(3) = wwc(3).laplacian(2) - vvc(3) ; 

	cout << "vdiff : " << endl ;
	vdiff.spectral_display() ; 
	
	cout << "Max of vdiff : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout <<  max(abs(vdiff(i))) << endl ; 
	}
	
	return EXIT_SUCCESS ; 
}
