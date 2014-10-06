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

char test_vector_df_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2014/10/06 15:12:56  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2008/08/27 08:54:28  jl_cornou
 * Added test for angular potential A
 *
 * Revision 1.6  2003/10/29 11:06:50  e_gourgoulhon
 * dec2_dzpuis() replaced by dec_dzpuis(2).
 *
 * Revision 1.5  2003/10/28 21:37:02  e_gourgoulhon
 * Change the name of the global char[] for code identification.
 *
 * Revision 1.4  2003/10/20 10:14:08  e_gourgoulhon
 * Added test with more general vectors.
 *
 * Revision 1.3  2003/10/19 20:05:07  e_gourgoulhon
 * Change of the argument list of Scalar::spectral_display
 * (cout now default).
 *
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
#include <cstdlib>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"


int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	int nz = 3 ; 	// Number of domains
	int nzm1 = nz - 1 ;  
	int nr = 17 ; 	// Number of collocation points in r in each domain
	int nt = 17 ; 	// Number of collocation points in theta in each domain
	int np = 32 ; 	// Number of collocation points in phi in each domain
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
	vvc.set(2) = 1 ; 
	vvc.set(3) = 0 ; 
	vvc.std_spectral_base() ; 

	Vector_divfree vvs = vvc ; 
	vvs.change_triad( map.get_bvect_spher() ) ; 
	
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
	
	cout << "Norme divergence vvc : " << norme( vvc.divergence(metc) ) << endl ; 
	cout << "Norme divergence vvs : " << norme( vvs.divergence(mets) ) << endl ; 
	cout << "Max divergence vvc : " << max( abs(vvc.divergence(metc)) ) << endl ; 
	cout << "Max divergence vvs : " << max( abs(vvs.divergence(mets)) ) << endl ; 
	arrete() ; 
		
	Vector_divfree vvs2(map, map.get_bvect_spher(), mets ) ;
	vvs2.set_vr_eta_mu(vvs(1), vvs.eta(), vvs.mu()) ; 
	Vector diff = vvs - vvs2 ; 
	cout << "diff : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << "Component " << i << " : " << endl ; 
		diff(i).spectral_display() ; 
		cout << "  norme: " << norme(diff(i)) << endl ; 
		arrete() ; 
	}
	

	cout << "========================================================" << endl ;
	cout << "                Test with a rotation vector" << endl ;
	cout << "                Cartesian comp.: V^i = (-y, x, 0)" << endl ; 
	cout << "========================================================" << endl ;

	vvc.set(1) = -y ; 
	vvc.set(2) = x ; 
	vvc.set(3) = 0 ; 
	vvc.annule_domain(nzm1) ; 
	vvc.std_spectral_base() ; 
	vvs.set_triad( map.get_bvect_cart() ) ;
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
	
	cout << "Norme divergence vvc : " << norme( vvc.divergence(metc) ) << endl ; 
	cout << "Norme divergence vvs : " << norme( vvs.divergence(mets) ) << endl ; 
	cout << "Max divergence vvc : " << max( abs(vvc.divergence(metc)) ) << endl ; 
	cout << "Max divergence vvs : " << max( abs(vvs.divergence(mets)) ) << endl ; 
	arrete() ; 
	
	vvs2.set_vr_eta_mu(vvs(1), vvs.eta(), vvs.mu()) ; 
	diff = vvs - vvs2 ; 
	cout << "diff : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << "Component " << i << " : " << endl ; 
		diff(i).spectral_display() ; 
		cout << "  norme: " << norme(diff(i)) << endl ; 
		arrete() ; 
	}
	

	cout << "========================================================" << endl ;
	cout << "                Test with a pretty general vector" << endl ;
	cout << "                           V = curl(A) " << endl ; 
	cout << "========================================================" << endl ;

	Vector aa(map, CON, map.get_bvect_cart()) ; 
	aa.set(1) = x + x*y - 3*z*z ; 
	aa.set(2) = z*z*x - 2 *y + 1 ; 
	aa.set(3) = z* ( 1 + x*x - y + x + x*y + z*z ) ; 
	aa.annule_domain(nzm1) ; 


 	Mtbl tced = sint*cosp / (r*r) + cost*cost / (r*r*r) ; 
 	aa.set(1).set_domain(nzm1) = tced(nzm1) ; 
 
 	tced =  1 / (r*r) + cost*cost*sint*sinp / (r*r*r) ; 
 	aa.set(2).set_domain(nzm1) = tced(nzm1) ; 
 
 	tced = cost / (r*r) + cost*sint*cosp / (r*r*r) ; 
 	aa.set(3).set_domain(nzm1) = tced(nzm1) ; 

	aa.std_spectral_base() ; 
	
	// cout << "aa : " << endl ; 
	// aa.spectral_display(1.e-14) ; 
	// arrete() ; 


	// Curl of aa:
	vvc.set(1) = aa(3).dsdy() - aa(2).dsdz() ; 
	vvc.set(2) = aa(1).dsdz() - aa(3).dsdx() ; 
	vvc.set(3) = aa(2).dsdx() - aa(1).dsdy() ; 
	vvc.dec_dzpuis(2) ; 
	
	vvs.set_triad( map.get_bvect_cart() ) ;
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
	
	cout << "Norme divergence vvc : " << norme( vvc.divergence(metc) ) << endl ; 
	cout << "Norme divergence vvs : " << norme( vvs.divergence(mets) ) << endl ; 
	cout << "Max divergence vvc : " << max( abs(vvc.divergence(metc)) ) << endl ; 
	cout << "Max divergence vvs : " << max( abs(vvs.divergence(mets)) ) << endl ; 

	Valeur vdiv = vvs.divergence(mets).get_spectral_va() ; 
	vdiv.coef() ; 
	cout << "Max cf vdiv : " << max(abs(*(vdiv.c_cf))) << endl ; 
	// vdiv.display_coef(1.e-12) ;

	arrete() ; 
	
	vvs2.set_vr_eta_mu(vvs(1), vvs.eta(), vvs.mu()) ; 
	diff = vvs - vvs2 ; 
	cout << "diff : " << endl ; 
	for (int i=1; i<=3; i++) {
		cout << "Component " << i << " : " << endl ; 
		diff(i).spectral_display() ; 
		cout << "  norme: " << norme(diff(i)) << endl ; 
		arrete() ; 
	}
	

	return EXIT_SUCCESS ; 
}
