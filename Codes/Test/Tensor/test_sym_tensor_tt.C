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
 * Revision 1.10  2003/11/06 10:33:41  e_gourgoulhon
 * Added tests for the mu part.
 *
 * Revision 1.9  2003/11/05 15:32:12  e_gourgoulhon
 * Use of new functions maxabs and diffrelmax acting on tensors.
 *
 * Revision 1.8  2003/11/04 23:04:43  e_gourgoulhon
 * First test of Sym_tensor_tt::set_rr_eta_mu.
 *
 * Revision 1.7  2003/11/04 15:01:16  e_gourgoulhon
 * Pretty general test for eta.
 *
 * Revision 1.6  2003/11/03 22:37:12  e_gourgoulhon
 * New version.
 *
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
	int nt = 7 ; 	// Number of collocation points in theta in each domain
	int np = 16 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi
	bool compact = true ; // external domain is compactified
  
	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double r_limits[] = {0., 2., 3., __infinity} ; 
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
	const Coord& phi = map.phi ; 

	
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

	cout << "Cartesian components : hhc : " << endl ;
	hhc.spectral_display() ; 
	arrete() ; 

	Tensor tmp = hhc ; 
	tmp.change_triad( map.get_bvect_spher() ) ; 
	Sym_tensor_trans hhs(map, map.get_bvect_spher(), mets ) ; 
	hhs = tmp ; 
	
	cout << "Spherical components : hhs : " << endl ;
	hhs.spectral_display() ; 
	arrete() ; 
	
	cout << "Maxabs divergence hhc : " << endl ; 
	maxabs( hhc.divergence(metc) ) ; 
	cout << "Maxabs divergence hhs : " << endl ; 
	maxabs( hhs.divergence(mets) ) ; 
		
	Scalar trace =  hhc.trace() ; 
	cout << "Maxabs Trace of hhc : " << endl ; 	
	maxabs( trace ) ; 	
		
	cout << "Maxabs Trace of hhs : " << endl ; 	
	maxabs( hhs.trace() ) ; 

	cout << "========================================================" << endl ;
	cout << "                Test with the tensor" << endl ;
	cout << "   Cart. comp. h^{ij} = d_i d_j Phi  with Lap(Phi) = 0 " << endl ; 
	cout << "========================================================" << endl ;
	
	arrete() ; 
	
	Scalar pot(map) ; 
	pot =  1 					// P_0^0
		   + x					// P_1^1 cos(p)
		   + y					// P_1^1 sin(p)
		   + (3*z*z - r*r) 		// P_2^0
		   + (x*x - y*y ) 		// P_2^2 cos(2p)
		   + x*y 	 			// P_2^2 sin(2p)
		   + x*(5*z*z-r*r)		// P_3^1 cos(p)
		   + y*(5*z*z-r*r)		// P_3^1 sin(p)
		   + x*(x*x-3*y*y)		// P_3^3 cos(3p)
		   + y*(y*y-3*x*x) ;	// P_3^3 sin(3p)
		   
	pot.annule_domain(nzm1) ; 

	Mtbl potced = 1 / r 		// P_0^0
		   + sint*cosp / pow(r,2)	// P_1^1 cos(p)
		   + sint*sinp / pow(r,2)	// P_1^1 sin(p)
		   + (3*cost*cost - 1) / pow(r,3)			// P_2^0
		   + sint*sint*cos(2*phi) / pow(r,3)		// P_2^2 cos(2p)
		   + sint*sint*sin(2*phi) / pow(r,3)	 	// P_2^2 sin(2p)
		   + sint*(15*cost*cost-3)*cosp / pow(r,4)		// P_3^1 cos(p)
		   + sint*(15*cost*cost-3)*sinp / pow(r,4)		// P_3^1 sin(p)
		   + pow(sint,3) * cos(3*phi) / pow(r,4)		// P_3^3 cos(3p)
		   + pow(sint,3) * sin(3*phi) / pow(r,4) ;		// P_3^3 sin(3p)
		   
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
	
	cout << "Maxabs divergence hhc : " << endl ; 
	maxabs( divc ) ; 
	
	arrete() ;
	
	Vector divs = hhs.divergence(mets) ; 
	cout << "Divergence of hhs : " << endl ; 
	divs.spectral_display() ; 
	
	cout << "Maxabs divergence hhs : " << endl ; 
	maxabs( divs ) ; 
		
	cout << "Maxabs of trace of hhc : " << endl ; 
	maxabs( hhc.trace() ) ; 
		
	cout << "Maxabs of trace of hhs : " << endl ; 
	maxabs( hhs.trace() ) ; 
		
	arrete() ; 
	
	Sym_tensor_tt  htt(map, map.get_bvect_spher(), mets ) ; 
	
	htt = hhs ; 
	
	cout << "Eta : " << endl ; 
	htt.eta().spectral_display() ; 

	cout << "Mu : " << endl ; 
	htt.mu().spectral_display() ; 
	
	
	arrete() ; 
	
	Sym_tensor_tt htt2(map, map.get_bvect_spher(), mets ) ; 
		
	htt2.set_rr_eta_mu(htt(1,1), htt.eta(), htt.mu()) ; 

	cout << "Relative difference (max) between htt2 and htt : " << endl ; 
	diffrelmax(htt2, htt) ; 

	cout << endl << "========================================================" << endl ;
	cout << "                Test with an ad hoc mu" << endl ;
	cout << "========================================================" << endl ;
	
	arrete() ; 
	
	Scalar mu0(map) ; 
	mu0 = z * pot.get_spectral_va() ;
		   
	mu0.annule_domain(nzm1) ; 

	Mtbl mu0ced = cost * potced  ;	
		   
	mu0.set_domain(nzm1) = mu0ced(nzm1) ; 
	mu0.std_spectral_base() ; 
	mu0.set_spectral_va().set_base_r(0, R_CHEBPIM_I) ; 
	mu0.set_spectral_va().set_base_t(T_COSSIN_CI) ; 
	

	cout << "mu0 : " << endl ; 
	mu0.spectral_display() ; 
	arrete() ; 

	Sym_tensor_tt htt3 = htt ; 
		
	htt3.set_rr_mu(htt(1,1), mu0) ; 

	cout << "Eta htt3 : " << endl ; 
	htt3.eta().spectral_display() ; 

	cout << "Mu htt3 : " << endl ; 
	htt3.mu().spectral_display() ; 
	
	arrete() ; 

	cout << "htt3 : " << endl ; 
	htt3.spectral_display() ; 
	arrete() ; 

	Sym_tensor_tt htt4(map, map.get_bvect_spher(), mets ) ; 
		
	htt4.set_rr_eta_mu(htt3(1,1), htt3.eta(), htt3.mu()) ; 

	cout << "Relative difference (max) between htt4 and htt3 : " << endl ; 
	diffrelmax(htt4, htt3) ; 

	arrete() ; 

	Sym_tensor_tt htt5 = htt3 ; 
	htt5.set(1,1) = htt3(1,1) ; // To force the deletion of eta and mu

	cout << "Relative difference (max) between htt5.eta and htt3.eta : " << endl ; 
	diffrelmax(htt5.eta(), htt3.eta(), cout) ; 
	cout << "Relative difference (max) between htt5.mu and htt3.mu : " << endl ; 
	diffrelmax(htt5.mu(), htt3.mu(), cout) ; 

		
	return EXIT_SUCCESS ; 
}
