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
 * Revision 1.8  2004/02/18 18:54:44  e_gourgoulhon
 * Method Tensor::scontract renamed Tensor::trace.
 *
 * Revision 1.7  2003/12/27 15:03:31  e_gourgoulhon
 * New tests (new index convention of covariant derivatives).
 *
 * Revision 1.6  2003/10/19 20:05:07  e_gourgoulhon
 * Change of the argument list of Scalar::spectral_display
 * (cout now default).
 *
 * Revision 1.5  2003/10/15 10:47:01  e_gourgoulhon
 * Reorganised the arrete()'s.
 *
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

	int arret = 0 ; 

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	int nz = 3 ; 	// Number of domains
	int nr = 9 ; 	// Number of collocation points in r in each domain
	int nt = 5 ; 	// Number of collocation points in theta in each domain
	int np = 12 ; 	// Number of collocation points in phi in each domain
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
	const Coord& y = map.y ; 
	const Coord& z = map.z ; 
	const Coord& r = map.r ; 
	const Coord& sint = map.sint ; 
	const Coord& sinp = map.sinp ; 

	Scalar uu(map) ; 
	Scalar tmp(map) ; 

	uu = x ; 	
	tmp = sint * sinp / r ;
        int nzm1 = nz-1 ; 
	uu.set_domain(nzm1) = tmp.domain(nzm1) ; // y/r^2 in the external domain
	
	uu.std_spectral_base() ;   // sets the standard spectral basis for 
									// expansion of a scalar field
									
	cout << "uu : "  << endl ; 
	uu.spectral_display() ; 
	arrete(arret) ; 
	
	// Gradient of the scalar field
	// ----------------------------
	
	Vector duuc = uu.derive_cov(metc)  ;  

	cout << "duuc : "  << endl ; 
        duuc.spectral_display() ; 

	arrete(arret) ; 

	Vector duus = uu.derive_cov(mets) ; 
	
	cout << "duus : "  << endl ; 
        duus.spectral_display() ; 

	arrete(arret) ; 
	
	// Test
	// ----
	
	Vector duus_c = duus ; 
	duus_c.change_triad( map.get_bvect_cart() ) ; 
	
	Vector diffc = duus_c - duuc ; 
	
	cout << "maxabs( diffc ) : " << maxabs(diffc) << endl ; 
	arrete() ; 
	
	
	Vector duuc_s = duuc ; 
	duuc_s.change_triad( map.get_bvect_spher() ) ; 
	
	Vector diffs = duuc_s - duus ; 

	cout << "maxabs(diffs) : " << maxabs(diffs) << endl ; 

	arrete() ; 
	
	cout << endl << 
	"================ TEST FOR A VECTOR FIELD =================\n" ; 
	
	Vector vvc(map, COV, map.get_bvect_cart()) ; 
	
	vvc.set(1) = x + x*x*y - y ; 
	vvc.set(2) = x*x + z*z * y ; 
	vvc.set(3) = z + x * z ; 
        vvc.annule_domain(nzm1) ; 
        
	vvc.std_spectral_base() ; 
        
        Tensor dvvc_ana(map, 2, COV, map.get_bvect_cart()) ;
        dvvc_ana.set(1,1) = 1 + 2*x*y ; 
        dvvc_ana.set(1,2) = x*x - 1 ; 
        dvvc_ana.set(1,3) = 0 ; 
        dvvc_ana.set(2,1) = 2*x ; 
        dvvc_ana.set(2,2) = z*z ; 
        dvvc_ana.set(2,3) = 2*z*y ; 
        dvvc_ana.set(3,1) = z ; 
        dvvc_ana.set(3,2) = 0; 
        dvvc_ana.set(3,3) = 1 + x ; 
        dvvc_ana.annule_domain(nzm1) ; 
        dvvc_ana.std_spectral_base() ;  
        
	Tensor dvvc = vvc.derive_cov(metc) ;
        
        Tensor diffvvc_ana = dvvc - dvvc_ana ;
        cout << "Maxabs(dvvc - dvvc_ana) : " << maxabs(diffvvc_ana) << endl ; 
        arrete() ; 
	
	Vector vvs = vvc ; 
	vvs.change_triad( map.get_bvect_spher() ) ; 
	
	Tensor dvvs = vvs.derive_cov(mets) ; 
	
	Tensor dvvs_c = dvvs ; 
	dvvs_c.change_triad( map.get_bvect_cart() ) ; 
	
	Tensor diffvvc = dvvs_c - dvvc ; 
	
        cout << "Maxabs(dvvs_c - dvvc) : \n" ;
        maxabs(diffvvc) ; 
        arrete() ; 


	Tensor dvvc_s = dvvc ; 
	dvvc_s.change_triad( map.get_bvect_spher() ) ; 

	Tensor diffvvs = dvvc_s - dvvs ; 
	
        cout << "Maxabs(dvvc_s - dvvs) : \n" ; 
        maxabs(diffvvs) ;
        arrete() ; 
	
        // Divergence
        // ----------
        Vector uvvc = vvc.up(0, metc) ; 
        Scalar divc = uvvc.divergence(metc) ; 
        Tensor udvvc = dvvc.up(0, metc) ; 
        Scalar diffdivc = divc - Scalar( udvvc.trace(0,1) ) ; 
        cout << "Error on the divergence (Cart.): " << endl ; 
        maxabs(diffdivc) ; 
        
        Vector uvvs = vvs.up(0, mets) ; 
        Scalar divs = uvvs.divergence(mets) ; 
        Tensor udvvs = dvvs.up(0, mets) ; 
        Scalar diffdivs = divs - Scalar( udvvs.trace(0,1) ) ; 
        cout << "Error on the divergence (spher.): " << endl ; 
        maxabs(diffdivs) ; 
        
        
	cout << endl << 
	"================ TEST FOR A RANK 2 TENSOR FIELD =================\n" ; 
	        
        Tensor ddvvc_ana(map, 3, COV, map.get_bvect_cart()) ;

        ddvvc_ana.set(1,1,1) = 2*y ; 
        ddvvc_ana.set(1,1,2) = 2*x ; 
        ddvvc_ana.set(1,1,3) = 0 ; 

        ddvvc_ana.set(1,2,1) = 2*x ; 
        ddvvc_ana.set(1,2,2) = 0 ; 
        ddvvc_ana.set(1,2,3) = 0 ; 

        ddvvc_ana.set(1,3,1) = 0 ; 
        ddvvc_ana.set(1,3,2) = 0 ; 
        ddvvc_ana.set(1,3,3) = 0 ; 

        ddvvc_ana.set(2,1,1) = 2 ; 
        ddvvc_ana.set(2,1,2) = 0 ; 
        ddvvc_ana.set(2,1,3) = 0 ; 

        ddvvc_ana.set(2,2,1) = 0 ; 
        ddvvc_ana.set(2,2,2) = 0 ; 
        ddvvc_ana.set(2,2,3) = 2*z ; 

        ddvvc_ana.set(2,3,1) = 0 ; 
        ddvvc_ana.set(2,3,2) = 2*z ; 
        ddvvc_ana.set(2,3,3) = 2*y ; 

        ddvvc_ana.set(3,1,1) = 0 ; 
        ddvvc_ana.set(3,1,2) = 0 ; 
        ddvvc_ana.set(3,1,3) = 1 ; 

        ddvvc_ana.set(3,2,1) = 0; 
        ddvvc_ana.set(3,2,2) = 0; 
        ddvvc_ana.set(3,2,3) = 0; 

        ddvvc_ana.set(3,3,1) = 1 ; 
        ddvvc_ana.set(3,3,2) = 0 ; 
        ddvvc_ana.set(3,3,3) = 0 ; 

        ddvvc_ana.annule_domain(nzm1) ; 
        
        Tensor ddvvc = dvvc.derive_cov(metc) ;         
        
        Tensor diffddvvc_ana = ddvvc - ddvvc_ana ;
        cout << "Maxabs(ddvvc - ddvvc_ana) : \n" ; 
        maxabs(diffddvvc_ana) ; 
        arrete() ; 
        
	Tensor ddvvs = dvvs.derive_cov(mets) ; 
	
//	Tensor ddvvs_c = ddvvs ; 
//	ddvvs_c.change_triad( map.get_bvect_cart() ) ; 
//	Tensor diffddvvc = ddvvs_c - ddvvc ; 
//        cout << "Maxabs(ddvvs_c - ddvvc) : \n" ;
//        maxabs(diffddvvc) ; 
//        arrete() ; 

//	Tensor ddvvc_s = ddvvc ; 
//	ddvvc_s.change_triad( map.get_bvect_spher() ) ; 
//	Tensor diffddvvs = ddvvc_s - ddvvs ; 
//        cout << "Maxabs(ddvvc_s - ddvvs) : \n" ; 
//        maxabs(diffddvvs) ;
//        arrete() ; 
	
        // Divergence
        // ----------
        
        Tensor u1dvvc = dvvc.up(1, metc) ; 
        Vector div_udvvc = u1dvvc.divergence(metc) ; 
        Tensor uddvvc = ddvvc.up(1, metc) ; 
        cout << "div_udvvc : \n " << endl ; 
        div_udvvc.spectral_display() ; 
        Vector diff_div_udvvc = div_udvvc - uddvvc.trace(1,2)  ; 
        cout << "Error on the divergence (Cart.): " << endl ; 
        maxabs(diff_div_udvvc) ; 
        
	return EXIT_SUCCESS ; 
}
