/*
 *  Code for testing the classes Sym_tensor, 
 *   in particular the transverse decomposion of a symmetric tensor.
 *
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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

char test_sym_tensor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/02/09 13:01:35  e_gourgoulhon
 * Added test of the TT decomposition.
 *
 * Revision 1.2  2003/12/10 10:18:24  e_gourgoulhon
 * First operational version.
 *
 * Revision 1.1  2003/11/27 16:02:24  e_gourgoulhon
 * First version of code test_sym_tensor
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
	int nr = 17 ; 	// Number of collocation points in r in each domain
	int nt = 9 ; 	// Number of collocation points in theta in each domain
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


	// Construction of a symmetric tensor field 
	// -----------------------------------------

	const Coord& x = map.x ; 
	const Coord& y = map.y ; 
	const Coord& z = map.z ; 
	const Coord& r = map.r ; 
	const Coord& cost = map.cost ; 
	const Coord& sint = map.sint ; 
	const Coord& cosp = map.cosp ; 
	const Coord& sinp = map.sinp ; 

	Sym_tensor hhc(map, CON, map.get_bvect_cart()) ; 
	
	
	Scalar tmp(map) ; 
		
    hhc.set(1,1) = 2 ; 
    tmp = 1 / (r*r*r*r)  ; 
	hhc.set(1,1).set_domain(nzm1) = tmp.domain(nzm1)  ; 
    hhc.set(1,1).set_dzpuis(4) ; 

    hhc.set(1,2) = x ;
    tmp = sint*cosp / (r*r) ; 
	hhc.set(1,2).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(1,2).set_dzpuis(4); 	 

	hhc.set(1,3) = y * z ; 
	tmp = sint*sinp*cost / (r*r*r) ;
	hhc.set(1,3).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(1,3).set_dzpuis(4); 	 

	hhc.set(2,2) = x*y ; 
	tmp = sint*sint*cosp*sinp / (r*r*r) ;
	hhc.set(2,2).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(2,2).set_dzpuis(4); 	 

	hhc.set(2,3) = z  ; 
	tmp = cost / (r*r) ;
	hhc.set(2,3).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(2,3).set_dzpuis(4); 	 

	hhc.set(3,3) = x ; 
	tmp = sint*cosp / (r*r) ;
	hhc.set(3,3).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(3,3).set_dzpuis(4); 	 




/*
	hhc.set(1,1) = 1. / (1 + pow(r,4)) ; ;
	tmp = 1. / (1. + 1./pow(r,4)) ;
	hhc.set(1,1).set_domain(nzm1) = tmp.domain(nzm1) ; 
	hhc.set(1,1).set_dzpuis(4); 	 

	hhc.set(1,2) = x / (1 + pow(r,5)) ; 
	tmp = sint*cosp / (1. + 1./pow(r,5)) ;
	hhc.set(1,2).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(1,2).set_dzpuis(4); 	 

	hhc.set(1,3) = y * z  / (1 + pow(r,6)) ; 
	tmp = sint*sinp*cost / (1. + 1./pow(r,6)) ;
	hhc.set(1,3).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(1,3).set_dzpuis(4); 	 

	hhc.set(2,2) = x*y / (1 + pow(r,6)) ; 
	tmp = sint*sint*cosp*sinp / (1. + 1./pow(r,6)) ;
	hhc.set(2,2).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(2,2).set_dzpuis(4); 	 

	hhc.set(2,3) = z / (1 + pow(r,5)) ; 
	tmp = cost / (1. + 1./pow(r,5)) ;
	hhc.set(2,3).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(2,3).set_dzpuis(4); 	 

	hhc.set(3,3) = x / (1 + pow(r,5)) - y*z*z / (1 + pow(r,7)) ; 
	tmp = sint*cosp / (1. + 1./pow(r,5))
		- sint*sinp*cost*cost / (1. + 1./pow(r,7)) ;
	hhc.set(3,3).set_domain(nzm1) = tmp.domain(nzm1)  ; 
	hhc.set(3,3).set_dzpuis(4); 	 
*/

	hhc.std_spectral_base() ; 

	Sym_tensor hhs = hhc ; 
	hhs.change_triad( map.get_bvect_spher() ) ; 
	
	cout << "Spherical components : hhs : " << endl ;
	hhs.spectral_display() ; 
	arrete() ; 
	
    // Transverse part
    // ---------------
    
	const Sym_tensor_trans& thhs = hhs.transverse(mets) ; 
	arrete() ; 

	cout << "Transverse part thhs : " << endl ;
	thhs.spectral_display() ; 
	arrete() ; 
	
    cout << "Max abs(trace) of the transverse part thhs : \n" ; 
    maxabs( thhs.trace() ) ; 
    
	cout << "Divergence of the transverse part thhs : " << endl ;
    Vector divt = thhs.divergence(mets) ; 
	// divt.spectral_display() ; 
	cout << "Max abs(divergence) of the transverse part thhs : " << endl ;
    maxabs( divt ) ; 
    
    // TT part
    // -------
    

    const Sym_tensor_tt& tthhs = thhs.tt_part() ; 
    
    cout << "TT part tthhs : \n" ; 
    // tthhs.spectral_display() ; 

    Scalar tr_tthhs = (tthhs.down(1, mets)).scontract(0,1) ;
	cout << "Max abs(trace) of the TT part tthhs : " << endl ;
    maxabs( tr_tthhs ) ; 
    
    Vector divtt = tthhs.divergence(mets) ; 
	cout << "Max abs(divergence) of the TT part tthhs : " << endl ;
    maxabs( divtt ) ; 
		
        
	return EXIT_SUCCESS ; 
}
