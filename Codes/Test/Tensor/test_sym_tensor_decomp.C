/*
 *  Code for testing the decomposition of a Sym_tensor in terms of eta/mu/W/X/T
 * and its transcription for the Poisson operator.
 *
 */

/*
 *   Copyright (c) 2005 Jerome Novak
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
 * Revision 1.4  2016/12/05 16:18:30  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2014/10/13 08:54:03  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2005/09/08 07:40:56  j_novak
 * Update of set_auxiliary arguments.
 *
 * Revision 1.1  2005/05/18 15:31:25  j_novak
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

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
	int nr = 33 ; 	// Number of collocation points in r in each domain
	int nt = 17 ; 	// Number of collocation points in theta in each domain
	int np = 32 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = SYM ; // no symmetry in phi
	bool compact = true ; // external domain is compactified
  
	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double R = 2. ;
	double r_limits[] = {0., 0.5*R, R, __infinity} ; 
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
	const Coord& r = map.r ; 

	Sym_tensor hhc(map, CON, map.get_bvect_cart()) ; 
	for (int i=1; i<=3; i++) 
	    for (int j=1; j<=i; j++)
		hhc.set(i,j) = 0 ;
	
	Scalar tmp(map) ;
	tmp = r*r*r*r*cos(x*y)*exp(-r*r)/(1+r*r) ;
	tmp.set_outer_boundary(nzm1, 0.)  ;
	tmp.annule(0,nz-2) ;
	hhc.set(1,1) = cos(x*y)*exp(-r*r)/(1+r*r) ;
	hhc.set(1,1).annule_domain(nzm1) ;
	hhc.set(1,1) += tmp ;
	for (int i=1; i<=3; i++)
	    for (int j=i; j<=3; j++)
		hhc.set(i,j).set_dzpuis(4) ;

	hhc.std_spectral_base() ;
	
	Sym_tensor hhs = hhc ; 
	hhs.change_triad( map.get_bvect_spher() ) ; 
    
	Scalar eta = hhs.eta() ; eta.div_r_dzpuis(4) ;
	Scalar mu = hhs.mu() ; mu.div_r_dzpuis(4) ;
	Scalar ww = hhs.www() ;
	Scalar xx = hhs.xxx() ;
	Scalar tt = hhs.ttt() ;

	// Test of the eta/mu/W/X decomposition
	//-------------------------------------
	Sym_tensor hh_new(map, CON, map.get_bvect_spher()) ;
	hh_new.set_auxiliary(hhs(1,1), eta, mu, ww, xx, tt) ;
	
	cout << "Test of the decomposition on eta/mu/W/X potentials : " << endl ;
	for (int i=1; i<=3; i++) 
	    for (int j=i; j<=3; j++) {
		Scalar raz = (hhs(i,j) - hh_new(i,j)) ;
		cout << "Component: " << i << " , " << j << ": " ;
		maxabs(raz) ;
		cout  << endl ;
	    }
	arrete() ;

    
	Sym_tensor thhs = hhs ;
	thhs.dec_dzpuis(4) ;

	Sym_tensor thhc = thhs ;
	thhc.change_triad(map.get_bvect_cart()) ;
	
	
	Sym_tensor sou_cc(map, CON, map.get_bvect_cart()) ; 
	for (int i=1; i<=3; i++)
	    for (int j=i; j<=3; j++) 
		sou_cc.set(i,j) = thhc(i,j).laplacian(4) ;

	Sym_tensor sou_ss = sou_cc ;
	sou_ss.change_triad(map.get_bvect_spher()) ;
	
	maxabs(thhc.divergence(metc), "Divergence cart. comp.") ;
	maxabs(thhs.divergence(mets), "Divergence spher. comp.") ;
	
	maxabs(sou_cc.divergence(metc), "Divergence cart. comp. for the source") ;
	maxabs(sou_ss.divergence(mets), "Divergence spher. comp. for the source") ;
	

	arrete() ;
	
	cout << "Test of the divergence formulas : " << endl ;
	Scalar dive(map) ;
	dive = thhs.eta().lapang() ;
	dive.div_r() ;
	dive += 2*thhs(1,1)  - thhs.ttt() ;
	dive.div_r() ;
	Scalar der = thhs(1,1).dsdr() ; der.dec_dzpuis(2) ;
	dive += der ;
	cout << "Comp. r : " << endl ;
	maxabs(dive) ;

	Scalar div_eta = thhs.eta().dsdr() ;
	div_eta.dec_dzpuis(2) ;
	tmp = 2*thhs.eta() ;
	tmp.div_r() ;
	div_eta += tmp + thhs.www().lapang() + 2*thhs.www() +0.5*thhs.ttt() ;
	cout << "Comp. eta : " << endl ;
	maxabs(div_eta.lapang()) ;


	Scalar div_mu = thhs.mu().dsdr() ;
	div_mu.dec_dzpuis(2) ;
	tmp = 2*thhs.mu() ;
	tmp.div_r() ;
	div_mu += tmp + thhs.xxx().lapang() + 2*thhs.xxx() ;
	cout << "Comp. mu : " << endl ;
	maxabs(div_mu.lapang()) ;
	
	arrete() ;
    
	// Test of the Poisson equations in terms of eta/mu/W/X
	//-----------------------------------------------------

	cout << "Test of Poisson equations in terms of eta/mu/W/X/T : " << endl ;

	Scalar trace = thhs(1,1) + thhs(2,2) + thhs(3,3) ;
	Scalar sou_trace = sou_ss(1,1) + sou_ss(2,2) + sou_ss(3,3) ;

	Scalar test(map) ;
	test = trace.laplacian(4) -sou_trace ;
	cout << "Eq. for trace : " << endl ;
	maxabs(test) ;
	cout << endl ;

	test = -4*thhs.eta().lapang() ;
	test.div_r_dzpuis(0) ;
	test -= 4*thhs(1,1)- 2*thhs.ttt() ; 
	test.div_r_dzpuis(2) ;
	test.div_r_dzpuis(4) ;
	test += thhs(1,1).laplacian(4) - sou_ss(1,1) ;
	cout << "Eq. for h^rr : " << endl ;
	maxabs(test) ;
	cout << endl ;

	test = 2*thhs(1,1) ;
	test.div_r_dzpuis(4) ;
	tmp = sou_ss.eta() ;
	tmp.inc_dzpuis() ;
	div_eta.div_r_dzpuis(4) ;
	test += thhs.eta().laplacian(4) - tmp -2*div_eta ;
	cout << "Eq. for eta : " << endl ;
	maxabs(test.lapang()) ;
	cout << endl ;

	tmp = sou_ss.mu() ;
	tmp.inc_dzpuis() ;
	div_mu.div_r_dzpuis(4) ;
	test = thhs.mu().laplacian(4) - tmp - 2*div_mu ; 
	cout << "Eq. for mu : " << endl ;
	maxabs(test) ;
	cout << endl ;

	test = 2*thhs.eta() ;
	test.div_r_dzpuis(0) ;
	test += 2*thhs.www() ;
	test.div_r_dzpuis(2) ; test.div_r_dzpuis(4) ;
	test += thhs.www().laplacian(4) - sou_ss.www() ;
	cout << "Eq. for W : " << endl ;
	maxabs(test) ;
	cout << endl ;

	test = 2*thhs.mu() ;
	test.div_r_dzpuis(0) ;
	test += 2*thhs.xxx() ;
	test.div_r_dzpuis(2) ; test.div_r_dzpuis(4) ;
	test += thhs.xxx().laplacian(4) - sou_ss.xxx() ;
	cout << "Eq. for X : " << endl ;
	maxabs(test) ;
	cout << endl ;

    return EXIT_SUCCESS ; 
}
