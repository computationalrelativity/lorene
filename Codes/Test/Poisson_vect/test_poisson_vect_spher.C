/*
 *  Code for testing the vector Poisson equation in spherical components.
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

char test_poisson_vect_spher_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/01/23 13:27:16  e_gourgoulhon
 * Scalar::set_val_inf --> Scalar::set_outer_boundary.
 *
 * Revision 1.4  2003/12/19 15:17:45  j_novak
 * *** empty log message ***
 *
 * Revision 1.3  2003/11/03 15:13:10  j_novak
 * Comparison with Philippe's Poisson solver.
 *
 * Revision 1.2  2003/10/29 11:05:44  e_gourgoulhon
 * Twice inc2_dzpuis() relaced by inc_dzpuis(4).
 *
 * Revision 1.1  2003/10/27 09:19:26  j_novak
 * Test file for vector Poisson equation in spherical components
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "metric.h"
#include "tenseur.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	const int nz = 3 ; 	// Number of domains
	int nr = 17 ; 	// Number of collocation points in r in each domain
	int nt = 5 ; 	// Number of collocation points in theta in each domain
	int np = 12 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi
    int nzm1 = nz - 1 ; 

	int nbr[] = {nr, 2*nr - 1, 2*nr - 1};
	int nbt[] = {nt, nt, nt} ;
	int nbp[] = {np, np, np} ;
	int tipe_r[] = {RARE, FIN, UNSURR} ;
  	assert( nz == 3 ) ;// since the above arrays are described in only 3 domains
  
	Mg3d mgrid(nz, nbr, tipe_r, nbt, symmetry_theta, nbp, symmetry_phi) ;
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double r_limits[] = {0., 1, 2., __infinity} ; 
  
	Map_af map(mgrid, r_limits) ; 
  	

	// Construction of flat metrics
	// ----------------------------

	Metric_flat mets(map, map.get_bvect_spher()) ; // spherical representation
	Metric_flat metc(map, map.get_bvect_cart()) ;  // Cartesian representation


	// Construction of a divergence free vector field 
	// ----------------------------------------------

	const Coord& x = map.x ; 
	const Coord& y = map.y ; 
	const Coord& z = map.z ; 
	const Coord& r = map.r ; 

	double lamda = 1./3. ; //for the vector equation Delta V + lamda grad(div V)

	Mtbl denom(mgrid) ;
	denom = r*r*r*r+ 1 ;
	Mtbl denom2 = denom*denom ;
	Mtbl denom3 = denom*denom2 ;
	
	Vector theo(map, CON, map.get_bvect_cart()) ; 
	theo.set(1) = x*x / denom ;
	theo.set(2) = 0 ;
	theo.set(3) =  0 ;
	theo.set(1).set_outer_boundary(nzm1, 0.) ;
	theo.annule_domain(nz-1) ;
	theo.std_spectral_base() ; 

	Vector vvc(map, CON, map.get_bvect_cart()) ;
	vvc.set(1) = 2*(1+lamda)/denom 
	  - ((36+20*lamda)*r*r + 8*lamda*x*x)*x*x/denom2
	  + (r*r+lamda*x*x)*32*r*r*r*r*x*x/denom3 ;
	vvc.set(1).set_outer_boundary(nzm1, 0.) ; ;
	vvc.set(2) = lamda*(-8*x*y*(r*r + x*x)/denom2 
			    + 32*r*r*r*r*x*x*x*y/denom3) ;
	vvc.set(2).set_outer_boundary(nzm1, 0.) ;
	vvc.set(3) = lamda*(-8*x*z*(r*r + x*x)/denom2 
			    + 32*r*r*r*r*x*x*x*z/denom3) ;
	vvc.set(3).set_outer_boundary(nzm1, 0.) ;
	vvc.annule_domain(nz-1) ;
	vvc.std_spectral_base() ;

  	Vector vvs = vvc ; 
  	vvs.change_triad( map.get_bvect_spher() ) ;
  	vvs.inc_dzpuis(4) ;
	Tenseur vvc_p(map, 1, CON, map.get_bvect_spher() ) ;
	vvc_p.set_etat_qcq() ;
	for (int i=0; i<3; i++) {
	  Cmp tmp_p(vvs(i+1)) ;
	  vvc_p.set(i) = tmp_p ;
	}
	vvc_p.change_triad(map.get_bvect_cart()) ;
	Tenseur vect_auxi (map, 1, CON, map.get_bvect_cart()) ;
	vect_auxi.set_etat_qcq() ;
	Tenseur scal_auxi (map) ;
	scal_auxi.set_etat_qcq() ;
	
	Tenseur resu_p(vvc_p.poisson_vect(lamda, vect_auxi, scal_auxi)) ;
	Tenseur theo_p(map, 1, CON, map.get_bvect_cart()) ;
	theo_p.set_etat_qcq() ;
	for (int i=0; i<3; i++) {
	  Cmp tmp_p(theo(i+1)) ;
	  theo_p.set(i) = tmp_p ;
	}

 	Vector resus = vvs.poisson(lamda) ;
 	Vector resu = resus ;
 	resu.change_triad(map.get_bvect_cart() ) ;

	cout << "Max of relative difference (in cartesian components): " << endl ; 
	for (int i=1; i<=3; i++) {
	  cout << "Component " << i << ": " << endl ;
      	  cout <<  "New version: " << diffrelmax(resu(i), theo(i)) << endl ; 
	  cout <<  "Grandclement et al.: " << 
	    diffrelmax(resu_p(i-1), theo_p(i-1)) << endl ; 
	  des_profile(resu_p(i-1) - theo_p(i-1), 0, 3, 1, 1) ;
	}

	return EXIT_SUCCESS ; 
}
