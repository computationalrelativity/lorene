/*
 *  Test of the resolution of the vector Poisson equation
 *
 */

/*
 *   Copyright (c) 2004  Jerome Novak
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

char test_poisson_vect_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/02/22 15:48:49  j_novak
 * New code for testing vector Poisson equation.
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

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	const int nz = 3 ; 	// Number of domains
	int nr = 33 ; 	// Number of collocation points in r in each domain
	int nt = 5 ; 	// Number of collocation points in theta in each domain
	int np = 6 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi

	int nbr[] = {nr, nr, nr};
	int nbt[] = {nt, nt, nt} ;
	int nbp[] = {np, np, np} ;
	int tipe_r[] = {RARE, FIN, UNSURR} ;
  	assert( nz == 3 ) ;// since the above arrays are described in only 3 domains
  
	Mg3d mgrid(nz, nbr, tipe_r, nbt, symmetry_theta, nbp, symmetry_phi) ;

	int nzm1 = nz - 1 ;
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	double R = 2.1 ;

	// Boundaries of each domains
	double r_limits[] = {0., 0.5*R, R, __infinity} ; 
  
	Map_af map(mgrid, r_limits) ; 
  	

	// Construction of flat metrics
	// ----------------------------

	Metric_flat mets(map, map.get_bvect_spher()) ; // spherical representation
	Metric_flat metc(map, map.get_bvect_cart()) ;  // Cartesian representation


	Scalar xx(map), yy(map), zz(map), rr(map) ;
	xx = map.x ;
	yy = map.y ;
	zz = map.z ;
	rr = map.r ;

	double lamda = 1. ; //for the vector equation Delta V + lamda grad(div V)
	int nn ;
	cout << "Enter n : (source decaying like 1/r^n+4, so n>0 is required)" 
	     << endl ;
	cin >> nn ;

	  Vector theo(map, CON, map.get_bvect_cart()) ; 
	  Vector theo_zec(map, CON, map.get_bvect_cart()) ; 
	  theo.set(1) = xx*rr*rr / (10 * (lamda+1) * pow(R, double(nn+5))) 
	    - double(nn + 5) *xx / 
	    (6 * (lamda + 1) * double(nn+3) * pow(R, double(nn+3))) ;
	  theo.set(1).annule_domain(nzm1) ;
	  theo_zec.set(1) = xx / ((lamda+1) * nn * (nn+3) * pow(rr, nn+3))
	    - (nn+5) * xx / ((lamda+1) * 15*nn * pow(R,double(nn)) * pow(rr,3)) ;
	  theo_zec.set(1).annule(0,nz-2) ;
	  theo.set(1) += theo_zec(1) ;
	  theo.set(1).set_outer_boundary(nz-1, 0.) ;
	  theo.set(2) = yy*rr*rr / (10 * (lamda+1) * pow(R, double(nn+5))) 
	    - double(nn + 5) *yy / (6 * (lamda + 1) * double(nn+3) * pow(R, double(nn+3))) ;
	  theo.set(2).annule_domain(nzm1) ;
	  theo_zec.set(2) = yy / ((lamda+1) * nn * (nn+3) * pow(rr, nn+3))
	    - (nn+5) * yy / ((lamda+1) * 15*nn * pow(R,double(nn)) * pow(rr,3)) ;
	  theo_zec.set(2).annule(0,nz-2) ;
	  theo.set(2) += theo_zec(2) ;
	  theo.set(2).set_outer_boundary(nz-1, 0.) ;
	  theo.set(3) = zz*rr*rr / (10 * (lamda+1) * pow(R, double(nn+5))) 
	    - double(nn + 5) *zz / (6 * (lamda + 1) * double(nn+3) * pow(R, double(nn+3))) ;
	  theo.set(3).annule_domain(nzm1) ;
	  theo_zec.set(3) = zz / ((lamda+1) * nn * (nn+3) * pow(rr, nn+3))
	    - (nn+5) * zz / ((lamda+1) * 15*nn * pow(R,double(nn)) * pow(rr,3)) ;
	  theo_zec.set(3).annule(0,nz-2) ;
	  theo.set(3) += theo_zec(3) ;
	  theo.set(3).set_outer_boundary(nz-1, 0.) ;
	  theo.std_spectral_base() ; 

	  Vector source(map, CON, map.get_bvect_cart()) ;
	  Vector source_zec(map, CON, map.get_bvect_cart()) ;
	  source.set(1) = xx / pow(R, double(nn+5)) ;
	  source.set(1).annule_domain(nzm1) ;
	  source_zec.set(1) = xx / pow(rr, nn+1) ;
	  source_zec.set(1).set_dzpuis(4) ;
	  source_zec.set(1).annule(0, nz-2) ;
	  source.set(1) += source_zec(1) ;
	  source.set(1).set_outer_boundary(nz-1, 0.) ;
	  source.set(2) = yy / pow(R, double(nn+5)) ;
	  source.set(2).annule_domain(nzm1) ;
	  source_zec.set(2) = yy / pow(rr, nn+1) ;
	  source_zec.set(2).set_dzpuis(4) ;
	  source_zec.set(2).annule(0, nz-2) ;
	  source.set(2) += source_zec(2) ;
	  source.set(2).set_outer_boundary(nz-1, 0.) ;
	  source.set(3) = zz / pow(R, double(nn+5)) ;
	  source.set(3).annule_domain(nzm1) ;
	  source_zec.set(3) = zz / pow(rr, nn+1) ;
	  source_zec.set(3).set_dzpuis(4) ;
	  source_zec.set(3).annule(0, nz-2) ;
	  source.set(3) += source_zec(3) ;
	  source.set(3).set_outer_boundary(nz-1, 0.) ;
	  source.std_spectral_base() ;

	  Vector source_s = source ; 
	  source_s.change_triad( map.get_bvect_spher() ) ;
	  Tenseur source_p(map, 1, CON, map.get_bvect_cart() ) ;
	  source_p.set_etat_qcq() ;
	  for (int i=0; i<3; i++) {
	    source_p.set(i) = Cmp(source(i+1)) ;
	  }
	  Tenseur vect_auxi (map, 1, CON, map.get_bvect_cart()) ;
	  vect_auxi.set_etat_qcq() ;
	  Tenseur scal_auxi (map) ;
	  scal_auxi.set_etat_qcq() ;
	
	  Tenseur resu_p(source_p.poisson_vect(lamda, vect_auxi, scal_auxi)) ;
	  Tenseur theo_p(map, 1, CON, map.get_bvect_cart()) ;
	  theo_p.set_etat_qcq() ;
	  for (int i=0; i<3; i++) {
	    theo_p.set(i) = Cmp(theo(i+1)) ;
	  }

	  Vector resus = source_s.poisson(lamda, mets) ;
	  Vector resu = resus ;
	  resu.change_triad(map.get_bvect_cart() ) ;

 	  cout << "Max of relative difference (in cartesian components): " << endl ; 
 	  for (int i=1; i<=3; i++) {
 	    cout << "Component " << i << ": " << endl ;
 	    cout <<  "New version: " << diffrelmax(resu(i), theo(i)) << endl ; 
 	    cout <<  "Grandclement et al.: " << 
 	      diffrelmax(resu_p(i-1), theo_p(i-1)) << endl ; 
 	  }

	return EXIT_SUCCESS ; 
}
