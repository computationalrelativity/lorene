/*
 *  Test of the resolution of the vector Poisson equation for a given l
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

char test_vpoisson_l_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/03/26 17:05:24  j_novak
 * Added new method n.3 using Tenseur::poisson_vect_oohara
 *
 * Revision 1.1  2004/03/26 15:35:46  j_novak
 * More vector Poisson testing
 *
 *
 * $Header$
 *
 */

// Lorene headers
#include "metric.h"
#include "cmp.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	const int nz = 3 ; 	// Number of domains
	int nr =33 ; 	// Number of collocation points in r in each domain
	int nt =17 ; 	// Number of collocation points in theta in each domain
	int np = 32 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi

	int nbr[] = {nr, nr, nr};
	int nbt[] =  {nt, nt, nt} ;
	int nbp[] = {np, np, np} ;
	int tipe_r[] = {RARE, FIN, UNSURR} ;
 
	Mg3d mgrid(nz, nbr, tipe_r, nbt, symmetry_theta, nbp, symmetry_phi) ;
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double r_limits[] = {0., 1.1, 2., __infinity} ; 
  
	Map_af map(mgrid, r_limits) ; 
	int nzm1 = nz - 1 ;
  	

	// Construction of flat metrics
	// ----------------------------

	const Metric_flat& mets = map.flat_met_spher() ; 

	Scalar r(map) ;
	r = map.r ; 
	Scalar xx(map) ;
	xx = map.x ;

	cout << "Entrer l: " << endl ;
	int lq; cin >> lq ;
	
	Scalar pot(map) ;
	pot = pow(xx,lq) / (1 + pow(r, 2*lq+2) ) ;
	pot.set_outer_boundary(nzm1,0.) ;
	pot.std_spectral_base() ;
	
	Vector sol(map, CON, map.get_bvect_cart()) ;
	sol.set(1) = pot.dsdy() ;
	sol.set(1).dec_dzpuis(2) ;
	sol.set(2) = 0 ;
	sol.set(3) = 0 ;

	sol.change_triad(map.get_bvect_spher()) ;

	Vector_divfree sol_df = sol.div_free(mets) ;

	Scalar poten = sol.potential(mets) ;

	Scalar diff = sol_df.divergence(mets) ;
	
	maxabs(diff, "Divergence de sol_df") ;

	double lambda = -1. ;

	Tensor grad = sol.derive_con(mets) ;
	Scalar div = sol.divergence(mets) ;
	
	Vector source = grad.divergence(mets) + lambda * div.derive_con(mets) ;

	source.inc_dzpuis() ;

	Vector_divfree sou_df = source.div_free(mets) ;

	Scalar pot_sou = source.potential(mets) ;

	diff = sou_df.divergence(mets) ;
	diff.dec_dzpuis() ;
	diff.dec_dzpuis(4) ;

	maxabs(diff, "Divergence de sou_df") ;

	Vector sol_num0 = source.poisson(lambda, 0) ;
	Vector sol_num1 = source.poisson(lambda, 1) ;
	Vector sol_num2 = source.poisson(lambda, 2) ;
	Vector sol_num3 = source.poisson(lambda, 3) ;

	cout << endl ;

	cout << "==================================================" << endl ;
	cout << " ---------    Erreur sur la solution    --------- " << endl ;
	cout << "==================================================" << endl ;


	maxabs(sol_num0 - sol, "Methode 0 (Comp. spheriques)") ;
	maxabs(sol_num1 - sol, "Methode 1 (Comp. spheriques)") ;
	maxabs(sol_num2 - sol, "Methode 2 (Comp. cartesiennes)") ;
	maxabs(sol_num3 - sol, "Methode 3 (Comp. cartesiennes)") ;


	return EXIT_SUCCESS ;
}
