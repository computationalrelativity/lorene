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

 


/*
 * $Id$
 * $Log$
 * Revision 1.8  2016/12/05 16:18:29  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:54:02  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2005/04/06 14:38:41  j_novak
 * Added method 6: block inversion with spherical components
 *
 * Revision 1.5  2004/12/29 12:25:51  j_novak
 * non-symmetric grid in phi.
 *
 * Revision 1.4  2004/07/27 09:40:34  j_novak
 * test of method 5.
 *
 * Revision 1.3  2004/05/07 15:33:24  j_novak
 * Treated the case where all components are null.
 *
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

using namespace Lorene ;

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	const int nz = 3 ; 	// Number of domains
	int nr =33 ; 	// Number of collocation points in r in each domain
	int nt =9 ; 	// Number of collocation points in theta in each domain
	int np = 16 ; 	// Number of collocation points in phi in each domain
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

	double lambda = 1./3. ;

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
	Vector sol_num4 = source.poisson(lambda, 4) ;
	Vector sol_num5 = source.poisson(lambda, 5) ;
	Vector sol_num6 = source.poisson(lambda, 6) ;

	cout << endl ;

	cout << "==================================================" << endl ;
	cout << " ---------    Erreur sur la solution    --------- " << endl ;
	cout << "==================================================" << endl ;


	maxabs(sol_num0 - sol, "Methode 0 (Comp. spheriques)") ;
	maxabs(sol_num1 - sol, "Methode 1 (Comp. spheriques)") ;
	maxabs(sol_num2 - sol, "Methode 2 (Comp. cartesiennes)") ;
	maxabs(sol_num3 - sol, "Methode 3 (Comp. cartesiennes)") ;
	maxabs(sol_num4 - sol, "Methode 4 (Comp. spheriques)") ;
	maxabs(sol_num5 - sol, "Methode 5 (Comp. spheriques)") ;
	maxabs(sol_num6 - sol, "Methode 6 (Comp. spheriques)") ;


	return EXIT_SUCCESS ;
}
