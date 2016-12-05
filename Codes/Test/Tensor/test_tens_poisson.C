/*
 *  Code for testing the tensorial Poisson equation.
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.7  2016/12/05 16:18:30  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.6  2014/10/13 08:54:03  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:12:56  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2004/02/18 18:58:36  e_gourgoulhon
 *  Method trace() renamed the_trace().
 *  Tensor::trace used instead of Tensor::scontract.
 *
 * Revision 1.3  2004/02/16 12:48:51  e_gourgoulhon
 * Methods Scalar::dsdx(),... do no longer have any argument.
 *
 * Revision 1.2  2004/01/04 21:02:31  e_gourgoulhon
 * Class Tensor_delta replaced by class Tensor_sym.
 *
 * Revision 1.1  2003/11/07 17:18:24  e_gourgoulhon
 * First version of test_tens_poisson
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
	cout << "                Test with the tensor" << endl ;
	cout << "   Cart. comp. h^{ij} = d_i d_j Phi  with Lap(Phi) = 0 " << endl ; 
	cout << "========================================================" << endl ;
		
	Scalar pot(map) ; 
	pot =  1 					// P_0^0
		   + x					// P_1^1 cos(p)
		   + y					// P_1^1 sin(p)
		   + (3*z*z - r*r) 	;	// P_2^0
/*		   + (x*x - y*y ) 		// P_2^2 cos(2p)
		   + x*y 	 			// P_2^2 sin(2p)
		   + x*(5*z*z-r*r)		// P_3^1 cos(p)
		   + y*(5*z*z-r*r)		// P_3^1 sin(p)
		   + x*(x*x-3*y*y)		// P_3^3 cos(3p)
		   + y*(y*y-3*x*x) ;	// P_3^3 sin(3p)
*/		   
	pot.annule_domain(nzm1) ; 

	Mtbl potced = 
			  (3*cost*cost - 1) / pow(r,3)			// P_2^0
		   + sint*sint*cos(2*phi) / pow(r,3)		// P_2^2 cos(2p)
		   + sint*sint*sin(2*phi) / pow(r,3)	 	// P_2^2 sin(2p)
		   + sint*(15*cost*cost-3)*cosp / pow(r,4)		// P_3^1 cos(p)
		   + sint*(15*cost*cost-3)*sinp / pow(r,4)		// P_3^1 sin(p)
		   + pow(sint,3) * cos(3*phi) / pow(r,4)		// P_3^3 cos(3p)
		   + pow(sint,3) * sin(3*phi) / pow(r,4) ;		// P_3^3 sin(3p)
		   
	//##	
	// pot.set_domain(nzm1) = potced(nzm1) ; 
	
	pot.std_spectral_base() ; 

	cout << "Potential : " << endl ; 
	pot.spectral_display() ; 
	arrete() ; 
	
	cout << "Laplacian of potential : " << endl ; 
	pot.laplacian().spectral_display() ; 
	arrete() ; 
	
	Sym_tensor_tt souc(map, map.get_bvect_cart(), metc ) ; 

	souc.set(1,1) = pot.dsdx().dsdx() ; 
	souc.set(1,2) = pot.dsdx().dsdy() ; 
	souc.set(1,3) = pot.dsdx().dsdz() ; 
	souc.set(2,2) = pot.dsdy().dsdy() ; 
	souc.set(2,3) = pot.dsdy().dsdz() ; 
	souc.set(3,3) = pot.dsdz().dsdz() ; 

	cout << "Cartesian components of the source : souc : " <<  endl ;
	souc.spectral_display() ; 

//	cout << "Maxabs divergence souc : " << endl ; 
//	maxabs( souc.divergence(metc) ) ; 
				
	cout << "Maxabs of trace of souc : " << endl ; 
	maxabs( souc.the_trace() ) ; 
		
	arrete() ; 

	Tensor tmp = souc ; 
	tmp.change_triad( map.get_bvect_spher() ) ; 
	Sym_tensor_tt sous(map, map.get_bvect_spher(), mets ) ; 
	sous = tmp ; 
	
	sous.inc_dzpuis(2) ; 

	cout << "Spherical components of the source : sous : " << endl ;
	sous.spectral_display() ; 
	
//	cout << "Maxabs divergence sous : " << endl ; 
//	maxabs( sous.divergence(mets) ) ; 
				
	cout << "Maxabs of trace of sous : " << endl ; 
	maxabs( sous.the_trace() ) ; 
		
	arrete() ; 
	
	// Resolution of Poisson equation
	// ------------------------------
	
	Sym_tensor_tt  htts = sous.poisson() ;  
	
	cout << "Solution htts (spherical components) : " << endl ;
	htts.spectral_display() ; 
	
	cout << "Maxabs divergence htts : " << endl ; 
	maxabs( htts.divergence(mets) ) ; 
		
	cout << "Maxabs of trace of htts : " << endl ; 
	maxabs( htts.the_trace() ) ; 
	
	arrete() ; 
	
	// Laplacian of the solution
	// -------------------------
	
	Tensor_sym dhtts0(map, CON, CON, COV, map.get_bvect_spher(), 0, 1)  ; 
        dhtts0 = htts.derive_cov(mets) ;
        
        Tensor dhtts = dhtts0.up(2, mets) ;
	dhtts.dec_dzpuis(2) ;  
	Tensor ddhtts = dhtts.derive_cov(mets) ;
	Sym_tensor lap_htts( ddhtts.trace(2,3) ) ; 
	
	tmp = htts ; 
	tmp.change_triad( map.get_bvect_cart() ) ;
	Sym_tensor_tt httc(map, map.get_bvect_cart(), metc ) ; 
	httc = tmp ; 
	
	cout << "Solution httc (Cartesian components) : " << endl ;
	httc.spectral_display() ; 
	
	cout << "Maxabs divergence httc : " << endl ; 
	maxabs( httc.divergence(metc) ) ; 
		
	cout << "Maxabs of trace of httc : " << endl ; 
	maxabs( httc.the_trace() ) ; 
	
	arrete() ; 
			
	Sym_tensor_tt lapc(map, map.get_bvect_cart(), metc ) ; 
	
	lapc.set(1,1) = httc(1,1).laplacian(2)  ; 
	lapc.set(1,2) = httc(1,2).laplacian(2)  ; 
	lapc.set(1,3) = httc(1,3).laplacian(2)  ; 
	lapc.set(2,2) = httc(2,2).laplacian(2)  ; 
	lapc.set(2,3) = httc(2,3).laplacian(2)  ; 
	lapc.set(3,3) = httc(3,3).laplacian(2)  ; 
	
	Sym_tensor diffc = lapc - souc ; 
	
	cout << "Maxabs of diffc : " << endl ; 
	maxabs( diffc ) ; 

	
	tmp = lapc ; 
	tmp.change_triad( map.get_bvect_spher() ) ; 
	Sym_tensor_tt laps(map, map.get_bvect_spher(), mets ) ; 
	laps = tmp ; 
	
	Sym_tensor diffs = laps - sous ; 
	cout << "diffs : " << endl ; 
	diffs.spectral_display() ; 
	cout << "Maxabs of diffs : " << endl ; 
	maxabs( diffs ) ; 

	arrete() ; 
	
	Sym_tensor difflap = lap_htts - laps ; 
	cout << "difflap : " << endl ; 
	difflap.spectral_display() ; 
	cout << "Maxabs of difflap : " << endl ; 
	maxabs( difflap ) ; 
	
		
	return EXIT_SUCCESS ; 
}
