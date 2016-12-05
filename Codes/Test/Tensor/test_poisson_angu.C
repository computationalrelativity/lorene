/*
 *  Code for testing the angular Poisson equation.
 *
 */

/*
 *   Copyright (c) 2003-2005 Eric Gourgoulhon & Jerome Novak
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
 * Revision 1.5  2014/10/06 15:12:55  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2005/05/18 15:31:25  j_novak
 * *** empty log message ***
 *
 * Revision 1.3  2005/04/04 21:36:31  e_gourgoulhon
 * Method Scalar::poisson_angu takes now a parameter
 * lambda, for the generalized angular Poisson equation
 * Lap_ang u + lambda u = source.
 * Test passed for lambda = 1 and lambda = 10.
 *
 * Revision 1.2  2003/10/19 20:05:07  e_gourgoulhon
 * Change of the argument list of Scalar::spectral_display
 * (cout now default).
 *
 * Revision 1.1  2003/10/15 21:15:25  e_gourgoulhon
 * Test for Scalar::poisson_angu().
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


using namespace Lorene ;

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	int nz = 3 ; 	// Number of domains
	int nr = 5 ; 	// Number of collocation points in r in each domain
	int nt = 5 ; 	// Number of collocation points in theta in each domain
	int np = 12 ; 	// Number of collocation points in phi in each domain
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

	uu = x + y + x*y + z*z*x - pow(z,4)*y*x ; 	
	tmp = sint * sinp / r ;
	uu.set_domain(nz-1) = tmp.domain(nz-1) ; // y/r^2 in the external domain
	
	uu.std_spectral_base() ;   // sets the standard spectral basis for 
									// expansion of a scalar field
									
	cout << "uu : " << uu << endl ; 
	uu.spectral_display() ; 
	arrete() ; 
	
        double lambda = 2. ;

	Scalar lap = uu.lapang() ;
        
        lap.set_spectral_va().ylm_i() ; 
        
        lap += lambda * uu ;  

	cout << "lap : " << endl ; 
	lap.spectral_display() ; 
	arrete() ; 
        

	Scalar uu1 = lap.poisson_angu(lambda) ; 
	
	cout << "uu1 : " << endl ; 
	uu1.spectral_display() ; 
	arrete() ; 

	Scalar diff = uu - uu1 ; 
	diff.set_spectral_va().ylm() ;
	diff.spectral_display("Diff en Ylm: ", 1.e-14) ;
	
	cout << "Norm of diff: " << norme(diff) << endl ; 	
	
	return EXIT_SUCCESS ; 
}
