/*****************************************************************************

            Explicit wave-equation solver in spherical symmetry

******************************************************************************/

#include <cmath>

// Lorene headers
#include "tensor.h"
#include "utilitaires.h"
#include "graphique.h"

using namespace Lorene ;

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	const int nz = 1 ; 	// Only a nucleus
	int nr; 
	cout << "Enter nr: " << endl ;
	cin >> nr ; // Number of collocation points in r
	int nt = 1 ; 	// Number of collocation points in theta 
	int np = 1 ; 	// Number of collocation points in phi
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = SYM ; // symmetry in phi

	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, false) ;
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	double Rlim = 2. ;

	// Boundaries of each domains
	double r_limits[] = {0., Rlim} ; 
  
	Map_af map(mgrid, r_limits) ; 
	const Coord& r = map.r ;

	// Setup of initial profile
	//-------------------------

	Scalar phim1(map) ;
	// Something should be done here...
	Scalar phi = phim1 ;
	double dt ;
	cout << "Enter dt: " << endl ;
	cin >> dt ;
	int iter = 0 ;
	int ndes = int(Rlim / dt) / 200 + 1; //Frequency for drawings

	// Definition of the "homogeneous solution
	//-----------------------------------------

	//-----------------------------------------
	//             Main loop
	//-----------------------------------------
	for (double tps=0.; tps<4.*Rlim; tps += dt) {
	    
	    Scalar phip1(map) ;
	    //Something should be done here...

 
	    // Drawing
	    //--------
	    if (iter%ndes == 0) 
		des_meridian(phip1, 0., Rlim, "\\gf", 1) ;
	    
	    // Preparation for next step
	    //--------------------------
	    iter++ ;
	    phim1 = phi ;
	    phi = phip1 ;
	}

	return EXIT_SUCCESS ; 
}
