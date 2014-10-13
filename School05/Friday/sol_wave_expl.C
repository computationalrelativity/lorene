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

	double sigma = 9;
	Scalar phim1(map) ;
  	phim1 = exp(-sigma*r*r) - exp(-sigma*Rlim*Rlim) ;
	phim1.std_spectral_base() ;
	Scalar phi = phim1 ;
	double dt ;
	cout << "Enter dt: " << endl ;
	cin >> dt ;
	int iter = 0 ;
	int ndes = int(Rlim / dt) / 200 + 1; //Frequency for drawings

	// Definition of the "homogeneous solution"
	//-----------------------------------------

	Scalar phi_hom(map) ;
	phi_hom.annule_hard() ;
	phi_hom.set_outer_boundary(nz-1, 1.) ;
	phi_hom.std_spectral_base() ;
	double d_hom = phi_hom.dsdr().val_grid_point(nz-1, 0, 0, nr-1) ;
	double a_bound = 3./(2.*dt) + 1./Rlim ; //For the boundary conditions
	double b_bound = 1. ;                   // here Sommerfeld outgoing
	double c_bound = 0. ;                   // condition

	//-----------------------------------------
	//             Main loop
	//-----------------------------------------
	for (double tps=0.; tps<4.*Rlim; tps += dt) {
	    
	    // Simple forward Euler explicit scheme
	    //-------------------------------------
	    Scalar phip1 = 2*phi - phim1 + dt*dt*phi.laplacian() ;
	    phip1.set_dzpuis(0) ;
	    
	    // Imposing boundary conditions
	    //-----------------------------
 	    c_bound = (4*phi.val_grid_point(nz-1, 0, 0, nr-1) 
 		       - phim1.val_grid_point(nz-1, 0, 0, nr-1) ) /(2*dt) ;
 	    double lambda = (c_bound - a_bound*phip1.val_grid_point(nz-1, 0, 0, nr-1) 
 		- b_bound*phip1.dsdr().val_grid_point(nz-1, 0, 0, nr-1)) /
 		( a_bound + b_bound*d_hom) ;
 	    phip1 += lambda*phi_hom ;

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
