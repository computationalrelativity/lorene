/*****************************************************************************

            Explicit wave-equation solver in spherical symmetry

******************************************************************************/

#include <math.h>

// Lorene headers
#include "tensor.h"
#include "utilitaires.h"
#include "graphique.h"

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

	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, false) ; // pas de domaine compactifié.. 
	
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
          phim1.set_dzpuis(0);
        double k=1. ; 
	phim1= exp (-k*r*r)-exp (- k *Rlim*Rlim);
        phim1.std_spectral_base ();
	Scalar phi = phim1 ;
         phi.std_spectral_base ();
         phi.set_dzpuis(0);
	double dt ;
	cout << "Enter dt: " << endl ;
	cin >> dt ;
	int iter = 0 ;
	int ndes = int(Rlim / dt) / 200 + 1; //Frequency for drawings

	// Definition of the "homogeneous solution  // on s'inspire ici de la solution; un dirac a la position r=Rlim, 0 ailleurs. Ca aide évidemment à fixer
// de manière simple les conditions en R, mais je ne vois pas pourquoi une telle fonction est solution de l'équation d'onde... ( surtout vis a vis de la 
// dernière équation, ajout de cette fonction homogène a la solution particulière pour obtenir la solution totale... 

 
    
	Scalar phi_h (map);
	     phi_h.std_spectral_base();
	  phi_h.annule_hard();
	   phi_h.set_dzpuis(0);
	phi_h.set_outer_boundary(nz-1, 1.);
 
        
	//-----------------------------------------



	//-----------------------------------------
	//             Main loop
	//-----------------------------------------
	for (double tps=0.; tps<4.*Rlim; tps += dt) {
	    
	    Scalar phip1(map) ;
      
 phip1.std_spectral_base ();
   
	    phip1= 2* phi-phim1 + dt*dt*phi.laplacian () ;
              phip1.set_dzpuis(0);
   

	Scalar lambda (map);  
	lambda.std_spectral_base();
        lambda= ((4*phi.val_grid_point(nz-1, 0, 0, nr-1)-phim1.val_grid_point(nz-1, 0, 0, nr-1))/(2*dt)-(3./(2.*dt)+1/Rlim)*phip1.val_grid_point(nz-1,0,0,nr-1)-phip1.dsdr().val_grid_point(nz-1, 0, 0, nr-1))/(( 3./(2.*dt)+ 1. /Rlim)*phi_h.val_grid_point(nz-1,0,0,nr-1)+ phi_h.dsdr().val_grid_point(nz-1, 0, 0, nr-1));
   
   
	phip1 += lambda*phi_h.val_grid_point(nz-1, 0, 0, nr-1); 
   
	    // Drawing
	    //--------
	    if (iter%ndes == 0){ 
		des_meridian(phip1, 0., Rlim, "\\gf", 1) ;
		;}
	    
	    // Preparation for next step
	    //--------------------------
	    iter++ ;
	    phim1 = phi ;
	    phi = phip1 ;
	}

	return EXIT_SUCCESS ; 
}



