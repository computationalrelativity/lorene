// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

using namespace Lorene ;

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  	const int nz = 4 ; 	// Number of domains
	int nr = 33 ; 	// Number of collocation points in r in each domain
	int nt = 17 ; 	// Number of collocation points in theta in each domain
	int np = 4 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = SYM ; // symmetric with respect to (x,y) -> (-x,-y)
	bool compact = true ;
	assert(nz > 1) ;

	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;

  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double Rlim = 4. ;
	Tbl r_limits(nz+1) ;
	r_limits.set_etat_qcq() ;
	for (int i=0; i<nz; i++)
	    r_limits.set(i) = double(i)/double(nz-1) *Rlim ;
	r_limits.set(nz) = __infinity ;

	Map_af map(mgrid, r_limits) ; 

	// Definition of coordinate fields 
	//--------------------------------
	const Coord& rr1 = map.r ;

	cout << "Here there's a mapping..." << endl ;
	cout << map ;
	

	return EXIT_SUCCESS ; 
}
