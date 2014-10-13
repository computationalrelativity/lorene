// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

using namespace Lorene ;

int main() {

	// Construction of a multi-grid (Mg3d)
	// -----------------------------------
  
	const int nz = 3 ; 	// Number of domains
	int nr = 33 ; 	// Number of collocation points in r in each domain
	int nt = 9 ; 	// Number of collocation points in theta in each domain
	int np = 8 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi
	bool compact = true ;
	assert(nz > 1) ;

	Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;
	cout << "Grid: " << endl << mgrid ;

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
	cout << "Mapping: " << endl << map ;

	// Definition of coordinate fields 
	//--------------------------------
	const Coord& r = map.r ;
	const Coord& x = map.x ;
	const Coord& y = map.y ;
	const Coord& sint = map.sint ;
	const Coord& cosp = map.cosp ;
	const Coord& sinp = map.sinp ;

	// Setup of a regular scalar field
	//--------------------------------
	Scalar phi(map) ;
	phi = x*exp(-r*r-y*y) ; 
  	phi.set_outer_boundary(nz-1, 0) ;
	phi.std_spectral_base() ;

	// Drawings
	//---------
 	des_coupe_z(phi, 0, 2) ;
 	des_meridian(phi, 0, 1.3*Rlim, "\\gf", 1) ;
	arrete() ;

	// Computation of the radial derivative
	//-------------------------------------
	Scalar dphidr = phi.dsdr() ;
	dphidr.dec_dzpuis(2) ;
 	des_coupe_z(dphidr, 0, 2) ;
 	des_meridian(dphidr, 0, 1.3*Rlim, "\\gf", 2) ;
	arrete() ;

	// Comparison with analytic formula
	//---------------------------------
	Scalar dphi_anal(map) ;
	dphi_anal = exp(-r*r-y*y)*(sint*cosp*(1 - 2*r*r) 
	    - 2*sint*sinp*x*y);
	dphi_anal.set_outer_boundary(nz-1, 0.) ;
	dphi_anal.set_spectral_base(dphidr.get_spectral_base()) ;

	maxabs(dphidr - dphi_anal, "Absolute error on the derivative of phi:") ;

	// Flat metric in spherical components
	//------------------------------------
	const Metric_flat& mets = map.flat_met_spher() ;

	// Vector field defined as the gradient of a scalar field
	//-------------------------------------------------------
	Vector grad_phi = phi.derive_con(mets) ;

	// Change to cartesian triad for drawing
	grad_phi.change_triad(map.get_bvect_cart()) ;
	des_coupe_vect_z(grad_phi, 0., -2, 0.5, 1) ;

	// Vector field defined by its cartesian components
	//-------------------------------------------------
	Vector v_cart(map, CON, map.get_bvect_cart()) ;
	v_cart.set(1) = y / (1. + r*r) ;
	v_cart.set(1).set_outer_boundary(nz-1, 0) ;
	v_cart.set(2) = -x*exp(-r*r) ;
	v_cart.set(2).set_outer_boundary(nz-1, 0) ;
	v_cart.set(3) = 0 ;
	v_cart.std_spectral_base() ;
	des_coupe_vect_z(v_cart, 0., -2, 0.5, 1) ;

	// Computation of the divergence w/ flat metric in cartesian components
	//---------------------------------------------------------------------
	Scalar dive = v_cart.divergence(map.flat_met_cart()) ;
	dive.dec_dzpuis(2) ;
	des_coupe_z(dive, 0, 2) ;

	return EXIT_SUCCESS ; 
}
