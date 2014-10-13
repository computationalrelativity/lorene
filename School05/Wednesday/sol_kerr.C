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
	const Coord& xx = map.x ;
	const Coord& yy = map.y ;
	const Coord& zz = map.z ;
	const Coord& rr1 = map.r ;
	Scalar x(map) ;
	x = xx ;
	Scalar y(map) ;
	y = yy ;
	Scalar z(map) ;
	z = zz ;
	Scalar r_rho(map) ;
	r_rho = rr1 ;

	// Flat metric and both standard triads
	//-------------------------------------
	const Metric_flat& mets = map.flat_met_spher() ;
	const Base_vect_spher& bspher = map.get_bvect_spher() ;
	const Base_vect_cart& bcart = map.get_bvect_cart() ; //Cartesian base vector

	// Parameters for the Boyer-Lindquist (Kerr) BH
	//---------------------------------------------
	double mass = 0.5 ;
	double kerr_a = 0.49999 ; 

	// Define rho radial coordinate
	//-----------------------------
	Scalar rho(map) ;
	rho = 0.5*(r_rho*r_rho - kerr_a*kerr_a) + 
	    sqrt( 0.25*pow(r_rho*r_rho - kerr_a*kerr_a , 2) 
		  + kerr_a*kerr_a*z*z ) ;
	rho = sqrt(rho) ;
	rho.std_spectral_base() ;
	rho.set_domain(0) = 1 ;

	// Define Ingoing null vector l_\mu
	//---------------------------------
	Scalar hhh(map) ;
	hhh = mass*pow(rho,3) / ( pow(rho,4) + kerr_a*kerr_a*z*z ) ;    
	hhh.set_outer_boundary(nz-1, 0) ;
	hhh.std_spectral_base() ;

	Vector lmu(map, COV, bcart) ; 
	lmu.set(1) = (rho*x + kerr_a*y) / (rho*rho + kerr_a*kerr_a) ;
	lmu.set(2) = (rho*y - kerr_a*x) / (rho*rho + kerr_a*kerr_a) ;
	lmu.set(3) = z / rho ;
	for (int i=0; i<3; i++)
	    lmu.set(i+1).set_outer_boundary(nz-1, 1) ;
	lmu.std_spectral_base() ;
	lmu.change_triad(bspher) ; // That's why we need at least 4 coefficients in phi!

	// Lapse 
	//------
	Scalar lapse(map) ;
	lapse = pow( 1. + 2.*hhh , -0.5 ) ;
	lapse.std_spectral_base() ;

	// 3-metric
	//---------
	const Metric gam(mets.cov() + 2*hhh*lmu*lmu) ;
   	
	//Shift
	//-----
	Vector beta = (2*hhh*lmu).up_down(gam) ;
   
	// Extrinsic curvature
	//--------------------
	Sym_tensor k_dd = 0.5*gam.cov().derive_lie(beta) / lapse ;   
	k_dd.dec_dzpuis(2) ;
	Sym_tensor k_uu = k_dd.up_down(gam) ;
	Scalar trK = k_dd.trace(gam) ; //trace of K_ij

	//------------------------------------------------
	//     Check of Einstein equations (3+1 form)
	//------------------------------------------------

	cout << "Checking the Einstein equations in 3+1 form: " << endl ;
	cout << "(results in the nucleus are not relevant!)" << endl ;

	// Hamiltonian constraint
	//-----------------------
	Scalar ham_constr = gam.ricci_scal() ;
	ham_constr.dec_dzpuis(3) ;
	ham_constr +=  trK*trK - contract(k_uu, 0, 1, k_dd, 0, 1) ;
	maxabs(ham_constr, "Hamiltonian constraint: ") ;
 
	// Momentum constraint
	//-------------------
	Vector mom_constr = k_uu.divergence(gam) - trK.derive_con(gam) ;
	mom_constr.dec_dzpuis(2) ;
	maxabs(mom_constr, "Momentum constraint: ") ;

	// Evolution equations
	//--------------------
	Sym_tensor evol_eq = lapse*gam.ricci() 
	    - lapse.derive_cov(gam).derive_cov(gam);
	evol_eq.dec_dzpuis() ;
	evol_eq += k_dd.derive_lie(beta) ;
	evol_eq.dec_dzpuis(2) ;
	evol_eq += lapse*(trK*k_dd - 2*contract(k_dd, 1, k_dd.up(0, gam), 0) ) ;
	maxabs(evol_eq, "Evolution equations: ") ;

	return EXIT_SUCCESS ; 
}
