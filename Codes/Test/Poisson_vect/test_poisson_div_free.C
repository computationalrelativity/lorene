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
	int nt = 5 ; 	// Number of collocation points in theta in each domain
	int np = 6 ; 	// Number of collocation points in phi in each domain
	int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
	int symmetry_phi = NONSYM ; // no symmetry in phi

	int nbr[] = {nr, nr , nr};
	int nbt[] = {nt, nt, nt} ;
	int nbp[] = {np, np, np} ;
	int tipe_r[] = {RARE, FIN, UNSURR} ;
  	assert( nz == 3 ) ;// since the above arrays are described in only 3 domains
	int nzm1 = nz - 1 ;
  
	Mg3d mgrid(nz, nbr, tipe_r, nbt, symmetry_theta, nbp, symmetry_phi) ;
	
  	// Construction of an affine mapping (Map_af)
  	// ------------------------------------------

	// Boundaries of each domains
	double r_limits[] = {0., 1., 2., __infinity} ; 
  
	Map_af map(mgrid, r_limits) ; 
  	

	// Construction of flat metrics
	// ----------------------------

	Metric_flat mets(map, map.get_bvect_spher()) ; // spherical representation
	Metric_flat metc(map, map.get_bvect_cart()) ;  // Cartesian representation


	// Construction of a divergence free vector field 
	// ----------------------------------------------

	const Coord& r = map.r ; 
	Scalar sint(map), cost(map) ;
	sint = map.sint ;
	cost = map.cost ; 

	Scalar comp_l2 = 3*cost*cost - 1 ;
	comp_l2.std_spectral_base() ;

	Scalar num(map) ;
	num = 3 - r*r*r*r ;
	num.std_spectral_base() ;

 	Base_val ma_base(nz) ;
 	ma_base.set_base_r(0, R_CHEBPIM_I) ;
 	for (int i=1; i<nz-1; i++) 
 	  ma_base.set_base_r(i, R_CHEB) ;
 	ma_base.set_base_r(nz-1, R_CHEBU) ;
 	ma_base.set_base_t(T_COSSIN_CP) ;
 	ma_base.set_base_p(P_COSSIN) ;

	Mtbl denom0(mgrid) ;
	denom0 = 1/(r*r*r*r) ;
	Scalar unsur4(map) ;
	unsur4 = 1 / (1 + r*r*r*r) ;
	unsur4.annule_domain(nz-1) ;
	unsur4.std_spectral_base() ;

	Scalar la_zec(map) ;
	la_zec = 1 ;
	la_zec.annule(0,nz-2) ;
	la_zec.set_dzpuis(4) ;
	Scalar tmp(map) ;
	tmp = 1 / (1 + denom0) ;
	tmp.annule(0, nz-2) ;
	la_zec *= tmp ;
	unsur4 += la_zec ;
	
	Scalar rtrois(map), denom1(map), denom2(map), denom3(map) ;
	rtrois = r*r*r ;
	rtrois.set_spectral_base(ma_base) ;
	denom1 = 1 / (r*r*r*r + 1) ;
	denom1.std_spectral_base() ;
	denom2 = unsur4 * denom1 ;
	denom3 = denom2*denom1 ;
	Scalar denom4 = denom3*denom1 ;

	Vector vvs1(map, CON, map.get_bvect_spher()) ;
	vvs1.set(1) =  denom1 *comp_l2 ;
 	vvs1.set(1).set_outer_boundary(nzm1, 0) ;
	vvs1.set(1).mult_r() ;
 	vvs1.set(1).set_outer_boundary(nzm1, 0) ;
	vvs1.set(2) = - num*denom1*denom1 ;
	vvs1.set(2).set_outer_boundary(nzm1, 0) ;
	vvs1.set(2).mult_rsint() ;
	vvs1.set(2).mult_cost() ;
	vvs1.set(3) = -denom1 ;
	vvs1.set(3).set_outer_boundary(nzm1, 0) ;
	vvs1.set(3).mult_rsint() ;
	vvs1.set(3).set_outer_boundary(nzm1, 0) ;
	vvs1.std_spectral_base() ;

	Vector source_s(map, CON, map.get_bvect_spher()) ;
	source_s.set(1) = -(32 * denom3 + 4 * denom2)*rtrois*comp_l2 ; 
	source_s.set(1).set_outer_boundary(nzm1, 0) ;
	source_s.set(2) = ( 384 * denom4 - 192 * denom3 - 12 * denom2)*rtrois ; 
	source_s.set(2).set_outer_boundary(nzm1, 0) ;
	source_s.set(2).set_spectral_base(ma_base) ;
	source_s.set(2).mult_sint() ;
	source_s.set(2).mult_cost() ;
	source_s.set(3) = (32*denom3 - 4 * denom2) * rtrois ;
	source_s.set(3).set_outer_boundary(nzm1, 0) ;
	source_s.std_spectral_base() ;
	source_s.set(3).set_spectral_base(ma_base) ;
	source_s.set(3).mult_sint() ;

	cout << "Divergence of source : " << endl ;
 	cout << "---------------------- " << endl << endl ;
	source_s.divergence(mets).spectral_display(0x0, 1.e-12) ;
 	int sd ; cin >> sd ;

	Scalar pot_source = source_s.potential(mets) ;
 	cout << "Source potential: " << endl ;
 	cout << "----------------- " << endl << endl ;
 	pot_source.spectral_display(0x0, 1.e-13) ;
 	cin >> sd ;

	Vector_divfree source_df = source_s.div_free(mets) ;

 	cout << "Difference with the div_free part of the source : " << endl;
 	cout << "------------------------------------------------- " << endl << endl ;
 	Vector v_diff = source_df - source_s ;
	v_diff.dec_dzpuis(4) ;
 	v_diff.spectral_display(0x0, 1.e-13) ;
 	cout << endl ;
 	cin >> sd ;

	Vector_divfree sol_df = source_df.poisson() ;
	Vector diff_df = sol_df - vvs1 ;
	cout << "Difference with analytic solution: " << endl ;
	cout << "---------------------------------- " << endl ;
	for (int i=1; i<4; i++) 
	  cout << max(abs(diff_df(i))) << endl ;
	cout << endl ;

	Vector sol = sol_df ;
	cout << "Divergence of the solution : " << endl ;
	cout << "---------------------------- " << endl ;
	Scalar dive = sol.divergence(mets) ;
	dive.dec_dzpuis(2) ;
	dive.spectral_display(0x0, 1.e-12) ;

	Scalar eta_theo = num*denom1*denom1*comp_l2 / 6 ;
	eta_theo.std_spectral_base() ;
	eta_theo.set_outer_boundary(nzm1,0) ;
	eta_theo.mult_r() ;

	Scalar eta = sol_df.eta() ;
	cout << "Difference in eta : " << endl ;
	cout << "------------------- " << endl ;
	cout << max(abs(eta- eta_theo)) << endl ;

	Scalar mu_theo(map) ;
	mu_theo = denom1 ;
	mu_theo.set_outer_boundary(nzm1,0) ;
	mu_theo.mult_r() ;
	mu_theo.mult_cost() ;

	Scalar mu = sol_df.mu() ;
	cout << "Diffrence in mu: " << endl  ;
	cout << "------------------- " << endl ;
	cout << max(abs(mu - mu_theo)) << endl ;

	return EXIT_SUCCESS ; 
}
