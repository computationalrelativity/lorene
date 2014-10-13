#include <math.h>

// Lorene headers
#include "metric.h"
#include "proto.h"
#include "graphique.h"
#include "diff.h"

double next_xi(double, double, double, double, double, int) ;

using namespace Lorene ;

int main() {

    // Construction of a multi-grid (Mg3d)
    // -----------------------------------
    
    const int nz = 1 ; 	// Only a nucleus
    int nr; 
    cout << "Enter nr: " << endl ;
    cin >> nr ; // Number of collocation points in r
    int nt = 5 ; 	// Number of collocation points in theta 
    int np = 6 ; 	// Number of collocation points in phi
    int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
    int symmetry_phi = SYM ; // symmetry in phi
    
    Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, false) ;
    
    // Construction of an affine mapping (Map_af)
    // ------------------------------------------
    
    double Rlim = 2. ;
    
    // Boundaries of each domains
    double r_limits[] = {0., Rlim} ; 
    
    Map_af map(mgrid, r_limits) ; 
    double alpha = map.get_alpha()[0] ;

    // Coordinate fields
    const Coord& r = map.r ;
    const Coord& x = map.x ;
    const Coord& y = map.y ;

    // Flat metric in spherical components
    //------------------------------------
    const Metric_flat& mets = map.flat_met_spher() ;
    
    // Initial data
    //-------------
    Scalar phim1(map) ;
    Scalar phi = phim1 ;

    // Time-step & boundary conditions
    //--------------------------------
    double dt ;
    cout << "Enter dt: " << endl ;
    cin >> dt ;
    int iter = 0 ;
    int ndes = int(Rlim / dt) / 200 + 1;

    // Some initializations
    //---------------------
    int l_q, m_q, base_r ;
    double emax = -1.e33 ;
    double emin = 1.e33 ;
    phi.set_spectral_va().ylm() ;
    const Base_val& base = phi.get_spectral_base() ;
    Mtbl_cf xi(mgrid.get_angu(), base) ;
    xi.annule_hard() ;
    Mtbl_cf xim1(mgrid.get_angu(), base) ;
    xim1.annule_hard() ;

    //----------------------------
    //         Main loop
    //----------------------------
    for (double tps=0.; tps<3.*Rlim; tps += dt) {
	
	// Crank-Nicholson 2nd-order scheme
	//---------------------------------
	Scalar source(map) ;
	phi.set_spectral_va().ylm() ;
	phim1.set_spectral_va().ylm() ;

	Mtbl_cf resu(mgrid, base) ;
	resu.annule_hard() ; //Initialized to zero
	    
    // Loop on theta, phi
    //-------------------
	for (int k=0; k<np; k++)
	    for (int j=0; j<nt; j++) 
		if (nullite_plm(j, nt, k, np, base) == 1) {
		    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
		    
		    // Recovering the elementary operators
		    //------------------------------------
		    Diff_sxdsdx dx(base_r, nr) ; const Matrice& mdx = dx.get_matrice() ;
		    Diff_dsdx2 dx2(base_r, nr) ; const Matrice& md2 = dx2.get_matrice() ;
		    Diff_sx2 sx2(base_r, nr) ; const Matrice& ms2 = sx2.get_matrice() ;
		    Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;
    
		    // Global operator (without BCs)
		    //------------------------------
		    Matrice ope(nr,nr) ;

		    // Boundary conditions on the last line
		    //-------------------------------------

		    Tbl rhs(nr) ; rhs.set_etat_qcq() ;
		    Tbl fi(nr) ; fi.set_etat_qcq() ;
		    Tbl fim1(nr) ; fim1.set_etat_qcq() ;
		    
		    for (int i=0; i<nr; i++) {
			rhs.set(i) = (*source.get_spectral_va().c_cf)(0, k, j, i) ;
			fi.set(i) = (*phi.get_spectral_va().c_cf)(nz-1, k, j, i) ;
			fim1.set(i) = (*phim1.get_spectral_va().c_cf)(nz-1, k, j, i) ;
		    }
		    double fi_1 = val1_dern_1d(0, fi, base_r) ;
		    double fim1_1 = val1_dern_1d(0, fim1, base_r) ;

		    // Enhanced outgoing BC (needs auxilliary equation to be solved)
		    //--------------------------------------------------------------
		    // To be done ....

		    // regularity conditions

		    // Solution of the system
		    //-----------------------
		    ope.set_lu() ;
		    Tbl sol = ope.inverse(rhs) ;
		    for (int i=0; i<nr; i++) 
			resu.set(0, k, j, i) = sol(i) ;
		}

	// Putting everything into place
	//------------------------------
	Scalar phip1(map) ;
	phip1.set_etat_qcq() ;
	phip1.set_spectral_va() = resu ;
	phip1.set_spectral_va().ylm_i() ;
	phi.set_spectral_va().ylm_i() ;
	phim1.set_spectral_va().ylm_i() ;

	// Energy calculation
	//-------------------
	Scalar dtphi = (phip1 - phim1)/(2*dt) ;
	double ener = (dtphi*dtphi 
		       + contract(phi.derive_con(mets), 0, 
				  phi.derive_cov(mets), 0)).integrale() ;
	
	emax = (ener > emax) ? ener : emax ;
	emin = (ener < emin) ? ener : emin ;
	
	cout << "Energy inside the grid: " << ener << endl ;
   
	// Drawing
	//--------
	if (iter%ndes == 0) {
	    des_meridian(phip1, 0., Rlim, "\\gf", 1) ;
	}
	    
	// Preparation for next step
	//--------------------------
	iter++ ;
	phim1 = phi ;
	phi = phip1 ;
    }
    cout << "Emin, Emax: " << emin << ", " << emax << endl ;
    cout << "Difference: " << (emax - emin) / emax << endl ;

    return EXIT_SUCCESS ; 
}

//----------------------------------------------------------------
// Auxilliary equation: wave-like equation on the boundary sphere
// 2nd-order Crank-Nicholson time-integration scheme
//----------------------------------------------------------------

// This might be usefull...
double next_xi(double xi, double xim1, double source, double Rlim, double dt,
	       int l_q) {

    double xip1 ; //\xi^{J+1}
    return xip1 ;

}
