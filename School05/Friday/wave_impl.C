/*****************************************************************************

    Implicit wave-equation solver in #D with Sommerfeld boundary conditions

******************************************************************************/

#include <cmath>

// Lorene headers
#include "metric.h"
#include "proto.h"
#include "graphique.h"
#include "diff.h"

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
    // Something should be done here ...
    Scalar phi = phim1 ;

    // Time-step & boundary conditions
    //--------------------------------
    double dt ;
    cout << "Enter dt: " << endl ;
    cin >> dt ;
    int iter = 0 ;
    int ndes = int(Rlim / dt) / 200 + 1;
    // Something should be done here ...

    // Some initializations
    //---------------------
    int l_q, m_q, base_r ;
    double emax = -1.e33 ;
    double emin = 1.e33 ;
    phi.set_spectral_va().ylm() ;
    const Base_val& base = phi.get_spectral_base() ;

    //----------------------------
    //         Main loop
    //----------------------------
    for (double tps=0.; tps<3.*Rlim; tps += dt) {

	Scalar source(map) ;
	// Something should be done here ...
	
	Mtbl_cf resu(mgrid, base) ;
	resu.annule_hard() ; //Initialized to zero
	    
    // Loop on theta, phi
    //-------------------
	for (int k=0; k<np; k++)
	    for (int j=0; j<nt; j++) 
		if (nullite_plm(j, nt, k, np, base) == 1) {
		    base.give_quant_numbers(0, k, j, m_q, l_q, base_r) ;
		    
		    // Recovering the elementary operators (Diff)
		    //------------------------------------
		    Diff_sxdsdx dx(base_r, nr) ; const Matrice& mdx = dx.get_matrice() ;
		    Diff_dsdx2 dx2(base_r, nr) ; const Matrice& md2 = dx2.get_matrice() ;
		    Diff_sx2 sx2(base_r, nr) ; const Matrice& ms2 = sx2.get_matrice() ;
		    Diff_id id(base_r, nr) ; const Matrice& mid = id.get_matrice() ;

		    Matrice ope(nr, nr) ;
	// Something should be done here ... to fill the matrix + BCs + regularity

		    Tbl rhs(nr) ; rhs.set_etat_qcq() ;
		    for (int i=0; i<nr; i++) {
			rhs.set(i) = (*source.get_spectral_va().c_cf)(0, k, j, i) ;
		    }
	// Something should be done here ... BCs still!


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
