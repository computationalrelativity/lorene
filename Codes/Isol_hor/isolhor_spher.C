
/*
 *  Main code for Isolated Horizon in spherical symmetry and flat
 *  conformal metric: to test kappa = const condition
 *
 */

/*
 *   Copyright (c) 2004-2005  Jose Luis Jaramillo
 *                            Francois limousin
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

 



// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cmath>
#include <cstring>

// Lorene headers
#include "tenseur.h"
#include "metric.h"
#include "evolution.h"
#include "param.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"
#include "time_slice.h"
#include "isol_hor.h"


using namespace Lorene ;

int main() {

    //======================================================================
    //      Construction and initialization of the various objects
    //======================================================================

    // Setup of a multi-domain grid (Lorene class Mg3d)
    // ------------------------------------------------

    int nz, nt, np, nr1, nrp1 ;

    ifstream fpar("par_hor.d") ;
    fpar.ignore(1000, '\n') ;
    fpar.ignore(1000, '\n') ;
    fpar >> nz; fpar.ignore(1000, '\n');
    fpar >> nt; fpar.ignore(1000, '\n');
    fpar >> np; fpar.ignore(1000, '\n');
    fpar >> nr1; fpar.ignore(1000, '\n');
    fpar >> nrp1; fpar.ignore(1000, '\n');


    int type_t = SYM ; // symmetry with respect to the equatorial plane
    int type_p = NONSYM ; // no symmetry in phi
 
    // Reading of the parameter file
    // ------------------------------

    int niter, bound_nn, bound_psi, bound_beta, solve_lapse, solve_psi ;
    int solve_shift ;
    double radius, relax, seuil, ang_vel, boost_x, boost_z, lim_nn ;
    fpar >> radius; fpar.ignore(1000, '\n');
    fpar >> relax; fpar.ignore(1000, '\n');
    fpar >> seuil; fpar.ignore(1000, '\n');
    fpar >> niter; fpar.ignore(1000, '\n');
    fpar >> ang_vel; fpar.ignore(1000, '\n');
    fpar >> boost_x ;
    fpar >> boost_z; fpar.ignore(1000, '\n');
    fpar >> bound_nn ;
    fpar >> lim_nn ;  fpar.ignore(1000, '\n');
    fpar >> bound_psi ;  fpar.ignore(1000, '\n');
    fpar >> bound_beta ;  fpar.ignore(1000, '\n');
    fpar >> solve_lapse ;  fpar.ignore(1000, '\n');
    fpar >> solve_psi ;  fpar.ignore(1000, '\n');
    fpar >> solve_shift ;  fpar.ignore(1000, '\n');
    

    int* nr_tab = new int[nz];
    int* nt_tab = new int[nz];
    int* np_tab = new int[nz];
    double* bornes = new double[nz+1];
    
    for (int l=0; l<nz; l++) {
      if (l==1) nr_tab[1] = nr1 ;
      else nr_tab[l] = nrp1 ;
      np_tab[l] = np ; 
      nt_tab[l] = nt ; 
      bornes[l] = pow(2., l-1) * radius ;
    }
    bornes[0] = 0. ;
    bornes[nz] = __infinity ; 

    // Type of r sampling : 
    int* type_r = new int[nz];
    type_r[0] = RARE ; 
    for (int l=1; l<nz-1; l++) {
      type_r[l] = FIN ; 
    }
    type_r[nz-1] = UNSURR ; 

    // Multi-domain grid construction:
    Mg3d mgrid(nz, nr_tab, type_r, nt_tab, type_t, np_tab, type_p) ;

    Map_af mp(mgrid, bornes) ;   // Mapping construction
  	
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = mp.r ;        // r field 
    const Coord& theta = mp.tet ;        // r field 
    const Coord& costt = mp.cost ;  // cos(theta) field
    const Coord& sintt = mp.sint ;  // sin(theta) field
    const Coord& cospp = mp.cosp ;  // cos(phi) field
    const Coord& sinpp = mp.sinp ;  // sin(phi) field

    Scalar rr(mp) ;
    rr = r ;
    rr.std_spectral_base() ;
 
    Scalar cos2t (mp) ;
    cos2t = cos(2.*theta) ;
    
    Scalar cost (mp) ;
    cost = costt ;
    
    Scalar cosp (mp) ;
    cosp = cospp ;
  
    Scalar sint (mp) ;
    sint = sintt ;
 
    Scalar sinp (mp) ;
    sinp = sinpp ;
 
    // Flat metric f
    // -------------

    const Metric_flat& ff = mp.flat_met_spher() ; 
    const Base_vect_spher& otriad = mp.get_bvect_spher() ;    

    // Working stuff
    // -------------
    
    Scalar tmp_scal(mp) ;
    Vector tmp_vect(mp, CON, otriad) ;
    Sym_tensor tmp_sym(mp, CON, otriad) ; 

    // Key function
    Scalar unsr(1./rr) ;    
    unsr.set_domain(0) = 1 ; // scalar set to 1 in the nucleus 
    unsr.std_spectral_base() ;

    Scalar expmr (exp (- 1/((rr-1)*(rr-1)))) ;
    expmr.std_spectral_base() ;

    Scalar expmrr (exp(-pow(rr-3,2.))) ;
    expmrr.std_spectral_base() ;
 

    // =================================================================
    // =================================================================

    // Physical Parameters
    //--------------------
    
    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(mp) ; 
    nn_init = 1. - 0.5*unsr ;
    nn_init.std_spectral_base() ;    // sets standard spectral bases

    // Set up of field Psi 
    // -------------------

    Scalar psi_init(mp) ; 
    psi_init =  1. ; //+ 0.1*unsr ;
    psi_init.std_spectral_base() ;    // sets standard spectral bases

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(mp, CON, mp.get_bvect_spher()) ; 
    beta_init.set(1) = 0.1*unsr*unsr ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;
    beta_init.annule_domain(0) ;
    beta_init.std_spectral_base() ;

    // TrK, TrK_point
    // --------------

    Scalar trK (mp) ;
    trK = 0*0.1*unsr*unsr ;
    trK.std_spectral_base() ;
    trK.inc_dzpuis(2) ;

    Scalar trK_point (mp) ;
    trK_point = 0. ;
    trK_point.std_spectral_base() ;
    trK_point.inc_dzpuis(2) ;
	
    // gamt, gamt_point
    // ----------------

    
    Sym_tensor gamt(mp, COV, mp.get_bvect_spher()) ;

    gamt = ff.cov() ;

    Metric met_gamt (gamt) ; 

    // Gamma-tilde_point
    //------------------
    
    Scalar khi(mp) ;
    Scalar mu(mp) ;

    khi = 0. ;
    khi.std_spectral_base() ;
    
    mu = 0. ;
    mu.std_spectral_base() ;

    Sym_tensor_tt hh_tmp (mp, otriad, ff) ;
    hh_tmp.set_khi_mu(khi, mu) ;

    Sym_tensor gamt_point(mp, CON, mp.get_bvect_spher()) ;
    gamt_point = hh_tmp ;
    gamt_point.inc_dzpuis(2) ;


    // =================================================================
    // =================================================================
    
    double mm, aaa, hh ;
    
    // Set up of extrinsic curvature
    // -----------------------------
    
    Metric met_gam(psi_init*psi_init*psi_init*psi_init*gamt) ;
    Sym_tensor kk_init (mp, CON, mp.get_bvect_spher()) ;

    int check ;
    check = 0 ;
    for (int k=0; k<np_tab[1]; k++)
	for (int j=0; j<nt_tab[1]; j++){
	    if (nn_init.val_grid_point(1, k, j , 0) < 1e-12){
		check = 1 ;
		break ;
	    }
	}
    if (check == 0)
	kk_init = - 0.5 * met_gam.con().derive_lie(beta_init) / nn_init ;
    else {
	Sym_tensor kk_temp (mp,  CON, mp.get_bvect_spher()) ;
	kk_temp = - 0.5 * met_gam.con().derive_lie(beta_init) ;

	Scalar nn_sxpun (division_xpun (Cmp(nn_init), 0)) ;
	nn_sxpun.set_domain(0) = 1. ;

	Scalar auxi (mp) ;
	for (int i=1 ; i<=3 ; i++)
	    for (int j=i ; j<=3 ; j++) {
		auxi = kk_temp(i, j) ;
		auxi.annule_domain(0) ;
		auxi = division_xpun (auxi, 0) ;
		kk_init.set(i,j) = auxi / nn_sxpun ;
	    }
    }

    Sym_tensor aa_init (mp, CON, mp.get_bvect_spher()) ;
    aa_init = psi_init*psi_init*psi_init*psi_init*kk_init 
	- 1./3. * trK * met_gamt.con() ;  
    
    
    //-------------------------------------
    //     Construction of the space-time
    //-------------------------------------

    Isol_hor isolhor(mp, nn_init, psi_init, beta_init, aa_init, met_gamt,
		     gamt_point, trK, trK_point, ff, 3) ;

   // In order to initialise isolhor.k_uu() and k_dd() at the good value
    Sym_tensor bidon (mp, CON, mp.get_bvect_spher()) ;
    bidon = isolhor.k_uu() ;
    Sym_tensor bidon2 (isolhor.k_dd()) ;
    Scalar bidon3 (isolhor.aa_quad()) ;
       
    
    //-----------------------------------------
    //          "Call to init_data.C" 
    //-----------------------------------------

    isolhor.set_omega(ang_vel) ;

    isolhor.init_data_spher(bound_nn, lim_nn, bound_psi, bound_beta, solve_lapse,
		      solve_psi, solve_shift, seuil, relax, niter) ;


    // Expansion
    //----------

    Scalar expansion = contract(isolhor.gam().radial_vect().derive_cov(isolhor.gam()), 0,1) 
      + contract(contract(isolhor.k_dd(), 0, isolhor.gam().radial_vect(), 0), 
	       0, isolhor.gam().radial_vect(), 0) 
      - isolhor.trk() ; 

    expansion.dec_dzpuis(2) ;


    double der_expansion_0 = expansion.derive_cov(isolhor.gam())(1)
      .val_point(1.00000001, 0, 0.) ;
    double der_expansion_1 = expansion.derive_cov(isolhor.gam())(1)
      .val_point(1.00000001, M_PI/4, 0.) ;
    double der_expansion_2 = expansion.derive_cov(isolhor.gam())(1)
      .val_point(1.00000001, M_PI/2, 0.) ;

    cout << "Radial derivative of the expansion at (1,0,0) = "
	 << der_expansion_0<<endl ;
    cout << "Radial derivative of the expansion at (1,Pi/4,0) = "
	 << der_expansion_1<<endl ;
    cout << "Radial derivative of the expansion at (1,Pi/2,0) = "
	 << der_expansion_2<<endl ;
    cout << "------------------------------------------------------"<<endl;
 
        
    //    des_meridian(expansion, 1, 1.2, "Expansion theta", 2) ;
    des_meridian(expansion, 1, 4., "Expansion theta", 1) ;
    //    des_meridian(expansion, 1, 10., "Expansion theta", 3) ;
    arrete() ;
    

    double max_exp = expansion.val_grid_point(1, 0, 0, 0) ;
    double min_exp = expansion.val_grid_point(1, 0, 0, 0) ;
    int nnp = mp.get_mg()->get_np(1) ;
    int nnt = mp.get_mg()->get_nt(1) ;
    for (int k=0 ; k<nnp ; k++)
      for (int j=0 ; j<nnt ; j++){
	if (expansion.val_grid_point(1, k, j, 0) > max_exp)
	  max_exp = expansion.val_grid_point(1, k, j, 0) ;
	if (expansion.val_grid_point(1, k, j, 0) < min_exp)
	  min_exp = expansion.val_grid_point(1, k, j, 0) ;
      }
    cout << "max_exp = " << max_exp << endl 
	 << "min_kss = " << min_exp << endl ;


    // Save in a file
    // --------------
    
    FILE* fresu = fopen("resu.d", "w") ;
    mgrid.sauve(fresu) ;
    mp.sauve(fresu) ;
    isolhor.sauve(fresu, true) ;
    fclose(fresu) ;     
    
    // Test of the constraints
    //------------------------

    /*    
    if (solve_shift == 1)
      isolhor.update_aa() ;
    else
	isolhor.aa_kerr_ww(mm, aaa) ;
    */
    cout<< "----------------------------------------" <<endl ;
    
    Tbl check_ham = isolhor.check_hamiltonian_constraint() ;
    Tbl check_mom = isolhor.check_momentum_constraint() ;

    cout<< "----------------------------------------" <<endl ;
 

    // Physical parameters of the Black Hole
    //--------------------------------------
    
    cout<< "------------------------------------------------" <<endl;
    cout<< "      Physical parameters of the Black Hole     " <<endl;
    cout<< "------------------------------------------------" <<endl;
    
    double rr_hor =  isolhor.radius_hor() ;
    cout<< "Radius of the horizon = " << rr_hor <<endl ;
    
    double jj_hor =  isolhor.ang_mom_hor() ;
    cout<< "Angular momentum of the horizon = " << jj_hor <<endl ; 

    double mm_hor = isolhor.mass_hor() ;
    cout<< "Mass of the horizon = " << mm_hor <<endl ;  

    double kappa_hor = isolhor.kappa_hor() ;
    cout<< "Surface gravity of the horizon = " << kappa_hor <<endl ; 

    double omega_hor = isolhor.omega_hor() ;
    cout<< "Orbital velocity of the horizon = " << omega_hor <<endl ; 


    // Physical parameters of the Bulk
    //--------------------------------

    cout.precision(8) ;
    cout<< endl;
    cout<< "------------------------------------------------" <<endl;
    cout<< "      Physical parameters of the Bulk           " <<endl;
    cout<< "------------------------------------------------" <<endl;
    
    double mm_adm = isolhor.adm_mass() ;
    cout << "ADM mass= " << mm_adm <<endl ;  

    double jj_adm = isolhor.ang_mom_adm() ;
    cout << "ADM angular momentum= " << jj_adm <<endl ;  

    double aa = jj_adm / mm_adm ;
    cout << "J / M (ADM) : " << aa << endl ;  

    double aasmm = aa / mm_adm ;
    cout << "J / M^2 : " << aasmm << endl ; 

    double epsa = isolhor.area_hor() / 
	(8*M_PI*(mm_adm*mm_adm + pow(pow(mm_adm, 4.) - jj_adm*jj_adm, 0.5))) ;
    cout << "epsilon A : " << epsa << endl ;
  
    double diff_mm = (mm_adm - mm_hor) / mm_adm ;
    cout << "diff mass : " << diff_mm << endl ;  

    double diff_jj = (jj_adm - jj_hor) / jj_adm ;
    cout << "diffangular momentum : " << diff_jj << endl ;  

    // Writing of all global quantities in a file
    ofstream resformat("resformat.d") ;
    resformat.precision(6) ;
    resformat << "# Ham. constr.    Mom. constr. " << endl ;
    resformat << max(check_ham) << "    " << max(check_mom) << endl ;
    resformat.precision(10) ;
    resformat << "# r_hor   J_hor   M_hor  kappa_hor  omega_hor" << endl ;
    resformat <<  rr_hor << "  " << jj_hor << "  " << mm_hor << "  " 
	      << kappa_hor << "  " << omega_hor << endl ;
    resformat << "# M_adm  J_adm  J/M   J/M2   Eps_a" << endl ;
    resformat << mm_adm << "  " << jj_adm << "  " << aa << "  " 
	      << aasmm << "  " << epsa << endl ;
    resformat << "# diff_mm    diff_jj" << endl ;
    resformat << diff_mm << "  " << diff_jj << endl ;

    resformat.close() ;

    //--------------------------------------
    //        Comparison
    //--------------------------------------

    cout<<"Tout va bien boudiou / Todo bien!!! (Viva Cai!)"<<endl ;

    return EXIT_SUCCESS ; 
}


    
    
