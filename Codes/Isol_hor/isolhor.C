/*
 *  Main code for Isolated Horizon in arbitrary gauge
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

char isolhor_C[] = "$Header$" ;

/* 
 * $Id$
 * $Log$
 * Revision 1.24  2005/04/02 15:52:05  f_limousin
 * New data member nz. Lichnerowicz choice for aquad. Delete function
 * compute_ww().
 *
 * Revision 1.23  2005/03/31 09:46:33  f_limousin
 * New functions compute_ww(..) and aa_kerr_ww().
 *
 * Revision 1.22  2005/03/30 12:08:33  f_limousin
 * Implementation of K^{ij} (Eq.(13) Of Sergio (2002)).
 *
 * Revision 1.21  2005/03/28 19:42:24  f_limousin
 * Implement the metric and A^{ij}A_{ij} of Sergio for pertubations
 * of Kerr black holes.
 *
 * Revision 1.20  2005/03/24 16:50:53  f_limousin
 * Add parameters solve_shift and solve_psi in par_isol.d and in function
 * init_dat(...). Implement Isolhor::kerr_perturb().
 *
 * Revision 1.19  2005/03/22 13:25:49  f_limousin
 * Small changes. The angular velocity and A^{ij} are computed
 * with a differnet sign.
 *
 * Revision 1.18  2005/03/15 19:28:11  f_limousin
 * Implement initial data for Kerr black hole in isotropic coordinates.
 *
 * Revision 1.17  2005/03/09 10:33:31  f_limousin
 * Delete functions init_data_b_neumann(...) and init_data_berlin(...)
 * --> New parameter solve_lapse in the function init_data(...).
 *
 * Revision 1.16  2005/03/03 10:18:57  f_limousin
 * Addition of the boost in x and z-direction.
 * The grid and the mapping are now saved in the output file.
 *
 * Revision 1.15  2004/12/31 12:30:07  f_limousin
 * Change the construction of an Isol_hor and the function sauve(FILE*, bool).
 *
 * Revision 1.14  2004/11/18 10:02:37  jl_jaramillo
 * gamt and gamt_point well constructed as tensor in spaherical
 * components
 *
 * Revision 1.13  2004/11/09 12:40:08  f_limousin
 * Add some printing
 *
 * Revision 1.12  2004/11/05 17:53:26  f_limousin
 * Perturbations of the conformal metric and its time derivative.
 *
 * Revision 1.8  2004/11/02 17:42:00  f_limousin
 * New method sauve(...) to save in a binary file.
 *
 * Revision 1.6  2004/10/29 15:41:02  jl_jaramillo
 * ADM angular momentum added
 *
 * Revision 1.5  2004/10/01 16:48:47  f_limousin
 * *** empty log message ***
 *
 * Revision 1.4  2004/09/28 15:57:45  f_limousin
 * Add the 2 lines  $Id and $Log to see the comments
 *
 *
 * Revision 1.1  2004/02/18 19:16:28  jl_jaramillo
 * First version: c'est loin d'etre pret tout ca !!!
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
    cos2t.std_spectral_base() ;
    
    Scalar cost (mp) ;
    cost = costt ;
    cost.std_spectral_base() ;
    
    Scalar cosp (mp) ;
    cosp = cospp ;
    cosp.std_spectral_base() ;
  
    Scalar sint (mp) ;
    sint = sintt ;
    sint.std_spectral_base() ;
 
    Scalar sinp (mp) ;
    sinp = sinpp ;
    sinp.std_spectral_base() ;
 
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
    psi_init =  1. + unsr ;
    psi_init.std_spectral_base() ;    // sets standard spectral bases

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(mp, CON, mp.get_bvect_spher()) ; 
    beta_init.set(1) = 0.001*unsr*unsr ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;
    beta_init.annule_domain(0) ;
    beta_init.std_spectral_base() ;

    // TrK, TrK_point
    // --------------

    Scalar trK (mp) ;
    trK = 0. ;//0.01*unsr*unsr ;
    trK.std_spectral_base() ;
    trK.inc_dzpuis(2) ;

    Scalar trK_point (mp) ;
    trK_point = 0. ;
    trK_point.std_spectral_base() ;
    trK_point.inc_dzpuis(2) ;
	
    // gamt, gamt_point
    // ----------------

    Scalar khi (mp) ;
    khi = 0*0.05*unsr*unsr*sint*sint ;
    khi.std_spectral_base() ;
    khi.annule_domain(0) ;
    
    Scalar mu (mp) ;
    mu = 0*0.05*unsr*unsr ;
    mu.std_spectral_base() ;
    mu.annule_domain(0) ;
    
    Sym_tensor_tt hh_tmp (mp, otriad, ff) ;
    hh_tmp.set_khi_mu(khi, mu) ;

    
    //Construction of a gamt
    //----------------------

    Sym_tensor gamt(mp, COV, mp.get_bvect_spher()) ;

    gamt = ff.cov() + hh_tmp.up_down(ff) ;

    // Determinant of gamma tilde is put to one 
    // ----------------------------------------

    Metric met_gamt_tmp (gamt) ;             
    Scalar det_ust = pow(met_gamt_tmp.determinant(), -1./3.) ;
    det_ust.std_spectral_base() ;
     
    gamt = gamt*det_ust ;
    Metric met_gamt (gamt) ; 

    // Gamma-tilde_point
    //------------------

    khi = 0. ;
    khi.std_spectral_base() ;
    
    mu = 0. ;
    mu.std_spectral_base() ;
    
    hh_tmp.set_khi_mu(khi, mu) ;

    Sym_tensor gamt_point(mp, CON, mp.get_bvect_spher()) ;
    gamt_point = hh_tmp ;
    gamt_point.inc_dzpuis(2) ;

    
    //--------------------------------------------------
    // Construction of Kerr Metric 
    //--------------------------------------------------

    Scalar a2(mp) ;
    Scalar b2(mp) ;

    // Parameters
    // -----------
    double mm = 2.5 ;
    double hh = 2. ;  // the radius is 1
    double aaa = pow (mm*mm-hh*hh, 0.5) ;
 
    a2 = 1. + 2.*mm/rr + (3.*mm*mm + aaa*aaa*cos2t)/(2.*rr*rr)
	+ (hh*hh*mm)/(2.*pow(rr, 3.)) + pow(hh,4.)/(16.*pow(rr,4.)) ;

    a2.std_spectral_base() ;
    a2.set_domain(0) = 1. ;

    b2 = ( pow(rr,8.) + 4.*mm*pow(rr,7.) + (7.*mm*mm + 
	   aaa*aaa*cost*cost)*pow(rr,6.) + mm*(7.*mm*mm+aaa*aaa)
	   *pow(rr,5.) + (4.*pow(mm,4.) + hh*hh*(3.*hh*hh/4.+aaa*aaa*sint
	   *sint)/2.)*pow(rr,4.) + hh*hh*mm*(2.*mm*mm-hh*hh/4.)
	   *pow(rr,3.) + pow(hh,4.)/16.*(7.*mm*mm + aaa*aaa*cost
	   *cost)*rr*rr + pow(hh,6.)*mm/16.*rr + pow(hh,8.)/256. ) 
	   / ( pow(rr,8.) + 2.*mm*pow(rr,7.) + (3.*mm*mm + aaa*aaa
           *cos2t)/2.*pow(rr,6.) + hh*hh*mm/2.*pow(rr,5.) 
	   + pow(hh,4.)/16.*pow(rr,4.)) ;

    b2.set_outer_boundary(nz-1 ,1.) ;
    b2.std_spectral_base() ;
    b2.set_domain(0) = 1. ;

    // Construction of the tilde metric
    // ---------------------------------

    Sym_tensor h_uu(mp, CON, mp.get_bvect_spher()) ;
    
    for (int i=1; i<=3; i++)
	for (int j=1; j<=i; j++){
	    if(i != j){
		h_uu.set(i,j) = 0. ;
	    }   
	}

    h_uu.set(1,1) = pow(b2/a2, 1./3.) - 1 ;
    h_uu.set(2,2) = pow(b2/a2, 1./3.) - 1 ;
    h_uu.set(3,3) = pow(a2/b2, 2./3.) - 1 ;
    h_uu.annule_domain(0) ;
    h_uu.std_spectral_base() ;

    Metric tgam (ff.con() + h_uu) ;
    gamt = tgam.cov() ;
    met_gamt = gamt ;

    cout << "norme de gamt" << endl << norme(gamt(1,1)) << endl << norme(gamt(2,1)) << endl << norme(gamt(3,1)) << endl << norme(gamt(2,2)) << endl << norme(gamt(3,2)) << endl << norme(gamt(3,3)) << endl ;

    ang_vel = 1.*aaa / (2*mm*(mm+pow(mm*mm-aaa*aaa, 0.5))) ;
    cout << "ang_vel = " << ang_vel << endl ;
    
    // Lapse function
    // --------------

    Scalar nnn (mp) ;
    nnn = ( pow(rr, 8) + 2*mm*pow(rr, 7) + (mm*mm+aaa*aaa*cost
					    *cost)*pow(rr, 6) 
	    - 0.5*hh*hh*mm*pow(rr, 5) - 0.5*hh*hh*(mm*mm+0.25*hh*hh
	    +aaa*aaa*cost*cost)*pow(rr,4) 
	    - pow(hh,4)*mm/8.*rr*rr*rr + pow(hh,4)/16.*(mm*mm +
	      aaa*aaa*cost*cost)*rr*rr + pow(hh,6)*mm/32.*rr
	    + pow(hh,8)/256.) / (pow(rr,8) + 4*mm*pow(rr,7) 
	    + (7*mm*mm+aaa*aaa*cost*cost)*pow(rr,6) 
	    + mm*(7*mm*mm+aaa*aaa)*pow(rr,5) + (4*mm*mm*mm*mm
	    + 0.5*hh*hh*(0.75*hh*hh+aaa*aaa*sint*sint))*pow(rr,4)
 	    + hh*hh*mm*(2*mm*mm-0.25*hh*hh)*pow(rr,3) 
	    + pow(hh,4)/16.*(7*mm*mm+aaa*aaa*cost*cost)*rr*rr 
	    + pow(hh,6)*mm/16.*rr + pow(hh,8)/256.) ;
			       
    nnn.std_spectral_base() ;
    nnn = pow(nnn, 0.5) ;
    nnn.set_outer_boundary(nz-1 ,1.) ;
    nnn.set_inner_boundary(1 ,0.) ;
    nnn.set_domain(0) = 1. ;
    nnn.std_spectral_base() ;
    
    nn_init = nnn ;

    // Shift vector 
    // ------------

    Scalar beta_phi (mp) ;
    beta_phi = 2*aaa*mm*unsr*unsr*sint/(a2*b2)
	*(1+mm*unsr+0.25*hh*hh*unsr*unsr) ;
    beta_phi.std_spectral_base() ;

    Vector beta_kerr (mp, CON, mp.get_bvect_spher()) ;
    beta_kerr.set(1) = 0. ;
    beta_kerr.set(2) = 0. ;
    beta_kerr.set(3) = beta_phi ;

    beta_kerr.std_spectral_base() ;
    beta_init = beta_kerr ;

    // Conformal factor Psi
    // ---------------------

    Scalar psi_kerr (pow(a2, 1./6.) * pow(b2,1./12.)) ;
    psi_kerr.std_spectral_base() ;
    psi_kerr.set_domain(0) = 1. ;
    psi_init = psi_kerr ;
    
    // End of the setup of Kerr metric
    
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
	kk_init = 0.5 * met_gam.con().derive_lie(beta_init) / nn_init ;
    else {
	Sym_tensor kk_temp (mp,  CON, mp.get_bvect_spher()) ;
	kk_temp = 0.5 * met_gam.con().derive_lie(beta_init) ;

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
    aa_init =  - psi_init*psi_init*psi_init*psi_init*kk_init 
	- 1./3. * trK * met_gamt.con() ;   // Voir le signe...

    
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
    
    //-------------------------------------------------------
    // Test of the formula for A^{ij}A_{ij} in Sergio's paper
    //-------------------------------------------------------
  
    // Put the metric tilde to the one in Sergio's paper
//    isolhor.met_kerr_perturb() ;
 
    isolhor.aa_kerr_ww(mm, aaa) ;

    // New initialisation of the metric quantities
    // --------------------------------------------

    psi_init = 1. ;
    psi_init.std_spectral_base() ;
    isolhor.set_psi_del_q(psi_init) ;
    
//    nn_init = 1. ;
//    nn_init.std_spectral_base() ;
 

    // Test of the constraints
    //------------------------

    cout<< "----------------------------------------" <<endl ;
    
    isolhor.check_hamiltonian_constraint() ;
    isolhor.check_momentum_constraint() ;

    cout<< "----------------------------------------" <<endl ;
    
    //-----------------------------------------
    //          "Call to init_data.C" 
    //-----------------------------------------

    isolhor.set_omega(ang_vel) ;
    isolhor.set_boost_x(boost_x) ;
    isolhor.set_boost_z(boost_z) ;

    isolhor.init_data(bound_nn, lim_nn, bound_psi, bound_beta, solve_lapse,
		      solve_psi, solve_shift, seuil, relax, niter) ;


    /*
    des_meridian(psi_kerr, 1.00000001, 4., "psi kerr", 12) ;
    des_meridian(isolhor.psi(), 1.00000001, 4., "psi", 13) ;
    des_meridian(isolhor.psi()-psi_kerr, 1.00000001, 4., "diff psi", 14) ;
    des_meridian(beta_phi, 1.00000001, 4., "beta(3) kerr", 15) ;
    des_meridian(isolhor.beta()(3), 1.00000001, 4., "beta(3)", 16) ;
    des_meridian(isolhor.beta()(3) - beta_phi, 1.00000001, 4., "diff beta(3)", 17) ;
    arrete() ;
    */

    // Save in a file
    // --------------
    
    FILE* fresu = fopen("resu.d", "w") ;
    mgrid.sauve(fresu) ;
    mp.sauve(fresu) ;
    isolhor.sauve(fresu, true) ;
    fclose(fresu) ;     
    
    // Test of the constraints
    //------------------------

//    isolhor.update_aa() ;

    cout<< "----------------------------------------" <<endl ;
    
    isolhor.check_hamiltonian_constraint() ;
    isolhor.check_momentum_constraint() ;

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


    //--------------------------------------
    //        Comparison
    //--------------------------------------

    cout<<"Tout va bien boudiou / Todo bien!!! (Viva Cai!)"<<endl ;

//    cout << isolhor << endl ;

    return EXIT_SUCCESS ; 
}


    
    
