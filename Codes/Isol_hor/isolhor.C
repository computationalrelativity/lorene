/*
 *  Main code for Isolated Horizon in arbitrary gauge
 *
 */

/*
 *   Copyright (c) 2004  Jose Luis Jaramillo
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
  
    int niter, bound_nn, bound_psi, bound_beta, solve_lapse ;
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

    Map_af map(mgrid, bornes) ;   // Mapping construction
  	
    // Denomination of various coordinates associated with the mapping 
    // ---------------------------------------------------------------

    const Coord& r = map.r ;        // r field 
    const Coord& costt = map.cost ;  // cos(theta) field
    const Coord& sintt = map.sint ;  // sin(theta) field
    const Coord& cospp = map.cosp ;  // cos(phi) field
    const Coord& sinpp = map.sinp ;  // sin(phi) field

    Scalar cost (map) ;
    cost = costt ;
    Scalar cosp (map) ;
    cosp = cospp ;
    Scalar sint (map) ;
    sint = sintt ;
    Scalar sinp (map) ;
    sinp = sinpp ;
 
    // Flat metric f
    // -------------

    const Metric_flat& ff = map.flat_met_spher() ; 
    const Base_vect_spher& otriad = map.get_bvect_spher() ;    

    // Working stuff
    // -------------
    
    Scalar tmp_scal(map) ;
    Vector tmp_vect(map, CON, otriad) ;
    Sym_tensor tmp_sym(map, CON, otriad) ; 

    // Key function
    Mtbl usr = 1 / r ;
    Scalar unsr(map) ;
    unsr = usr ;
    
    unsr.set_domain(0) = 1 ; // scalar set to 1 in the nucleus 
    unsr.std_spectral_base() ;

    Mtbl expr = exp (- 1/((r-1)*(r-1))) ;
    Scalar expmr(map) ;
    expmr = expr ;
    expmr.std_spectral_base() ;
//    cout << "expmr" << endl << expmr << endl ;
//    des_meridian(expmr, 0, 5, "expunsr", 0) ;
//    arrete() ;

    Mtbl exprr = exp(-pow(r-3,2.)) ;
    Scalar expmrr(map) ;
    expmrr = exprr ;
 
    // Physical Parameters
    //--------------------
    
    // Set up of lapse function N
    // --------------------------
    
    Scalar nn_init(map) ; 
    nn_init = 1. - 0.5*unsr ;
    nn_init.std_spectral_base() ;    // sets standard spectral bases

    // Set up of field Psi 
    // -------------------

    Scalar psi_init(map) ; 
    psi_init =  1. + unsr ;
    psi_init.std_spectral_base() ;    // sets standard spectral bases

    // Set up of shift vector beta
    // ---------------------------    

    Vector beta_init(map, CON, map.get_bvect_spher()) ; 
    beta_init.set(1) = 0.001*unsr*unsr ;
    beta_init.set(2) = 0. ;
    beta_init.set(3) = 0. ;
    beta_init.annule_domain(0) ;
    beta_init.std_spectral_base() ;

    // TrK, TrK_point
    // --------------

    Scalar trK (map) ;
    trK = 0. ;//0.01*unsr*unsr ;
    trK.std_spectral_base() ;
    trK.inc_dzpuis(2) ;

    Scalar trK_point (map) ;
    trK_point = 0. ;
    trK_point.std_spectral_base() ;
    trK_point.inc_dzpuis(2) ;
	
    // gamt, gamt_point
    // ----------------

    Scalar khi (map) ;
    khi = 0*0.05*unsr*unsr*sint*sint ;
    khi.std_spectral_base() ;
    khi.annule_domain(0) ;
    
    Scalar mu (map) ;
    mu = 0*0.05*unsr*unsr ;
    mu.std_spectral_base() ;
    mu.annule_domain(0) ;
    
    Sym_tensor_tt hh_tmp (map, otriad, ff) ;
    hh_tmp.set_khi_mu(khi, mu) ;

    
    //Construction of a gamt
    //----------------------

    Sym_tensor gamt(map, COV, map.get_bvect_spher()) ;


    gamt = ff.cov() + hh_tmp.up_down(ff) ;
    cout << "norme de gamt" << endl << norme(gamt(1,1)) << endl << norme(gamt(2,1)) << endl << norme(gamt(3,1)) << endl << norme(gamt(2,2)) << endl << norme(gamt(3,2)) << endl << norme(gamt(3,3)) << endl ;

/*
    // Construction de la metrique de Kerr
    
    Scalar a2(map) ;
    Scalar b2(map) ;
    double hh = 2. ;
    double aaa = 0.7 ;
    double mm ;
    mm = pow(hh*hh+aaa*aaa, 0.5) ;

    const Coord& rr = map.r ;
    const Coord& theta = map.tet ;

    a2 = 1. + 2.*mm/rr + (3.*mm*mm + aaa*aaa*cos(2.*theta))/(2.*rr*rr)
	+ (hh*hh*mm)/(2.*pow(rr, 3.)) + pow(hh,4.)/(16.*pow(rr,4.)) ;

    a2.std_spectral_base() ;

    b2 = ( pow(rr,8.) + 4.*mm*pow(rr,7.) + (7.*mm*mm + 
	   aaa*aaa*cos(theta)*cos(theta))*pow(rr,6.) + mm*(7.*mm*mm+aaa*aaa)
	   *pow(rr,5.) + (4.*pow(mm,4.) + hh*hh*(3.*hh*hh/4.+aaa*aaa*sin(theta)
	   *sin(theta))/2.)*pow(rr,4.) + hh*hh*mm*(2.*mm*mm-hh*hh/4.)
	   *pow(rr,3.) + pow(hh,4.)/16.*(7.*mm*mm + aaa*aaa*cos(theta)
	   *cos(theta))*rr*rr + pow(hh,6.)*mm/16.*rr + pow(hh,8.)/256. ) 
	   / ( pow(rr,8.) + 2.*mm*pow(rr,7.) + (3.*mm*mm + aaa*aaa
           *cos(2.*theta))/2.*pow(rr,6.) + hh*hh*mm/2.*pow(rr,5.) 
	   + pow(hh,4.)/16.*pow(rr,4.)) ;

    b2.set_outer_boundary(nz-1 ,1.) ;
    b2.std_spectral_base() ;

    Sym_tensor h_uu(map, CON, map.get_bvect_spher()) ;
    
    for (int i=1; i<=3; i++)
	for (int j=1; j<=i; j++){
	    if(i != j){
		h_uu.set(i,j) = 0 ;
	    }   
	}

    h_uu.set(1,1) = pow(b2/a2, 1./3.) - 1 ;
    h_uu.set(2,2) = pow(b2/a2, 1./3.) - 1 ;
    h_uu.set(3,3) = pow(a2/b2, 2./3.) - 1 ;
    h_uu.annule_domain(0) ;
    h_uu.std_spectral_base() ;

    gamt = ff.cov() + h_uu.up_down(ff) ;

    ang_vel = aaa / (2*mm*(mm+pow(mm*mm-aaa*aaa, 0.5))) ;
    cout << "ang_vel = " << ang_vel << endl ;
    

    Scalar nnn (map) ;
    nnn = ( pow(rr, 8) + 2*mm*pow(rr, 7) + (mm*mm+aaa*aaa*cos(theta)
					    *cos(theta))*pow(rr, 6) 
	    - 0.5*hh*hh*mm*pow(rr, 5) - 0.5*hh*hh*(mm*mm+0.25*hh*hh
	    +aaa*aaa*cos(theta)*cos(theta))*pow(rr,4) 
	    - pow(hh,4)*mm/8.*rr*rr*rr + pow(hh,4)/16.*(mm*mm +
	      aaa*aaa*cos(theta)*cos(theta))*rr*rr + pow(hh,6)*mm/32.*rr
	    + pow(hh,8)/256.) / (pow(rr,8) + 4*mm*pow(rr,7) 
	    + (7*mm*mm+aaa*aaa*cos(theta)*cos(theta))*pow(rr,6) 
	    + mm*(7*mm*mm+aaa*aaa)*pow(rr,5) + (4*mm*mm*mm*mm
	    + 0.5*hh*hh*(0.75*hh*hh+aaa*aaa*sin(theta)*sin(theta)))*pow(rr,4)
 	    + hh*hh*mm*(2*mm*mm-0.25*hh*hh)*pow(rr,3) 
	    + pow(hh,4)/16.*(7*mm*mm+aaa*aaa*cos(theta)*cos(theta))*rr*rr 
			+ pow(hh,6)*mm/16.*rr + pow(hh,8)/256.) ;
			       
    nnn.std_spectral_base() ;
    nnn = pow(nnn, 0.5) ;
    nnn.set_outer_boundary(nz-1 ,1.) ;
    nnn.std_spectral_base() ;
    
    Scalar beta_phi (map) ;
    beta_phi = 2*aaa*mm*unsr*unsr*unsr/(a2*b2)*(1+mm*unsr+0.25*hh*hh*unsr) ;
    beta_phi.std_spectral_base() ;
    
    Vector beta_kerr (map, CON, map.get_bvect_spher()) ;
    beta_kerr.set(1) = 0. ;
    beta_kerr.set(2) = 0. ;
    beta_kerr.set(3) = beta_phi ;
    beta_kerr.std_spectral_base() ;

    Scalar psi_kerr (pow(a2, 1./6.) * pow(b2,1./12.)) ;
    psi_kerr.std_spectral_base() ;

*/
    // Construction de Kerr Shild
/*
    const Coord& rr = map.r ;
    const Coord& theta = map.tet ;
    double aaa = 0. ;
    double mm = 0.5 ;
    Mtbl rho2 = rr*rr + aaa*aaa*cos(theta)*cos(theta) ;

    Scalar nnn(map) ;
    nnn = 1. + 2.*mm*rr/rho2 ;
    nnn.set_outer_boundary(nz-1 ,1.) ;
    nnn.std_spectral_base() ;
    nnn = pow(nnn, -0.5) ;
    nnn.std_spectral_base() ;
//    nnn.annule_domain(0) ;

    Scalar beta_r (map) ;
    beta_r = 1. / (rho2/(2*mm*rr) + 1. ) ;
    beta_r.set_outer_boundary(nz-1 ,0.) ;
    beta_r.std_spectral_base() ;
   

    Vector beta_kerr (map, CON, map.get_bvect_spher()) ;
    beta_kerr.set(1) = beta_r ;
    beta_kerr.set(2) = 0. ;
    beta_kerr.set(3) = 0. ;
    beta_kerr.std_spectral_base() ;
    beta_kerr.annule_domain(0) ;

    Sym_tensor gamma(map, COV, map.get_bvect_spher()) ;
    gamma.set(1,1) = 1. + 2.*mm*rr/rho2 ;
    gamma.set(1,1).set_outer_boundary(nz-1 ,1.) ;
    gamma.set(2,1) = 0. ;
    gamma.set(3,1) = - (1+2.*mm*rr/rho2)*aaa*sin(theta)/rr ;
    gamma.set(3,1).set_outer_boundary(nz-1 ,0.) ;
    gamma.set(2,2) = rho2/(rr*rr) ;
    gamma.set(2,2).set_outer_boundary(nz-1 ,1.) ;
    gamma.set(3,2) = 0. ;
    gamma.set(3,3) = (rr*rr + aaa*aaa + 2.*mm*rr/rho2*aaa*aaa*sin(theta)
		      *sin(theta))/(rr*rr) ;
    gamma.set(3,3).set_outer_boundary(nz-1 ,1.) ;
    gamma.std_spectral_base() ;

    Metric met_gamma (gamma) ;     
    Scalar psi4 (map) ;
    psi4 = pow(met_gamma.determinant(), 1./3.) ;
    psi4.std_spectral_base() ;

    gamt = gamma / psi4 ;
    trK = (beta_kerr.ope_killing(met_gamma)).trace(0, 1, met_gamma) ;
    cout << "trK" << endl << norme(trK) << endl ;
    
*/    



    Metric met_gamt_tmp (gamt) ;     
        
    Scalar det_ust = pow(met_gamt_tmp.determinant(), -1./3.) ;
    det_ust.std_spectral_base() ;
    
    gamt = gamt*det_ust ;
    Metric met_gamt (gamt) ; 

    cout << "norme de gamt" << endl << norme(gamt(1,1)) << endl << norme(gamt(2,1)) << endl << norme(gamt(3,1)) << endl << norme(gamt(2,2)) << endl << norme(gamt(3,2)) << endl << norme(gamt(3,3)) << endl ;


    // Gamma-tilde_point
    //------------------
    khi = 0. ;
    khi.std_spectral_base() ;
    
    mu = 0. ;
    mu.std_spectral_base() ;
    
    hh_tmp.set_khi_mu(khi, mu) ;

    Sym_tensor gamt_point(map, CON, map.get_bvect_spher()) ;
    gamt_point = hh_tmp ;
    gamt_point.inc_dzpuis(2) ;

    // Set up of extrinsic curvature
    // -----------------------------
    
    Metric met_gam(psi_init*psi_init*psi_init*psi_init*gamt) ;
    Sym_tensor kk_init (map,  CON, map.get_bvect_spher()) ;
    for (int i=1; i<=3; i++) {
      for (int j=1; j<=i; j++) {
	kk_init.set(i,j) = 0.5 * (beta_init.derive_con(met_gam)(i,j)  
	  + beta_init.derive_con(met_gam)(j,i)) ;
      }
    }
    kk_init = kk_init/nn_init ;

    Sym_tensor aa_init (map, CON, map.get_bvect_spher()) ;
    aa_init = psi_init*psi_init*psi_init*psi_init*kk_init
	- 1./3. * trK * met_gamt.con() ;
   

    //-------------------------------------
    //     Construction of the space-time
    //-------------------------------------

    Isol_hor isolhor(map, nn_init, psi_init, beta_init, aa_init, met_gamt,
		     gamt_point, trK, trK_point, ff, 3) ;
 
  
    //-----------------------------------------
    //          "Call to init_data.C" 
    //-----------------------------------------

    isolhor.set_omega(ang_vel) ;
    isolhor.set_boost_x(boost_x) ;
    isolhor.set_boost_z(boost_z) ;

    isolhor.init_data(bound_nn, lim_nn, bound_psi, bound_beta, solve_lapse,
			  seuil, relax, niter) ;
    
    // Save in a file
    // --------------
    
    FILE* fresu = fopen("resu.d", "w") ;
    mgrid.sauve(fresu) ;
    map.sauve(fresu) ;
    isolhor.sauve(fresu, true) ;
    fclose(fresu) ;     
    
    // Test of the constraints
    //------------------------

    cout<< "----------------------------------------" <<endl ;
    
    isolhor.check_hamiltonian_constraint() ;
    isolhor.check_momentum_constraint() ;

    cout<< "----------------------------------------" <<endl ;
 

/*
    // Comparaison avec Kerr en isotropique

    des_profile(isolhor.nn(), 1.00001, 10, M_PI/2., 0., "isolhor.nn()") ;
    des_profile(nnn, 1.00001, 10, M_PI/2., 0., "nnn pour Kerr") ;
    des_profile(nnn-isolhor.nn(), 1.00001, 10, M_PI/2., 0., "diff nn") ;
    
    des_profile(isolhor.psi(), 1.00001, 10, M_PI/2., 0., "isolhor.psi()") ;
    des_profile(psi_kerr, 1.00001, 10, M_PI/2., 0., "psi pour Kerr") ;
    des_profile(psi_kerr-isolhor.psi(), 1.00001, 10, M_PI/2., 0., "diff psi") ;

    des_profile(isolhor.beta()(3), 1.00001, 10, M_PI/2., 0., "isolhor.beta()(3)") ;
    des_profile(beta_phi, 1.00001, 10, M_PI/2., 0., "beta_phi pour Kerr") ;
    des_profile(beta_phi-isolhor.beta()(3), 1.00001, 10, M_PI/2., 0., "diff beta") ;
*/

/*
    // Comparaison avec Kerr en isotropique

    des_profile(isolhor.nn(), 1.00001, 10, M_PI/2., 0., "isolhor.nn()") ;
    des_profile(nnn, 1.00001, 10, M_PI/2., 0., "nnn pour Kerr") ;
    des_profile(nnn-isolhor.nn(), 1.00001, 10, M_PI/2., 0., "diff nn") ;
    
    des_profile(isolhor.psi4(), 1.00001, 10, M_PI/2., 0., "isolhor.psi4()") ;
    des_profile(psi4, 1.00001, 10, M_PI/2., 0., "psi4 pour Kerr") ;
    des_profile(psi4-isolhor.psi4(), 1.00001, 10, M_PI/2., 0., "diff psi4") ;

    des_profile(isolhor.beta()(1), 1.00001, 10, M_PI/2., 0., "isolhor.beta()(1)") ;
    des_profile(beta_r, 1.00001, 10, M_PI/2., 0., "beta_r pour Kerr") ;
    des_profile(beta_r-isolhor.beta()(1), 1.00001, 10, M_PI/2., 0., "diff beta") ;
*/

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
    cout << "aa_adm : " << aa << endl ;  

    double aasmm = aa / mm_adm ;
    cout << "aa / M : " << aasmm << endl ; 

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


    
    
