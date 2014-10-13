/*
 * Reads a binary black hole configuration 
 *
 */

/*
 *   Copyright (c) 2005 Francois Limousin
 *                      Jose Luis Jaramillo
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

char lit_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2014/10/13 08:53:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:09:42  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2005/11/16 15:33:39  f_limousin
 * Adaptation of the code for spherical components
 *
 * Revision 1.5  2005/11/16 14:27:20  f_limousin
 * Output of boundary conditions
 *
 * Revision 1.4  2005/09/24 08:28:31  f_limousin
 * Implementation of Smarr formula. Computation of dt_psi/psi.
 *
 * Revision 1.3  2005/06/09 16:17:21  f_limousin
 * Many different changes.
 *
 * Revision 1.2  2005/03/04 09:42:25  f_limousin
 * New construction of the object Bin_hor.
 *
 * Revision 1.1  2005/03/03 13:51:56  f_limousin
 * First version
 *
 * 
 * $Header$
 *
 */

//standard
#include <cstdlib>
#include <cmath>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "param.h"
#include "tenseur.h"
#include "metric.h"
#include "isol_hor.h"
#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"


using namespace Lorene ;

int main(int argc, char** argv) {

  using namespace Unites ;
  
    if (argc <2) {
	cout <<" Passer nom du ficher en arguments SVP !" << endl ;
	abort() ;
    }
    
    char* name_fich = argv[1] ;
  
//    char* name_fich = "bin.dat" ;


    // Construction of the binary
    // --------------------------

    int depth = 3 ;
    FILE* fich = fopen(name_fich, "r") ;    
    Mg3d grid (fich) ;
    Map_af map_un (grid, fich) ;
    Map_af map_deux (grid, fich) ;
    Bin_hor bin (map_un, map_deux, fich, true, depth) ;
    fclose(fich) ;
        
    // Inititialisation of fields :
    // ---------------------------- 

    bin.set(1).n_comp (bin(2)) ;
    bin.set(1).psi_comp (bin(2)) ;
    bin.set(2).n_comp (bin(1)) ;
    bin.set(2).psi_comp (bin(1)) ;
    bin.decouple() ;
    bin.extrinsic_curvature() ;

    cout << "axi_break = " << bin(1).axi_break() << endl ;
//    abort() ;

    // Verification of boundary conditions
    // -----------------------------------

    Metric_flat fff (bin(1).get_mp().flat_met_spher()) ;

    // lapse
    cout << "Lapse boundary condition" << endl ;
    for(int j=0; j<bin(1).get_mp().get_mg()->get_nt(1); j++)
	for(int k=0; k<bin(1).get_mp().get_mg()->get_np(1); k++){
	    cout << bin(1).nn().val_grid_point(1, k, j, 0) << " " 
		 << (0.3 - bin(1).nn()).val_grid_point(1, k, j, 0)  << endl ;
	}
     arrete() ;

    // Psi
    cout << "Psi boundary condition" << endl ;
    Scalar bound_temp1 (contract(bin(1).k_dd(), 0, 1, bin(1).tradial_vect_hor()
 		   * bin(1).tradial_vect_hor(), 0, 1) / bin(1).psi()) ;
    Scalar bound_temp2 ( contract(bin(1).k_dd(), 0, 1, bin(1).tradial_vect_hor() * bin(1).tradial_vect_hor(), 0, 1) / bin(1).psi() + bin(1).psi() * bin(1).tradial_vect_hor().divergence(fff) + 4 * contract(bin(1).tradial_vect_hor() * bin(1).psi().derive_cov(fff), 0, 1) ) ;
    
    for(int j=0; j<bin(1).get_mp().get_mg()->get_nt(1); j++)
	for(int k=0; k<bin(1).get_mp().get_mg()->get_np(1); k++){
	    cout << bound_temp1.val_grid_point(1, k, j, 0) 
		 << " " << bound_temp2.val_grid_point(1, k, j, 0) 
		 << endl ;  
	}
     arrete() ;

    // Shift 
    
    Scalar bb (bin(1).get_mp()) ;
    Vector shiftt (bin(1).beta()) ;
    shiftt.change_triad(bin(1).get_mp().get_bvect_cart()) ;

    Scalar xa (bin(1).get_mp()) ;
    xa = bin(1).get_mp().xa ;
    xa.annule_domain(bin(1).get_mp().get_mg()->get_nzone()-1) ;
    Scalar ya (bin(1).get_mp()) ;
    ya = bin(1).get_mp().ya ;
    ya.annule_domain(bin(1).get_mp().get_mg()->get_nzone()-1) ;
    
    double orientation = bin(1).get_mp().get_rot_phi() ;
    assert ((orientation == 0) || (orientation == M_PI)) ;
    int aligne = (orientation == 0) ? 1 : -1 ;

    shiftt.set(1) = shiftt(1) - aligne * bin(1).get_omega() * ya ;
    shiftt.set(2) = shiftt(2) + aligne * bin(1).get_omega() * xa ;
    
    Vector normal (bin(1).radial_vect_hor()) ;
    normal.change_triad(bin(1).get_mp().get_bvect_cart()) ;
    bb = contract(shiftt, 0, normal.up_down(bin(1).gam()), 0) ;
    

    cout << "Shift boundary condition" << endl ;
    for(int j=0; j<bin(1).get_mp().get_mg()->get_nt(1); j++)
	for(int k=0; k<bin(1).get_mp().get_mg()->get_np(1); k++){
	  cout << bb.val_grid_point(1, k, j, 0) << " " << bin(1).nn().val_grid_point(1, k, j, 0) << " " << (bb - bin(1).nn()).val_grid_point(1, k, j, 0)  << endl ;
	}
     arrete() ;
	    

    // Calculation of global quantities
    // --------------------------------

    double distance = map_un.get_ori_x() - map_deux.get_ori_x() ; 
    double beta = distance/bin(1).get_radius() ;
    double omega = bin.get_omega() ;
    double adm = bin.adm_mass() ;
    double komar = bin.komar_mass() ;
    double moment_inf = bin.ang_mom_adm() ;
    double moment_hor = bin.ang_mom_hor() ;
    double mass_area = sqrt(bin(1).area_hor()/16/M_PI) + 
	sqrt(bin(2).area_hor()/16/M_PI) ;
        
    cout << "Beta              : " << beta << endl ;
    cout << "Omega             : " << omega << endl ;
    cout << "ADM mass          : " << adm << endl ;
    cout << "Komar mass        : " << komar << endl ;
    cout << "Mass area         : " << mass_area << endl ;
    cout << "ADM ang. mom.     : " << moment_inf << endl ;
    cout << "horizon ang.mom.  : " << moment_hor << endl ;
    
    // Verification of \partial_t \Psi << Omega
    // -----------------------------------------
    
    const Metric& flat1 (map_un.flat_met_spher()) ;
    const Metric& flat2 (map_deux.flat_met_spher()) ;
    
    Vector omdsdp (map_un, CON, map_un.get_bvect_cart()) ;
    Scalar yya (map_un) ;
    yya = map_un.ya ;
    Scalar xxa (map_un) ;
    xxa = map_un.xa ;
    
    double om = bin.get_omega() ;

    if (fabs(map_un.get_rot_phi()) < 1e-10){ 
	omdsdp.set(1) = - om * yya ;
	omdsdp.set(2) = om * xxa ;
	omdsdp.set(3).annule_hard() ;
    }
    else{
	omdsdp.set(1) = om * yya ;
	omdsdp.set(2) = - om * xxa ;
	omdsdp.set(3).annule_hard() ;
    }
    omdsdp.annule_domain(map_un.get_mg()->get_nzone()-1) ;

    omdsdp.set(1).set_spectral_va()
	.set_base(*(map_un.get_mg()->std_base_vect_cart()[0])) ;
    omdsdp.set(2).set_spectral_va()
	.set_base(*(map_un.get_mg()->std_base_vect_cart()[1])) ;
    omdsdp.set(3).set_spectral_va()
	.set_base(*(map_un.get_mg()->std_base_vect_cart()[2])) ;
 
    omdsdp.change_triad(map_un.get_bvect_spher()) ;
    Vector shift (bin(1).beta() + omdsdp) ;

    Scalar dt_psi (contract(shift, 0, bin(1).psi().derive_cov(flat1), 0) + 
		   bin(1).psi()*bin(1).beta().divergence(flat1)/6.) ;
		       
    cout << "Max(dt_psi/psi) : " << endl << max(dt_psi/bin(1).psi()) << endl ;
    cout << "Norme(dt_psi/psi) : " << endl << norme(dt_psi/bin(1).psi())/
      (map_un.get_mg()->get_nr(1)*map_un.get_mg()->get_nt(1)*
       map_un.get_mg()->get_np(1)) << endl ;
    cout << " omega = " << bin(1).get_omega() << endl ;
 
    // Verification of Smarr :
    // -----------------------
        
    Vector integrand_un (map_un, COV, map_un.get_bvect_spher()) ;
    integrand_un = bin(1).nn().derive_cov(flat1)*pow(bin(1).psi(), 2)
	- bin(1).nn()*contract(bin(1).k_dd(), 1,
		 bin(1).gam().radial_vect(), 0)*pow(bin(1).psi(), 2) ;
    integrand_un.std_spectral_base() ;
    //integrand_un.change_triad(map_un.get_bvect_spher()) ;

    Vector integrand_deux (map_deux, COV, map_deux.get_bvect_spher()) ;
    integrand_deux = bin(2).nn().derive_cov(flat2)*pow(bin(2).psi(), 2)
	- bin(2).nn()*contract(bin(2).k_dd(), 1,
		      bin(2).gam().radial_vect(), 0)*pow(bin(2).psi(), 2) ;
    integrand_deux.std_spectral_base() ;
    //integrand_deux.change_triad(map_deux.get_bvect_spher()) ;

    double horizon = map_un.integrale_surface(integrand_un(1), 
					      bin(1).get_radius())+
	map_deux.integrale_surface(integrand_deux(1), bin(2).get_radius()) ;
	
    horizon /= 4*M_PI ;

    double j_test = (komar - horizon) / 2. / bin.get_omega() ;
    
    cout.precision(10) ;
    cout << "------------------------------------------" << endl ;
    cout << "Difference between the two J : " << fabs(moment_inf - moment_hor)
	/ moment_inf << endl ;
    cout << "Difference between ADM and Smarr  : " 
	 << fabs(moment_inf - j_test) / moment_inf << endl ;
    cout << "Difference between horizon and Smarr : " 
	 << fabs(moment_hor - j_test) / moment_hor << endl ;

    cout << "------------------------------------------" << endl ;
    cout << "Difference Komar-ADM : " << fabs(komar-adm)/fabs(adm) << endl ;
    cout << "Comparison to Kepler    : " << 4*moment_inf * pow(omega, 1./3.)
	/ pow(adm, 5./3.)
	<< endl ;
    cout << "------------------------------------------" << endl ;
    cout << "ADM mass      : " << adm/ggrav/msol << " solar masses"<< endl ;
    cout << "Frequence       : " << omega/2/M_PI*f_unit << " Hz" << endl ;

    cout <<"--------------------------------------------------------" << endl ;

    
    // Definition of the surface
    // -------------------------

    Cmp surface_un (map_un) ;
    surface_un = pow(map_un.r, 2.)-pow(bin(1).get_radius(), 2.) ;
    surface_un.annule(grid.get_nzone()-1) ;
    surface_un.std_base_scal() ;
    
    Cmp surface_deux (map_deux) ;
    surface_deux = pow(map_deux.r, 2.)-pow(bin(2).get_radius(), 2.) ;
    surface_deux.annule(grid.get_nzone()-1) ;
    surface_deux.std_base_scal() ;
    
    // ---------------------
    // Some drawings
    // ---------------------


    // Filter
    // -------

    double ta = 18.5 ;
    Scalar filtre_un (map_un) ;
    int zex = grid.get_nzone()-1 ;
    filtre_un = 1. + 1e-15 ;
    double alpha = map_un.get_alpha()[zex] ;
    double rext = 1/(-2*alpha) ;
    int nr = grid.get_nr(zex) ;
    int np = grid.get_np(zex) ;
    int nt = grid.get_nt(zex) ;
    
    double uu, coloc ;
    for (int i=0 ; i<nr-1 ; i++) {
	coloc = -cos(M_PI*i/(nr-1)) ;
	uu = alpha*(coloc-1) ;
	if (uu <= 1./2/rext)
	    for (int j=0 ; j<nt ; j++)
		for (int k=0 ; k<np ; k++)
		    filtre_un.set_grid_point(zex, k, j, i) = 
		    0.5*(cos(M_PI*2*rext*(uu-1./2/rext))+1) ;
    }
    for (int j=0 ; j<nt ; j++)
	for (int k=0 ; k<np ; k++)
	    filtre_un.set_grid_point(zex, k, j, nr-1) = 0 ;
    filtre_un.std_spectral_base() ;
    
    Scalar filtre_deux (map_deux) ;
    filtre_deux.set_etat_qcq() ;
    filtre_deux.set_spectral_va() = filtre_un.get_spectral_va() ;
    
    // Shift
    // ------

    Vector shift_un (bin(1).beta_auto()) ;
    shift_un.change_triad(map_un.get_bvect_cart()) ;

    Scalar xa_un (map_un) ;
    xa_un = omega/2*map_un.xa ;
    Scalar ya_un (map_un) ;
    ya_un = omega/2*map_un.ya ;
    shift_un.set(1) = shift_un(1)-ya_un ;
    shift_un.set(2) = shift_un(2)+xa_un ;
     for (int i=1 ; i<=3 ; i++) {
	shift_un.set(i).annule_domain(0) ;
	shift_un.set(i)= filtre_un*shift_un(i) ;
	shift_un.set(i).set_outer_boundary(zex, 0) ;
	}
	
    Vector shift_deux (bin(2).beta_auto()) ;
    shift_deux.change_triad(map_deux.get_bvect_cart()) ;
    shift_deux.change_triad(map_un.get_bvect_cart()) ;
    Scalar xa_deux (map_deux) ;
    xa_deux = omega/2*map_deux.xa ;
    Scalar ya_deux (map_deux) ;
    ya_deux = omega/2*map_deux.ya ;
    shift_deux.set(1) = shift_deux(1)-ya_deux ;
    shift_deux.set(2) = shift_deux(2)+xa_deux ;
     for (int i=1 ; i<=3 ; i++) {
	shift_deux.set(i).annule_domain(0) ;
	shift_deux.set(i)= filtre_deux*shift_deux(i) ;
	shift_deux.set(i).set_outer_boundary(zex, 0) ;
	}
	
    shift_deux.std_spectral_base() ;

    Tenseur beta_un (map_un, 1, CON, map_un.get_bvect_cart()) ;
    beta_un.set_etat_qcq() ;
    Cmp beta1_x (shift_un(1)) ;
    Cmp beta1_y (shift_un(2)) ;
    Cmp beta1_z (shift_un(3)) ;
    beta_un.set(0) = beta1_x ;
    beta_un.set(1) = beta1_y ;
    beta_un.set(2) = beta1_z ;
    beta_un.set_std_base() ;

    Tenseur beta_deux (map_deux, 1, CON, map_un.get_bvect_cart()) ;
    beta_deux.set_etat_qcq() ;
    Cmp beta2_x (shift_deux(1)) ;
    Cmp beta2_y (shift_deux(2)) ;
    Cmp beta2_z (shift_deux(3)) ;
    beta_deux.set(0) = beta2_x ;
    beta_deux.set(1) = beta2_y ;
    beta_deux.set(2) = beta2_z ;
    beta_deux.set_std_base() ;

    des_vect_bin_z (beta_un, beta_deux, 0, 200, 1, -ta, ta, -ta, ta, 
    "Shift vector (Z=0)", &surface_un, &surface_deux, false, 12, 12) ;
    

    // Lapse 
    // -------

    ta = 18.5 ;
    Cmp dessin_un (bin(1).n_auto()) ;
    dessin_un.annule(0) ;
    
    Cmp dessin_deux (bin(2).n_auto()) ;
    dessin_deux.annule(0) ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "Lapse function (Z=0)", &surface_un, &surface_deux, 
	false, 15, 300, 300) ;
     
    // Psi
    // -----

    dessin_un = bin(1).psi_auto() ;
    dessin_un.annule(0) ;
    
    dessin_deux = bin(2).psi_auto() ;
    dessin_deux.annule(0) ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "Conformal factor (Z=0)", &surface_un, &surface_deux, 
	false, 15, 300, 300) ;
    

    // Extrinsic curvature
    // --------------------

    Sym_tensor aa_auto_un (bin(1).aa_auto()) ;
    Sym_tensor aa_auto_deux (bin(2).aa_auto()) ;

    ta = 18.5 ;
    dessin_un = aa_auto_un(1, 1) ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = aa_auto_deux(1, 1) ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
 
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\urr\\d (Z=0)", &surface_un, &surface_deux, false
	, 20, 300, 300) ;
    
    dessin_un = aa_auto_un(1, 3) ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = aa_auto_deux(1, 3) ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\urp\\d (Z=0)", &surface_un, &surface_deux, false, 
	20, 300, 300) ;
    
    // In cartesian ccordinates.

    aa_auto_un.change_triad(bin(1).get_mp().get_bvect_cart()) ;
    aa_auto_deux.change_triad(bin(2).get_mp().get_bvect_cart()) ;

    dessin_un = aa_auto_un(1, 1) ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = aa_auto_deux(1, 1) ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uXX\\d (Z=0)", &surface_un, &surface_deux, false
	, 20, 300, 300) ;
 

    dessin_un = aa_auto_un(2, 1) ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = aa_auto_deux(2, 1) ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uXY\\d (Z=0)", &surface_un, &surface_deux, false
	, 20, 300, 300) ;

 
    dessin_un = aa_auto_un(2, 2) ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = aa_auto_deux(2, 2) ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uYY\\d (Z=0)", &surface_un, &surface_deux, false
	, 20, 300, 300) ;    


    dessin_un = aa_auto_un(3, 3) ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = aa_auto_deux(3, 3) ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uZZ\\d (Z=0)", &surface_un, &surface_deux, false
	, 20, 300, 300) ;    

    return 1; 
}
