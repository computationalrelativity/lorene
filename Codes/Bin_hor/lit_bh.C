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
#include <stdlib.h>
#include <math.h>

// LORENE
#include "type_parite.h"
#include "nbr_spx.h"
#include "proto.h"
#include "param.h"
#include "coord.h"
#include "scalar.h"
#include "cmp.h"
#include "tensor.h"
#include "tenseur.h"
#include "isol_hor.h"
#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"


int main(int argc, char** argv) {

  using namespace Unites ;
    if (argc <2) {
	cout <<" Passer nom du ficher en arguments SVP !" << endl ;
	abort() ;
    }
    
    char* name_fich = argv[1] ;

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
    
    // Calculation of global quantities
    // --------------------------------

    double distance = map_un.get_ori_x() - map_deux.get_ori_x() ; 
    double beta = distance/bin(1).get_radius() ;
    double omega = bin.get_omega() ;
    double adm = bin.adm_mass() ;
    double komar = bin.komar_mass() ;
    double moment_inf = bin.ang_mom_adm() ;
    double moment_hor = bin(1).ang_mom_hor() + bin(2).ang_mom_hor() ;
    //   double distance_propre = bin.distance_propre() ;
    double mass_area = sqrt(bin(1).area_hor()/16/M_PI) + 
	sqrt(bin(2).area_hor()/16/M_PI) ;
        
    cout << "Beta              : " << beta << endl ;
    cout << "Omega             : " << omega << endl ;
    cout << "ADM mass          : " << adm << endl ;
    cout << "Komar mass        : " << komar << endl ;
    cout << "Mass area         : " << mass_area << endl ;
    cout << "ADM ang. mom.     : " << moment_inf << endl ;
    cout << "horizon ang.mom.  : " << moment_hor << endl ;
//    cout << "Distance propre : " << distance_propre << endl ;
    
    // Verification of Smarr :
    // -----------------------

    Scalar integrand_un (bin(1).nn().dsdr()*
	pow(bin(1).psi(), 2)) ;
    integrand_un.std_spectral_base() ;
    integrand_un.raccord(1) ;
    Scalar integrand_deux (bin(2).nn().dsdr()*
	pow(bin(2).psi(), 2)) ;
    integrand_deux.std_spectral_base() ;
    integrand_deux.raccord(1) ;
    
    double horizon = map_un.integrale_surface(integrand_un, 
					      bin(1).get_radius())+
	map_deux.integrale_surface(integrand_deux, bin(2).get_radius()) ;
	
    horizon /= 4*M_PI ;

    double j_test = (komar - horizon) / 2 / bin.get_omega() ;
    
    cout.precision(10) ;
    cout << "------------------------------------------" << endl ;
    cout << "Difference between the two J : " << fabs(moment_inf - moment_hor)
	/ moment_inf << endl ;
    cout << "Difference between ADM and Smarr  : " 
	 << fabs(moment_inf - j_test) / moment_inf << endl ;
    cout << "Difference between horizon and Smarr : " 
	 << fabs(moment_hor - j_test) / moment_inf << endl ;

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
    
    
    // Some drawings
    // -------------

    double ta = 15 ;
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
    
    ta = 13.5 ;
    Cmp dessin_un (bin(1).nn()) ;
    dessin_un.annule(0) ;
    
    Cmp dessin_deux (bin(2).nn()) ;
    dessin_deux.annule(0) ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "Lapse function (Z=0)", &surface_un, &surface_deux, 
	false, 15, 300, 300) ;
    
    dessin_un = bin(1).psi_auto() ;
    dessin_un.annule(0) ;
    
    dessin_deux = bin(2).psi_auto() ;
    dessin_deux.annule(0) ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "Conformal factor (Z=0)", &surface_un, &surface_deux, 
	false, 15, 300, 300) ;
	
    ta = 18.5 ;
    dessin_un = bin(1).aa_auto()(1, 1) ;
    dessin_un.std_base_scal() ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = bin(2).aa_auto()(1, 1) ;
    dessin_deux.std_base_scal() ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uXX\\d (Z=0)", &surface_un, &surface_deux, false
	, 15, 300, 300) ;
    
    dessin_un = bin(1).aa_auto()(2, 1) ;
    dessin_un.std_base_scal() ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = bin(2).aa_auto()(2, 1) ;
    dessin_deux.std_base_scal() ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uXY\\d (Z=0)", &surface_un, &surface_deux, false, 
	15, 300, 300) ;
    
    dessin_un = bin(1).aa_auto()(2, 2) ;
    dessin_un.std_base_scal() ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = bin(2).aa_auto()(2, 2) ;
    dessin_deux.std_base_scal() ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uYY\\d (Z=0)", &surface_un, &surface_deux, 
	false, 15, 300, 300) ;
    

    return 1; 
}
