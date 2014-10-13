/*
 * Reads a binary black hole configuration 
 *
 */

/*
 *   Copyright (c) 2002 Philippe Grandclement
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

char lit_holes_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.7  2014/10/13 08:53:52  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:09:41  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2005/08/25 12:14:11  p_grandclement
 * Addition of a new method to solve the scalar Poisson equation, based on a multi-domain Tau-method
 *
 * Revision 1.4  2004/03/25 12:35:35  j_novak
 * now using namespace Unites
 *
 * Revision 1.3  2003/11/14 08:14:52  p_grandclement
 * Correction of the codes for binary black holes in circular orbit
 *
 * Revision 1.2  2003/04/09 15:05:09  e_gourgoulhon
 * Added copyright and log messages.
 *
 *
 * Revision 1.1
 * date: 2003/03/31 16:08:23;  author: p_grandclement;  state: Exp;
 * Add codes for the binary black holes computation :
 *       - init_bh.C : computes the Misner-Lindquist solution (used as initial guess for the coal code)
 *       - param_statiques.d : example of parameter file for init_bh
 *       - coal_bh.C : computation of the binary
 *       - param_coal.d : example of the paramater file for coal_bh
 *       - lit_holes_bin.C : reads a configuration.
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
#include "coord.h"
#include "cmp.h"
#include "tenseur.h"
#include "bhole.h"
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

    FILE* fich = fopen(name_fich, "r") ;
    
    Mg3d grille (fich) ;
    Map_af map_un (grille, fich) ;
    Map_af map_deux (grille, fich) ;
    Bhole hole_un (map_un, fich) ;
    Bhole hole_deux (map_deux, fich) ;
    fclose(fich) ;
    
    assert (hole_un.get_omega() == hole_deux.get_omega()) ;
    
    Bhole_binaire systeme (map_un, map_deux) ;
    systeme.set(1) = hole_un ;
    systeme.set(2) = hole_deux ;
    systeme.set_omega(hole_un.get_omega()) ;
        
    // On initialise les grandeurs derivees :
    systeme.set(1).fait_n_comp (systeme(2)) ;
    systeme.set(1).fait_psi_comp (systeme(2)) ;
    systeme.set(2).fait_n_comp (systeme(1)) ;
    systeme.set(2).fait_psi_comp (systeme(1)) ;
    systeme.fait_decouple() ;
    systeme.fait_tkij() ;
    
    double distance = map_un.get_ori_x() - map_deux.get_ori_x() ;
    
    double beta = distance/systeme(1).get_rayon() ;
    double omega = systeme.get_omega() ;
    double adm = systeme.adm_systeme() ;
    double komar = systeme.komar_systeme() ;
    double moment_inf = systeme.moment_systeme_inf() ;
    double moment_hor = systeme.moment_systeme_hor() ;
    double distance_propre = systeme.distance_propre() ;
    double aire = sqrt(systeme(1).area()/16/M_PI) + 
    sqrt(systeme(2).area()/16/M_PI) ;
    
    cout << "Beta            : " << beta << endl ;
    cout << "Omega           : " << omega << endl ;
    cout << "Masse ADM       : " << adm << endl ;
    cout << "Masse Komar     : " << komar << endl ;
    cout << "Masse aire      : " << aire << endl ;
    cout << "Moment infiny   : " << moment_inf << endl ;
    cout << "Moment horizon  : " << moment_hor << endl ;
    cout << "Distance propre : " << distance_propre << endl ;
    
    // On verifie smarr :
    Cmp integrant_un (systeme(1).get_n_tot()().dsdr()*
	pow(systeme(1).get_psi_tot()(), 2)) ;
    integrant_un.std_base_scal() ;
    integrant_un.raccord(1) ;
    Cmp integrant_deux (systeme(2).get_n_tot()().dsdr()*
	pow(systeme(2).get_psi_tot()(), 2)) ;
    integrant_deux.std_base_scal() ;
    integrant_deux.raccord(1) ;
    
    double horizon = map_un.integrale_surface(integrant_un, systeme(1).get_rayon())+
	map_deux.integrale_surface(integrant_deux, systeme(2).get_rayon()) ;
	
    horizon /= 4*M_PI ;

    double j_test = (komar-horizon)/2/systeme.get_omega() ;
    
    cout.precision(10) ;
    cout << "------------------------------------------" << endl ;
    cout << "Ecart entre les deux J : " << fabs(moment_inf-moment_hor)/moment_inf << endl ;
    cout << "Ecart infini et Smarr  : " << fabs(moment_inf-j_test)/moment_inf << endl ;
    cout << "Ecart horizon et Smarr : " << fabs(moment_hor-j_test)/moment_inf << endl ;

    cout << "------------------------------------------" << endl ;
    cout << "Ecart Komar-ADM : " << fabs(komar-adm)/fabs(adm) << endl ;
    cout << "Ecart Kepler    : " << 4*moment_inf*pow(omega, 1./3.)/pow(adm, 5./3.)
	<< endl ;
    cout << "------------------------------------------" << endl ;
    cout << "Masse ADM       : " << adm/ggrav/msol << "Masses solaires" << endl ;
    cout << "Frequence       : " << omega/2/M_PI*f_unit << " Hz" << endl ;

    cout <<"---------------------------------------------------------" << endl ;
    double misner = adm_serie(systeme(1).get_rayon(), distance, 1e-15) ;
    cout << "Ecart ADM-Misner  : " << fabs(misner-adm)/misner << endl ;
    cout << "----------------------------------------------------------" << endl ;
    
    cout << beta << " " << omega << " " << adm << " " << moment_inf << " " <<
    aire << endl ;
    
    
      // PLOTS //
    // Les Cmp pour annuler :
    Cmp surface_un (map_un) ;
    surface_un = pow(map_un.r, 2.)-pow(systeme(1).get_rayon(), 2.) ;
    surface_un.annule(grille.get_nzone()-1) ;
    surface_un.std_base_scal() ;
    
    Cmp surface_deux (map_deux) ;
    surface_deux = pow(map_deux.r, 2.)-pow(systeme(2).get_rayon(), 2.) ;
    surface_deux.annule(grille.get_nzone()-1) ;
    surface_deux.std_base_scal() ;
    
    
    // Des dessins, des dessins, des dessins
    double ta = 15 ;
    Cmp filtre_un (map_un) ;
    int zex = grille.get_nzone()-1 ;
    filtre_un = 1 ;
    double alpha = map_un.get_alpha()[zex] ;
    double rext = 1/(-2*alpha) ;
    int nr = grille.get_nr(zex) ;
    int np = grille.get_np(zex) ;
    int nt = grille.get_nt(zex) ;
    
    double uu, coloc ;
    for (int i=0 ; i<nr-1 ; i++) {
	coloc = -cos(M_PI*i/(nr-1)) ;
	uu = alpha*(coloc-1) ;
	if (uu <= 1./2/rext)
	    for (int j=0 ; j<nt ; j++)
		for (int k=0 ; k<np ; k++)
		    filtre_un.set(zex, k, j, i) = 
		    0.5*(cos(M_PI*2*rext*(uu-1./2/rext))+1) ;
    }
    for (int j=0 ; j<nt ; j++)
	for (int k=0 ; k<np ; k++)
	    filtre_un.set(zex, k, j, nr-1) = 0 ;
    filtre_un.std_base_scal() ;
    
    Cmp filtre_deux (map_deux) ;
    filtre_deux.set_etat_qcq() ;
    filtre_deux.va = filtre_un.va ;
    
    Tenseur shift_un (systeme(1).get_shift_auto()) ;
    Cmp xa_un (map_un) ;
    xa_un = omega/2*map_un.xa ;
    Cmp ya_un (map_un) ;
    ya_un = omega/2*map_un.ya ;
    shift_un.set(0) = shift_un(0)-ya_un ;
    shift_un.set(1) = shift_un(1)+xa_un ;
     for (int i=0 ; i<3 ; i++) {
	shift_un.set(i).annule(0) ;
	shift_un.set(i)= filtre_un*shift_un(i) ;
	shift_un.set(i).set_val_inf(0) ;
	}
	
    shift_un.set_std_base() ;
    
    Tenseur shift_deux (map_deux, 1, CON, map_un.get_bvect_cart()) ;
    shift_deux.set_etat_qcq() ;
    shift_deux.set(0) = -systeme(1).get_shift_auto()(0) ;
    shift_deux.set(1) = -systeme(1).get_shift_auto()(1) ;
    shift_deux.set(2) = systeme(1).get_shift_auto()(2) ;
    Cmp xa_deux (map_deux) ;
    xa_deux = omega/2*map_deux.xa ;
    Cmp ya_deux (map_deux) ;
    ya_deux = omega/2*map_deux.ya ;
    shift_deux.set(0) = shift_deux(0)-ya_deux ;
    shift_deux.set(1) = shift_deux(1)+xa_deux ;
    for (int i=0 ; i<3 ; i++) {
	shift_deux.set(i).annule(0) ;
	shift_deux.set(i) = filtre_deux*shift_deux(i) ;
	shift_deux.set(i).set_val_inf(0) ;
	}
    shift_deux.set_std_base() ;
    
    des_vect_bin_z (shift_un, shift_deux, 0, 200, 1, -ta, ta, -ta, ta, 
    "Shift vector (Z=0)", &surface_un, &surface_deux, false, 12, 12) ;
    
    ta = 13.5 ;
    Cmp dessin_un (systeme(1).get_n_auto()()) ;
    dessin_un.annule(0) ;
    
    Cmp dessin_deux (systeme(2).get_n_auto()()) ;
    dessin_deux.annule(0) ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "Lapse function (Z=0)", &surface_un, &surface_deux, 
	false, 15, 300, 300) ;
    
    dessin_un = systeme(1).get_psi_auto()() ;
    dessin_un.annule(0) ;
    
    dessin_deux = systeme(2).get_psi_auto()() ;
    dessin_deux.annule(0) ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "Conformal factor (Z=0)", &surface_un, &surface_deux, 
	false, 15, 300, 300) ;
	
    ta = 18.5 ;
    dessin_un = systeme(1).get_tkij_auto()(0, 0) ;
    dessin_un.std_base_scal() ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = systeme(2).get_tkij_auto()(0, 0) ;
    dessin_deux.std_base_scal() ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uXX\\d (Z=0)", &surface_un, &surface_deux, false
	, 15, 300, 300) ;
    
    dessin_un = systeme(1).get_tkij_auto()(1, 0) ;
    dessin_un.std_base_scal() ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = systeme(2).get_tkij_auto()(1, 0) ;
    dessin_deux.std_base_scal() ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uXY\\d (Z=0)", &surface_un, &surface_deux, false, 
	15, 300, 300) ;
    
    dessin_un = systeme(1).get_tkij_auto()(1, 1) ;
    dessin_un.std_base_scal() ;
    dessin_un.annule(0) ;
    dessin_un.dec2_dzpuis() ;
    
    dessin_deux = systeme(2).get_tkij_auto()(1, 1) ;
    dessin_deux.std_base_scal() ;
    dessin_deux.annule(0) ;
    dessin_deux.dec2_dzpuis() ;
    
    des_coupe_bin_z (dessin_un, dessin_deux, 0, 
	-ta, ta, -ta, ta, "A\\uYY\\d (Z=0)", &surface_un, &surface_deux, 
	false, 15, 300, 300) ;

    return 1; 
}
