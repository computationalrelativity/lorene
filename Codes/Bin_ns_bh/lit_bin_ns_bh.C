/*
 * Code for reading a equilibrium configuration of a NS-BH binary system.
 *
 */

/*
 *   Copyright (c) 2005  Philippe Grandclement
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

char lit_bin_ns_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2008/09/26 08:44:04  p_grandclement
 * Mixted binaries with non vanishing spin
 *
 * Revision 1.6  2006/09/06 11:52:46  p_grandclement
 * Update of the Bin_ns_bh codes
 *
 * Revision 1.4  2006/04/25 07:22:00  p_grandclement
 * Various changes for the NS_BH project
 *
 * Revision 1.2  2005/10/18 13:12:34  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.1  2005/08/29 15:10:19  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 * *
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
#include "cmp.h"
#include "tenseur.h"
#include "bhole.h"
#include "bin_ns_bh.h"
#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"
#include "eos.h"


int main(int argc, char** argv) {

  using namespace Unites ;
    if (argc <2) {
	cout <<" Passer nom du ficher en arguments SVP !" << endl ;
	abort() ;
    }
    
    char* name_fich = argv[1] ;

    FILE* fich = fopen(name_fich, "r") ;
    Mg3d mg_ns(fich) ;
    Map_et mp_ns(mg_ns, fich) ;
    Eos* peos = Eos::eos_from_file(fich) ;

    Mg3d mg_bh(fich) ;
    Map_af mp_bh(mg_bh, fich) ;

    Bin_ns_bh bibi(mp_ns, *peos, mp_bh, fich) ;    
    fclose(fich) ;
    
   // On initialise les grandeurs derivees 
    bibi.set_ns().update_metric( bibi.get_bh() ) ;
    bibi.set_bh().fait_n_comp( bibi.get_ns() ) ;
    bibi.set_bh().fait_psi_comp( bibi.get_ns() ) ;
    bibi.set_bh().fait_taij_auto( ) ;
    bibi.set_ns().update_metric_der_comp( bibi.get_bh() ) ;
    bibi.set_bh().update_metric (bibi.get_ns()) ;
    bibi.fait_tkij() ;
    
   // Initialisation of hydro quantities for NS
    // -----------------------------------------
    bibi.set_ns().equation_of_state() ;
    bibi.set_ns().kinematics(bibi.get_omega(), bibi.get_x_axe()) ;
    bibi.set_ns().fait_d_psi() ;
    bibi.set_ns().hydro_euler() ;
  
    double omega = bibi.get_omega() ;
    double omega_local = bibi.get_bh().get_omega_local() ;
    double adm = bibi.adm_systeme() ;
    double adm_vol = bibi.adm_systeme_volume() ;
    double komar = bibi.komar_systeme() ;
    double moment = bibi.moment_systeme_inf() ;
    double moment_hor = bibi.moment_systeme_hor() ;
    double masse_smarr = bibi.smarr() ;
    double area_bh = bibi.get_bh().area() ;
    double M_bh = sqrt(area_bh/16/M_PI) ; 
    double Mb_ns = bibi.get_ns().mass_b() ;
    double Mg_ns = bibi.get_ns().mass_g() ;
    
    // Calcul de Ki 
    int nzet = bibi.get_ns().get_nzet() ;
    Cmp dent (bibi.get_ns().get_ent()().dsdr()) ;
    int nr = mg_ns.get_nr(nzet-1) ;
    int nt = mg_ns.get_nt(nzet-1) ;
    double ki = dent(nzet-1, 0, nt-1, nr-1) / dent(nzet-1, 0,0, nr-1) ;

    cout.precision(10) ;
    cout << "Omega           : " << omega*f_unit << " Hz" << endl ;
    cout << "Local omega     : " << omega_local*f_unit << " Hz" << endl ;
    cout << "Masse ADM       : " << adm/ggrav/msol << " solar mass" << endl ;
    cout << "Masse ADM vol.  : " << adm_vol/ggrav/msol << " solar mass" << endl ;
    cout << "Masse Komar     : " << komar/ggrav/msol << " solar mass" << endl ;   
    cout << "Masse  Smarr    : " << masse_smarr/ggrav/msol << " solar mass" << endl ;
    cout << "Moment inf      : " << moment << " (Lorene units)" << endl ;
    cout << "Moment hor      : " << moment_hor << " (Lorene units)" << endl ;
    cout << "Regularization  : " << bibi.get_bh().get_regul() << endl ;
    cout << "BH area mass    : " << M_bh/msol/ggrav << " solar mass" <<  endl ;
    cout << "NS baryon mass  : " << Mb_ns/msol << " solar mass" << endl ;
    cout << "NS grav. mass   : " << Mg_ns/msol << " solar mass" << endl ;
    cout << "Deformation     : " << ki << endl ;
    
    double x_bh = bibi.get_bh().get_mp().get_ori_x() ;
    double x_ns = bibi.get_ns().get_mp().get_ori_x() ;
    double centre = (x_ns+x_bh)/2 ;
    double taille = 1.3*(fabs (x_bh) + fabs(x_ns)) ; 
  // PLOTS //
    // Les Cmp pour annuler :
    Cmp surface_bh (mp_bh) ;
    surface_bh = pow(mp_bh.r, 2.)-pow(bibi.get_bh().get_rayon(), 2.) ;
    surface_bh.annule(mg_bh.get_nzone()-1) ;
    surface_bh.std_base_scal() ;
    
    Cmp surface_ns (bibi.get_ns().get_ent()()) ;
    
    des_coupe_bin_z (bibi.get_ns().get_n_auto()(), bibi.get_bh().get_n_auto()(), 0, x_bh-1, x_bh+1, -1, 1, "Lapse function (Z=0)", 
    	&surface_ns, &surface_bh, 
	false, 15, 300, 300) ;
	
	/*
	taille -= 3 ;
       Cmp p1 (bibi.get_ns().get_taij_auto()(0,0)) ;
        p1.dec2_dzpuis() ;
        Cmp p2 (bibi.get_bh().get_taij_auto()(0,0)) ;
        p2.dec2_dzpuis() ;
	des_coupe_bin_z (p1, p2, 0, centre-taille, centre+taille, -taille, taille, "A\\uXX\\d (Z=0)", &surface_ns, &surface_bh, false, 15, 300, 300) ;
	
    
    des_coupe_bin_z (bibi.get_ns().get_n_auto()(), bibi.get_bh().get_n_auto()(), 0, centre-taille, centre+taille, -taille, taille, "Lapse") ;
    des_coupe_bin_z (bibi.get_ns().get_confpsi_auto()(), bibi.get_bh().get_psi_auto()(), 0, centre-taille, centre+taille, -taille, taille, "Psi") ;
    Tenseur u_euler (bibi.get_ns().get_u_euler()) ;
    u_euler.annule(1, bibi.get_ns().get_mp().get_mg()->get_nzone()-1) ;
    des_coupe_vect_z (u_euler, 0, -1, 0.5, 1, "U euler") ; 
 
    des_coupe_z (bibi.get_ns().get_ent()(), 0, x_ns-4, x_ns+4, -4, 4, "Enthalpie") ; 
    des_vect_bin_z (bibi.get_ns().get_shift_auto(),bibi.get_bh().get_shift_auto(), 0, 500, 0.5, centre-taille, centre+taille, 
                              -taille, taille, "Shift") ;
			      
    for (int i=0 ; i<2 ; i++)
      for (int j=0 ; j<2 ; j++) {
	Cmp p1 (bibi.get_ns().get_taij_auto()(i,j)) ;
        p1.dec2_dzpuis() ;
        Cmp p2 (bibi.get_bh().get_taij_auto()(i,j)) ;
        p2.dec2_dzpuis() ;
	des_coupe_bin_z (p1, p2, 0, centre-taille, centre+taille, -taille, taille, "A_ij") ;
      }
   */
   
    return EXIT_SUCCESS; 
}
