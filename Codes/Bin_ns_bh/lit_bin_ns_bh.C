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
    double adm = bibi.adm_systeme() ;
    double adm_vol = bibi.adm_systeme_volume() ;
    double komar = bibi.komar_systeme() ;
    double moment = bibi.moment_systeme_inf() ;
    double moment_hor = bibi.moment_systeme_hor() ;
    double masse_smarr = bibi.smarr() ;
    double linear = bibi.linear_momentum_systeme_inf()(1) ;
    double area_bh = bibi.get_bh().area() ;
    double M_bh = sqrt(area_bh/16/M_PI) ; 
    double Mb_ns = bibi.get_ns().mass_b() ;
    double Mg_ns = bibi.get_ns().mass_g() ;
    double d_bh = bibi.distance_propre_axe_bh() ;
    double d_ns = bibi.distance_propre_axe_ns() ;
    
    cout << "Omega           : " << omega << endl ;  
    cout << "Masse ADM       : " << adm/ggrav/msol << " solar mass" << endl ;
    cout << "Masse ADM vol.  : " << adm_vol/ggrav/msol << " solar mass" << endl ;
    cout << "Masse Komar     : " << komar/ggrav/msol << " solar mass" << endl ;   
    cout << "Masse  smarr    : " << masse_smarr/ggrav/msol << " solar mass" << endl ;
    cout << "Moment inf      : " << moment << endl ;
    cout << "Moment hor      : " << moment_hor << endl ;
    cout << "Quantite move   : " << linear << endl ;
    cout << "Regularisation  : " << bibi.get_bh().get_regul() << endl ;
    cout << "BH area mass    : " << M_bh/msol/ggrav << " solar mass" <<  endl ;
    cout << "NS baryon mass  : " << Mb_ns/msol << " solar mass" << endl ;
    cout << "NS grav. mass   : " << Mg_ns/msol << " solar mass" << endl ;
    cout << "Distance BH     : " << d_bh << endl ;
    cout << "Distance NS     : " << d_ns << endl ;
    
    double x_bh = bibi.get_bh().get_mp().get_ori_x() ;
    double x_ns = bibi.get_ns().get_mp().get_ori_x() ;
    double regul = bibi.get_bh().get_regul() ;
    cout << "Coord bh         : " << fabs(x_bh) << endl ;
    cout << "Coord ns         : " << fabs(x_ns) << endl ;
    cout << "Regularisation   : " << regul << endl ;
    
    double mtot =  M_bh + Mb_ns*ggrav ; 
    cout << omega*mtot << " " << (adm - mtot)/mtot << endl ;
    
    double centre = 0 ;
    double taille = fabs (x_bh) + fabs(x_ns) ;
    cout << mg_ns << endl ;
    cout << mg_bh << endl ; 
 
    des_map_et (mp_ns, 0) ;
    des_map_et (mp_ns, 1) ;

    des_coupe_bin_z (bibi.get_ns().get_n_auto()(), bibi.get_bh().get_n_auto()(), 0, centre-taille, centre+taille, -taille, taille, "Lapse") ;
    des_coupe_bin_z (bibi.get_ns().get_confpsi_auto()(), bibi.get_bh().get_psi_auto()(), 0, centre-taille, centre+taille, -taille, taille, "Psi") ;
    Tenseur u_euler (bibi.get_ns().get_u_euler()) ;
    u_euler.annule(1, bibi.get_ns().get_mp().get_mg()->get_nzone()-1) ;
    des_coupe_vect_z (u_euler, 0, -1, 0.5, 1, "U euler") ; 
 
    des_coupe_z (bibi.get_ns().get_ent()(), 0, x_ns-4, x_ns+4, -4, 4, "Enthalpie") ;
    des_vect_bin_z (bibi.get_ns().get_shift_auto(),bibi.get_bh().get_shift_auto(), 0, 1000, 0.5, centre-taille, centre+taille, 
                              -taille, taille, "Shift") ;
			      
    return EXIT_SUCCESS; 
}
