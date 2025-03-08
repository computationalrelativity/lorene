/*
 * Main code for computation of equilibrium configuration of a NS-BH binary system.
 *
 */

/*
 *   Copyright (c) 2004  Philippe Grandclement, Keisuke Taniguchi,
 *              Eric Gourgoulhon
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

 

/*
 * $Id$
 * $Log$
 * Revision 1.15  2016/12/05 16:18:22  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.14  2014/10/13 08:53:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.13  2014/10/06 15:09:42  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.12  2008/09/26 08:44:04  p_grandclement
 * Mixted binaries with non vanishing spin
 *
 * Revision 1.8  2006/06/01 12:47:54  p_grandclement
 * update of the Bin_ns_bh project
 *
 * Revision 1.7  2006/04/25 07:22:00  p_grandclement
 * Various changes for the NS_BH project
 *
 * Revision 1.6  2005/10/18 13:12:34  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.4  2004/03/25 12:35:35  j_novak
 * now using namespace Unites
 *
 * Revision 1.3  2003/11/29 07:01:25  k_taniguchi
 * Large modification. However, the code has not worked yet.
 *
 * Revision 1.2  2002/12/19 14:58:20  e_gourgoulhon
 * First non-empty version.
 *
 * Revision 1.1  2002/12/18 10:33:10  e_gourgoulhon
 * Computations of NS - BH binaries
 *
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cmath>

// Lorene headers
#include "bhole.h"
#include "bin_ns_bh.h"
#include "nbr_spx.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"
#include "graphique.h"


using namespace Lorene ;

int main(int argc, char** argv) {

     using namespace Unites ;
       
    //Lecture du fichier de parametres :
     if (argc <3) {
	cout <<" Passer nom des fichiers en arguments SVP !" << endl ;
	abort() ;    } 
     //------------------------------------------------------------------
     //	    Parameters of the computation
     //------------------------------------------------------------------
    char blabla[120] ;
    double distance, precis, relax, search, m1, m2, scale_ome_local, spin, mirr, lapse_hori ;
    int itemax_equil, itemax_mp_et ;

    char* name_fich = argv[1] ;
    ifstream fpar(name_fich) ;
    fpar >> distance ; fpar.getline(blabla, 120) ;
    fpar >> m1 ; fpar >> m2 ; fpar.getline(blabla, 120) ;
    fpar >> search ; fpar.getline(blabla, 120) ;
    fpar >> precis ; fpar.getline(blabla, 120) ;
    fpar >> relax ; fpar.getline(blabla, 120) ;
    fpar >> itemax_equil ; fpar.getline(blabla, 120) ;
    fpar >> itemax_mp_et ; fpar.getline(blabla, 120) ;    
    fpar >> spin ; fpar.getline(blabla, 120) ;
    fpar >> scale_ome_local ; fpar.getline(blabla, 120) ;
    fpar >> lapse_hori ; fpar.getline(blabla, 120) ;
    fpar.close() ;
    
    //------------------------------------------------------------------
    //	    Read of the initial conditions
    //------------------------------------------------------------------
    name_fich = argv[2] ;
    FILE* fich = fopen(name_fich, "r") ;
    Mg3d mg_ns(fich) ;
    Map_et mp_ns(mg_ns, fich) ;
    Eos* peos = Eos::eos_from_file(fich) ;
    Mg3d mg_bh(fich) ;
    Map_af mp_bh(mg_bh, fich) ;
    cout << "Make an object of Bin_ns_bh" << endl ;
    Bin_ns_bh bin(mp_ns, *peos, mp_bh, fich) ;
    fclose(fich) ;

     
    //------------------------------------------------------------------
    //	    Update of the initial conditions
    //------------------------------------------------------------------
    cout << "Update the initial conditions for NS" << endl ;
    bin.set_ns().update_metric(bin.get_bh()) ;
    bin.set_bh().fait_n_comp (bin.get_ns())  ;
    bin.set_bh().fait_psi_comp (bin.get_ns()) ;
    bin.set_bh().fait_taij_auto( ) ;
    
    cout << "Update the initial conditions for BH" << endl ;
    // Initialisation of gradients of companion potentials
    // ---------------------------------------------------
    bin.set_ns().update_metric_der_comp (bin.get_bh()) ;
    bin.set_bh().update_metric (bin.get_ns()) ;
    bin.fait_tkij(0, lapse_hori) ;
   
    // Initialisation of hydro quantities for NS
    // -----------------------------------------
    bin.set_ns().equation_of_state() ;
    bin.set_ns().hydro_euler() ;
    bin.analytical_omega() ;
    bin.set_ns().kinematics (bin.get_omega(), bin.get_x_axe()) ;
    bin.set_ns().fait_d_psi() ;
    
    // Masses in good units :
    m1 *= ggrav*msol ;
    m2 *= msol ;  
    spin *= m1*m1 ;
    mirr = sqrt(0.5*(m1*m1+sqrt(m1*m1*m1*m1-spin*spin))) ;
  

    double ent_c_init = bin.get_ns().get_ent()()(0,0,0,0) ;
    bin.coal (precis, relax, itemax_equil, itemax_mp_et, ent_c_init, search, distance, m1, m2, spin, 
										scale_ome_local, 1, 0, lapse_hori) ;

    // On sauve
    char name[20] ;
    sprintf(name, "bin.dat") ;
    FILE* fresu = fopen(name, "w") ; 
    mg_ns.sauve(fresu) ;
    mp_ns.sauve(fresu) ;
    peos->sauve(fresu) ;
    mg_bh.sauve(fresu) ;
    mp_bh.sauve(fresu) ;
    bin.sauve(fresu) ;
    fclose(fresu) ;
    
    return 0 ;
}
