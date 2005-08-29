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

char coal_ns_bh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2005/08/29 15:10:19  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
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
#include <stdlib.h>
#include <math.h>

// Lorene headers
#include "bhole.h"
#include "bin_ns_bh.h"
#include "nbr_spx.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"
#include "graphique.h"

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
    double precis, relax, search_ome, m1, m2 ;
    int nbr_ome, itemax_equil, itemax_mp_et ;

    char* name_fich = argv[1] ;
    //char* name_fich = "par_coal.d" ;
    ifstream fpar(name_fich) ;
    fpar >> m1 ; fpar >> m2 ; fpar.getline(blabla, 120) ;
    fpar >> search_ome ; fpar.getline(blabla, 120) ;
    fpar >> precis ; fpar.getline(blabla, 120) ;
    fpar >> relax ; fpar.getline(blabla, 120) ;
    fpar >> nbr_ome ; fpar.getline(blabla, 120) ;
    fpar >> itemax_equil ; fpar.getline(blabla, 120) ;
    fpar >> itemax_mp_et ; fpar.getline(blabla, 120) ;
    fpar.close() ;
    
    //------------------------------------------------------------------
    //	    Read of the initial conditions
    //------------------------------------------------------------------
    name_fich = argv[2] ;
    //name_fich = "statiques.dat" ;
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
    bin.fait_tkij() ;
   
    // Initialisation of hydro quantities for NS
    // -----------------------------------------
    bin.set_ns().equation_of_state() ;
    bin.set_ns().kinematics (bin.get_omega(), bin.get_x_axe()) ;
    bin.set_ns().fait_d_psi() ;
    bin.set_ns().hydro_euler() ;
    
    bin.coal (precis, relax, itemax_equil, itemax_mp_et, nbr_ome, search_ome, m1, m2, 1) ;

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
