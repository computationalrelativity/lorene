/*
 * Main code for computation of equilibrium configuration of a NS-BH binary system.
 *
 */

/*
 *   Copyright (c) 2002  Philippe Grandclement, Keisuke Taniguchi,
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
#include "bin_ns_bh.h"
#include "nbr_spx.h"
#include "eos.h"
#include "utilitaires.h"

int  main(){


    // Identification of all the subroutines called by the code :

    // system("ident coal_ns_bh > identif.d") ;

    // For the display :
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;
    char display_normal[] = "x[0m" ; display_normal[0] = 27 ;

    #include "unites.h"
    // To avoid some compilation warnings
    if (display_bold == 0x0) {
	cout << qpig << f_unit << mevpfm3 << km << endl ;
    }

    //------------------------------------------------------------------
    //	    Parameters of the computation
    //------------------------------------------------------------------

    char blabla[120], nomini[80] ;
    int mermax, mermax_eqb, prompt, graph, fmer_stop, fmer_save, mermax_poisson ;
    int mermax_potvit, mer_masse, fmer_upd_met, ind_rel_met ;
    double seuil, relax_poisson, relax_potvit, relax, aexp_masse ;
    double mbar_voulue, mirr_voulue, fact_separ, relax_met, relax_omeg ;
    double fact_omeg_min, fact_omeg_max, thres_adapt ;

    ifstream fpar("par_coal.d") ;
    fpar.getline(blabla, 120) ;
    fpar.getline(blabla, 120) ;
    fpar.getline(nomini, 80) ;
    fpar >> fact_separ ; fpar.getline(blabla, 120);
    fpar >> mbar_voulue ; fpar.getline(blabla, 120) ;
    fpar >> mirr_voulue ; fpar.getline(blabla, 120) ;
    mbar_voulue *= msol ;
    mirr_voulue *= msol ;
    fpar.getline(blabla, 120) ;
    fpar >> mermax ; fpar.getline(blabla, 120) ;
    fpar >> relax ; fpar.getline(blabla, 120) ;
    fpar >> mermax_eqb ; fpar.getline(blabla, 120) ;
    fpar >> prompt ; fpar.getline(blabla, 120) ;
    fpar >> graph ; fpar.getline(blabla, 120) ;
    fpar >> seuil ; fpar.getline(blabla, 120) ;
    fpar >> fmer_stop ; fpar.getline(blabla, 120) ;
    fpar >> fmer_save ; fpar.getline(blabla, 120) ;
    fpar >> mermax_poisson ; fpar.getline(blabla, 120) ;
    fpar >> relax_poisson ; fpar.getline(blabla, 120) ;
    fpar >> mermax_potvit ; fpar.getline(blabla, 120) ;
    fpar >> relax_potvit ; fpar.getline(blabla, 120) ;
    fpar >> mer_masse ; fpar.getline(blabla, 120) ;
    fpar >> aexp_masse ; fpar.getline(blabla, 120) ;
    fpar >> fmer_upd_met ; fpar.getline(blabla, 120);
    fpar >> ind_rel_met ; fpar.getline(blabla, 120);
    fpar >> relax_met ; fpar.getline(blabla, 120);
    if (ind_rel_met == 0) relax_met = 1. ;
    fpar >> relax_omeg ; fpar.getline(blabla, 120);
    fpar >> fact_omeg_min ; fpar.getline(blabla, 120);
    fpar >> fact_omeg_max ; fpar.getline(blabla, 120);
    fpar >> thres_adapt ; fpar.getline(blabla, 120);
    fpar.close() ;

    cout << endl
	 << "==========================================================" << endl
	 << "                    Physical parameters                   " << endl
	 << "=========================================================="
	 << endl ;
    cout << endl << endl ;
    cout << "File containing the initial conditions : " << nomini << endl ;
    cout << "Factor by which the initial separation will be multiplied : "
	 << fact_separ << endl ;
    if ( abs(mer_masse) < mermax ) {
	cout << "Baryon mass required for the NS [M_sol] : "
	     << mbar_voulue / msol << endl ;
	cout << "Irreducible mass required for the BH  [M_sol] : "
	     << mirr_voulue / msol << endl ;
    }
    cout << endl
	 << "==========================================================" << endl
	 << "              Parameters of the computation               " << endl
	 << "=========================================================="
	 << endl ;
    cout << "Maximum number of steps in the main iteration : "
	 << mermax << endl ;
    cout << "Relaxation factor in the main iteration  : "
	 << relax << endl ;
    cout << "Maximum number of steps in Etoile_bin::equilibrium : "
	 << mermax_eqb << endl ;
    cout << "Threshold on the enthalpy relative change for ending the computation : "
	 << seuil << endl ;
    cout << "Step interval between safeguards of the whole configuration  : "
	 << fmer_save << endl ;
    cout << "Maximum number of steps in Map_et::poisson : "
	 << mermax_poisson << endl ;
    cout << "Relaxation factor in Map_et::poisson : "
	 << relax_poisson << endl ;
    cout << "Maximum number of steps in Map_radial::poisson_compact : "
	 << mermax_potvit << endl ;
    cout << "Relaxation factor in Map_radial::poisson_compact : "
	 << relax_potvit << endl ;
    cout << "Step from which the baryon mass is forced to converge : "
	 << mer_masse << endl ;
    cout << "Exponent for the increase factor of the central enthalpy : "
	 << aexp_masse << endl ;
    cout << "Step interval between metric updates : "
	 << fmer_upd_met << endl ;
    if (ind_rel_met == 1) {
	cout << "Relaxation factor of the metric : "
	     << relax_met << endl ;
    }
    else {
	cout << "No relaxation on the metric" << endl ;
    }
    cout << "Relaxation factor on Omega (orbital angular velocity) : "
	 << relax_omeg << endl ;
    cout << "Relative low bound in the omega search :  "
	 << fact_omeg_min << endl ;
    cout << "Relative high bound in the omega search : "
	 << fact_omeg_max << endl ;
    cout <<
    "Threshold on |dH/dr|_eq / |dH/dr|_pole for the adaptation of the mapping for the NS"
    << endl << thres_adapt << endl ;

    arrete(prompt) ;

    //------------------------------------------------------------------
    //	    Read of the initial conditions
    //------------------------------------------------------------------

    FILE* fich = fopen(nomini, "r") ;

    int mer_ini ;
    fread(&mer_ini, sizeof(int), 1, fich) ;

    Mg3d mg_ns(fich) ;
    Map_et mp_ns(mg_ns, fich) ;
    Eos* peos = Eos::eos_from_file(fich) ;

    Mg3d mg_bh(fich) ;
    Map_af mp_bh(mg_bh, fich) ;

    Bin_ns_bh bibi(mp_ns, *peos, mp_bh, fich) ;

    fclose(fich) ;


    //------------------------------------------------------------------
    //	    Modification of the separation between the two stars
    //------------------------------------------------------------------

        double ori_x = ((bibi.get_ns()).get_mp()).get_ori_x() ;
        ori_x *= fact_separ ;
        ((bibi.set_ns()).set_mp()).set_ori(ori_x, 0., 0.) ;

        ori_x = ((bibi.get_bh()).get_mp()).get_ori_x() ;
        ori_x *= fact_separ ;
        ((bibi.set_bh()).set_mp()).set_ori(ori_x, 0., 0.) ;


    //------------------------------------------------------------------
    //	    Update of the initial conditions
    //------------------------------------------------------------------

    // Initialisation of A^2, N, etc...
    // ---------------------------------
    (bibi.set_ns()).update_metric( bibi.get_bh() ) ;

    // Initialisation of gradients of companion potentials
    // ---------------------------------------------------

    (bibi.set_ns()).update_metric_der_comp( bibi.get_bh() ) ;

        bibi.set_ns().equation_of_state() ;
	cout << "Initial data: " << bibi << endl ; 

}
