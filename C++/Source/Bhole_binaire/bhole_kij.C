/*
 *   Copyright (c) 2001 Philippe Grandclement
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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


char bhole_kij_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/05/27 07:10:22  p_grandclement
 * Correction of some shadowed variables
 *
 * Revision 1.2  2002/10/16 14:36:33  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.6  2001/05/07  09:12:02  phil
 * *** empty log message ***
 *
 * Revision 1.5  2001/04/30  09:22:52  phil
 * cas omega = 0
 *
 * Revision 1.4  2001/04/27  11:44:01  phil
 * correction devant assurer la symetrie entre les deux trous
 *
 * Revision 1.3  2001/04/26  12:04:44  phil
 * *** empty log message ***
 *
 * Revision 1.2  2001/04/25  15:54:30  phil
 * corrections diverses rien de bien mechant
 *
 * Revision 1.1  2001/04/25  15:10:23  phil
 * Initial revision
 *
 *
 * $Header$
 *
 */


//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "nbr_spx.h"
#include "tenseur.h"
#include "bhole.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

//calcul de kij total. (la regularisation ayant ete faite)
void Bhole_binaire::fait_tkij () {
    
    hole1.tkij_tot.set_etat_qcq() ;
    hole2.tkij_tot.set_etat_qcq() ;
    
    // On construit a_ij a partir du shift ...
    // taij tot doit etre nul sur les deux horizons.
    hole1.fait_taij_auto () ;
    hole2.fait_taij_auto () ;
    
    // On trouve les trucs du compagnon
    hole1.taij_comp.set_etat_qcq() ;
    hole2.taij_comp.set_etat_qcq() ;
    
    Tenseur sans_dz_un (hole1.taij_auto) ;
    sans_dz_un.dec2_dzpuis() ;
    Tenseur sans_dz_deux (hole2.taij_auto) ;
    sans_dz_deux.dec2_dzpuis() ;
    
    // ON DOIT VERIFIER SI LES DEUX Aij sont alignes ou non !
    // Les bases des deux vecteurs doivent etre alignees ou non alignees :
    double orientation_un = hole1.taij_auto.get_mp()->get_rot_phi() ;
    assert ((orientation_un==0) || (orientation_un==M_PI)) ;
    double orientation_deux = hole2.taij_auto.get_mp()->get_rot_phi() ;
    assert ((orientation_deux==0) || (orientation_deux==M_PI)) ;
    int same_orient = (orientation_un == orientation_deux) ? 1 : -1 ;
    
    // Les importations :
    if (hole2.taij_auto.get_etat() == ETATZERO)
	hole1.taij_comp.set_etat_zero() ;
    else {
	hole1.taij_comp.set(0, 0).import_asymy(sans_dz_deux(0, 0)) ;
	hole1.taij_comp.set(0, 1).import_symy(sans_dz_deux(0, 1)) ;
	hole1.taij_comp.set(0, 2).import_asymy(same_orient*sans_dz_deux(0, 2)) ;
	hole1.taij_comp.set(1, 1).import_asymy(sans_dz_deux(1, 1)) ;
	hole1.taij_comp.set(1, 2).import_symy(same_orient*sans_dz_deux(1, 2)) ;
	hole1.taij_comp.set(2, 2).import_asymy(sans_dz_deux(2, 2)) ;
    }
    
     if (hole1.taij_auto.get_etat() == ETATZERO)
	hole2.taij_comp.set_etat_zero() ;
    else {
	hole2.taij_comp.set(0, 0).import_asymy(sans_dz_un(0, 0)) ;
	hole2.taij_comp.set(0, 1).import_symy(sans_dz_un(0, 1)) ;
	hole2.taij_comp.set(0, 2).import_asymy(same_orient*sans_dz_un(0, 2)) ;
	hole2.taij_comp.set(1, 1).import_asymy(sans_dz_un(1, 1)) ;
	hole2.taij_comp.set(1, 2).import_symy(same_orient*sans_dz_un(1, 2)) ;
	hole2.taij_comp.set(2, 2).import_asymy(sans_dz_un(2, 2)) ;
    }
    
    hole1.taij_comp.set_std_base() ;
    hole2.taij_comp.set_std_base() ;
    hole1.taij_comp.inc2_dzpuis() ;
    hole2.taij_comp.inc2_dzpuis() ;
    
    // Et enfin les trucs totaux...
    hole1.taij_tot = hole1.taij_auto + hole1.taij_comp ;
    hole2.taij_tot = hole2.taij_auto + hole2.taij_comp ;
    
    if ((hole1.taij_tot.get_etat() == ETATZERO) && 
	(hole2.taij_tot.get_etat() == ETATZERO)) {
	    
	    hole1.tkij_tot.set_etat_zero() ;
	    hole1.tkij_auto.set_etat_zero() ;
	    hole2.tkij_tot.set_etat_zero() ;
	    hole2.tkij_auto.set_etat_zero() ;
	}
    else {
    int nz_un = hole1.mp.get_mg()->get_nzone() ;
    int nz_deux = hole2.mp.get_mg()->get_nzone() ;
    
    Cmp ntot_un (hole1.n_tot()) ;
    ntot_un = division_xpun (ntot_un, 0) ;
    ntot_un.raccord(1) ;
    
    Cmp ntot_deux (hole2.n_tot()) ;
    ntot_deux = division_xpun (ntot_deux, 0) ;
    ntot_deux.raccord(1) ;
    
    // Boucle sur les composantes :
    for (int lig = 0 ; lig<3 ; lig++)
	for (int col = lig ; col<3 ; col++) {
	    
	    // Le sens d orientation
	    int ind = 1 ;
	    if (lig !=2)
		ind *= -1 ;
	    if (col != 2)
		ind *= -1 ;
	    if (same_orient == 1)
		ind = 1 ;
	    
	    // Pres de H1 :
	    Cmp auxi_un (hole1.taij_tot(lig, col)/2.) ;
	    auxi_un.dec2_dzpuis() ;
	    auxi_un = division_xpun (auxi_un, 0) ;
	    auxi_un = auxi_un / ntot_un ;
	    auxi_un.raccord(1) ;
	    
	    // Pres de H2 :
	    Cmp auxi_deux (hole2.taij_tot(lig, col)/2.) ;
	    auxi_deux.dec2_dzpuis() ;
	    auxi_deux = division_xpun (auxi_deux, 0) ;
	    auxi_deux = auxi_deux / ntot_deux ;
	    auxi_deux.raccord(1) ;
	    
	    // copie :
	    Cmp copie_un (hole1.taij_tot(lig, col)) ;
	    copie_un.dec2_dzpuis() ;
	    
	    Cmp copie_deux (hole2.taij_tot(lig, col)) ;
	    copie_deux.dec2_dzpuis() ;
	    
	    // Double les rayons limites :
	    double lim_un = hole1.mp.get_alpha()[1] + hole1.mp.get_beta()[1] ;
	    double lim_deux = hole2.mp.get_alpha()[1] + hole2.mp.get_beta()[1] ;
	    
	    Mtbl xabs_un (hole1.mp.xa) ;
	    Mtbl yabs_un (hole1.mp.ya) ;
	    Mtbl zabs_un (hole1.mp.za) ;
	    
	    Mtbl xabs_deux (hole2.mp.xa) ;
	    Mtbl yabs_deux (hole2.mp.ya) ;
	    Mtbl zabs_deux (hole2.mp.za) ;
	    
	    double xabs, yabs, zabs, air, theta, phi ;
	    
	    // On boucle sur les autres zones :
	    for (int l=2 ; l<nz_un ; l++) {
		
		int nr = hole1.mp.get_mg()->get_nr (l) ;
		
		if (l==nz_un-1)
		    nr -- ;
		
		int np = hole1.mp.get_mg()->get_np (l) ;
		int nt = hole1.mp.get_mg()->get_nt (l) ;
		
		for (int k=0 ; k<np ; k++)
		    for (int j=0 ; j<nt ; j++)
			for (int i=0 ; i<nr ; i++) {
			    
			    xabs = xabs_un (l, k, j, i) ;
			    yabs = yabs_un (l, k, j, i) ;
			    zabs = zabs_un (l, k, j, i) ;
			    
			    // les coordonnees du point en 2 :
			    hole2.mp.convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
			
			    if (air >= lim_deux)
				// On est loin du trou 2 : pas de pb :
				auxi_un.set(l, k, j, i) = 
		copie_un(l, k, j, i) / ntot_un (l, k, j, i)/2. ;
			    else 
			    // On est pres du trou deux :
				auxi_un.set(l, k, j, i) = 
		ind * auxi_deux.val_point (air, theta, phi) ;
			    }
			    
		// Cas infini :
		if (l==nz_un-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    auxi_un.set(nz_un-1, k, j, nr) = 0 ;
	    }
	    
	    // Le second trou :
	    for (int l=2 ; l<nz_deux ; l++) {
		
		int nr = hole2.mp.get_mg()->get_nr (l) ;
		
		if (l==nz_deux-1)
		    nr -- ;
		
		int np = hole2.mp.get_mg()->get_np (l) ;
		int nt = hole2.mp.get_mg()->get_nt (l) ;
		
		for (int k=0 ; k<np ; k++)
		    for (int j=0 ; j<nt ; j++)
			for (int i=0 ; i<nr ; i++) {
			    
			    xabs = xabs_deux (l, k, j, i) ;
			    yabs = yabs_deux (l, k, j, i) ;
			    zabs = zabs_deux (l, k, j, i) ;
			    
			    // les coordonnees du point en 1 :
			    hole1.mp.convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
			
			    if (air >= lim_un)
				// On est loin du trou 1 : pas de pb :
				auxi_deux.set(l, k, j, i) = 
		copie_deux(l, k, j, i) / ntot_deux (l, k, j, i) /2.;
			    else 
			    // On est pres du trou deux :
				auxi_deux.set(l, k, j, i) = 
		ind * auxi_un.val_point (air, theta, phi) ;
			    }
		// Cas infini :
		if (l==nz_deux-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    auxi_un.set(nz_deux-1, k, j, nr) = 0 ;
		}
	    
	    auxi_un.inc2_dzpuis() ;
	    auxi_deux.inc2_dzpuis() ;
	    
	    hole1.tkij_tot.set(lig, col) = auxi_un ;
	    hole2.tkij_tot.set(lig, col) = auxi_deux ;
	}
    
    hole1.tkij_auto.set_etat_qcq() ;
    hole2.tkij_auto.set_etat_qcq() ;
    
    for (int lig=0 ; lig<3 ; lig++)
	for (int col=lig ; col<3 ; col++) {
	    hole1.tkij_auto.set(lig, col) = hole1.tkij_tot(lig, col)*
		hole1.decouple ;
	    hole2.tkij_auto.set(lig, col) = hole2.tkij_tot(lig, col)*
		hole2.decouple ;
	}
    }
}

void Bhole_binaire::fait_decouple () {
    
    int nz_un = hole1.mp.get_mg()->get_nzone() ;
    int nz_deux = hole2.mp.get_mg()->get_nzone() ;
    
    // On determine R_limite (pour le moment en tout cas...) :
    double distance = hole1.mp.get_ori_x() - hole2.mp.get_ori_x() ;
    double lim_un = distance/2. ;
    double lim_deux = distance/2. ;
    double int_un = distance/6. ;
    double int_deux = distance/6. ;
    
    // Les fonctions de base
    Cmp fonction_f_un (hole1.mp) ;
    fonction_f_un = 0.5*pow(
	cos((hole1.mp.r-int_un)*M_PI/2./(lim_un-int_un)), 2.)+0.5 ;
    fonction_f_un.std_base_scal();
    
    Cmp fonction_g_un (hole1.mp) ;
    fonction_g_un = 0.5*pow
	(sin((hole1.mp.r-int_un)*M_PI/2./(lim_un-int_un)), 2.) ;
    fonction_g_un.std_base_scal();
    
    Cmp fonction_f_deux (hole2.mp) ;
    fonction_f_deux = 0.5*pow(
	cos((hole2.mp.r-int_deux)*M_PI/2./(lim_deux-int_deux)), 2.)+0.5 ;
    fonction_f_deux.std_base_scal();
    
    Cmp fonction_g_deux (hole2.mp) ;
    fonction_g_deux = 0.5*pow(
	sin((hole2.mp.r-int_deux)*M_PI/2./(lim_un-int_deux)), 2.) ;
    fonction_g_deux.std_base_scal();
    
    // Les fonctions totales :
    Cmp decouple_un (hole1.mp) ;
    decouple_un.allocate_all() ;
    Cmp decouple_deux (hole2.mp) ;
    decouple_deux.allocate_all() ;
    
    Mtbl xabs_un (hole1.mp.xa) ;
    Mtbl yabs_un (hole1.mp.ya) ;
    Mtbl zabs_un (hole1.mp.za) ;
	    
    Mtbl xabs_deux (hole2.mp.xa) ;
    Mtbl yabs_deux (hole2.mp.ya) ;
    Mtbl zabs_deux (hole2.mp.za) ;
	    
    double xabs, yabs, zabs, air_un, air_deux, theta, phi ;
	    
    // On boucle sur les autres zones :
    for (int l=0 ; l<nz_un ; l++) {
	int nr = hole1.mp.get_mg()->get_nr (l) ;
		
	if (l==nz_un-1)
	    nr -- ;
		
	int np = hole1.mp.get_mg()->get_np (l) ;
	int nt = hole1.mp.get_mg()->get_nt (l) ;
		
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
			    
		    xabs = xabs_un (l, k, j, i) ;
		    yabs = yabs_un (l, k, j, i) ;
		    zabs = zabs_un (l, k, j, i) ;
			    
		    // les coordonnees du point :
		    hole1.mp.convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    hole2.mp.convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;
		
		    if (air_un <= lim_un)
			if (air_un < int_un)
			    decouple_un.set(l, k, j, i) = 1 ;
			else
			// pres du trou un :
			decouple_un.set(l, k, j, i) = 
			    fonction_f_un (l, k, j, i) ;
		    else 
			if (air_deux <= lim_deux)
			    if (air_deux < int_deux)
				decouple_un.set(l, k, j, i) = 0 ;
			    else
			// On est pres du trou deux :
				decouple_un.set(l, k, j, i) = 
		fonction_g_deux.val_point (air_deux, theta, phi) ;
		
			else
			    // On est loin des deux trous :
			    decouple_un.set(l, k, j, i) = 0.5 ;
		}
			    
		// Cas infini :
		if (l==nz_un-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    decouple_un.set(nz_un-1, k, j, nr) = 0.5 ;
	    }
    
    for (int l=0 ; l<nz_deux ; l++) {
	int nr = hole2.mp.get_mg()->get_nr (l) ;
		
	if (l==nz_deux-1)
	    nr -- ;
		
	int np = hole2.mp.get_mg()->get_np (l) ;
	int nt = hole2.mp.get_mg()->get_nt (l) ;
		
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++) {
			    
		    xabs = xabs_deux (l, k, j, i) ;
		    yabs = yabs_deux (l, k, j, i) ;
		    zabs = zabs_deux (l, k, j, i) ;
			    
		    // les coordonnees du point  :
		    hole1.mp.convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    hole2.mp.convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;
		    
		    if (air_deux <= lim_deux)
			if (air_deux < int_deux)
			    decouple_deux.set(l, k, j, i) = 1 ;
			else
			// pres du trou deux :
			decouple_deux.set(l, k, j, i) = 
			    fonction_f_deux (l, k, j, i) ;
		    else 
			if (air_un <= lim_un)
			    if (air_un < int_un)
				decouple_deux.set(l, k, j, i) = 0 ;
			    else
			// On est pres du trou un :
				decouple_deux.set(l, k, j, i) = 
		fonction_g_un.val_point (air_un, theta, phi) ;
		
			else
			    // On est loin des deux trous :
			    decouple_deux.set(l, k, j, i) = 0.5 ;
		}
			    
		// Cas infini :
		if (l==nz_deux-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    decouple_deux.set(nz_un-1, k, j, nr) = 0.5 ;
   }
   
   hole1.decouple = decouple_un ;
   hole2.decouple = decouple_deux ;
}
