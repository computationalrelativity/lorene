/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose Luis Jaramillo                       
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


char binhor_kij_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/12/31 15:41:54  f_limousin
 * Correction of an error
 *
 * Revision 1.1  2004/12/29 16:12:03  f_limousin
 * First version
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
#include "tensor.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

void Bin_hor::extrinsic_curvature () {
    
    // Computation of A^{ij}_auto
    Sym_tensor aa_auto_un (hole1.mp, CON, hole1.mp.get_bvect_spher()) ;
    Sym_tensor aa_auto_deux (hole2.mp, CON, hole2.mp.get_bvect_spher()) ;
    
    aa_auto_un = ( hole1.beta_auto().ope_killing_conf(hole1.tgam()) + 
	hole1.gamt_point*hole1.decouple ) / (2.* hole1.nn()) ;            
    aa_auto_deux = ( hole2.beta_auto().ope_killing_conf(hole2.tgam()) + 
    hole2.gamt_point*hole2.decouple ) / (2.* hole2.nn()) ;            

    double ttime = hole1.the_time[hole1.jtime] ;
    hole1.aa_auto_evol.update(aa_auto_un, hole1.jtime, ttime) ;
    hole2.aa_auto_evol.update(aa_auto_deux, hole2.jtime, ttime) ;

    // Computation of A^{ij}_comp

    aa_auto_un.dec_dzpuis(2) ;
    aa_auto_deux.dec_dzpuis(2) ;

    Sym_tensor aa_comp_un (hole1.mp, CON, hole1.mp.get_bvect_cart()) ;
    aa_comp_un.set_etat_qcq() ;
    Sym_tensor aa_comp_deux (hole2.mp, CON, hole2.mp.get_bvect_cart()) ;
    aa_comp_deux.set_etat_qcq() ;
    
    aa_auto_deux.change_triad(hole2.mp.get_bvect_cart()) ;
    aa_auto_deux.change_triad(hole1.mp.get_bvect_cart()) ;
    assert(*(aa_auto_deux.get_triad()) == *(aa_comp_un.get_triad())) ;
    // importations :
    aa_comp_un.set(1, 1).import_asymy(aa_auto_deux(1, 1)) ;
    aa_comp_un.set(1, 2).import_symy(aa_auto_deux(1, 2)) ;
    aa_comp_un.set(1, 3).import_asymy(aa_auto_deux(1, 3)) ;
    aa_comp_un.set(2, 2).import_asymy(aa_auto_deux(2, 2)) ;
    aa_comp_un.set(2, 3).import_symy(aa_auto_deux(2, 3)) ;
    aa_comp_un.set(3, 3).import_asymy(aa_auto_deux(3, 3)) ;

    aa_comp_un.std_spectral_base() ;
    aa_comp_un.inc_dzpuis(2) ;
    aa_comp_un.change_triad(hole1.mp.get_bvect_spher()) ;

    aa_auto_un.change_triad(hole1.mp.get_bvect_cart()) ;
    aa_auto_un.change_triad(hole2.mp.get_bvect_cart()) ;
    assert(*(aa_auto_un.get_triad()) == *(aa_comp_deux.get_triad())) ;
    // importations :
    aa_comp_deux.set(1, 1).import_asymy(aa_auto_un(1, 1)) ;
    aa_comp_deux.set(1, 2).import_symy(aa_auto_un(1, 2)) ;
    aa_comp_deux.set(1, 3).import_asymy(aa_auto_un(1, 3)) ;
    aa_comp_deux.set(2, 2).import_asymy(aa_auto_un(2, 2)) ;
    aa_comp_deux.set(2, 3).import_symy(aa_auto_un(2, 3)) ;
    aa_comp_deux.set(3, 3).import_asymy(aa_auto_un(3, 3)) ;

    aa_comp_deux.std_spectral_base() ;
    aa_comp_deux.inc_dzpuis(2) ;
    aa_comp_deux.change_triad(hole2.mp.get_bvect_spher()) ;

    hole1.aa_comp_evol.update(aa_comp_un, hole1.jtime, ttime) ;
    hole2.aa_comp_evol.update(aa_comp_deux, hole2.jtime, ttime) ;

    
    // Computation of A^{ij}_ total
    hole1.aa_evol.update(hole1.aa_auto() + hole1.aa_comp(), 
			 hole1.jtime, ttime) ;
    hole2.aa_evol.update(hole2.aa_auto() + hole2.aa_comp(), 
			 hole2.jtime, ttime) ;


}

void Bin_hor::decouple () {
    
    int nz_un = hole1.mp.get_mg()->get_nzone() ;
    int nz_deux = hole2.mp.get_mg()->get_nzone() ;
    
    // We determine R_limite :
    double distance = hole1.mp.get_ori_x() - hole2.mp.get_ori_x() ;
    double lim_un = distance/2. ;
    double lim_deux = distance/2. ;
    double int_un = distance/6. ;
    double int_deux = distance/6. ;
    
    // The functions used.
    Scalar fonction_f_un (hole1.mp) ;
    fonction_f_un = 0.5*pow(
	cos((hole1.mp.r-int_un)*M_PI/2./(lim_un-int_un)), 2.)+0.5 ;
    fonction_f_un.std_spectral_base();
    
    Scalar fonction_g_un (hole1.mp) ;
    fonction_g_un = 0.5*pow
	(sin((hole1.mp.r-int_un)*M_PI/2./(lim_un-int_un)), 2.) ;
    fonction_g_un.std_spectral_base();
    
    Scalar fonction_f_deux (hole2.mp) ;
    fonction_f_deux = 0.5*pow(
	cos((hole2.mp.r-int_deux)*M_PI/2./(lim_deux-int_deux)), 2.)+0.5 ;
    fonction_f_deux.std_spectral_base();
    
    Scalar fonction_g_deux (hole2.mp) ;
    fonction_g_deux = 0.5*pow(
	sin((hole2.mp.r-int_deux)*M_PI/2./(lim_un-int_deux)), 2.) ;
    fonction_g_deux.std_spectral_base();
    
    // The functions total :
    Scalar decouple_un (hole1.mp) ;
    decouple_un.allocate_all() ;
    Scalar decouple_deux (hole2.mp) ;
    decouple_deux.allocate_all() ;
    
    Mtbl xabs_un (hole1.mp.xa) ;
    Mtbl yabs_un (hole1.mp.ya) ;
    Mtbl zabs_un (hole1.mp.za) ;
	    
    Mtbl xabs_deux (hole2.mp.xa) ;
    Mtbl yabs_deux (hole2.mp.ya) ;
    Mtbl zabs_deux (hole2.mp.za) ;
	    
    double xabs, yabs, zabs, air_un, air_deux, theta, phi ;
	    
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
			    
		    // Coordinates of the point
		    hole1.mp.convert_absolute 
			(xabs, yabs, zabs, air_un, theta, phi) ;
		    hole2.mp.convert_absolute 
			(xabs, yabs, zabs, air_deux, theta, phi) ;
		
		    if (air_un <= lim_un)
			if (air_un < int_un)
			    decouple_un.set_grid_point(l, k, j, i) = 1 ;
			else
			// Close to hole 1 :
			decouple_un.set_grid_point(l, k, j, i) = 
			    fonction_f_un.val_grid_point(l, k, j, i) ;
		    else 
			if (air_deux <= lim_deux)
			    if (air_deux < int_deux)
				decouple_un.set_grid_point(l, k, j, i) = 0 ;
			    else
			// Close to hole 2 :
				decouple_un.set_grid_point(l, k, j, i) = 
		fonction_g_deux.val_point (air_deux, theta, phi) ;
		
			else
			    // Far from both holes :
			    decouple_un.set_grid_point(l, k, j, i) = 0.5 ;
		}
			    
		// Case infinity :
		if (l==nz_un-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			    decouple_un.set_grid_point(nz_un-1, k, j, nr)=0.5 ;
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
			    decouple_deux.set_grid_point(l, k, j, i) = 1 ;
			else
			// pres du trou deux :
			decouple_deux.set_grid_point(l, k, j, i) = 
			    fonction_f_deux.val_grid_point(l, k, j, i) ;
		    else 
			if (air_un <= lim_un)
			    if (air_un < int_un)
				decouple_deux.set_grid_point(l, k, j, i) = 0 ;
			    else
			// On est pres du trou un :
				decouple_deux.set_grid_point(l, k, j, i) = 
			 fonction_g_un.val_point (air_un, theta, phi) ;
		
			else
			    // On est loin des deux trous :
			    decouple_deux.set_grid_point(l, k, j, i) = 0.5 ;
		}
			    
		// Cas infini :
		if (l==nz_deux-1)
		    for (int k=0 ; k<np ; k++)
			for (int j=0 ; j<nt ; j++)
			 decouple_deux.set_grid_point(nz_un-1, k, j, nr)=0.5 ;
   }
   
   hole1.decouple = decouple_un ;
   hole2.decouple = decouple_deux ;
}
