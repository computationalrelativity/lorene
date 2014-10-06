/*
 *   Copyright (c) 2005 Francois Limousin
 *                      Jose Luis Jaramillo
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


char single_regul_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2014/10/06 15:13:11  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.1  2007/04/13 15:28:35  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 *
 * $Header$
 *
 */


//Standard
#include <cstdlib>
#include <cmath>

//Lorene
#include "isol_hor.h"
#include "nbr_spx.h"
#include "tensor.h"

double Single_hor::regularisation (const Vector& shift_auto_temp, 
				 const Vector& shift_comp_temp, double om) {
    
    Vector shift_auto(shift_auto_temp) ;
    shift_auto.change_triad(shift_auto.get_mp().get_bvect_cart()) ;
    Vector shift_comp(shift_comp_temp) ;
    shift_comp.change_triad(shift_comp.get_mp().get_bvect_cart()) ;
    Vector shift_old (shift_auto) ;
    
    double orientation = shift_auto.get_mp().get_rot_phi() ;
    assert ((orientation==0) || (orientation == M_PI)) ;
    double orientation_autre = shift_comp.get_mp().get_rot_phi() ;
    assert ((orientation_autre==0) || (orientation_autre == M_PI)) ;
    
    int alignes = (orientation == orientation_autre) ? 1 : -1 ;
    
    int np = shift_auto.get_mp().get_mg()->get_np(1) ;
    int nt = shift_auto.get_mp().get_mg()->get_nt(1) ;
    int nr = shift_auto.get_mp().get_mg()->get_nr(1) ;
    
    // Minimisation of the derivative of the shift on r
    Vector shift_tot (shift_auto.get_mp(), CON, *shift_auto.get_triad()) ;
    shift_tot.set(1).import(alignes*shift_comp(1)) ;
    shift_tot.set(2).import(alignes*shift_comp(2)) ;
    shift_tot.set(3).import(shift_comp(3)) ;

    shift_tot = shift_tot + shift_auto ;
 
    double indic = (orientation == 0) ? 1 : -1 ;
    
    Vector tbi (shift_tot) ;
    if (om != 0) {
	for (int i=1 ; i<=3 ; i++) {
	    tbi.set(i).set_spectral_va().coef_i() ;
	    tbi.set(i).set_spectral_va().set_etat_c_qcq() ;
	    }
	    
	tbi.set(1) = *shift_tot(1).get_spectral_va().c - indic *om * shift_tot.get_mp().ya ;
	tbi.set(2) = *shift_tot(2).get_spectral_va().c + indic *om * shift_tot.get_mp().xa ;
	tbi.std_spectral_base() ;
	tbi.set(1).annule_domain(nz-1) ;
	tbi.set(2).annule_domain(nz-1) ;
    }
      
    Vector derive_r (shift_auto.get_mp(), CON, *shift_auto.get_triad()) ;
    for (int i=1 ; i<=3 ; i++) 
	derive_r.set(i) = tbi(i).dsdr() ;
	
    
    // We substract a function in order that Kij is regular
    
    Valeur val_hor (shift_auto.get_mp().get_mg()) ;
    Valeur fonction_radiale (shift_auto.get_mp().get_mg()) ;
    Scalar enleve (shift_auto.get_mp()) ;
  
    double erreur = 0 ;
    for (int comp=1 ; comp<=3 ; comp++) {
	    val_hor.annule_hard() ; 
	    for (int k=0 ; k<np ; k++)
		for (int j=0 ; j<nt ; j++)
		    for (int i=0 ; i<nr ; i++)
		    val_hor.set(1, k, j, i) = derive_r(comp).
			val_grid_point(1, k, j, 0) ;
			     
	    double r_0 = shift_auto.get_mp().val_r (1, -1, 0, 0) ;
	    double r_1 = shift_auto.get_mp().val_r (1, 1, 0, 0) ;
	    
	    fonction_radiale = pow(r_1-shift_auto.get_mp().r, 3.)*
		    (shift_auto.get_mp().r-r_0)/pow(r_1-r_0, 3.) ;
	    fonction_radiale.annule(0) ;
	    fonction_radiale.annule(2, nz-1) ;
	      
	    enleve = fonction_radiale * val_hor ;
	    enleve.set_spectral_va().set_base (shift_auto(comp).
					       get_spectral_va().get_base()) ;
	    
	    if (norme(enleve)(1) != 0)
		shift_auto.set(comp) = shift_auto(comp) - enleve ;
	    if (norme(shift_auto(comp))(1) > 1e-5) {
		Tbl diff (diffrelmax (shift_auto(comp), shift_old(comp))) ;
		if (erreur < diff(1))
		    erreur = diff(1) ;
	    }
    }

    shift_auto.change_triad(shift_auto.get_mp().get_bvect_spher()) ;

    beta_auto = shift_auto ; 

    return erreur ;
}


// Regularisation if only one black hole :
double Single_hor::regularise_one () {

    Vector shift (beta) ;
    
    shift.change_triad(mp.get_bvect_cart()) ;
    // Vector B (without boost and rotation)
    Vector tbi (shift) ;
 
   for (int i=1 ; i<=3 ; i++) {
	tbi.set(i).set_spectral_va().coef_i() ;
	tbi.set(i).set_spectral_va().set_etat_c_qcq() ;
    }
	
    for (int i=1 ; i<=3 ; i++)
	shift(i).get_spectral_va().coef_i() ;
	
    tbi.set(1) = *shift(1).get_spectral_va().c - omega*mp.y ;
    tbi.set(2) = *shift(2).get_spectral_va().c + omega*mp.x ;
    if (shift(3).get_etat() !=  ETATZERO)
	tbi.set(3) = *shift(3).get_spectral_va().c ;
    else 
	tbi.set(3) = 0. ;
    tbi.std_spectral_base() ;
    
    // We only need values at the horizon
    tbi.set(1).annule_domain(mp.get_mg()->get_nzone()-1) ;
    tbi.set(2).annule_domain(mp.get_mg()->get_nzone()-1) ;
      
    Vector derive_r (mp, CON, mp.get_bvect_cart()) ;
    derive_r.set_etat_qcq() ;
    for (int i=1 ; i<=3 ; i++)
	derive_r.set(i) = tbi(i).dsdr() ;
    
    Valeur val_hor (mp.get_mg()) ;
    Valeur fonction_radiale (mp.get_mg()) ;
    Scalar enleve (mp) ;
    
    double erreur = 0 ;
    int np = mp.get_mg()->get_np(1) ;
    int nt = mp.get_mg()->get_nt(1) ;
    int nr = mp.get_mg()->get_nr(1) ;
    
    double r_0 = mp.val_r(1, -1, 0, 0) ;
    double r_1 = mp.val_r(1, 1, 0, 0) ;
    
    for (int comp=1 ; comp<=3 ; comp++) {
	val_hor.annule_hard() ;
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		for (int i=0 ; i<nr ; i++)
		    val_hor.set(1, k, j, i) = derive_r(comp).val_grid_point(1, k, j, 0) ;
	
	fonction_radiale = pow(r_1-mp.r, 3.)* (mp.r-r_0)/pow(r_1-r_0, 3.) ;
	fonction_radiale.annule(0) ;
	fonction_radiale.annule(2, nz-1) ;
    
	enleve = fonction_radiale*val_hor ;
	enleve.set_spectral_va().base = shift(comp).get_spectral_va().base ;
	
	Scalar copie (shift(comp)) ;
	shift.set(comp) = shift(comp)-enleve ;
	shift.std_spectral_base() ;

	assert (shift(comp).check_dzpuis(0)) ;
	
	// Intensity of the correction (if nonzero)
	    Tbl norm (norme(shift(comp))) ;
	    if (norm(1) > 1e-5) {
		Tbl diff (diffrelmax (copie, shift(comp))) ;
		if (erreur<diff(1)) 
		    erreur = diff(1) ;
	    }
	}
    
    shift.change_triad(mp.get_bvect_spher()) ;
    beta = shift ;

    return erreur ;
}
