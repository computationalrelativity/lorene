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


char bhole_glob_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.5  2001/06/28  12:11:19  eric
 * Correction d'un memory leak dans Bhole_binaire::moment_systeme_inf()
 *   (ajout de delete [] bornes).
 *
 * Revision 2.4  2001/05/07  12:24:19  phil
 * correction de calcul de J
 *
 * Revision 2.3  2001/05/07  09:12:09  phil
 * *** empty log message ***
 *
 * Revision 2.2  2001/03/22  10:49:47  phil
 * *** empty log message ***
 *
 * Revision 2.1  2001/03/22  10:44:20  phil
 * *** empty log message ***
 *
 * Revision 2.0  2001/03/01  08:18:13  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */



//standard
#include <stdlib.h>
#include <math.h>
#include <fstream.h>

// Lorene
#include "nbr_spx.h"
#include "map.h"
#include "valeur.h"
#include "cmp.h"
#include "tenseur.h"
#include "bhole.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

double Bhole_binaire::adm_systeme() const {
    Cmp der_un (hole1.psi_auto().dsdr()) ;
    Cmp der_deux (hole2.psi_auto().dsdr()) ;
    
    double masse = hole1.mp.integrale_surface_infini(der_un) + 
	hole2.mp.integrale_surface_infini(der_deux) ;
    
    masse /= -2*M_PI ;
    return masse ;
}

double Bhole_binaire::komar_systeme() const {
    Cmp der_un (hole1.n_auto().dsdr()) ;
    Cmp der_deux (hole2.n_auto().dsdr()) ;
    
    double masse = hole1.mp.integrale_surface_infini(der_un) + 
	hole2.mp.integrale_surface_infini(der_deux) ;
    
    masse /= 4*M_PI ;
    return masse ;
}

double Bhole_binaire::moment_systeme_hor() const {
    
    if (omega == 0)
	return 0 ;
    else {
	// Alignes ou non ?
	double orientation_un = hole1.mp.get_rot_phi() ;
	assert ((orientation_un==0) || (orientation_un==M_PI)) ;
	double orientation_deux = hole2.mp.get_rot_phi() ;
	assert ((orientation_deux==0) || (orientation_deux==M_PI)) ;
	int same_orient = (orientation_un == orientation_deux) ? 1 : -1 ;
	
	// Integrale sur le premiere horizon :
	Cmp xa_un (hole1.mp) ;
	xa_un = hole1.mp.xa ;
	xa_un.std_base_scal() ;
	
	Cmp ya_un (hole1.mp) ;
	ya_un = hole1.mp.ya ;
	ya_un.std_base_scal() ;
	
	Tenseur vecteur_un (hole1.mp, 1, CON, hole1.mp.get_bvect_cart()) ;
	vecteur_un.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_un.set(i) = (-ya_un*hole1.tkij_tot(0, i)+
		xa_un * hole1.tkij_tot(1, i)) ;
	vecteur_un.set_std_base() ;
	vecteur_un.annule(hole1.mp.get_mg()->get_nzone()-1) ;
	vecteur_un.change_triad (hole1.mp.get_bvect_spher()) ;
	
	Cmp integrant_un (pow(hole1.psi_tot(), 6)*vecteur_un(0)) ;
	integrant_un.std_base_scal() ;
	double moment_un = hole1.mp.integrale_surface
	    (integrant_un, hole1.rayon)/8/M_PI ;
	
	//Integrale sur le second horizon :
	Cmp xa_deux (hole2.mp) ;
	xa_deux = hole2.mp.xa ;
	xa_deux.std_base_scal() ;
	
	Cmp ya_deux (hole2.mp) ;
	ya_deux = hole2.mp.ya ;
	ya_deux.std_base_scal() ;
	
	Tenseur vecteur_deux (hole2.mp, 1, CON, hole2.mp.get_bvect_cart()) ;
	vecteur_deux.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_deux.set(i) = -ya_deux*hole2.tkij_tot(0, i)+
		xa_deux * hole2.tkij_tot(1, i) ;
	vecteur_deux.set_std_base() ;
	vecteur_deux.annule(hole2.mp.get_mg()->get_nzone()-1) ;
	vecteur_deux.change_triad (hole2.mp.get_bvect_spher()) ;
	
	Cmp integrant_deux (pow(hole2.psi_tot(), 6)*vecteur_deux(0)) ;
	integrant_deux.std_base_scal() ;
	double moment_deux = hole2.mp.integrale_surface
	    (integrant_deux, hole2.rayon)/8/M_PI ;
	
	return moment_un+same_orient*moment_deux ;
	}
}

double Bhole_binaire::moment_systeme_inf() {
    
    if (omega == 0)
	return 0 ;
    else {
	// Alignes ou non ?
	double orientation_un = hole1.mp.get_rot_phi() ;
	assert ((orientation_un==0) || (orientation_un==M_PI)) ;
	double orientation_deux = hole2.mp.get_rot_phi() ;
	assert ((orientation_deux==0) || (orientation_deux==M_PI)) ;
	int same_orient = (orientation_un == orientation_deux) ? 1 : -1 ;

	// On construit une grille et un mapping auxiliaire :
	int nzones = hole1.mp.get_mg()->get_nzone() ;
	double* bornes = new double [nzones+1] ;
	double courant = (hole1.mp.get_ori_x()-hole2.mp.get_ori_x())+1 ;
	for (int i=nzones-1 ; i>0 ; i--) {
	    bornes[i] = courant ;
	    courant /= 2. ;
	    }
	bornes[0] = 0 ;
	bornes[nzones] = __infinity ;
	
	Map_af mapping (*hole1.mp.get_mg(), bornes) ;
	
	delete [] bornes ; 
	
	// On construit k_total dessus :
	Tenseur_sym k_total (mapping, 2, CON, mapping.get_bvect_cart()) ;
	k_total.set_etat_qcq() ;
	
	Tenseur shift_un (mapping, 1, CON, mapping.get_bvect_cart()) ;
	shift_un.set_etat_qcq() ;
	
	Tenseur shift_deux (mapping, 1, CON, mapping.get_bvect_cart()) ;
	shift_deux.set_etat_qcq() ;
	
	shift_un.set(0).import_asymy(hole1.shift_auto(0)) ;
	shift_un.set(1).import_symy(hole1.shift_auto(1)) ;
	shift_un.set(2).import_asymy(hole1.shift_auto(2)) ;
	
	shift_deux.set(0).import_asymy(same_orient*hole2.shift_auto(0)) ;
	shift_deux.set(1).import_symy(same_orient*hole2.shift_auto(1)) ;
	shift_deux.set(2).import_asymy(hole2.shift_auto(2)) ;
	
	Tenseur shift_tot (shift_un+shift_deux) ;
	shift_tot.set_std_base() ;
	shift_tot.annule(0, nzones-2) ;
	// On enleve les residus
	shift_tot.inc2_dzpuis() ;
	shift_tot.dec2_dzpuis() ;
	
	Tenseur grad (shift_tot.gradient()) ;
	Tenseur trace (grad(0, 0)+grad(1, 1)+grad(2, 2)) ;
	for (int i=0 ; i<3 ; i++) {
	    k_total.set(i, i) = grad(i, i)-trace()/3. ;
	    for (int j=i+1 ; j<3 ; j++)
		k_total.set(i, j) = (grad(i, j)+grad(j, i))/2. ;
	    }	
	
	for (int lig=0 ; lig<3 ; lig++)
	   for (int col=lig ; col<3 ; col++)
		k_total.set(lig, col).mult_r_zec() ;
	
	Tenseur vecteur_un (mapping, 1, CON, mapping.get_bvect_cart()) ;
	vecteur_un.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_un.set(i) = k_total(0, i) ;
	vecteur_un.change_triad (mapping.get_bvect_spher()) ;
	Cmp integrant_un (vecteur_un(0)) ;
	
	Tenseur vecteur_deux (mapping, 1, CON, mapping.get_bvect_cart()) ;
	vecteur_deux.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    vecteur_deux.set(i) = k_total(1, i) ;
	vecteur_deux.change_triad (mapping.get_bvect_spher()) ;
	Cmp integrant_deux (vecteur_deux(0)) ;
	
	// On multiplie par y et x :
	integrant_un.va = integrant_un.va.mult_st() ;
	integrant_un.va = integrant_un.va.mult_sp() ;
	
	integrant_deux.va = integrant_deux.va.mult_st() ;
	integrant_deux.va = integrant_deux.va.mult_cp() ;
	
	double moment = mapping.integrale_surface_infini (-integrant_un+integrant_deux) ;
	
	moment /= 8*M_PI ;
	return moment ;
	}
}

double Bhole_binaire::distance_propre(const int nr) const {
    
    // On determine les rayons coordonnes des points limites de l'integrale :
    double x_un = hole1.mp.get_ori_x() - hole1.rayon ;
    double x_deux = hole2.mp.get_ori_x() + hole2.rayon ;
    
    // Les coefficients du changement de variable :
    double pente = 2./(x_un-x_deux) ;
    double constante = - (x_un+x_deux)/(x_un-x_deux) ;
    
    
    double ksi ; // variable d'integration.
    double xabs ; // x reel.
    double air_un, theta_un, phi_un ; // coordonnee spheriques 1
    double air_deux, theta_deux, phi_deux ; // coordonnee spheriques 2
    
    double* coloc = new double[nr] ;
    double* coef = new double[nr] ;
    int* deg = new int[3] ;
    deg[0] = 1 ; deg[1] = 1 ; deg[2] = nr ;
    
    for (int i=0 ; i<nr ; i++) {
	ksi = -cos (M_PI*i/(nr-1)) ;
	xabs = (ksi-constante)/pente ;
	
	hole1.mp.convert_absolute (xabs, 0, 0, air_un, theta_un, phi_un) ;
	hole2.mp.convert_absolute (xabs, 0, 0, air_deux, theta_deux, phi_deux) ;
	
	coloc[i] = pow (hole1.psi_auto().val_point (air_un, theta_un, phi_un) +
		   hole2.psi_auto().val_point (air_deux, theta_deux, phi_deux), 2.) ;
    }
    
    // On prend les coefficients de la fonction
    cfrcheb(deg, deg, coloc, deg, coef) ;
    
    // On integre
    double* som = new double[nr] ;
    som[0] = 2 ;
    for (int i=2 ; i<nr ; i+=2)
	som[i] = 1./(i+1)-1./(i-1) ;
    for (int i=1 ; i<nr ; i+=2)
	som[i] = 0 ;
    
    double res = 0 ;
    for (int i=0 ; i<nr ; i++)
	res += som[i]*coef[i] ;
	
    res /= pente ;
    
    delete [] deg ;
    delete [] coef ;
    delete [] coloc ;
    delete [] som ;
    
    return res ;
}
