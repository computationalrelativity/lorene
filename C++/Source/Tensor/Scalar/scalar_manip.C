/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *   
 *   Copyright (c) 2000-2001 Philippe Grandclement (Cmp version)
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


char scalar_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.2  2003/10/08 14:24:09  j_novak
 * replaced mult_r_zec with mult_r_ced
 *
 * Revision 1.1  2003/09/25 09:33:36  j_novak
 * Added methods for integral calculation and various manipulations
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "tensor.h"
#include "proto.h"
#include "utilitaires.h"
/*
 * Annule les n derniers coefficients en r dans la derniere zone
 */
 
void Scalar::filtre (int n) {
    
    assert (etat != ETATNONDEF) ;
    if ( (etat == ETATZERO) || (etat == ETATUN) )
	return ;
    
    int nz = mp->get_mg()->get_nzone() ;
    int np = mp->get_mg()->get_np(nz-1) ;
    int nt = mp->get_mg()->get_nt(nz-1) ;
    int nr = mp->get_mg()->get_nr(nz-1) ;
    
    del_deriv() ;
    
    va.coef() ;
    va.set_etat_cf_qcq() ;
    
    for (int k=0 ; k<np+1 ; k++)
	if (k!=1)
	    for (int j=0 ; j<nt ; j++)
		for (int i=nr-1 ; i>nr-1-n ; i--)
		    va.c_cf->set(nz-1, k, j, i) = 0 ;
}

/*
 * Annule les n derniers coefficients en phi dans zone nz
 */
 
void Scalar::filtre_phi (int n, int nz) {
    assert (etat != ETATNONDEF) ;
    if ( (etat == ETATZERO) || (etat == ETATUN) )
	return ;
    
    del_deriv() ;
    
    va.coef() ;
    va.set_etat_cf_qcq() ;
    int np = mp->get_mg()->get_np(nz) ;
    int nt = mp->get_mg()->get_nt(nz) ;
    int nr = mp->get_mg()->get_nr(nz) ;
    
    for (int k=np+1-n ; k<np+1 ; k++)
	for (int j=0 ; j<nt ; j++)
	    for (int i=0 ; i<nr ; i++)
		va.c_cf->set(nz, k, j, i) = 0 ;
}

/*
 * Fixe la valeur a l'infini (si la derniere zone est compactifiee) 
 * d'un Scalar a val
 * Utile quand on a affaire a des nan0x10000000
 */

void Scalar::set_val_inf (double val) {
    
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO)
	if (val == 0)
	    return ;
	else
	    annule_hard() ;
    
    del_deriv() ;
    
    int nz = mp->get_mg()->get_nzone() ;
    
    // On verifie la compactification
    assert (mp->get_mg()->get_type_r(nz-1) == UNSURR) ;
    
    int nr = mp->get_mg()->get_nr(nz-1) ;
    int nt = mp->get_mg()->get_nt(nz-1) ;
    int np = mp->get_mg()->get_np(nz-1) ;
    
    va.coef_i() ;
    va.set_etat_c_qcq() ;
    
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    va.set(nz-1, k, j, nr-1) = val ;
}

/*
 * Fixe la valeur d'un Scalar a val, sur la frontiere interne de la coquille zone.
 * Utile quand on a affaire a des nan0x10000000
 */

void Scalar::set_val_hor (double val, int zone) {
    
    assert (etat != ETATNONDEF) ;
     if (etat == ETATZERO)
	if (val == 0)
	    return ;
	else
	    annule_hard() ;
    
    assert ((zone>0) && (zone < mp->get_mg()->get_nzone())) ;
    del_deriv() ;
    
    int nt = mp->get_mg()->get_nt(zone) ;
    int np = mp->get_mg()->get_np(zone) ;
    
    va.coef_i() ;
    va.set_etat_c_qcq() ;
    
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    va.set(1, k, j, 0) = val ;
}

/*
 * Permet de fixer la decroissance du cmp a l infini en viurant les 
 * termes en 1/r^n
 */
void Scalar::fixe_decroissance (int puis) {
    
    if (puis<dzpuis)
	return ;
    else {
	
	int nbre = puis-dzpuis ;
	
	// le confort avant tout ! (c'est bien le confort ...)
	int nz = mp->get_mg()->get_nzone() ;
	int np = mp->get_mg()->get_np(nz-1) ;
	int nt = mp->get_mg()->get_nt(nz-1) ;
	int nr = mp->get_mg()->get_nr(nz-1) ;
	
	const Map_af* map  = dynamic_cast<const Map_af*>(mp) ;
	if (map == 0x0) {
	    cout << "Le mapping doit etre affine" << endl ;
	    abort() ;
	}
	
	double alpha = map->get_alpha()[nz-1] ;
	
	Scalar courant (*this) ;
	
	va.coef() ;
	va.set_etat_cf_qcq() ;
	
	for (int conte=0 ; conte<nbre ; conte++) {
	    
	    int base_r = courant.va.base.get_base_r(nz-1) ;
	    
	    courant.va.coef() ;
	    
	    // On calcul les coefficients de 1/r^conte
	    double* coloc = new double [nr] ;
	    int * deg = new int[3] ;
	    deg[0] = 1 ; 
	    deg[1] = 1 ;
	    deg[2] = nr ;
		    
	    for (int i=0 ; i<nr ; i++)
		coloc[i] =pow(alpha, double(conte))*
		    pow(-1-cos(M_PI*i/(nr-1)), double(conte)) ;
		    
	    cfrcheb(deg, deg, coloc, deg, coloc) ;
	    
	    for (int k=0 ; k<np+1 ; k++)
		if (k != 1)
		for (int j=0 ; j<nt ; j++) {
		    
		    // On doit determiner le coefficient du truc courant :
		    double* coef = new double [nr] ;
		    double* auxi = new double[1] ;
		    for (int i=0 ; i<nr ; i++)
			coef[i] = (*courant.va.c_cf)(nz-1, k, j, i) ;
		    switch (base_r) {
			case R_CHEBU :
			som_r_chebu (coef, nr, 1, 1, 1, auxi) ;
			break ;
		    default :
			som_r_pas_prevu (coef, nr, 1, 1, 1, auxi) ;
			break ;
		    }
		    
		    // On modifie le cmp courant :
		    courant.va.coef() ;
		    courant.va.set_etat_cf_qcq() ;
		    courant.va.c_cf->set(nz-1, k, j, 0) -= *auxi ;  
			
		    for (int i=0 ; i<nr ; i++)
		    	this->va.c_cf->set(nz-1, k, j, i) -= *auxi * coloc[i] ;

			  
		    delete [] coef ;
		    delete [] auxi ;
		}
	    delete [] coloc ;
	    delete [] deg ;
	    
	    courant.mult_r_ced() ;
	}
    }
}
