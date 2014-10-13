/*
 *   Copyright (c) 2000-2001 Philippe Grandclement
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


char cmp_manip_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2014/10/13 08:52:47  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.5  2014/10/06 15:13:03  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.4  2008/08/19 06:41:59  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.3  2003/10/23 09:41:27  p_grandclement
 * small modif of set_val_hor (one can work at the origin now)
 *
 * Revision 1.2  2003/10/03 15:58:44  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2001/05/25  09:29:58  phil
 * ajout de filtre_phi
 *
 * Revision 2.1  2001/02/12  18:08:51  phil
 * ajout de fixe_decroissance
 *
 * Revision 2.0  2000/10/19  09:23:37  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//standard
#include <cstdlib>
#include <cmath>

// Lorene
#include "cmp.h"
#include "proto.h"

/*
 * Annule les n derniers coefficients en r dans la derniere zone
 */
 
namespace Lorene {
void Cmp::filtre (int n) {
    
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO)
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
 
void Cmp::filtre_phi (int n, int nz) {
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO)
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
 * d'un Cmp a val
 * Utile quand on a affaire a des nan0x10000000
 */

void Cmp::set_val_inf (double val) {
    
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO) {
	if (val == 0)
	    return ;
	else
	    annule_hard() ;
    }
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
 * Fixe la valeur d'un Cmp a val, sur la frontiere interne de la coquille zone.
 * Utile quand on a affaire a des nan0x10000000
 */

void Cmp::set_val_hor (double val, int zone) {
    
    assert (etat != ETATNONDEF) ;
    if (etat == ETATZERO) {
	if (val == 0)
	    return ;
	else
	    annule_hard() ;
    }
    assert (zone < mp->get_mg()->get_nzone()) ;
    del_deriv() ;
    
    int nt = mp->get_mg()->get_nt(zone) ;
    int np = mp->get_mg()->get_np(zone) ;
    
    va.coef_i() ;
    va.set_etat_c_qcq() ;
    
    for (int k=0 ; k<np ; k++)
	for (int j=0 ; j<nt ; j++)
	    va.set(zone, k, j, 0) = val ;
}

/*
 * Permet de fixer la decroissance du cmp a l infini en viurant les 
 * termes en 1/r^n
 */
void Cmp::fixe_decroissance (int puis) {
    
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
	
	Cmp courant (*this) ;
	
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
	    
	    courant.mult_r_zec() ;
	}
    }
}
}
