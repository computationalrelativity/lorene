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


 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:17:49  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:52:48  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2003/10/03 15:58:45  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2001/10/10  13:53:27  eric
 * Modif Joachim: sqrt(2) --> sqrt(double(2))
 *
 * Revision 2.1  2001/04/02  12:16:39  phil
 * *** empty log message ***
 *
 * Revision 2.0  2001/03/30  13:37:32  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */



//standard
#include <cstdlib>
#include <cmath>

// LORENE
#include "matrice.h"
#include "cmp.h"
#include "proto.h"


// Calcul des Cnp
namespace Lorene {
int cnp (int n, int p) {
    
    assert (p<=n) ;
    
    if ((p==0) || (p==n))
	return 1 ;
    else {
    int fact_un = 1 ;
    for (int conte=n ; conte >n-p ; conte --)
	fact_un *= conte ;
    
    int fact_deux = 1 ;
    for (int conte = 1 ; conte<p+1 ; conte++)
	fact_deux *= conte ;
    
    return int(fact_un/fact_deux) ;
    }
}

// Fait le raccord dans la zec ...
// Suppose (pour le moment, le meme nbre de points sur les angles ...)
// et que la zone precedente est une coquille

void Cmp::raccord_externe (int power, int nbre, int lmax) {

     va.coef() ;
     va.ylm() ;
     
     Base_val base_devel (va.base) ;
     int base_r, m_quant, l_quant ;

     // Confort :
     int zone = mp->get_mg()->get_nzone()-2 ;
     int nt = mp->get_mg()->get_nt(zone) ; 
     int np = mp->get_mg()->get_np(zone) ;
     int nr = mp->get_mg()->get_nr(zone) ;
     
     // Le mapping doit etre affine :
    const Map_af* map  = dynamic_cast<const Map_af*>(mp) ;
    if (map == 0x0) {
	cout << "Le mapping doit etre affine" << endl ;
	abort() ;
    }
     
     // Mappinhg en r
     double alpha = map->get_alpha()[zone] ;
     double beta = map->get_beta()[zone] ;
     
     // Mapping en 1/r 
     double new_alpha = -alpha/(beta*beta-alpha*alpha) ;
     double new_beta = beta/(beta*beta-alpha*alpha) ;
     
    // Mapping dans la zec :
    double alpha_zec = map->get_alpha()[zone+1] ;
    
    // Maintenant on construit les matrices de passage :
    // Celle de ksi a T
    Matrice tksi (nbre, nbre) ;
    tksi.set_etat_qcq() ;
    
    // Premier polynome
    tksi.set(0, 0) = sqrt(double(2)) ;
    for (int i=1 ; i<nbre ; i++)
	tksi.set(0, i) = 0 ;
    
    //Second polynome
    tksi.set(1, 0) = 0 ;
    tksi.set(1, 1) = sqrt(double(2)) ;
    for (int i=2 ; i<nbre ; i++)
	tksi.set(1, i) = 0 ;
	
    // On recurre :
    for (int lig=2 ; lig<nbre ; lig++) {
	tksi.set(lig, 0) = -tksi(lig-2, 0) ;
	for (int col=1 ; col<nbre ; col++)
	    tksi.set(lig, col) = 2*tksi(lig-1, col-1)-tksi(lig-2, col) ;
	}
    
    // Celle de u/new_alpha a ksi :
    Matrice ksiu (nbre, nbre) ;
    ksiu.set_etat_qcq() ;
    
    for (int lig=0 ; lig<nbre ; lig++) {
	for (int col=0 ; col<=lig ; col++)
	    ksiu.set(lig, col) = cnp(lig, col)*
		pow(-new_beta/new_alpha, lig-col) ;
	for (int col = lig+1 ; col<nbre ; col++)
	    ksiu.set(lig, col) = 0 ;
	}
    
    // La matrice totale :
    Matrice tu (nbre, nbre) ;
    tu.set_etat_qcq() ;
    double somme ;
    for (int lig=0 ; lig<nbre ; lig++)
	for (int col=0 ; col<nbre ; col++) {
	    somme = 0 ;
	    for (int m=0 ; m<nbre ; m++)
		somme += tksi(lig, m)*ksiu(m, col) ;
	    tu.set(lig, col) = somme ;
	}
     
    // On calcul les coefficients de u^n dans la zec
    Tbl coef_u (nbre+lmax, nr) ;
    coef_u.set_etat_qcq() ;
    int* dege = new int [3] ;
    dege[0] = 1 ; dege[1] = 1 ; dege[2] = nr ;
    double* ti = new double [nr] ;
    
    for (int puiss=0 ; puiss<nbre+lmax ; puiss++) {
	for (int i=0 ; i<nr ; i++)
	    ti[i] = pow(-cos(M_PI*i/(nr-1))-1, puiss) ;
	cfrcheb (dege, dege, ti, dege, ti) ;
	for (int i=0 ; i<nr ; i++)
	    coef_u.set(puiss, i) = ti[i] ;
    }
     
     // Avant d entrer dans la boucle :
     dege[2] = nbre ;
     double *coloc = new double[nbre] ;
     double *auxi = new double [1] ;
     
     Tbl coef_zec (np+2, nt,  nr) ;
     coef_zec.annule_hard() ;
     
     // Boucle sur les harmoniques :
     
     for (int k=0 ; k<np+2 ; k++)
	for (int j=0 ; j<nt ; j++)
	    if (nullite_plm (j, nt, k, np, base_devel)==1) {
	    donne_lm (zone+2, zone+1, j, k, base_devel, m_quant, 
		l_quant, base_r) ;
	    if (l_quant <= lmax) {
		
    // On bosse :
    // On recupere les valeus aux points de colocation en 1/r :
    double ksi, air ;
    for (int i=0 ; i<nbre ; i++) {
	ksi = -cos(M_PI*i/(nbre-1)) ;
	air = 1./(new_alpha*ksi+new_beta) ;
	ksi = (air-beta)/alpha ;
	for (int m=0 ; m<nr ; m++)
	    ti[m] = (*va.c_cf)(zone, k, j, m) ;
	som_r_cheb (ti, nr, 1, 1, ksi, auxi) ;
	coloc[i] = auxi[0]/ 
	    pow (-new_alpha*cos(M_PI*i/(nbre-1))+new_beta, power+l_quant);
	}
	
    cfrcheb (dege, dege, coloc, dege, coloc) ;
    
    Tbl expansion (nbre) ;
    expansion.set_etat_qcq() ;
    for (int i=0 ; i<nbre ; i++) {
	somme = 0 ;
	for (int m=0 ; m<nbre ; m++)
	    somme += coloc[m]*tu(m, i) ;
	expansion.set(i) = somme ;
    }
    
    for (int i=0 ; i<nr ; i++) {
	somme = 0 ;
	for (int m=0 ; m<nbre ; m++)
	    somme += coef_u(m+l_quant, i)*expansion(m)*
			pow(alpha_zec, m+l_quant)/
			pow(new_alpha, m) ;
	coef_zec.set(k, j, i) = somme ;
    }
    }
    }
    
    va.set_etat_cf_qcq() ;
    va.c_cf->set_etat_qcq() ;
    va.c_cf->t[zone+1]->set_etat_qcq() ;
    
    for (int k=0 ; k<np+2 ; k++)
	for (int j=0 ; j<nt ; j++)
	    for (int i=0 ; i<nr ; i++)
		va.c_cf->set(zone+1, k, j, i) = coef_zec(k, j, i) ;
    
    set_dzpuis(power) ;
    va.ylm_i() ;
    
    delete[] auxi ;
    delete [] dege ;
    delete [] ti ;
    delete [] coloc ;
}
}
