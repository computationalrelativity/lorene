/*
 * Plots the spectral coefficients of a Valeur.
 *
 * (see file graphique.h for the documentation).
 *
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char des_coef_valeur_C[] = "$Header$" ;



/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.5  2000/02/25  10:28:02  eric
 * Suppression des appels a Mtbl_cf::nettoie().
 *
 * Revision 1.4  1999/12/20  14:27:21  eric
 * Amelioration des legendes.
 *
 * Revision 1.3  1999/12/20  10:57:33  eric
 * Ajout des arguments device, newgraph, nxpage et nypage.
 *
 * Revision 1.2  1999/12/10  12:30:44  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/10  12:14:28  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header C
#include <stdlib.h>
#include <string.h>

// Header Lorene
#include "valeur.h"
#include "graphique.h"

			//-------------------------//
			//	xi coefficients	   //
			//-------------------------//

void des_coef_xi(const Valeur& uu, int l, int k, int j, double pzero, 
		 char* nomy, char* title, char* device, 
	         int newgraph, int nxpage, int nypage) {

    assert(uu.get_etat() != ETATNONDEF) ; 	
    uu.coef() ; 
    
    int nr = uu.get_mg()->get_nr(l) ; 
    
    double* cf = new double[nr] ; 
    
    // Are all the coefficients zero ?
    // -------------------------------
    if (uu.get_etat() == ETATZERO) {
	for (int i=0; i<nr; i++) {
	    cf[i] = 0 ; 
	}
    }
    else{
	assert(uu.get_etat() == ETATQCQ) ;
	for (int i=0; i<nr; i++) {
	    cf[i] = (*(uu.c_cf))(l, k, j, i) ; 
	}
    }

    char* nomx = "i" ; 

    char title1[80] ;
    char nomy1[80] ;
    char bslash[2] = {92,  '\0'} ;  // 92 is the ASCII code for the backslash 
				    // character
    char nom_l[3] ; 
    sprintf(nom_l, "%d", l) ; 
    char nom_k[4] ; 
    sprintf(nom_k, "%d", k) ; 
    char nom_j[4] ; 
    sprintf(nom_j, "%d", j) ; 

    if (title == 0x0) {
	strcpy(title1, bslash ) ; 
	strcat(title1, "gc coef. for k=" ) ; 
	strcat(title1, nom_k) ; 
	strcat(title1, ", j=" ) ; 
	strcat(title1, nom_j) ; 
	strcat(title1, " (domain " ) ;
	strcat(title1, nom_l) ; 
	strcat(title1, ")" ) ;  
    }
    else{
	strncpy(title1, title, 80) ; 
    }

    if (nomy == 0x0) {
	strcpy(nomy1, "log| c" ) ;
	strcat(nomy1, bslash ) ; 
	strcat(nomy1, "d" ) ;
	strcat(nomy1, nom_k ) ;
	strcat(nomy1, "," ) ;
	strcat(nomy1, nom_j ) ;
	strcat(nomy1, "," ) ;
	strcat(nomy1, "i" ) ;
	strcat(nomy1, bslash ) ; 
	strcat(nomy1, "u |" ) ;	
    }
    else{
	strncpy(nomy1, nomy, 80) ; 	
    }
    
    des_coef(cf, nr, pzero, nomx, nomy1, title1, device, newgraph, 
	     nxpage, nypage) ;    
    
    delete [] cf ; 
    
} 

			//------------------------------//
			//	theta coefficients	//
			//------------------------------//

void des_coef_theta(const Valeur& uu, int l, int k, int i, double pzero, 
		 char* nomy, char* title, char* device, 
	         int newgraph, int nxpage, int nypage) {

    assert(uu.get_etat() != ETATNONDEF) ; 	
    uu.coef() ; 
    
    int nt = uu.get_mg()->get_nt(l) ; 
    
    double* cf = new double[nt] ; 
    
    // Are all the coefficients zero ?
    // -------------------------------
    if (uu.get_etat() == ETATZERO) {
	for (int j=0; j<nt; j++) {
	    cf[j] = 0 ; 
	}
    }
    else{
	assert(uu.get_etat() == ETATQCQ) ;
	for (int j=0; j<nt; j++) {
	    cf[j] = (*(uu.c_cf))(l, k, j, i) ; 
	}
    }

    char* nomx = "j" ; 

    char title1[80] ;
    char nomy1[80] ;
    char bslash[2] = {92,  '\0'} ;  // 92 is the ASCII code for the backslash 
				    // character
    char nom_l[3] ; 
    sprintf(nom_l, "%d", l) ; 
    char nom_k[4] ; 
    sprintf(nom_k, "%d", k) ; 
    char nom_i[4] ; 
    sprintf(nom_i, "%d", i) ; 

    if (title == 0x0) {
	strcpy(title1, bslash ) ; 
	strcat(title1, "gh coef. for k=" ) ; 
	strcat(title1, nom_k) ; 
	strcat(title1, ", i=" ) ; 
	strcat(title1, nom_i) ; 
	strcat(title1, " (domain " ) ;
	strcat(title1, nom_l) ; 
	strcat(title1, ")" ) ;  
    }
    else{
	strncpy(title1, title, 80) ; 
    }

    if (nomy == 0x0) {
	strcpy(nomy1, "log| c" ) ;
	strcat(nomy1, bslash ) ; 
	strcat(nomy1, "d" ) ;
	strcat(nomy1, nom_k ) ;
	strcat(nomy1, ",j," ) ;
	strcat(nomy1, nom_i ) ;
	strcat(nomy1, bslash ) ; 
	strcat(nomy1, "u |" ) ;	
    }
    else{
	strncpy(nomy1, nomy, 80) ; 	
    }
    
    des_coef(cf, nt, pzero, nomx, nomy1, title1, device, newgraph, 
	     nxpage, nypage) ;    
    
    delete [] cf ; 
    
} 


			//------------------------------//
			//	phi coefficients	//
			//------------------------------//

void des_coef_phi(const Valeur& uu, int l, int j, int i, double pzero, 
		 char* nomy, char* title, char* device, 
	         int newgraph, int nxpage, int nypage) {

    assert(uu.get_etat() != ETATNONDEF) ; 	
    uu.coef() ; 
    
    int np = uu.get_mg()->get_np(l) + 2 ; 
    
    double* cf = new double[np] ; 
    
    // Are all the coefficients zero ?
    // -------------------------------
    if (uu.get_etat() == ETATZERO) {
	for (int k=0; k<np; k++) {
	    cf[k] = 0 ; 
	}
    }
    else{
	assert(uu.get_etat() == ETATQCQ) ;
	for (int k=0; k<np; k++) {
	    cf[k] = (*(uu.c_cf))(l, k, j, i) ; 
	}
    }

    char* nomx = "k" ; 

    char title1[80] ;
    char nomy1[80] ;
    char bslash[2] = {92,  '\0'} ;  // 92 is the ASCII code for the backslash 
				    // character
    char nom_l[3] ; 
    sprintf(nom_l, "%d", l) ; 
    char nom_j[4] ; 
    sprintf(nom_j, "%d", j) ; 
    char nom_i[4] ; 
    sprintf(nom_i, "%d", i) ; 

    if (title == 0x0) {
	strcpy(title1, bslash ) ; 
	strcat(title1, "gf coef. for j=" ) ; 
	strcat(title1, nom_j) ; 
	strcat(title1, ", i=" ) ; 
	strcat(title1, nom_i) ; 
	strcat(title1, " (domain " ) ;
	strcat(title1, nom_l) ; 
	strcat(title1, ")" ) ;  
    }
    else{
	strncpy(title1, title, 80) ; 
    }

    if (nomy == 0x0) {
	strcpy(nomy1, "log| c" ) ;
	strcat(nomy1, bslash ) ; 
	strcat(nomy1, "dk," ) ;
	strcat(nomy1, nom_j ) ;
	strcat(nomy1, "," ) ;
	strcat(nomy1, nom_i ) ;
	strcat(nomy1, bslash ) ; 
	strcat(nomy1, "u |" ) ;	
    }
    else{
	strncpy(nomy1, nomy, 80) ; 	
    }
    
    des_coef(cf, np, pzero, nomx, nomy1, title1, device, newgraph, 
	     nxpage, nypage) ;    
    
    delete [] cf ; 
    
} 
