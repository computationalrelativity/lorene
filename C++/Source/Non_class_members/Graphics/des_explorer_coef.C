/*
 *   Copyright (c) 2001 Eric Gourgoulhon
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


char des_explorer_coef_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2001/05/22  13:32:09  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C
#include <strings.h>
#include <math.h>

// Headers Lorene
#include "cmp.h"


void des_explorer_coef(const Cmp& uu, int lz, const char* filename) {
    
    // Calcul des coef
    uu.va.coef() ;
    
    // Extraction de cf_rho dans le tbl t
    Mtbl_cf* mt = uu.va.c_cf ;
    Tbl* t = (mt->t)[lz] ;   // La couche No l0
    
    // Extraction du nombre de points
    int nr = t->get_dim(0) ;
    int nt = t->get_dim(1) ;
    int np = t->get_dim(2) ;	// comprend le plus 2

    // Ouverture du fichier
    int taille = nr*nt*np ;
    float* pr = (float*)(malloc( taille * sizeof(float) ) ) ;

    int i ;
    // etat = ZERO ?
    if (t->get_etat() == ETATZERO) {
	for (i=0 ; i<taille ; i++) {
	    pr[i] = 0 ;
	}
    } else {
	for (i=0 ; i<taille ; i++) {
	    float x = fabs( (t->t)[i] ) ;
	    if (x < 1.e-14) x = 1e-14 ;
	    pr[i] = 14 + log10(x) ;
	}
    }
    
    FILE* fd = fopen( filename, "w" ) ;
    fwrite (&np, sizeof(int), 1, fd);
    fwrite (&nt, sizeof(int), 1, fd);
    fwrite (&nr, sizeof(int), 1, fd);
    
    fwrite (pr, sizeof(float), taille, fd);
    
    // Ecriture des bornes
    float zero = 0 ;
    float npx ;
    npx = np - 1 ;
    fwrite (&zero, sizeof(float), 1, fd) ; 
    fwrite (&npx, sizeof(float), 1, fd) ;
    npx = nt -1 ;
    fwrite (&zero, sizeof(float), 1, fd) ; 
    fwrite (&npx, sizeof(float), 1, fd) ;
    npx = nr -1 ;
    fwrite (&zero, sizeof(float), 1, fd) ; 
    fwrite (&npx, sizeof(float), 1, fd) ;
    
    fclose( fd ) ;
    
    // menage
    free(pr) ;
}
    
