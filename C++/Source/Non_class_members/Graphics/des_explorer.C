/*
 *  Sets of routines producing outputs of different Lorene objects for
 *  visualization by Iris Explorer.
 */

/*
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char des_explorer_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2000/02/11  09:56:34  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header C
#include <string.h>

// Headers Lorene
#include "etoile.h"
#include "binaire.h"
#include "graphique.h"

void des_explorer(const Etoile& star, const char* name) {
    
    // Number of domains to visualize  ---> nzm1 
    int nz = (star.get_mp()).get_mg()->get_nzone() ; 
    int nzm1 = nz - 1 ;
    
    char nom_fich[80] ;
    char nom_l[3] ; 

    for (int l=0; l<nzm1; l++) {
	
	strncpy(nom_fich, name, 80) ; 
	sprintf(nom_l, "%d", l) ; 
	strcat(nom_fich, nom_l) ; 
	strcat(nom_fich, ".d") ; 

	des_explorer(star.get_nbar()(), l, nom_fich) ; 	
	
    }
    
}

void des_explorer(const Binaire& bibi, const char* name) {
    
    char nom_fich[80] ;
    char num_star[2] ; 

    for (int i=1; i<=2; i++) {
	
	sprintf(num_star, "%d", i) ; 

	strncpy(nom_fich, name, 80) ; 
	strcat(nom_fich, "_st") ; 
	strcat(nom_fich, num_star) ; 
	strcat(nom_fich, "_dm") ; 

	des_explorer(bibi(i), nom_fich) ; 	
	
    }
    
}

