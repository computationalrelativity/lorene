/*
 *  Plots the coefficients of the functions F_l(theta', phi') and
 *   G_l(theta',  phi') which define a mapping of the class Map_et.
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.5  2016/12/05 16:18:06  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:22  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:16:05  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 2.4  2000/11/20  21:43:57  eric
 * Correction erreur affichage k=3 --> k=4 pour G.
 *
 * Revision 1.3  2000/11/14  15:12:27  eric
 * Traitement du cas np=1
 *
 * Revision 1.2  1999/12/20  14:48:42  eric
 * *** empty log message ***
 *
 * Revision 1.1  1999/12/20  14:27:36  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header C
#include <cstdio>
#include <cstring>

// Headers Lorene
#include "map.h"
#include "graphique.h"

namespace Lorene {
void des_map_et(const Map_et& mp, int lz) {
    
    double pzero = 1.e-14 ; 
    
    char nom_zone[3] ; 
    sprintf(nom_zone, "%d", lz) ; 
    
    char title[50] ; 
    char nomy[50] ;
    char nom_k[4] ; 
    
    char bslash[2] = {92,  '\0'} ;  // 92 is the ASCII code for the backslash 
				    // character

    strcpy(nomy, "log(abs( c" ) ; 
    strcat(nomy, bslash) ; 
    strcat(nomy, "dkj" ) ; 
    strcat(nomy, bslash) ; 
    strcat(nomy, "u ))" ) ;     
     
    const char* device = 0x0 ; 
    int newgraph = 1 ;	// to open the graphic device
    int nxpage = 2 ; 
    int nypage = 2 ;     

    int np = mp.get_mg()->get_np(lz) ; 

    int k ; 
    if ( (lz == 0) && (np > 1) ) {
	k = 2 ;
    }
    else {
	k = 0 ; 
    }
    sprintf(nom_k, "%d", k) ; 

    strcpy(title, " ") ; 
    strcat(title, "Theta coef. of F for k=" ) ;
    strcat(title, nom_k) ; 
    strcat(title, " (domain ") ;
    strcat(title, nom_zone) ; 
    strcat(title, ")") ;  

    des_coef_theta(mp.get_ff(), lz, k, 0, pzero, nomy, title, device, 
	           newgraph, nxpage, nypage) ;

    k = 0 ; 
    sprintf(nom_k, "%d", k) ; 
    strcpy(title, " ") ; 
    strcat(title, "Theta coef. of G for k=" ) ;
    strcat(title, nom_k) ; 
    strcat(title, " (domain ") ;
    strcat(title, nom_zone) ; 
    strcat(title, ")") ;  

    newgraph = 0 ;	// graphic device already opened
    if (np == 1) newgraph = 2 ; 
    des_coef_theta(mp.get_gg(), lz, k, 0, pzero, nomy, title, device, 
	           newgraph, nxpage, nypage) ;
    
    if (np > 1) {

	k = (lz == 0) ? 3 : 2 ; 
        
	sprintf(nom_k, "%d", k) ; 

	strcpy(title, " ") ; 
	strcat(title, "Theta coef. of F for k=" ) ;
	strcat(title, nom_k) ; 
	strcat(title, " (domain ") ;
	strcat(title, nom_zone) ; 
	strcat(title, ")") ;  
	des_coef_theta(mp.get_ff(), lz, k, 0, pzero, nomy, title, device, 
	           newgraph, nxpage, nypage) ;
    
    
	k = (lz == 0) ? 4 : 2 ; 
	
	sprintf(nom_k, "%d", k) ;

	strcpy(title, " ") ; 
	strcat(title, "Theta coef. of G for k=" ) ;
	strcat(title, nom_k) ; 
	strcat(title, " (domain ") ;
	strcat(title, nom_zone) ; 
	strcat(title, ")") ;  

	newgraph = 2 ;	// closes the graphic device
	des_coef_theta(mp.get_gg(), lz, k, 0, pzero, nomy, title, device, 
	           newgraph, nxpage, nypage) ;
    
	int j = 0 ; 

	strcpy(title, " ") ; 
	strcat(title, "Phi coef. of F for j=0 (domain ") ;
	strcat(title, nom_zone) ; 
	strcat(title, ")") ;  

	newgraph = 1 ; 
	nxpage = 2 ; 
	nypage = 1 ; 
	des_coef_phi(mp.get_ff(), lz, j, 0, pzero, nomy, title, device, 
	           newgraph, nxpage, nypage) ;
    
	strcpy(title, " ") ; 
	strcat(title, "Phi coef. of G for j=0 (domain ") ;
	strcat(title, nom_zone) ; 
	strcat(title, ")") ;  

	newgraph = 2 ; 
	des_coef_phi(mp.get_gg(), lz, j, 0, pzero, nomy, title, device, 
	           newgraph, nxpage, nypage) ;
    
    }   
    
}
}
