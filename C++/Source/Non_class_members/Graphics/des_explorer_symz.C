/*
 * Prepares a file for an Explorer visualisation with both hemispheres
 *
 */

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


char des_explorer_symz_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2001/03/07  10:47:54  eric
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


		    //----------------------------//
		    //	    Cmp version		  //
		    //----------------------------//
		    

void des_explorer_symz(const Cmp& uu, int lz, const char* filename){
    
    const Map& mp = *(uu.get_mp()) ; 
    
    const Mg3d& mg = *(mp.get_mg()) ; 
    
    assert( mg.get_type_r(lz) != UNSURR ) ; 
    assert( mg.get_type_p() == NONSYM ) ; 
    
    int nr = mg.get_nr(lz) ; 
    int nt = mg.get_nt(lz) ; 
    int nt2 = 2*nt - 1  ;	// theta in [0,pi] (full range)
    int np = mg.get_np(lz) ; 
    int np1 = np + 1 ;		// + 1 to close the grid
    
    int ntot = nr*nt2*np1 ; 
    
    // Absolute Cartesian coordinates
    Mtbl xa = mp.xa ; 
    Mtbl ya = mp.ya ; 
    Mtbl za = mp.za ; 

    float* xx = new float[ntot] ; 
    float* yy = new float[ntot] ; 
    float* zz = new float[ntot] ; 
    
    float* fuu = new float[ntot] ; 
    
    float* pxx = xx ; 
    float* pyy = yy ; 
    float* pzz = zz ; 
    float* pfuu = fuu ; 

    for (int k=0; k<np; k++) {
	// theta in [0,pi/2] : 
	// ------------------
	for (int j=0; j<nt; j++) {
	    for (int i=0; i<nr; i++) {
		*pxx = xa(lz, k, j, i) ; 
		*pyy = ya(lz, k, j, i) ; 
		*pzz = za(lz, k, j, i) ; 
		*pfuu = uu(lz, k, j, i) ;
		pxx++ ; 
		pyy++ ; 
		pzz++ ; 
		pfuu++ ; 
	    }
	}

	// theta in ]pi/2,pi] : 
	// ------------------
	for (int j=nt; j<nt2; j++) {
	    int j_sym = nt2 - 1 - j ;	// point symmetric / plane z=0
	    for (int i=0; i<nr; i++) {
		*pxx = xa(lz, k, j_sym, i) ; 
		*pyy = ya(lz, k, j_sym, i) ; 
		*pzz = - za(lz, k, j_sym, i) ; 
		*pfuu = uu(lz, k, j_sym, i) ;
		pxx++ ; 
		pyy++ ; 
		pzz++ ; 
		pfuu++ ; 
	    }
	}


    }

    // Last point in phi (2 Pi) = first point (0) (closed grid) : 
    
    // theta in [0,pi/2] : 
    // ------------------
    for (int j=0; j<nt; j++) {
	for (int i=0; i<nr; i++) {
	    *pxx = xa(lz, 0, j, i) ; 
	    *pyy = ya(lz, 0, j, i) ; 
	    *pzz = za(lz, 0, j, i) ; 
	    *pfuu = uu(lz, 0, j, i) ;
	    pxx++ ; 
	    pyy++ ; 
	    pzz++ ; 
	    pfuu++ ; 
	}
    }
    
    // theta in ]pi/2,pi] : 
    // ------------------
    for (int j=nt; j<nt2; j++) {
	int j_sym = nt2 - 1 - j ;	// point symmetric / plane z=0
	for (int i=0; i<nr; i++) {
	    *pxx = xa(lz, 0, j_sym, i) ; 
	    *pyy = ya(lz, 0, j_sym, i) ; 
	    *pzz = - za(lz, 0, j_sym, i) ; 
	    *pfuu = uu(lz, 0, j_sym, i) ;
	    pxx++ ; 
	    pyy++ ; 
	    pzz++ ; 
	    pfuu++ ; 
	}
    }

    
    FILE* fich = fopen(filename, "w" ) ;
    fwrite (&nr, sizeof(int), 1, fich);
    fwrite (&nt2, sizeof(int), 1, fich);
    fwrite (&np1, sizeof(int), 1, fich);
    for (int i=0 ; i<ntot ; i++ ) {
	fwrite (&(xx[i]), sizeof(float), 1, fich);
	fwrite (&(yy[i]), sizeof(float), 1, fich);
	fwrite (&(zz[i]), sizeof(float), 1, fich);
    }
    fwrite (fuu, sizeof(float), ntot, fich);
    fclose(fich) ;
   
    
    delete [] xx ; 
    delete [] yy ; 
    delete [] zz ; 
    delete [] fuu ; 
    
}

		    //----------------------------//
		    //	    Etoile version	  //
		    //----------------------------//
		    

void des_explorer_symz(const Etoile& star, const char* name) {
    
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

	des_explorer_symz(star.get_nbar()(), l, nom_fich) ; 	
	
    }
    
}

		    //----------------------------//
		    //	    Binaire version	  //
		    //----------------------------//
		    

void des_explorer_symz(const Binaire& bibi, const char* name) {
    
    char nom_fich[80] ;
    char num_star[2] ; 

    for (int i=1; i<=2; i++) {
	
	sprintf(num_star, "%d", i) ; 

	strncpy(nom_fich, name, 80) ; 
	strcat(nom_fich, "_st") ; 
	strcat(nom_fich, num_star) ; 
	strcat(nom_fich, "_dm") ; 

	des_explorer_symz(bibi(i), nom_fich) ; 	
	
    }
    
}
