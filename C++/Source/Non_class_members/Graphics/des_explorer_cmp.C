/*
 * Prepares a file for an Explorer visualisation of a Cmp in a given domain
 *
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


char des_explorer_cmp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:29  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.0  2000/02/10  20:21:19  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "cmp.h"

void des_explorer(const Cmp& uu, int lz, const char* filename){
    
    const Map& mp = *(uu.get_mp()) ; 
    
    const Mg3d& mg = *(mp.get_mg()) ; 
    
    assert( mg.get_type_r(lz) != UNSURR ) ; 
    assert( mg.get_type_p() == NONSYM ) ; 
    
    int nr = mg.get_nr(lz) ; 
    int nt = mg.get_nt(lz) ; 
    int np = mg.get_np(lz) ; 
    int np1 = np + 1 ;		// + 1 to close the grid
    
    int ntot = nr*nt*np1 ; 
    
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
    }

    // Last point in phi (2 Pi) = first point (0) (closed grid) : 
    
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

    
    FILE* fich = fopen(filename, "w" ) ;
    fwrite (&nr, sizeof(int), 1, fich);
    fwrite (&nt, sizeof(int), 1, fich);
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
