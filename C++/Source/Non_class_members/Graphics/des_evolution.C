/*
 *  Plot of time evolution
 *
 *    (see file graphique.h for documentation).
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon. 
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char des_evolution_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/05/11 20:09:47  e_gourgoulhon
 * Corrected bug when j_min != 0.
 * Added version of des_evol for plot on the whole Evolution's time range.
 *
 * Revision 1.1  2004/02/17 22:16:08  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */


// Lorene headers
#include "graphique.h"
#include "evolution.h"

// Plot on  the whole time range
//------------------------------

void des_evol(const Evolution<double>& uu, const char* nomy, 
    const char* title, int ngraph, bool closeit, bool show_time, 
    const char* nomx) {
    
    int jmin = uu.j_min() ; 
    int jmax = uu.j_max() ; 

    des_evol(uu, jmin, jmax, nomy, title, ngraph, closeit,
             show_time, nomx) ; 
}


// Plot within a specified time range
//------------------------------------

void des_evol(const Evolution<double>& uu, int j_min, int j_max, 
    const char* nomy, const char* title, int ngraph, bool closeit, 
    bool show_time, const char* nomx) {

    int npt = j_max - j_min + 1 ; 

    float* uutab = new float[npt] ;	    // Values of uu at the npt points
    float* xtab = new float[npt] ;	    // Values of t at the npt points
    
    for (int j=j_min; j<=j_max; j++) {
	uutab[j-j_min] = uu[j] ; 
    }

    if (show_time) {
        for (int j=j_min; j<=j_max; j++) {
            xtab[j-j_min] = uu.get_time(j) ; 
        }
    }
    else{
        for (int j=j_min; j<=j_max; j++) {
            xtab[j-j_min] = j ; 
        }
    }
        
    if (nomx == 0x0) nomx = (show_time) ? "t" : "j" ;

    if (nomy == 0x0) nomy = "" ;

    if (title == 0x0) title = "" ;
   
    des_profile_mult(uutab, 1, npt, xtab, nomx, nomy, title, 0x0, ngraph,
        closeit) ; 
    
    delete [] uutab ; 
    
}

