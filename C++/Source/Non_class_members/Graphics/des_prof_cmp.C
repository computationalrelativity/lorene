/*
 * Draws the profile of a {\tt Cmp} along some radial axis determined by
 *  a fixed value of $(\theta, \phi)$.
 */

/*
 *   Copyright (c) 1999-2004 Eric Gourgoulhon
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


char des_prof_cmp_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.6  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.5  2004/03/25 10:29:25  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2004/02/12 16:20:36  e_gourgoulhon
 * Functions des_profile for Scalar's are now in the new file
 *  des_prof_scalar.C
 *
 * Revision 1.3  2004/02/04 14:28:14  p_grandclement
 * Ajout de la version Scalar de des_profile
 *
 * Revision 1.2  2003/06/03 09:59:35  e_gourgoulhon
 * Added a new des_profile with scale and nomx in arguments
 *
 * Revision 1.1.1.1  2001/11/20 15:19:29  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  1999/12/09  16:38:31  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header C
#include <math.h>


// Header Lorene
#include "cmp.h"
#include "graphique.h"
#include "unites.h"

//******************************************************************************

void des_profile(const Cmp& uu, double r_min, double r_max, 
		     double theta, double phi, const char* nomy, const char* title) {
		
using namespace Unites ;  

    const int npt = 400 ;   // Number of points along the axis
    
    float uutab[npt] ;	    // Value of uu at the npt points
    
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    for (int i=0; i<npt; i++) {
    
	double r = hr * i + r_min ; 
	
	uutab[i] = float(uu.val_point(r, theta, phi)) ; 
	
    }
    
    float xmin = float(r_min / km) ;
    float xmax = float(r_max / km) ;
    
    const char* nomx = "r [km]" ; 
    
    if (title == 0x0) {
	title = "" ;
    }

    if (nomy == 0x0) {
	nomy = "" ;
    }
    
    
    des_profile(uutab, npt, xmin, xmax, nomx, nomy, title) ; 
    
} 

//******************************************************************************

void des_profile(const Cmp& uu, double r_min, double r_max, double scale,
		     double theta, double phi, const char* nomx, const char* nomy, const char* title) {
		

    const int npt = 400 ;   // Number of points along the axis
    
    float uutab[npt] ;	    // Value of uu at the npt points
    
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    for (int i=0; i<npt; i++) {
    
	double r = hr * i + r_min ; 
	
	uutab[i] = float(uu.val_point(r, theta, phi)) ; 
	
    }
    
    float xmin = float(r_min * scale) ;
    float xmax = float(r_max * scale) ;
    
    
    if (title == 0x0) {
	title = "" ;
    }

    if (nomx == 0x0) {
	nomx = "" ;
    }
    
    if (nomy == 0x0) {
	nomy = "" ;
    }
    
    
    des_profile(uutab, npt, xmin, xmax, nomx, nomy, title) ; 
    
} 


