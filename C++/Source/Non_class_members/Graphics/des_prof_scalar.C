/*
 * Draws the profile of a {\tt Scalar} along some radial axis determined by
 *  a fixed value of $(\theta, \phi)$.
 */

/*
 *   Copyright (c) 2004 Eric Gourgoulhon & Philippe Grandclement
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


char des_prof_scalar_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/04/05 14:42:02  e_gourgoulhon
 * Added functions des_meridian.
 *
 * Revision 1.4  2004/02/17 22:18:00  e_gourgoulhon
 * Changed prototype of des_profile_mult (added radial_scale, theta and
 * phi can now vary from one profile to the other, etc...)
 *
 * Revision 1.3  2004/02/16 13:23:33  e_gourgoulhon
 * Function des_profile_mult: added delete [] uutab at the end.
 *
 * Revision 1.2  2004/02/15 21:57:45  e_gourgoulhon
 * des_profile_mult: changed argument Scalar* to Scalar**.
 *
 * Revision 1.1  2004/02/12 16:21:28  e_gourgoulhon
 * Functions des_profile for Scalar's transfered from file des_prof_cmp.C.
 * Added new function des_profile_mult.
 *
 *
 * $Header$
 *
 */

// Header C
#include <math.h>


// Header Lorene
#include "scalar.h"
#include "graphique.h"

//******************************************************************************


// VERSION SCALAR SANS UNITES 

void des_profile(const Scalar& uu, double r_min, double r_max, 
		     double theta, double phi, char* nomy, char* title) {
  

    const int npt = 400 ;   // Number of points along the axis
    
    float uutab[npt] ;	    // Value of uu at the npt points
    
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    for (int i=0; i<npt; i++) {
    
	double r = hr * i + r_min ; 
	
	uutab[i] = uu.val_point(r, theta, phi) ; 
	
    }
    
    float xmin = r_min ;
    float xmax = r_max  ;
    
    char* nomx = "r" ; 
    
    if (title == 0x0) {
	title = "" ;
    }

    if (nomy == 0x0) {
	nomy = "" ;
    }
    
    
    des_profile(uutab, npt, xmin, xmax, nomx, nomy, title) ; 
    
} 

//******************************************************************************

void des_profile(const Scalar& uu, double r_min, double r_max, double scale,
		     double theta, double phi, char* nomx, char* nomy, char* title) {
		

    const int npt = 400 ;   // Number of points along the axis
    
    float uutab[npt] ;	    // Value of uu at the npt points
    
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    for (int i=0; i<npt; i++) {
    
	double r = hr * i + r_min ; 
	
	uutab[i] = uu.val_point(r, theta, phi) ; 
	
    }
    
    float xmin = r_min * scale ;
    float xmax = r_max * scale ;
    
    
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


//******************************************************************************

void des_profile_mult(const Scalar** uu, int nprof, double r_min, double r_max, 
	const double* theta, const double* phi, double radial_scale,
        bool closeit, const char* nomy, const char* title, int ngraph,
        const char* nomx, const int* line_style) {
		

    const int npt = 400 ;   // Number of points along the axis
    double rr[npt] ; 
    
    float* uutab = new float[npt*nprof] ; // Value of uu at the npt points
    					  // for each of the nprof profiles
    
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    for (int i=0; i<npt; i++) {
        rr[i] = hr * i + r_min ; 
    }
	
    
    for (int j=0; j<nprof; j++) {
	
        const Scalar& vv = *(uu[j]) ; 
		
        for (int i=0; i<npt; i++) {
            uutab[j*npt+i] = vv.val_point(rr[i], theta[j], phi[j]) ; 
        }
    }

    
    float xmin = radial_scale * r_min ;
    float xmax = radial_scale * r_max ;
    
    if (nomx == 0x0) nomx = "r" ;

    if (nomy == 0x0) nomy = "" ;

    if (title == 0x0) title = "" ;
   
    
    des_profile_mult(uutab, nprof, npt, xmin, xmax, nomx, nomy, title, 
                     line_style, ngraph, closeit) ; 
                     
      
    delete [] uutab ; 
    
} 

//******************************************************************************

void des_meridian(const Scalar& uu, double r_min, double r_max,
                  const char* nomy, int ngraph) {

        const Scalar* des[] = {&uu, &uu, &uu, &uu, &uu} ; 
        double phi1[] = {0., 0., 0., 0.25*M_PI, 0.25*M_PI} ; 
        double theta1[] = {0., 0.25*M_PI, 0.5*M_PI, 0., 0.25*M_PI} ;
         
        des_profile_mult(des, 5, r_min, r_max, theta1, phi1, 1., false, 
            nomy, 
            "phi=0: th=0, pi/4, pi/2, phi=pi/4: th=0, pi/4",
            ngraph) ;
        
}

//******************************************************************************


void des_meridian(const Sym_tensor& hh, double r_min, double r_max,
                  const char* name) {
    
    char nomy[80] ;
    
    int k = 0 ; 
    for (int i=1; i<=3; i++) {
        for (int j=i; j<=3; j++) {

                char nom_i[3] ; 
                sprintf(nom_i, "%d", i) ; 
                char nom_j[3] ; 
                sprintf(nom_j, "%d", j) ; 
                strncpy(nomy, name, 40) ; 
                strcat(nomy, "  comp. ") ; 
	            strcat(nomy, nom_i) ; 
	            strcat(nomy, nom_j) ; 
    
                des_meridian(hh(i,j), r_min, r_max, nomy, 50+k) ; 
                k++ ; 
                                
        }
    }              

}


