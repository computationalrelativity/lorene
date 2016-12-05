/*
 * Draws the profile of a {\tt Scalar} along some radial axis determined by
 *  a fixed value of $(\theta, \phi)$.
 */

/*
 *   Copyright (c) 2004-2005 Eric Gourgoulhon & Philippe Grandclement
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
 * Revision 1.13  2016/12/05 16:18:06  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.12  2014/10/13 08:53:22  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.11  2014/10/06 15:16:05  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.10  2012/01/17 10:35:40  j_penner
 * added point plot
 *
 * Revision 1.9  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
 * Revision 1.8  2005/03/25 19:57:04  e_gourgoulhon
 * Added plot of domain boundaries (new argument draw_bound).
 *
 * Revision 1.7  2004/05/17 19:47:25  e_gourgoulhon
 *  -- Function des_profile_mult(const Scalar**,...): added argument
 *     device.
 *  -- Functions des_meridian: added arguments device and closeit.
 *
 * Revision 1.6  2004/04/06 07:47:29  j_novak
 * Added a #include "string.h"
 *
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
#include <cmath>
#include <cstring>

// Header Lorene
#include "scalar.h"
#include "graphique.h"

#include <vector>
namespace Lorene {
//******************************************************************************
// VERSION SCALAR SANS UNITES 

void des_profile(const Scalar& uu, double r_min, double r_max, 
		     double theta, double phi, const char* nomy, const char* title,
                     bool draw_bound) {
  

    const int npt = 400 ;   // Number of points along the axis
    
    float uutab[npt] ;	    // Value of uu at the npt points
    
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    for (int i=0; i<npt; i++) {
    
	double r = hr * i + r_min ; 
	
	uutab[i] = float(uu.val_point(r, theta, phi)) ; 
	
    }
    
    float xmin = float(r_min) ;
    float xmax = float(r_max)  ;
    
    const char* nomx = "r" ; 
    
    if (title == 0x0) {
	title = "" ;
    }

    if (nomy == 0x0) {
	nomy = "" ;
    }
    
    // Preparations for the drawing of boundaries
    // ------------------------------------------
    const Map& mp = uu.get_mp() ; 
    int nz = mp.get_mg()->get_nzone() ;         
    int l_max = (mp.get_mg()->get_type_r(nz-1) == UNSURR) ? nz-2 : nz-1 ; 
    
    float* xbound = new float[l_max+1] ; 
    int nbound = 0 ; 

    if (draw_bound) {
        const double xi_max = 1. ; 
        for (int l=0; l<=l_max; l++) {
    
            double rb = mp.val_r(l, xi_max, theta, phi) ; 
        
            if ((rb >= r_min) && (rb <= r_max)) {
                xbound[nbound] = float(rb) ; 
                nbound++ ;    
            }
        }
    }
    
    des_profile(uutab, npt, xmin, xmax, nomx, nomy, title, 0x0, 
                nbound, xbound) ; 
    
    delete [] xbound ; 
    
} 

//******************************************************************************

void des_profile(const Scalar& uu, double r_min, double r_max, double scale,
		     double theta, double phi, const char* nomx, const char* nomy, 
                     const char* title, bool draw_bound) {
		

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
    
    // Preparations for the drawing of boundaries
    // ------------------------------------------
    const Map& mp = uu.get_mp() ; 
    int nz = mp.get_mg()->get_nzone() ;         
    int l_max = (mp.get_mg()->get_type_r(nz-1) == UNSURR) ? nz-2 : nz-1 ; 
    
    float* xbound = new float[l_max+1] ; 
    int nbound = 0 ; 

    if (draw_bound) {
        const double xi_max = 1. ; 
        for (int l=0; l<=l_max; l++) {
    
            double rb = mp.val_r(l, xi_max, theta, phi) ; 
        
            if ((rb >= r_min) && (rb <= r_max)) {
                xbound[nbound] = float(rb) ; 
                nbound++ ;    
            }
        }
    }

    // Call to the low level routine
    // -----------------------------
    des_profile(uutab, npt, xmin, xmax, nomx, nomy, title, 0x0, 
                nbound, xbound) ; 
    
    delete [] xbound ; 
    
} 


//******************************************************************************

void des_profile_mult(const Scalar** uu, int nprof, double r_min, double r_max, 
	const double* theta, const double* phi, double radial_scale,
        bool closeit, const char* nomy, const char* title, int ngraph,
        const char* nomx, const int* line_style, const char* device,
        bool draw_bound) {
		
    // Special case of no graphical output:
    if (device != 0x0) {
        if ((device[0] == '/') && (device[1] == 'n')) return ; 
    }

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
            uutab[j*npt+i] = float(vv.val_point(rr[i], theta[j], phi[j])) ; 
        }
    }

    
    float xmin = float(radial_scale * r_min) ;
    float xmax = float(radial_scale * r_max) ;
    
    if (nomx == 0x0) nomx = "r" ;

    if (nomy == 0x0) nomy = "" ;

    if (title == 0x0) title = "" ;

    // Preparations for the drawing of boundaries
    // ------------------------------------------
    
    int nbound_max = 100 * nprof ; 
    float* xbound = new float[nbound_max] ; 
    int nbound = 0 ; 

    if (draw_bound) {
        const double xi_max = 1. ; 
        for (int j=0; j<nprof; j++) {
    
            const Map& mp = uu[j]->get_mp() ; 
            int nz = mp.get_mg()->get_nzone() ;         
            int l_max = (mp.get_mg()->get_type_r(nz-1) == UNSURR) ? nz-2 : nz-1 ; 
        
            for (int l=0; l<=l_max; l++) {
            
                double rb = mp.val_r(l, xi_max, theta[j], phi[j]) ; 
        
                if ((rb >= r_min) && (rb <= r_max)) {
                    xbound[nbound] = float(rb * radial_scale) ; 
                    nbound++ ; 
                    if (nbound > nbound_max-1) {
                        cout << "des_profile_mult : nbound too large !" << endl ; 
                        abort() ; 
                    }   
                }
            }
        }
    }

    // Call to the low level routine
    // -----------------------------
    
    des_profile_mult(uutab, nprof, npt, xmin, xmax, nomx, nomy, title, 
                     line_style, ngraph, closeit, device, nbound, xbound) ; 
                     
      
    delete [] uutab ; 
    delete [] xbound ; 
    
} 

//******************************************************************************

void des_meridian(const Scalar& uu, double r_min, double r_max,
                  const char* nomy, int ngraph, const char* device,
                  bool closeit, bool draw_bound) {

    // Special case of no graphical output:
    if (device != 0x0) {
        if ((device[0] == '/') && (device[1] == 'n')) return ; 
    }

    const Scalar* des[] = {&uu, &uu, &uu, &uu, &uu} ; 
    double phi1[] = {0., 0., 0., 0.25*M_PI, 0.25*M_PI} ; 
    double theta1[] = {0., 0.25*M_PI, 0.5*M_PI, 0., 0.25*M_PI} ;
         
    des_profile_mult(des, 5, r_min, r_max, theta1, phi1, 1., closeit, 
            nomy, 
            "phi=0: th=0, pi/4, pi/2, phi=pi/4: th=0, pi/4",
            ngraph, 0x0, 0x0, device, draw_bound) ;
        
}

//******************************************************************************


void des_meridian(const Sym_tensor& hh, double r_min, double r_max,
                  const char* name, int ngraph0, const char* device,
                  bool closeit) {
    
    // Special case of no graphical output:
    if (device != 0x0) {
        if ((device[0] == '/') && (device[1] == 'n')) return ; 
    }

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
    
                des_meridian(hh(i,j), r_min, r_max, nomy, ngraph0+k, device,
                             closeit) ; 
                k++ ; 
                                
        }
    }              

}


//******************************************************************************
// VERSION SCALAR SANS UNITES 

void des_points(const Scalar& uu, 
		     double theta, double phi, const char* nomy, const char* title,
                     bool draw_bound) {
  
    const Map& mp = uu.get_mp() ; 
    int nz = mp.get_mg()->get_nzone() ;         
    int nt = mp.get_mg()->get_nt(nz-1) ; 
    int np = mp.get_mg()->get_np(nz-1) ;

//    const int npt = *(uu.get_mp().get_mg())->get_nzone() ;   // Number of points along the axis
    

    int npt=0;

    for(int ii = 0; ii<nz; ii++)
        npt += (uu.get_mp().get_mg())->get_nr(ii) ;

    float *uutab = new float[npt] ;	    // define a dynamic array
    float *xtab = new float[npt] ;	    // define a dynamic array
   
    Mtbl r = *(mp.r.c);

    for(int ii = 0; ii<nz; ii++){
	int nr = (uu.get_mp().get_mg())->get_nr(ii) ; 
	for(int ij=0; ij<nr; ij++){
	uutab[ii*nr+ij] = float(uu.val_grid_point(ii,np-1,nt-1,ij)) ; 
	xtab[ii*nr+ij] = float(r(ii,np-1,nt-1,ij)) ; 
	}
    }
    
    float xmin = float(totalmin(r)) ;
    float xmax = float(totalmax(r)) ;
    
    const char* nomx = "r" ; 
    
    if (title == 0x0) {
	title = "" ;
    }

    if (nomy == 0x0) {
	nomy = "" ;
    }
    
    // Preparations for the drawing of boundaries
    // ------------------------------------------
    int l_max = (mp.get_mg()->get_type_r(nz-1) == UNSURR) ? nz-2 : nz-1 ; 
    
    float* xbound = new float[l_max+1] ; 
    int nbound = 0 ; 

    if (draw_bound) {
        const double xi_max = 1. ; 
        for (int l=0; l<=l_max; l++) {
    
            double rb = mp.val_r(l, xi_max, theta, phi) ; 
        
            if ((rb >= xmin) && (rb <= xmax)) {
                xbound[nbound] = float(rb) ; 
                nbound++ ;    
            }
        }
    }
    
    des_profile(uutab, npt, xtab, nomx, nomy, title, 0x0, 
                nbound, xbound) ; 
    
    delete [] xbound ; 
    
} 

//******************************************************************************

void des_points(const Scalar& uu, double scale,
		     double theta, double phi, const char* nomx, const char* nomy, 
                     const char* title, bool draw_bound) {
		
    const Map& mp = uu.get_mp() ; 
    int nz = mp.get_mg()->get_nzone() ;         
    int nt = mp.get_mg()->get_nt(nz-1) ; 
    int np = mp.get_mg()->get_np(nz-1) ;

//    const int npt = *(uu.get_mp().get_mg())->get_nzone() ;   // Number of points along the axis
    

    int npt=0;

    for(int ii = 0; ii<nz; ii++)
        npt += (uu.get_mp().get_mg())->get_nr(ii) ;

    float *uutab = new float[npt] ;	    // define a dynamic array
    float *xtab = new float[npt] ;	    // define a dynamic array
   
    Mtbl r = *(mp.r.c);
   
    for(int ii = 0; ii<nz; ii++){
	int nr = (uu.get_mp().get_mg())->get_nr(ii) ; 
	for(int ij=0; ij<nr; ij++){
	uutab[ii*nr+ij] = float(uu.val_grid_point(ii,np-1,nt-1,ij)) ; 
	xtab[ii*nr+ij] = float(r(ii,np-1,nt-1,ij)) ; 
	}
    }
    
    float xmin = float(totalmin(r) * scale) ;
    float xmax = float(totalmax(r) * scale) ;
    
    
    if (title == 0x0) {
	title = "" ;
    }

    if (nomx == 0x0) {
	nomx = "" ;
    }
    
    if (nomy == 0x0) {
	nomy = "" ;
    }
    
    // Preparations for the drawing of boundaries
    // ------------------------------------------
    int l_max = (mp.get_mg()->get_type_r(nz-1) == UNSURR) ? nz-2 : nz-1 ; 
    
    float* xbound = new float[l_max+1] ; 
    int nbound = 0 ; 

    if (draw_bound) {
        const double xi_max = 1. ; 
        for (int l=0; l<=l_max; l++) {
    
            double rb = mp.val_r(l, xi_max, theta, phi) ; 
        
            if ((rb >= xmin/scale) && (rb <= xmax/scale)) {
                xbound[nbound] = float(rb) ; 
                nbound++ ;    
            }
        }
    }

    // Call to the low level routine
    // -----------------------------
    des_profile(uutab, npt, xtab, nomx, nomy, title, 0x0, 
                nbound, xbound) ; 
    
    delete [] xbound ; 
    
}
}
