/*
 * Draws the profile of a {\tt Scalar} along radial axis directions, with a
 * log-scale in the y-axis.
 */

/*
 *   Copyright (c) 2022 Jerome Novak
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
 * Revision 1.1  2022/02/10 16:33:53  j_novak
 * New drawing functions with logscale in y : des_profile_log, des_meridian_log.
 *
 *
 * $Header$
 *
 */

// Header Lorene
#include "scalar.h"
#include "graphique.h"

namespace Lorene {


  //*****************************************************************************
  
  void des_profile_log(const Scalar& uu, double r_min, double r_max, 
		       double theta, double phi, double pzero, const char* nomy,
		       const char* title, bool draw_bound) {
    
    const int npt = 400 ;   // Number of points along the axis
    
    float uutab[npt] ;	    // Value of uu at the npt points
    
    double hr = (r_max - r_min) / double(npt-1) ;

    for (int i=0; i<npt; i++) {
      
      double r = hr * i + r_min ; 

      double tmp_val = fabs(uu.val_point(r, theta, phi)) ;
      if (tmp_val < pzero) tmp_val = pzero ;
      uutab[i] = float(log10(tmp_val)) ; 
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
                nbound, xbound, true) ; 
    
    delete [] xbound ; 
    
  } 


  //******************************************************************************
  
  void des_prof_mult_log(const Scalar** uu, int nprof, double r_min, double r_max, 
			 const double* theta, const double* phi, double pzero,
			 double radial_scale, bool closeit, const char* nomy,
			 const char* title, int ngraph, const char* nomx,
			 const int* line_style, const char* device,
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
	double tmp_val = fabs(vv.val_point(rr[i], theta[j], phi[j]))  ;
	if (tmp_val < pzero) tmp_val = pzero ;
	uutab[j*npt+i] = float(log10(tmp_val)) ;
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
                     line_style, ngraph, closeit, device, nbound, xbound, true) ; 
    
    
    delete [] uutab ; 
    delete [] xbound ; 
    
  } 

  //******************************************************************************

  void des_meridian_log(const Scalar& uu, double r_min, double r_max,
			const char* nomy, int ngraph, double pzero,
			const char* device,
			bool closeit, bool draw_bound) {
    
    // Special case of no graphical output:
    if (device != 0x0) {
      if ((device[0] == '/') && (device[1] == 'n')) return ; 
    }
    
    const Scalar* des[] = {&uu, &uu, &uu, &uu, &uu} ; 
    double phi1[] = {0., 0., 0., 0.25*M_PI, 0.25*M_PI} ; 
    double theta1[] = {0., 0.25*M_PI, 0.5*M_PI, 0., 0.25*M_PI} ;
    
    des_prof_mult_log(des, 5, r_min, r_max, theta1, phi1, pzero, 1., closeit, 
		     nomy,  "phi=0: th=0, pi/4, pi/2, phi=pi/4: th=0, pi/4",
		     ngraph, 0x0, 0x0, device, draw_bound) ;
        
  }

} // end of namespace Lorene
