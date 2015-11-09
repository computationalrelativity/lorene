/*
 *  Code for reading scalar-hair BH from a file
 *  Scalar-hair BH are the KBHsSH described in Herdeiro & Radu (2015)
 *  CQG, 32, 144001 
 */

/*
 *   Copyright (c) 2015 Frederic Vincent, Eric Gourgoulhon
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

char scalarbh_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2015/11/09 16:01:42  f_vincent
 * Updated scalarBH code
 *
 * Revision 1.2  2015/11/05 17:32:11  f_vincent
 * Updated code for class scalarBH.
 *
 * Revision 1.1  2015/10/22 09:20:02  f_vincent
 * New code scalarbh
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>
using namespace std ;

// Lorene headers
#include "compobj.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "proto.h"
#include "graphique.h"

using namespace Lorene ;

int main() {

  // Parameters of the computation
  // -----------------------------
    
  ifstream fpar("par_scalarbh.d") ;
  if ( !fpar.good() ) {
    cerr << "Problem in opening the file par_scalarbh.d ! " << endl ;
    abort() ;
  }
    
  char file_name[256] ; 
  fpar.getline(file_name, 256) ;
  cout << "File to be read: " << file_name << endl ; 

  double mass ; // M
  fpar >> mass ; fpar.ignore(1000,'\n') ;
    
  int graphic_out ; // flag for graphical outputs
  fpar >> graphic_out ; fpar.ignore(1000,'\n') ; 

  int nr ; // Number of collocation points in r in each domain
  fpar >> nr; fpar.ignore(1000,'\n') ;

  int nt ; // Number of collocation points in theta in each domain
  fpar >> nt; fpar.ignore(1000,'\n') ;

  int np ; // Number of collocation points in phi in each domain
  fpar >> np; fpar.ignore(1000,'\n') ;

  int nz ; // Number of domains
  fpar >> nz ; fpar.ignore(1000,'\n') ;
  int nzm1 = nz - 1 ; // Index of outermost domain

  fpar.ignore(1000,'\n') ; // skip title
  double* r_limits = new double[nz+1];  // inner boundaries of each domain in units of M      
  for (int l=0; l<nz; l++) 
    {
      fpar >> r_limits[l]; 
    }
  r_limits[nz] = __infinity ;
    
  fpar.close();


  // Setup of a multi-domain grid (Lorene class Mg3d)
  // ------------------------------------------------
  
  int symmetry_theta = SYM ; // symmetry with respect to the equatorial plane
  int symmetry_phi = SYM ; // symmetry with respect to phi --> phi + pi
  bool compact = true ; // external domain is compactified

  Mg3d mgrid(nz, nr, nt, np, symmetry_theta, symmetry_phi, compact) ;

  //cout << mgrid << endl ; 

  
  // Setup of an affine mapping : grid --> physical space (Lorene class Map_af)
  // --------------------------------------------------------------------------
  
  Map_af map(mgrid, r_limits) ;

  // Construction of the ScalarBH object:
  // ----------------------------------

  ScalarBH object(map, file_name) ; 

  //object.update_metric() ; // TO BE ADDED WHEN WRITTEN

  cout.precision(15) ; 
  //cout << endl << "******* object ******** " << endl ;  
  //cout << object << endl ; 

  /*cout << "Lapse= " << endl;
  cout << object.get_nn() << endl;

  cout << "F0= " << endl;
  cout << object.get_ff0() << endl;*/

  cout << "betap far= " << object.get_beta()(3).val_point(1e5,M_PI/2.,0.) << endl;
  //cout << "lapse far= " << object.get_nn().val_point(10000,0.157,0.) << endl;
    
  // Drawings    
  if (graphic_out == 1) 
    {
      double r_max = 4;//1.5*map.val_r(nzm1,-1.,0.,0.) ; 
      //des_meridian(object.get_sfield(), 0, r_max, "Phi", 1) ; 

      //des_meridian(object.get_nn(), 0, r_max, "N", 1) ; 
	
      //des_meridian(object.get_gamma().cov()(1,1), 0.58, 0.64, "gamma_11", 2) ; 
      //des_meridian(object.get_gamma().cov()(2,2), 0.58, 0.64, "gamma_22", 3) ; 
      //des_meridian(object.get_gamma().cov()(3,3), 0., r_max, "gamma_33", 4) ; 

      des_meridian(object.get_beta()(3), 0., r_max, "Nphi", 5) ; 
      des_meridian(object.get_beta()(3), 1e10, 1e12, "Nphi", 6) ; 
      //des_meridian(object.get_beta()(3), 1.62, 1.63, "Nphi", 6) ; 
	
      //des_meridian(object.get_kk()(1,3), 0, r_max, "K_(r)(ph)", 7) ; 
      //des_meridian(object.get_kk()(2,3), 0, r_max, "K_(th)(ph)", 8) ; 
      arrete() ; 
    }
    
  //----------------------
  // Output file for GYOTO
  //----------------------
  object.gyoto_data("gyoto_scalarBH.d") ;    
    
  return EXIT_SUCCESS ; 
}













