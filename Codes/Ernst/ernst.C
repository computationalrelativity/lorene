/*
 * Simple code for solving the Ernst equation for Kerr boundary data with Lorene
 *
 * 18.09.2002
 *
 */

/*
 *   Copyright (c) 2002  Jörg Frauendiener
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

char name_of_this_file_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/28 15:48:29  j_frauendiener
 * Multiplication with cos(phi) added in comparison of results
 *
 * Revision 1.2  2002/10/17 07:56:49  j_frauendiener
 * Initial revision
 *
 * $Header$
 *
 */

// C++ headers
#include <iostream.h>
#include <fstream.h>

// C headers
#include <stdlib.h>

// Lorene headers
#include "cmp.h"
#include "nbr_spx.h"
#include "graphique.h"

const double alpha = 1.2;
const double M = 1.0;

const double c = cos(alpha);
const double c2 = cos(2*alpha);
const double s = sin(alpha);
const double cc = c*c;
const double Mc = M*c;

const double mu = 0.5;

//=============================================================

int main() {

  // Read grid parameters from a file
  // --------------------------------

  ifstream parfile("ernstpar.d") ; 
  char blabla[80] ; 
  int nz, nr, nt, np ;
  int MaxIt;
  
  parfile >> nz ; parfile.getline(blabla, 80) ;  
  parfile >> nr ; parfile.getline(blabla, 80) ; 
  parfile >> nt ; parfile.getline(blabla, 80) ; 
  parfile >> np ; parfile.getline(blabla, 80) ;
  parfile >> MaxIt ; parfile.getline(blabla, 80) ;  
  parfile.close() ;


  // Construction of a multi-grid (Mg3d)
  // -----------------------------------
  
  Mg3d mgrid(nz, nr, nt, np, SYM, NONSYM, true) ;

  cout << "Mult_grid : " << mgrid << endl ; 

  // Construction of an affine mapping (Map_af)
  // ------------------------------------------

  double* r_limits = new double[nz+1] ; 
  assert( nz == 3 ) ; 
  r_limits[0] = 0 ; 
  r_limits[1] = 1 ; 
  r_limits[2] = 2 ; 
  r_limits[3] = __infinity ; 
  
  Map_af map(mgrid, r_limits) ; 
  
  // Construction of a scalar field (Cmp)
  // ------------------------------------

  const Coord& z = map.z ; 
  const Coord& r = map.r ; 
  const Coord& cost = map.cost ; 
  //  const Coord& sint = map.sint ; 
  
  Cmp U(map) ; 
  Cmp V(map) ;
  Cmp Usource(map);
  Cmp Vsource(map);

  // Auxiliary fields
  //------------------------

  Cmp Rm(map);
  Cmp Rp(map);
  Cmp F(map);
  Cmp I(map);
  
  Rp = sqrt(Mc*Mc + r*r + 2*Mc*z);
  Rm = sqrt(Mc*Mc + r*r - 2*Mc*z);

  F =  Rp*Rp + Rm*Rm +2*c2*Rp*Rm + 4*M*cc*(Rp + Rm) + 4*Mc*Mc; 
  I = -4.0*M*s*c*(Rp-Rm)/F;
  F = (Rp*Rp + Rm*Rm +2*c2*Rp*Rm - 4*Mc*Mc)/F - 1.0;
  for (int i=0;i<np;i++)
    for(int j=0; j<nt;j++)
      F.set(nz-1,i,j,nr-1) = 0.0;
  F.std_base_scal();		// F contains the exact solution
  

  // Auxiliary fields for boundary data
  //-----------------------------------------------------
  
  Valeur rminus( mgrid.get_angu() ) ; 
  Valeur rplus( mgrid.get_angu() ) ; 

  Valeur N( mgrid.get_angu() ) ; 
  Valeur bcU( mgrid.get_angu() ) ; 
  Valeur bcV( mgrid.get_angu() ) ;// constructs a field on the angular grid
  
  Map_af mpa( *(mgrid.get_angu()), r_limits ) ; 

  rplus  = sqrt(Mc*Mc + mpa.r*mpa.r + 2*Mc*mpa.z);
  rminus = sqrt(Mc*Mc + mpa.r*mpa.r - 2*Mc*mpa.z);

  N = rplus*rplus + rminus*rminus +2*c2*rplus*rminus + 4*M*cc*(rplus + rminus) + 4*Mc*Mc;

  
  // boundary values
  
  bcU = (rplus*rplus + rminus*rminus +2*c2*rplus*rminus - 4*Mc*Mc)/N - 1.0;
  bcU.set(0) = bcU(1) ;
  
  bcV = 4*Mc*s*(rminus-rplus)/N;
  bcV.set(0) = bcV(1) ;

  bcU.std_base_scal() ;
  bcV.std_base_scal() ;// sets standard spectral bases 
  bcV.set_base_t(T_COSSIN_CI) ;
  
    
  // initial values
  U = 1.0;
  U.std_base_scal() ;

  V = cost/pow(r,2);
  
  
  V.std_base_scal() ;
  V.va.set_base_t(T_COSSIN_CI) ;


  // No stopping criterion provided for the moment
  
  for (int iter=0; iter < MaxIt;iter++)
    {
      
      cout << "Iteration " << iter << endl;
      
      Usource = (U.dsdr()*U.dsdr() + U.srdsdt()*U.srdsdt()
		 - V.dsdr()*V.dsdr() - V.srdsdt()*V.srdsdt())/(1.0+U);
      Vsource = 2.0*(V.dsdr()*U.dsdr()+V.srdsdt()*U.srdsdt())/(1.0+U) ;
      Usource.set(0) = 1.0;
      Vsource.set(0) = 1.0;

      Cmp UV = Usource.poisson_dirichlet(bcU, 0) ;
      UV.set(0) = 1.0;

      Tbl diff =  norme(abs(U-UV));
      cout << diff(1) << "   " << diff(2) << endl;
      U = (1-mu)*U + mu*UV;
      
      UV = Vsource.poisson_dirichlet(bcV, 0) ;
      UV.set(0) = 0.0;
      V = (1-mu)*V + mu*UV;
      
    }

  des_profile(U-F,1.0001,20.0,M_PI/2,0);
  des_profile(V-I,1.0001,20.0,0,M_PI/2);
  
//   des_coupe_x(U, 0., 2, "field U at x=0") ; 
//   des_coupe_y(U, 0., 2, "field U at y=0") ; 
//   des_coupe_z(U, 0., 2, "field U at z=0") ; 

//   des_coupe_x(V, 0., 2, "field V at x=0") ; 
//   des_coupe_y(V, 0., 2, "field V at y=0") ; 
//   des_coupe_z(V, 0., 2, "field V at z=0") ; 

  return EXIT_SUCCESS ; 

}
