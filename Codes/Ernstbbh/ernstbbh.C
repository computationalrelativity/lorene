/*
 * Simple code for solving the Ernst equation for Kerr 
 * boundary data in rotating coordinates with Lorene. 
 * The boundary data are given on  a sphere, the focus is
 * on the treatement of the light cylinder
 *
 * 28.01.03
 *
 */

/*
 *   Copyright (c) 2002  E. Gourgoulhon, C. Klein
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
 *
 * $Header: /cvsroot/Lorene/Codes/Ernst/ernstbbh.C
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>

// Lorene headers
#include "cmp.h"
#include "nbr_spx.h"
#include "graphique.h"
#include "utilitaires.h"






//=============================================================

int main() {

  // Read grid parameters from a file
  // --------------------------------

  ifstream parfile("ernstpar.d") ; 
  char blabla[80] ; 
  int nz, nr, nt, np ;
  double mu, alpha, M, tol;
  int MaxIt;
  
  parfile >> nz ; parfile.getline(blabla, 80) ;  
  parfile >> nr ; parfile.getline(blabla, 80) ; 
  parfile >> nt ; parfile.getline(blabla, 80) ; 
  parfile >> np ; parfile.getline(blabla, 80) ;
  parfile >> mu ; parfile.getline(blabla, 80) ;
  parfile >> MaxIt ; parfile.getline(blabla, 80) ; 
  parfile >> tol ; parfile.getline(blabla, 80) ;
  parfile >> M ; parfile.getline(blabla, 80) ;
  parfile >> alpha ; parfile.getline(blabla, 80) ;
  parfile.close() ;

  const double c = cos(alpha);
// const double c2 = cos(2*alpha);
const double s = sin(alpha);
const double cc = c*c;
const double Mc = M*c;
const double J = M*M*s; // angular momentum
const double Om = 0.5*tan(0.5*alpha)/M;    // angular velocity of the horizon                    

cout << "mass: " << Mc <<"\t" << J <<"\t" << Om  << endl;


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
  const Coord& sint = map.sint ; 
  
  Cmp U(map) ; 
  Cmp V(map) ;
  Cmp Usource(map);
  Cmp Vsource(map);

  // Auxiliary fields for the Kerr solution
  //------------------------

  Cmp F(map); // real part of the Ernst potential
  Cmp B(map); // imaginary part of the Ernst potential
  Cmp A(map); // metric function -g_tphi
  Cmp Rp(map);
  Cmp Rm(map);
  
  
  Rp = sqrt(Mc*Mc + r*r + 2*Mc*z);
  Rm = sqrt(Mc*Mc + r*r - 2*Mc*z);
  Cmp X = 0.5*(Rp+Rm)/Mc;
  Cmp Y = 0.5*(Rp-Rm)/Mc;
  
  
  
  Cmp N = (c*X+1)*(c*X+1)+s*s*Y*Y;
  F = (cc*X*X+s*s*Y*Y-1)/N;
  B = -2*s*Y/N;
  A = 2*M*s*(1-Y*Y)*(1+c*X)/N;


  for (int i=0;i<np;i++) // set at infinity
  {
    for(int j=0; j<nt;j++)
    {
	F.set(nz-1,i,j,nr-1) = 1.0;
	B.set(nz-1,i,j,nr-1) = 0.0;
	A.set(nz-1,i,j,nr-1) = 0.0;
    }
  }
    
  F.std_base_scal();		
  B.std_base_scal();
  B.va.set_base_t(T_COSSIN_CI) ;
  A.std_base_scal();
  

// rotating coordinates and functions vanishing at infinity
  Cmp Frot(map);
  Cmp Brrot(map);
  Cmp R(map);
  Cmp Ct(map);
  Cmp St(map);
  Cmp At(map); //A_theta/r/rho
  Cmp Ft(map); //F_theta/r/rho
  
  
  R = r;
  Ct = cost;
  St = sint;
  Cmp Omrsint(map);
  Omrsint = Om*R*St;
  Omrsint.std_base_scal();
//   Omrsint.mult_rsint();
  Cmp H = (A*A-R*St*R*St)/F+R*St*R*St+2*M*R*St*St+2*M*M*St*St;
  Frot = F-1+2*Om*A+Om*Om*H;
  At = 2*X*(1+c*X)*(1+c*X)*(1+c*X)+s*s*X*(1+c*X)*(1+Y*Y)+s*s*c*Y*Y*(1-Y*Y);
  At = At*2*Mc*M*s*Y/Rp/Rm/N/N/R;   
  Ft = (c-s*s*X)*(1+c*X)*(1+c*X)-c*s*s*Y*Y;
  Ft = Ft*2*Mc*Y/Rp/Rm/N/N/R;  

  
  Cmp H1 = F*F+2*Om*A*F+Om*Om*A*A+Omrsint*Omrsint;
//   At = A;
//   At.srdsdt();     

// At.div_rsint();
 //   arrete();
  Cmp H2 = (H1*A-2*Omrsint*R*St*(F+Om*A))/F;
//   Cmp Ft = F;
//   Ft.srdsdt();
//   Ft.div_rsint();
  Cmp H3 = 2*Om*Ct*(F+Om*A);  
  Brrot = (H1*At-H2*Ft/F-H3)/F+2*Om*Ct;
  
  
  
   for (int i=0;i<np;i++)
  {
    for(int j=0; j<nt;j++)
    {
	Frot.set(nz-1,i,j,nr-1) = 0.0;
 	Brrot.set(nz-1,i,j,nr-1) = 0.0;
    }
  }
  Frot.std_base_scal();		
  Brrot.std_base_scal();  
  Brrot.va.set_base_t(T_COSSIN_CI) ;
//    des_profile(A,1.0001,20.0,M_PI/2,0);
//    des_profile(At,1.0001,20.0,M_PI/2,0);
// des_coupe_z(Brrot,0,nz-1);
// des_coupe_x(Brrot,0,nz-1);


  // Boundary Values
  //-----------------------------------------------------
  
    

  Valeur bcU( mgrid.get_angu());
  Valeur bcV( mgrid.get_angu());
  bcU.annule_hard();
  bcV.annule_hard();  
  
  for (int m=0; m < nt; m++)
    {
      for (int l=0; l < np; l++)
	{
	  bcU.set(1,l,m,0) = Frot(1,l,m,nr-2);
	  bcV.set(1,l,m,0) = Brrot(1,l,m,nr-2);	  
	}
    }
  
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
  
//   for (int iter=0; iter < MaxIt;iter++)
//     {
//       
//       cout << "Iteration " << iter << endl;
      
      Cmp Urot = U+1-Om*Om*St*St*(R*R+2*M*R+2*M*M); //F with correct asymptotic behavior
      Cmp Usource1 = 4*Om*Om*(1+2*M/R*(1-St*St)+M*M/R/R*(2-3*St*St));
      Usource = (U.dsdr()-2*Om*Om*St*St*(M+R)*(M+R))*(U.dsdr()-2*Om*Om*St*St*(M+R)*(M+R))
		+(U.srdsdt()-2*Om*Om*St*Ct*(R+2*M+2*M*M/R))*(U.srdsdt()-2*Om*Om*St*Ct*(R+2*M+2*M*M/R))
		- (V.dsdr()-2*Om*Ct)*(V.dsdr()-2*Om*Ct)
		- (V.srdsdt()+2*Om*St+6*Om*Om*M*M*s*St*St*St/R-4*Om*M*St/R)
		*(V.srdsdt()+2*Om*St+6*Om*Om*M*M*s*St*St*St/R-4*Om*M*St/R);
      Usource = Usource/Urot+Usource1;
	
      Cmp Vsource1 = 8*Om*M*Ct/R/R-24*Om*Om*M*M*s*Ct*St*St/R/R;
      Vsource = (V.dsdr()-2*Om*Ct)*(U.dsdr()-2*Om*Om*St*St*(R+M))
      +(V.srdsdt()+2*Om*St+6*Om*Om*M*M*s*St*St*St/R-4*Om*M*St/R)
      *(U.srdsdt()-2*Om*Om*St*Ct*(R+2*M+2*M*M/R));
      Vsource = Vsource/Urot+Vsource1;
      
      Usource.set(0) = 1.0;
      Vsource.set(0) = 1.0;
      for (int i=0;i<np;i++) // set at infinity
  {
    for(int j=0; j<nt;j++)
    {
	Usource.set(nz-1,i,j,nr-1) = 0.0;
	Vsource.set(nz-1,i,j,nr-1) = 0.0;
    }
  }


      Cmp UV = Usource.poisson_dirichlet(bcU, 0) ;
      UV.set(0) = 1.0;

      Tbl diff =  norme(abs(U-UV));
      cout << diff(1) << "   " << diff(2) << endl;
      U = (1-mu)*U + mu*UV;
      
      UV = Vsource.poisson_neumann(bcV, 0) ;
      UV.set(0) = 0.0;
      V = (1-mu)*V + mu*UV;
      
//     }

//   des_profile(U-F,1.0001,20.0,M_PI/2,0);
//   des_profile(V-B,1.0001,20.0,0,M_PI/2);
  
  
//   des_coupe_x(U, 0., 2, "field U at x=0") ; 
//   des_coupe_y(U, 0., 2, "field U at y=0") ; 
//   des_coupe_z(U, 0., 2, "field U at z=0") ; 

//   des_coupe_x(V, 0., 2, "field V at x=0") ; 
//   des_coupe_y(V, 0., 2, "field V at y=0") ; 
//   des_coupe_z(V, 0., 2, "field V at z=0") ; 

  return EXIT_SUCCESS ; 

}
