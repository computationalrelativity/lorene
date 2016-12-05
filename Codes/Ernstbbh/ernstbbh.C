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

 

/*
 *
 * $Header: /cvsroot/Lorene/Codes/Ernst/ernstbbh.C
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <cstdlib>

// Lorene headers
#include "cmp.h"
#include "nbr_spx.h"
#include "graphique.h"
#include "utilitaires.h"






//=============================================================

using namespace Lorene ;

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
  
  Mg3d mgrid(nz, nr, nt, np, SYM, SYM, true) ;

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
//   const Coord& cost = map.cost ; 
//   const Coord& sint = map.sint ; 
  
  Cmp U(map) ; 
  Cmp V(map) ;
  Cmp Usource(map); // RHS for the real part of the Ernst equation
  Cmp Vsource(map); // RHS for the imaginary part of the Ernst equation

  // Auxiliary fields for the Kerr solution
  //------------------------

  Cmp F(map); // real part of the Ernst potential
  Cmp B(map); // imaginary part of the Ernst potential
  Cmp A(map); // metric function -g_tphi
  Cmp Rp(map);
  Cmp Rm(map);
  
  
  Rp = sqrt(Mc*Mc + r*r + 2*Mc*z);
  Rm = sqrt(Mc*Mc + r*r - 2*Mc*z);
  Cmp X = 0.5*(Rp+Rm)/Mc;//parabolic coordinates
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
   
  Cmp Un(map); // set F=1 in the nucleus to avoid problems due to the ergosphere
  Un = 1;
  Un.annule(1,2);
  F.annule(0);
  F = F +Un;
  F.std_base_scal();		
  B.std_base_scal();
//   B.va.set_base_t(T_COSSIN_CI) ;
  B.va.set_base_t(T_COS_I);
  A.std_base_scal();
  

// rotating coordinates and functions vanishing at infinity
  Cmp Frot(map); // F in corotating coordinates
  Cmp Brrot(map); // B_r in corotating coordinates
  Cmp Btrot(map); // B_theta in corotating coordinates
  Cmp At(map); //A_theta/r/sint
  Cmp Ft(map); //F_theta/r/sint
  Cmp Ar(map); //A_r
  Cmp Fr(map); //F_r
  Cmp Omsint(map);// Omega*sin(theta)
  Omsint = Om;
  Omsint.std_base_scal();
  Omsint.va = Omsint.va.mult_st();

  
  Cmp Ftemp = (F+Om*A);
  Cmp Help2 = Ftemp;
  Ftemp.div_r();
  Frot = (Ftemp*Ftemp-Omsint*Omsint)/F;
  Cmp Rminus(map);// Rminus = 1/r
  Rminus = 1;
  Rminus.std_base_scal();
  Rminus.div_r();
  
// Ftilde:  F in corotating coordinates minus diverging and constant terms
  Cmp Ftilde = Frot - Rminus*Rminus+Omsint*Omsint*(1+2*M*Rminus+2*M*M*Rminus*Rminus); 
  Ftilde.std_base_scal();
  Ftilde.mult_r();
  Ftilde.mult_r();
  Ftilde.inc2_dzpuis();
  
  
  
  At.set_etat_qcq();
  Ft.set_etat_qcq();
  Ar.set_etat_qcq();
  Fr.set_etat_qcq();
  At.va = A.va.dsdt().ssint();
  Ft.va = F.va.dsdt().ssint();
  Ar = A.dsdr();
  Fr = F.dsdr();
  Ar.dec2_dzpuis();
  Fr.dec2_dzpuis();
  Cmp Help1 = Ftemp*Ftemp+Omsint*Omsint;
  Cmp Omsint2(map);// Omega*sin(theta)^2
  Omsint2.set_etat_qcq();
  Omsint2.va = Omsint.va.mult_st();
  Cmp Omcost(map);// Omega*cos(theta)
  Omcost = Om;
  Omcost.std_base_scal();
  Omcost.va = Omcost.va.mult_ct();
  Brrot = Help1*At+(2*Omsint2*Help2-Help1*A)*Ft/F-2*Omcost*Help2;
  Brrot = Brrot/F;
  Cmp Brtilde = Brrot+2*Omcost; //Brot_r minus diverging and constant terms
  Btrot = Help1*Ar+(2*Omsint2*Help2-Help1*A)*Fr/F;
  Btrot.va = Btrot.va.ssint();
  Btrot.mult_r();
  Btrot.inc_dzpuis();
  Btrot = -Btrot+2*Omsint*Help2;
  Btrot = Btrot/F;
  // Brot_theta minus diverging and constant terms
  Cmp Bttilde = Btrot-2*Omsint-6*M*M*s*Omsint*Omsint*Omsint/Om*Rminus+4*M*Omsint*Rminus;
  
  

  
  
  


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
	  bcU.set(1,l,m,0) = Ftilde(1,l,m,0);
	  bcV.set(1,l,m,0) = Brtilde(1,l,m,0);	  
	}
    }
 
  bcU.set_base(Ftilde.va.base);
  bcV.set_base(Brtilde.va.base);

  // initial values
//   U = 1;
  U = Ftilde;
  U.std_base_scal() ;

  V = U;
  V.annule_hard();
  V.std_base_scal() ;
//   V.va.set_base_t(T_COSSIN_CI) ;
  V.va.set_base_t(T_COS_I);

  // No stopping criterion provided for the moment
  
//   for (int iter=0; iter < MaxIt;iter++)
//     {
//       
//       cout << "Iteration " << iter << endl;
      
      Cmp Ur = U;
      Ur.div_r();
      Ur.div_r();
      Cmp Urot = Ur+Rminus*Rminus-Omsint*Omsint*(1+2*M*Rminus+2*M*M*Rminus*Rminus); //F with correct asymptotic behavior
      
      Cmp Usource1 = 4*Omsint*Omsint*(2*M*Rminus+3*M*M*Rminus*Rminus)-4*Om*Om*(1+2*M*Rminus+2*M*M*Rminus*Rminus);
      Cmp LHS = Frot*(Ftilde.laplacien(0)+Usource1); // test Kerr
     
      Cmp Urr = U.dsdr();
      Urr.div_r();
      Cmp Utr = U.srdsdt();
      Utr.div_r();
      Cmp Vsr = V.dsdr();
      Cmp Vst = V.srdsdt();
      Vsr.dec2_dzpuis();
      Vst.dec2_dzpuis();
      Urr.dec2_dzpuis();
      Utr.dec2_dzpuis();
      
//       Test of Kerr
      Vsr = Brtilde;
      Vst = Bttilde;
      
      
      Cmp Ustemp = -(Vsr-2*Omcost)*(Vsr-2*Omcost)
		- (Vst+2*Omsint+6*Omsint*Omsint*Omsint*M*M*s*Rminus/Om-4*Omsint*M*Rminus)
		*(Vst+2*Omsint+6*Omsint*Omsint*Omsint*M*M*s*Rminus/Om-4*Omsint*M*Rminus);
      Ustemp.div_r();
      Ustemp.div_r();
      Usource = (Urr-2*Omsint*Omsint*(1+M*Rminus))*(Urr-2*Omsint*Omsint*(1+M*Rminus))
		+(Utr-2*Omsint*Omcost*(1+2*M*Rminus+2*M*M*Rminus*Rminus))*(Utr-2*Omsint*Omcost*(1+2*M*Rminus+2*M*M*Rminus*Rminus))
 		+Ustemp;
      Cmp test = Usource;
      
//       cout << test << endl;
//       arrete();
      Usource = Usource/Urot-Usource1;
       
   
	
      Cmp Vsource1 = (-8*Omcost*M+24*Omcost*Omsint*Omsint*M*M*s/Om)*Rminus*Rminus;
      Vsource = (Vsr-2*Omcost)*(Urr-2*Omsint*Omsint*(1+M*Rminus))
      +(Vst+2*Omsint+6*Omsint*Omsint*Omsint*M*M*s*Rminus/Om-4*Omsint*M*Rminus)
      *(Utr-2*Omsint*Omcost*(1+2*M*Rminus+2*M*M*Rminus*Rminus));
      Vsource.div_r();
      Vsource = 2*Vsource/Urot-Vsource1;
      
      Usource.set(0) = 1.0;
      Vsource.set(0) = 1.0;
      test.set(0) = 1.0;
      LHS.set(0) = 1.0;
//       for (int i=0;i<np;i++) // set at infinity
//   {
//     for(int j=0; j<nt;j++)
//     {
// 	test.set(nz-1,i,j,nr-1) = 0.0;
// 	Usource.set(nz-1,i,j,nr-1) = 0.0;
// 	Vsource.set(nz-1,i,j,nr-1) = 0.0;
//     }
//   }      
//    Valeur ** Coeff = Usource.asymptot(3);
//    cout <<*Coeff [1];
   des_profile(test-LHS,1.0001,20.0,1.2,0);
//    des_profile(LHS,1.0001,20.0,1.2,0);

      Usource.inc2_dzpuis();
      Usource.inc_dzpuis();
  

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
