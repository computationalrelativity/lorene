/*
 * Code for solving the statinary axisymmetric Einstein equations for given boundary data
 *
 */

/*
 *   Copyright (c) 2002  Jörg Frauendiener, Christian Klein
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
 * Revision 1.1  2002/10/17 07:56:49  j_frauendiener
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


const int MaxIt = 200;
const double shift = 0.1;	// used to shift the cos-term


//=============================================================

int main() {

  // Read grid parameters from a file
  // --------------------------------

  ifstream parfile("par.d") ; 
  char blabla[80] ;
  char datafilename[80];
  int nz, nr, np ;
  double mu;			// relaxation parameter
  double tol;			// tolerance between successive iterates  
  
  parfile >> nz ; parfile.getline(blabla, 80) ;  
  parfile >> nr ; parfile.getline(blabla, 80) ; 
  parfile >> np ; parfile.getline(blabla, 80) ;
  parfile >> mu ; parfile.getline(blabla, 80) ;
  parfile >> tol ; parfile.getline(blabla, 80) ;  
  parfile.getline(datafilename,80);
  parfile.close() ;


  // read the size in theta direction and the values of the radii
  // ----------------------------------------------------------------------------------------

  cout << "Reading data from file " << datafilename << endl;
  
  ifstream datafile(datafilename);
  int nt;
  double r1, r2;

  datafile >> nt       ; datafile.getline(blabla, 80) ;
  datafile >> r1 >> r2 ; datafile.getline(blabla, 80) ;
  
  
  // Construction of a multi-grid (Mg3d)
  // -----------------------------------
  
  Mg3d mgrid(nz, nr, nt, np, SYM, NONSYM, true) ;

  cout << "Mult_grid : " << mgrid << endl ; 

  // Construction of an affine mapping (Map_af)
  // ------------------------------------------

  double* r_limits = new double[nz+1] ; 
  assert( nz == 3 ) ; 
  r_limits[0] = 0 ; 
  r_limits[1] = r1 ; 
  r_limits[2] = r2 ; 
  r_limits[3] = __infinity ; 
  
  Map_af map(mgrid, r_limits) ; 
  
  // Construction of a scalar field (Cmp)
  // ------------------------------------

  const Coord& rr = map.r;
  const Coord& sint = map.sint;
  const Coord& pp= map.phi;
  

  Cmp rho2(map);
  rho2 = rr*rr*sint*sint;

  Cmp cos_phi(map);
  cos_phi = cos(pp+shift);
  cos_phi.std_base_scal();

  Cmp x(map) ;
  x = rr*sint;
  x = x*cos_phi;
  
  
  rho2.std_base_scal();
  x.std_base_scal();
  
  // construct fields on the angular grid
  // for boundary conditions   and comparisons

  Valeur bcU( mgrid.get_angu() ) ;  bcU.annule_hard(); 
  Valeur bcV( mgrid.get_angu() ) ;  bcV.annule_hard();
  
  Map_af mpa( *(mgrid.get_angu()), r_limits ) ;

  // read data from file per domain
  // domain 0: boundary condition at r=r1
  //  format of datafile must be changed to read things in one loop...

  for (int m=0; m < nt; m++)
    {
      double u, v;
      datafile >> u >> v ; datafile.getline(blabla, 80) ;
      for (int l=0; l < np; l++)
	{
	  
	  bcU.set(0,l,m,0) = u;
	  bcV.set(0,l,m,0) = v;
	}
    }
  
  // domain 1: values at r=r2  
  for (int m=0; m < nt; m++)
    {
      double u, v;
      datafile >> u >> v ; datafile.getline(blabla, 80) ;
      for (int l=0; l < np; l++)
	{
	  
	  bcU.set(1,l,m,0) = u;
	  bcV.set(1,l,m,0) = v;
	}
    }  

  bcV = bcV*cos(mpa.phi+shift);	// now bcV contains the bc for V*rho*cos(phi+shift)
  bcU.std_base_scal() ;
  bcV.std_base_scal() ;

  //  fields and initial data
  //---------------------

  Cmp U(map) ;			// log of the lapse
  Cmp V(map) ;			// component of the shift 
  Cmp Usource(map);		// rhs of U-eqn
  Cmp Vsource(map);		// rhs of V-eqn
  Cmp N4(map);			// exp(-4U)
  
  U.annule_hard() ;
  V = U;

  U.std_base_scal() ;
  V.std_base_scal() ;
  Usource.std_base_scal() ;	// sets the bases for spectral expansions
  Vsource.std_base_scal() ;	// to be the standard ones



  // solution of equations
  //---------------------------------  
  cout << "Starting iteration, shooting for tolerance " << tol << endl;
  
  int iter = 0;
  while (iter< MaxIt)
    {
      iter++;
      
      N4 = exp(-4.0*U);
      N4.std_base_scal();

      Usource = 0.5*rho2*N4*(V.dsdr()*V.dsdr() + V.srdsdt()*V.srdsdt());
      Vsource = 4.0*x*(V.dsdr()*U.dsdr() + V.srdsdt()*U.srdsdt()) ;
      Vsource.set(0) = 1.0;
      Usource.set(0) = 1.0;

      // force the rhs to be zero at infinity
      for (int i=0;i<np;i++)
	for(int j=0; j<nt;j++) {
	  Usource.set(nz-1,i,j,nr-1) = 0.0;
	  Vsource.set(nz-1,i,j,nr-1) = 0.0;
	}
  
      
      Cmp U1 = Usource.poisson_dirichlet(bcU, 0) ;
      Cmp V1 = Vsource.poisson_dirichlet(bcV, 0) ;
      U1.set(0) = 0.5;

      V1.div_rsint();

      V1 = V1/cos_phi;
      V1.set(0) = 0.5;

      Tbl diff =  norme(abs(U-U1)) + norme(abs(V-V1));
      
      //      cout << "Iteration " << iter 
      //	   << ":  Difference:  "
      //   << diff(1) << "\t" << diff(2) << endl;

      U = (1-mu)*U + mu*U1;
      V = (1-mu)*V + mu*V1;

      if ((diff(1) < tol) && (diff(2) < tol)) break;
    }

  // no convergence to specified tolerance
  if (iter >= MaxIt)
    {
      cout << "Tolerance could not be achieved after "<< iter << " Iterations."<< endl;
      cout << "Stop." << endl;
      return EXIT_SUCCESS;
    }
  
  /*

  Valeur Uexact( mgrid.get_angu());
  Valeur Vexact( mgrid.get_angu());
  Uexact.annule_hard();
  Vexact.annule_hard();  
  
  for (int m=0; m < nt; m++)
    {
      for (int l=0; l < np; l++)
	{
	  Uexact.set(1,l,m,0) = U(1,l,m,nr-1);
	  Vexact.set(1,l,m,0) = V(1,l,m,nr-1);	  
	}
    }
  */

  V.mult_rsint();
  cout.setf(ios::scientific,ios::floatfield);
  cout.precision(7);
  
  for (int j=0; j<nt; j++)
    {
      cout << "U: " << U(1,0,j,nr-1) << "\t" << bcU(1,0,j,0) << "\t"
	   << U(1,0,j,nr-1) - bcU(1,0,j,0) << endl;
      cout << "V: " << V(1,0,j,nr-1) << "\t" << bcV(1,0,j,0) << "\t"
	   << V(1,0,j,nr-1) - bcV(1,0,j,0) << endl;
      cout << "------------------------------------------------------------------------------" << endl;
      
    }

  //  des_profile(U,1.01,20.0,M_PI/2,M_PI/2);
  //  des_profile(V,1.01,20.0,M_PI/2,M_PI/2);
 
  return EXIT_SUCCESS ; 

}
