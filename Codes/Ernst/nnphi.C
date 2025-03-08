/*
 * Code for solving the stationary axisymmetric Einstein equations for given boundary data
 *
 */

/*
 *   Copyright (c) 2002  J�rg Frauendiener, Christian Klein
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
 * $Id$
 * $Log$
 * Revision 1.8  2016/12/05 16:18:24  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:53:56  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:09:44  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2003/02/23 20:54:17  c_klein
 * new output format
 *
 * Revision 1.4  2003/01/09 11:07:50  j_novak
 * headcpp.h is now compliant with C++ norm.
 * The include files have been ordered, as well as the local_settings_linux
 *
 * Revision 1.3  2002/12/09 16:07:25  j_frauendiener
 * binary output added
 *
 * Revision 1.2  2002/10/28 15:48:29  j_frauendiener
 * Multiplication with cos(phi) added in comparison of results
 *
 * Revision 1.1  2002/10/17 07:56:49  j_frauendiener
 * Initial revision
 *
 * $Header$
 *
 */


// C++ headers
#include "headcpp.h"


// C headers
#include <cstdlib>

// Lorene headers
#include "cmp.h"
#include "nbr_spx.h"
#include "graphique.h"


const int MaxIt = 500;
const double shift = 0.1;	// used to shift the cos-term


//=============================================================

using namespace Lorene ;

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
      double u1, v1, u2, v2 ;
      datafile >> u1 >> v1 >> u2 >> v2; datafile.getline(blabla, 80) ;
      for (int l=0; l < np; l++)
	{
	  bcU.set(0,l,m,0) = u1;
	  bcV.set(0,l,m,0) = v1;
	  bcU.set(1,l,m,0) = u2;
	  bcV.set(1,l,m,0) = v2;
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
      Vsource.set(0) = 0.5;
      Usource.set(0) = 0.5;

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
      
      cout << "Iteration " << iter 
      	   << ":  Difference:  "
	   << diff(1) << "\t" << diff(2) << endl;

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
  bcV = bcV/cos(mpa.phi+shift);

     ofstream out("out.dat");
     out.setf(ios::scientific,ios::floatfield);
     out.precision(7);

//   double a[nt][4];

   for (int j=0; j<nt; j++)
    {
            cout <<  U(1,0,j,nr-1) << "\t" << bcU(1,0,j,0) << "\t"
      	   << U(1,0,j,nr-1) - bcU(1,0,j,0) << endl;
            cout <<  V(1,0,j,nr-1) << "\t" << bcV(1,0,j,0) << "\t"
      	   << V(1,0,j,nr-1) - bcV(1,0,j,0) << endl;
            cout << "------------------------------------------------------------------------------" << endl;
	     out << bcU(1,0,j,0)  - U(1,0,j,nr-1) << endl;
//       a[j][0] = U(1,0,j,nr-1);    
//       a[j][1] = V(1,0,j,nr-1);    
//       a[j][2] = bcU(1,0,j,0);    
//       a[j][3] = bcV(1,0,j,0);    
    }

 

     
    out.close();  
  return EXIT_SUCCESS ; 

}
