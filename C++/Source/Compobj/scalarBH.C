/*
 *  Methods of the class ScalarBH
 *
 *    (see file compobj.h for documentation).
 *
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

char scalarBH_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2015/10/22 09:18:36  f_vincent
 * New class ScalarBH
 *
 *
 * $Header$
 *
 */


// C headers
#include <cassert>

// Lorene headers
#include "compobj.h"
#include "unites.h"
#include "nbr_spx.h"

//--------------//
// Constructors //
//--------------//

// Standard constructor
// --------------------
namespace Lorene {
  ScalarBH::ScalarBH(Map& mpi, const char* file_name) :
    Compobj(mpi),
    ff0(mpi),
    ff1(mpi),
    ff2(mpi),
    ww(mpi),
    sfield(mpi)
  {

    ifstream file(file_name) ; 
    if ( !file.good() ) {
      cerr << "Problem in opening the file " << file_name << endl ;
      abort() ;
    }
  
    const Mg3d* mg = mp.get_mg() ; 
    double rHor, rH2 ;
    int nrfile, nthetafile;
    file >> nrfile >> nthetafile ;
    file >> rHor ;
    rH2 = rHor*rHor ;

    cout << "nr, ntheta from file = " << nrfile << " " << nthetafile << endl;
    cout << "rHor from file = " << rHor << endl;

    double* rfile = new double[nrfile * nthetafile] ;
    double* thetafile = new double[nthetafile * nthetafile] ;
    double* f0file = new double[nrfile * nthetafile] ;
    double* f1file = new double[nrfile * nthetafile] ;
    double* f2file = new double[nrfile * nthetafile] ;
    double* wwfile = new double[nrfile * nthetafile] ;
    double* sfieldfile = new double[nrfile * nthetafile] ;

    cout << "Reading metric data... ";
    for (int ii=0;ii<nrfile*nthetafile;ii++){
      // there are empty lines in Carlos file, but it doesn't seem to be a pb
      file >> rfile[ii] ; 
      file >> thetafile[ii] ;
      file >> f0file[ii] ;
      file >> f1file[ii] ;
      file >> f2file[ii] ;
      file >> wwfile[ii] ;
      file >> sfieldfile[ii] ;
      //cout << ii << " " << rfile[ii] << " " << thetafile[ii] << " " << f0file[ii] << " " << f1file[ii] << " " << f2file[ii] << " " << wwfile[ii] << " " << sfieldfile[ii] << endl;
    }  
    cout << "done" << endl;
    file.close() ; 


    int nz = mg->get_nzone() ; 
    cout << "nz : " << nz << endl ; 
    cout << "Initializing metric scalars... ";
    ff0.allocate_all() ; // Memory allocation for F_0
    ff1.allocate_all() ; // Memory allocation for F_1
    ff2.allocate_all() ; // Memory allocation for F_2
    ww.allocate_all() ; // Memory allocation for W
    sfield.allocate_all() ; // Memory allocation for scalar field phi
    cout << "done." << endl;

    Mtbl rr(mp.r); 
    int l_min = 0; // this should be 0 for a spacetime without horizon, 1 with
    for (int l=l_min; l<nz; l++) {
      cout << "l = " << l << endl ; 
      int nr = mg->get_nr(l) ; 
      int nt = mg->get_nt(l) ; 
      int np = mg->get_np(l) ;
      cout << "Starting loop k j i" << endl;
      for (int k=0; k<np; k++){
	for (int j=0; j<nt; j++){
	  for (int i=0; i<nr; i++){
	    cout << "Calling r at grid point" << endl;
	    double r0 = rr(l, k, j, i);
	    cout << "Few computation" << endl;
	    double x0 = sqrt(r0*r0 - rH2);
	    double xx0 = x0 / (x0+1);
	    double f0 = 0; // here interpolation!
	    //ff0.set_spectral_va().set(l,0,j,i) = f0;
	    cout << "affecting value to metric coef" << endl;
	    ff0.set_grid_point(l,k,j,i) = f0 ;
	    cout << "the end" << endl;
	  }
	} 
      }
    }
   
  //bbb.std_spectral_base() ;

    // Deleting arrays
    delete[] rfile;
    delete[] thetafile;
    delete[] f0file;
    delete[] f1file;
    delete[] f2file;
    delete[] wwfile;
    delete[] sfieldfile;
   
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
  }

  // Copy constructor
  // --------------------
  ScalarBH::ScalarBH(const ScalarBH& other) :
    Compobj(other),
    ff0(other.ff0),
    ff1(other.ff0),
    ff2(other.ff0),
    ww(other.ff0),
    sfield(other.ff0)
  {
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
  }


  // Constructor from a file
  // -----------------------
  ScalarBH::ScalarBH(Map& mpi, FILE* ) :
    Compobj(mpi),
    ff0(mpi),
    ff1(mpi),
    ff2(mpi),
    ww(mpi),
    sfield(mpi)
  {
    // Pointers of derived quantities initialized to zero : 
    set_der_0x0() ;
    
    // Read of the saved fields:
    // ------------------------

  }

  //------------//
  // Destructor //
  //------------//

  ScalarBH::~ScalarBH(){

    del_deriv() ; 

  }


  //----------------------------------//
  // Management of derived quantities //
  //----------------------------------//

  void ScalarBH::del_deriv() const {

    Compobj::del_deriv() ; 


    ScalarBH::set_der_0x0() ; 
  }			    


  void ScalarBH::set_der_0x0() const {
 	 
  }			    

  //--------------//
  //  Assignment  //
  //--------------//

  // Assignment to another ScalarBH
  // --------------------------------
  void ScalarBH::operator=(const ScalarBH& other) {

    // Assignment of quantities common to all the derived classes of Compobj
    Compobj::operator=(other) ;	    
    
    del_deriv() ;  // Deletes all derived quantities
  }	

  //--------------//
  //	  Outputs   //
  //--------------//

  // Save in a file
  // --------------
  void ScalarBH::sauve(FILE* ) const {

    
  }

  // Printing
  // --------

  ostream& ScalarBH::operator>>(ostream& ost) const {

    using namespace Unites ;
	
    Compobj::operator>>(ost) ; 
    
    ost << endl << "Black hole with scalar hair (class ScalarBH) " << endl ; 
    //    ost << description1 << endl ; 
    //    ost << description2 << endl ; 
   
    return ost ; 
      
  }

  //-------------------------//
  //	Computational methods  //
  //-------------------------//
			    
  // Updates the extrinsic curvature
  // -------------------------------

  //void ScalarBH::extrinsic_curvature() {

    // FV: commenting out October 2015 to compile

    // // Special treatment for axisymmetric case:
    
    // if ( (mp.get_mg())->get_np(0) == 1) {
    
    //   // What follows is valid only for a mapping of class Map_radial :   
    //   assert( dynamic_cast<const Map_radial*>(&mp) != 0x0 ) ;
        
    //   Scalar tmp = krphi ; 
    //   tmp.mult_sint() ;  // multiplication by sin(theta)
    //   kk.set(1,3) = tmp ; 
        
    //   kk.set(2,3) = 0 ; 

    //   kk.set(1,1) = 0 ; 
    //   kk.set(1,2) = 0 ; 
    //   kk.set(2,2) = 0 ; 
    //   kk.set(3,3) = 0 ; 
    // }
    // else {

    //   // General case:

    //   Compobj::extrinsic_curvature() ; 
    // }
    
    // // Computation of A^2 K_{ij} K^{ij}
    // // --------------------------------
        
    // ak_car = 2 * ( kk(1,3)*kk(1,3) +  kk(2,3)*kk(2,3) ) / b_car ;
    
    // del_deriv() ; 

  //}
}
