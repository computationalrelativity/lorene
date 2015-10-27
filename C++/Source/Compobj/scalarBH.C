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
 * Revision 1.2  2015/10/27 10:53:23  f_vincent
 * Updated class scalarBH
 *
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

    if (rHor<0. || nrfile<0. || nthetafile<0.){
      cerr << "In scalarBH::scalarBH(Map,char*): "
	   << "Bad parameters rHor, nrfile, nthetafile" << endl;
      abort();
    }

    cout << "nr, ntheta from file = " << nrfile << " " << nthetafile << endl;
    cout << "rHor from file = " << rHor << endl;

    double* Xfile = new double[nrfile * nthetafile] ;
    double* thetafile = new double[nrfile * nthetafile] ;
    double* f0file = new double[nrfile * nthetafile] ;
    double* f1file = new double[nrfile * nthetafile] ;
    double* f2file = new double[nrfile * nthetafile] ;
    double* wwfile = new double[nrfile * nthetafile] ;
    double* sfieldfile = new double[nrfile * nthetafile] ;

    cout << "Reading metric data... ";
    for (int ii=0;ii<nrfile*nthetafile;ii++){
      // there are empty lines in Carlos file, but it doesn't seem to be a pb
      file >> Xfile[ii] ; 
      file >> thetafile[ii] ;
      file >> f0file[ii] ;
      file >> f1file[ii] ;
      file >> f2file[ii] ;
      file >> wwfile[ii] ;
      file >> sfieldfile[ii] ;
      //cout << ii << " " << Xfile[ii] << " " << thetafile[ii] << " " << f0file[ii] << " " << f1file[ii] << " " << f2file[ii] << " " << wwfile[ii] << " " << sfieldfile[ii] << endl;
    }  
    cout << "done" << endl;
    file.close() ; 

    double Xbefmax = Xfile[nrfile-2];

    int nz = mg->get_nzone() ; 
    //cout << "nz : " << nz << endl ; 
    ff0.allocate_all() ; // Memory allocation for F_0
    ff1.allocate_all() ; // Memory allocation for F_1
    ff2.allocate_all() ; // Memory allocation for F_2
    ww.allocate_all() ; // Memory allocation for W
    sfield.allocate_all() ; // Memory allocation for scalar field phi

    double delta_theta = 0.1*fabs(thetafile[nrfile] - thetafile[0]); // small wrt the theta step in Carlos grid

    cout << "Starting interpolating Lorene grid... " ;
    Mtbl rr(mp.r); 
    Mtbl theta(mp.tet); 
    int l_min; // this should be 0 for a spacetime without horizon, 1 with
    if (rHor>0.){
      l_min = 1; 
    }else{
      l_min = 0; 
    }
    for (int l=l_min; l<nz; l++) {
      int nr = mg->get_nr(l) ; 
      int nt = mg->get_nt(l) ; 
      int np = mg->get_np(l) ;
      //cout << "Starting loop k j i" << endl;
      for (int k=0; k<np; k++){
	for (int j=0; j<nt; j++){
	  double th0 = theta(l, k, j, 0);
	  for (int i=0; i<nr; i++){
	    double r0 = rr(l, k, j, i);
	    double x0, xx0;
	    if (r0 < __infinity){
	      x0 = sqrt(r0*r0 - rH2);
	      xx0 = x0 / (x0+1);
	    }else{
	      xx0 = 1.;
	    }
	    
	    //cout << "Loene radial stuff= " << r0 << " " << x0 << " " << xx0 << endl;
	    //cout << "Lorene points= " << xx0 << " " << th0 << endl;
	      
	    int ith=0;
	    int ithbis=0;
	    double thc = thetafile[ith];
	    while (fabs(th0 - thc) > delta_theta){
	      ith += nrfile;
	      ithbis++;
	      thc = thetafile[ith];
	    }
	    int ir=ith;
	    double xc = Xfile[ir];
	    if (xc > 0.){
	      cerr << "In scalarBH::ScalarBH(): "
		"r should be zero here" << endl;
	      abort();
	    }
	    int doradinterp=1;
	    if (xx0 == 0.){
	      doradinterp=0;
	    }else if(xx0 == 1.){
	      ir += nrfile-1;
	      xc = Xfile[ir];
	      doradinterp=0;
	    }else if(xx0 > Xbefmax && xx0< 1.){
	      ir+=nrfile-2;
	      xc = Xfile[ir];
	    }else{
	      //cout << "Indice= " << ithbis << " " << ir << " " << xc << endl;

	      while (xx0 > xc){
		ir ++;
		xc = Xfile[ir];
		//cout << "xc, xx0 in loop= " << xc << " " << xx0 << endl;
	      }	      
	    }

	    double f0interp, f1interp, f2interp, winterp, sfieldinterp;
	    if (doradinterp){
	      double xcinf = Xfile[ir-1], xcsup = Xfile[ir+1];
	      if (fabs(thc-th0)>delta_theta){
		cerr << "scalarBH::ScalarBH(): theta problem in grid" << endl;
		cerr << "Theta info: " << thc << " " << th0 << endl;
		abort();
	      }
	      if (xx0 <= Xbefmax &&
		  (xx0 < xcinf || xx0 > xc || xx0 > xcsup)){
		cerr << "scalarBH::ScalarBH(): rad problem in grid" << endl;
		cerr << "Radial info: " << xcinf << " " << xx0 << " " 
		     << xc << " " << xcsup << endl;
		abort();
	      }else if (xx0 > Xbefmax &&
			(xx0 < xcinf || xx0 < xc || xx0 > xcsup)){
		cerr << "scalarBH::ScalarBH(): special rad "
		  "problem in grid" << endl;
		cerr << "Radial info: " << xcinf << " " << xx0 << " " 
		     << xc << " " << xcsup << endl;
		abort();
	      }
	      // Radial polynomials
	      double polyrinf = (xx0-xc)*(xx0-xcsup)/((xcinf-xc)*(xcinf-xcsup)),
		polyrmid = (xx0-xcinf)*(xx0-xcsup)/((xc-xcinf)*(xc-xcsup)),
		polyrsup = (xx0-xcinf)*(xx0-xc)/((xcsup-xcinf)*(xcsup-xc));
	      // Grid values of all Scalars
	      double f0inf = f0file[ir-1], f0mid = f0file[ir], 
		f0sup = f0file[ir+1], f1inf = f1file[ir-1], f1mid = f1file[ir], 
		f1sup = f1file[ir+1], f2inf = f2file[ir-1], f2mid = f2file[ir], 
		f2sup = f2file[ir+1], winf = wwfile[ir-1], wmid = wwfile[ir], 
		wsup = wwfile[ir+1], sfinf = sfieldfile[ir-1], 
		sfmid = sfieldfile[ir], sfsup = sfieldfile[ir+1];
	      
	      /*cout << "Interpolating" << endl;
	      cout << "Carlos points= " << xcinf << " " << xc << " " << xcsup << endl;
	      cout << "f0= " << f0inf << " " << f0mid << " " << f0sup << endl;*/
	      
	      // Interpolate Scalars
	      f0interp = f0inf*polyrinf + f0mid*polyrmid + f0sup*polyrsup;
	      f1interp = f1inf*polyrinf + f1mid*polyrmid + f1sup*polyrsup;
	      f2interp = f2inf*polyrinf + f2mid*polyrmid + f2sup*polyrsup;
	      winterp = winf*polyrinf + wmid*polyrmid + wsup*polyrsup;
	      sfieldinterp = sfinf*polyrinf + sfmid*polyrmid + sfsup*polyrsup;
	      
	      //cout << "f0interp= " << f0interp << endl;
	    }else{
	      
	      /*cout << "Not interpolating" << endl;
	      cout << "Carlos point= " << xc << endl;
	      cout << "f0= " << f0file[ir] << endl;*/
	      // No interpolation at grid ends
	      f0interp = f0file[ir] ;
	      f1interp = f1file[ir] ;
	      f2interp = f2file[ir] ;
	      winterp = wwfile[ir] ;
	      sfieldinterp = sfieldfile[ir] ;
	    }
	    //cout << "theta stuff= " << thc << " should be equal to " << th0 << endl;
	    //cout << "X stuff= " << Xfile[ir-1] << " " << xx0 << " " << xc << " " << Xfile[ir+1] << endl;

	    ff0.set_grid_point(l,k,j,i) = f0interp ;
	    ff1.set_grid_point(l,k,j,i) = f1interp ;
	    ff2.set_grid_point(l,k,j,i) = f2interp ;
	    ww.set_grid_point(l,k,j,i) = winterp ;
	    sfield.set_grid_point(l,k,j,i) = sfieldinterp ;
	  }
	} 
      }
    }

    // Deleting arrays useless now
    delete[] Xfile;
    delete[] thetafile;
    delete[] f0file;
    delete[] f1file;
    delete[] f2file;
    delete[] wwfile;
    delete[] sfieldfile;

    ff0.std_spectral_base() ;
    ff1.std_spectral_base() ;
    ff2.std_spectral_base() ;
    ww.std_spectral_base() ;
    sfield.std_spectral_base() ;

    cout << "done." << endl;

    // At this point the Scalar ff0, ff1, ff2, ww, sfield
    // are initialized on the Lorene grid to proper interpolated values

    // Next step: from the Lorene metric
   
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
