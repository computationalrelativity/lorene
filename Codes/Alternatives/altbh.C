/*
 *  Code for reading alternative BH spacetime from a file
 */

/*
 *   Copyright (c) 2013 Odele Straub, Eric Gourgoulhon
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

char altBH_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/02/12 16:44:54  o_straub
 * New output : calculates radii for a range of spin and stores them in file
 *
 * Revision 1.3  2014/01/14 20:54:47  e_gourgoulhon
 * Comparison with Kerr_QI; better outputs
 *
 * Revision 1.2  2013/04/17 13:03:29  e_gourgoulhon
 * New member krphi
 *
 * Revision 1.1  2013/04/16 15:29:33  e_gourgoulhon
 * New code for reading Enrico's data
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


int main() {

    // Parameters of the computation
    // -----------------------------
    
    ifstream fpar("par_altbh.d") ;
    if ( !fpar.good() ) {
        cerr << "Problem in opening the file par_altbh.d ! " << endl ;
        abort() ;
    }
    
    char file_name[256] ; 
    fpar.getline(file_name, 256) ;
    cout << "File to be read: " << file_name << endl ; 

    double mass ; // M
    fpar >> mass ; fpar.ignore(1000,'\n') ;
    
    double a_ov_m ; 
    fpar >> a_ov_m ; fpar.ignore(1000,'\n') ;

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
    
   // cout << "M = " << mass << ",  a/M = " << a_ov_m << endl ; 
    double r_hor = double(0.5)*mass*sqrt(double(1)-a_ov_m*a_ov_m) ; 

   // cout << "Value of coordinate r at the event horizon : " << r_hor  << endl ; 
   // cout << "r_limits : " ; 
   // for (int l=0; l<nz+1; l++) {
   //     cout << r_limits[l] << "  " ; 
   //  }
   //  cout << endl ; 
    
   //  arrete() ; 


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

    // Construction of the AltBH_QI object:
    // ----------------------------------

    AltBH_QI bh(map, file_name, a_ov_m) ; 
    Kerr_QI bh2(map, 1., a_ov_m) ; 

    bh.update_metric() ; 
    bh2.update_metric() ; 


    cout.precision(15) ; 
    cout << endl << "******* bh ******** " << endl ;  
    cout << bh << endl ; 
    cout << endl << "******* bh2 ******** " << endl ; 
    cout << bh2 << endl ; 


	// SOME TESTS
	//    Scalar diffn = bh.get_nn() - bh2.get_nn();
	//    cout << "max(N - N(2)) : " << max(diffn) << endl ; 
	//    Scalar tmp = bh.get_nn().dsdr() ; 
	//    tmp.annule_domain(nzm1) ;
	//    tmp = tmp.dsdr() ; 
	//    Scalar tmp2 = bh2.get_nn().dsdr() ; 
	//    tmp2.annule_domain(nzm1) ;
	//    tmp2 = tmp2.dsdr() ; 
	//    Scalar diffddn = tmp - tmp2;
	//    cout << "max(d2N - d2N(2)) : " << max(diffddn) << endl ; 
	//    Scalar diffa2 = bh.get_a_car() - bh2.get_a_car();
	//    cout << "max(A^2 - A^2(2)) : " << max(diffa2) << endl ; 
	//    Scalar diffb2 = bh.get_b_car() - bh2.get_b_car();
	//    cout << "max(B^2 - B^2(2)) : " << max(diffb2) << endl ; 
	//    Scalar diffnp = bh.get_nphi() - bh2.get_nphi();
	//    cout << "max(Nphi - Nphi(2)) : " << max(diffnp) << endl ; 
	//    arrete() ;
	    



	// ISCO from Eq. (21) of Bardeen, Press & Teukolsky, ApJ 178, 347 (1972):
	double third = double(1)/double(3) ;
	double z1 = 1 + pow(1-a_ov_m*a_ov_m, third)*( pow(1+a_ov_m, third) +
	                                              pow(1-a_ov_m, third) ) ;
	double z2 = sqrt(3*a_ov_m*a_ov_m + z1*z1) ;

	double R_isco     =  3 + z2 - sqrt((3-z1)*(3+z1+2*z2))  ; 
	double R_isco_QI  = bh2.r_isco(0)/mass + mass*(1-a_ov_m*a_ov_m)/(4*bh2.r_isco(0)) + 1 ;
	double R_isco_alt = bh.r_isco(0)/mass  + mass*(1-a_ov_m*a_ov_m)/(4*bh.r_isco(0))  + 1 ;


	cout.precision(15) ; 
	cout << "Analytic value of R_ISCO (Kerr):         " << R_isco     << " M " << endl ; 
	cout << "Numerical value of R_ISCO (Kerr_QI):     " << R_isco_QI  << " M " << endl ; 
	cout << "Numerical value of R_ISCO (alternative): " << R_isco_alt << " M " << endl ;  



	// R_mb analogue to Eq. (19) of Bardeen, Press & Teukolsky, ApJ 178, 347 (1972):
	double R_mb     = 2 - a_ov_m + 2 * sqrt(1 - a_ov_m) ;
	double R_mb_QI  = bh2.r_mb(0)/mass + mass*(1 - a_ov_m*a_ov_m)/(4*bh2.r_mb(0)) + 1 ;
	double R_mb_alt = bh.r_mb(0)/mass  + mass*(1 - a_ov_m*a_ov_m)/(4*bh.r_mb(0))  + 1 ;

	cout.precision(15) ; 
	cout << "Analytical value of R_mb (Kerr):       " << R_mb     << " M " << endl ;  
	cout << "Numerical value of R_mb (Kerr_QI):     " << R_mb_QI  << " M " << endl ;
	cout << "Numerical value of R_mb (alternative): " << R_mb_alt << " M " << endl ;




	// Drawings    
	    if (graphic_out == 1) 
	    {
	       double r_max = 1.5*map.val_r(nzm1,-1.,0.,0.) ; 
	       des_meridian(bh.get_nn(), 0, r_max, "N", 1) ; 
	       //des_meridian(diffn , 0, r_max, "N - N(2)", 2) ; 
	       //des_meridian(diffddn , 0, r_max, "diff d2N", 3) ; 

	       des_meridian(bh.get_nphi(), 0, r_max, "Nphi", 4) ; 
	       des_meridian(bh.get_nphi().dsdr(), 0, r_max, "dNphi/dr", 5) ; 

	       des_meridian(bh.get_gamma().cov()(1,1), 0, r_max, "gamma_11", 6) ; 
	       des_meridian(bh.get_gamma().cov()(3,3), 0, r_max, "gamma_33", 7) ; 
	 
	       des_meridian(bh.get_kk()(1,3), 0.4, 0.6, "K_(r)(ph)", 8) ; 
	       des_meridian(bh.get_krphi(), 0.4, 0.6, "K_(r)(ph)/sin(theta) from file", 9) ; 
	       des_meridian(bh.get_kk()(2,3), 0, r_max, "K_(th)(ph)", 10) ; 
	       //arrete() ; 
	     }



	//----------------------
	// Output file for GYOTO
	//----------------------
	    bh.gyoto_data("gyoto_altBH_QI.d") ;
	    


	// a small loop to get Risco and Rmb for a range of spins
	// and store the results in a file
	    for (a_ov_m=0.0001; a_ov_m<0.3; a_ov_m+=0.01)
	    {  

	    AltBH_QI bh(map, file_name, a_ov_m) ; 
	    
	    bh.update_metric() ; 
	    
	    // ISCO from Eq. (21) of Bardeen, Press & Teukolsky, ApJ 178, 347 (1972):
	    double third = double(1)/double(3) ;
	    double z1 = 1 + pow(1-a_ov_m*a_ov_m, third)*( pow(1+a_ov_m, third) +
	                                                  pow(1-a_ov_m, third) ) ;
	    double z2 = sqrt(3*a_ov_m*a_ov_m + z1*z1) ;

	    double R_isco     =  3 + z2 - sqrt((3-z1)*(3+z1+2*z2))  ; 
	    double R_isco_alt = bh.r_isco(0)/mass  + mass*(1-a_ov_m*a_ov_m)/(4*bh.r_isco(0))  + 1 ;

	    // R_mb analogue to Eq. (19) of Bardeen, Press & Teukolsky, ApJ 178, 347 (1972):
	    double R_mb     = 2 - a_ov_m + 2 * sqrt(1 - a_ov_m) ;
	    double R_mb_alt = bh.r_mb(0)/mass  + mass*(1 - a_ov_m*a_ov_m)/(4*bh.r_mb(0))  + 1 ;

	     
	    //-------------------------------
	    // Output files for Risco and Rmb
	    //-------------------------------
	    freopen("radbound.txt", "a", stdout) ; // "a" for append (during loop)
	    cout.precision(15) ;
	    cout << a_ov_m << " " << R_mb << " " << R_mb_alt << " " << endl ;
	    
	    freopen("radisco.txt", "a", stdout) ;  
	    cout.precision(15) ;
	    cout << a_ov_m << " " << R_isco << " " << R_isco_alt << " " << endl ;

	    freopen("dummy.txt", "w", stdout) ;  // empty log file (needed for clean stdout - not sure why)
	    cout << " " << endl ;
	    
	    }

	    fclose (stdout) ;

	return EXIT_SUCCESS ; 

}













