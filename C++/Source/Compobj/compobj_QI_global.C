/*
 * Method of class Compobj_QI to compute the location of the ISCO
 *
 * (see file compobj.h for documentation).
 *
 */

/*
 *   Copyright (c) 2012 Odele Straub, Claire Some, Eric Gourgoulhon
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


char compobj_QI_global_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2013/11/13 11:20:01  j_novak
 * Minor correction to compile with older versions of g++
 *
 * Revision 1.5  2013/07/25 19:44:11  o_straub
 * calculation of the marginally bound radius
 *
 * Revision 1.4  2013/04/04 15:32:32  e_gourgoulhon
 * Better computation of the ISCO
 *
 * Revision 1.3  2012/11/22 16:04:51  c_some
 * Minor modifications
 *
 * Revision 1.2  2012/11/21 14:53:45  c_some
 * corrected mom_euler
 *
 * Revision 1.1  2012/11/16 16:14:11  c_some
 * New class Compobj_QI
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "compobj.h"
#include "param.h"
#include "utilitaires.h"

double funct_compobj_QI_isco(double, const Param&) ; 
double funct_compobj_QI_rmb(double, const Param&) ;


			//----------------------------//
			//	Angular momentum      //
			//----------------------------//

double Compobj_QI::angu_mom() const {

    if (p_angu_mom == 0x0) {    // a new computation is required
	
		cerr << "Compobj_QI::angu_mom() : not implemented yet !" << endl ; //## provisory
		abort() ; 
    }
    
    return *p_angu_mom ; 

}



//=============================================================================
//		r_isco()
//=============================================================================

double Compobj_QI::r_isco(int lmin, ostream* ost) const {

    if (p_r_isco == 0x0) {    // a new computation is required

    // First and second derivatives of metric functions
    // ------------------------------------------------

    int nzm1 = mp.get_mg()->get_nzone() - 1 ;
    Scalar dnphi = nphi.dsdr() ;
    dnphi.annule_domain(nzm1) ;
    Scalar ddnphi = dnphi.dsdr() ;    	// d^2/dr^2 N^phi

    Scalar tmp = nn.dsdr() ;
    tmp.annule_domain(nzm1) ;
    Scalar ddnnn = tmp.dsdr() ; 		// d^2/dr^2 N

    tmp = bbb.dsdr() ;
    tmp.annule_domain(nzm1) ;
    Scalar ddbbb = tmp.dsdr() ; 		// d^2/dr^2 B

    // Constructing the velocity of a particle corotating with the star
    // ----------------------------------------------------------------

    Scalar bdlog = bbb.dsdr() / bbb ;
    Scalar ndlog = nn.dsdr() / nn ;
    Scalar bsn = bbb / nn ;

    Scalar r(mp) ;
    r = mp.r ;

    Scalar r2= r*r ;

    bdlog.annule_domain(nzm1) ;
    ndlog.annule_domain(nzm1) ;
    bsn.annule_domain(nzm1) ;
    r2.annule_domain(nzm1) ;

    // ucor_plus - the velocity of corotating particle on the circular orbit
    Scalar ucor_plus = (r2 * bsn * dnphi +
        sqrt ( r2 * r2 * bsn *bsn * dnphi * dnphi +
		4 * r2 * bdlog * ndlog + 4 * r * ndlog) ) /
		2 / (1 + r * bdlog ) ;

    Scalar factor_u2 = r2 * (2 * ddbbb / bbb - 2 * bdlog * bdlog +
    					  4 * bdlog * ndlog ) +
       2 * r2 * r2 * bsn * bsn * dnphi * dnphi +
       4 * r * ( ndlog - bdlog ) - 6 ;

    Scalar factor_u1 = 2 * r * r2 * bsn * ( 2 * ( ndlog - bdlog ) *
       				dnphi - ddnphi ) ;

    Scalar factor_u0 = - r2 * ( 2 * ddnnn / nn - 2 * ndlog * ndlog +
					 4 * bdlog * ndlog ) ;

    // Scalar field the zero of which will give the marginally stable orbit
    Scalar orbit = factor_u2 * ucor_plus * ucor_plus +
				factor_u1 * ucor_plus + factor_u0 ;
    orbit.std_spectral_base() ;
    
    // Search for the zero
    // -------------------

    double r_ms, theta_ms, phi_ms, xi_ms, 
    	   xi_min = -1, xi_max = 1; 
  
    int l_ms = lmin, l ;
    bool find_status = false ; 

    const Valeur& vorbit = orbit.get_spectral_va() ;
        
    // Preliminary location of the zero
    // of the orbit function with an error = 0.01
    theta_ms = M_PI / 2. ; // pi/2
    phi_ms = 0. ;

    for(l = lmin; l <= nzm1; l++) { 

    xi_min = -1. ;
    xi_max = 1. ;

    double resloc_old ;
    double xi_f = xi_min;
    
    double resloc = vorbit.val_point(l, xi_min, theta_ms, phi_ms) ;
	
	for (int iloc=0; iloc<200; iloc++) {
		xi_f = xi_f + 0.01 ;
     	resloc_old = resloc ;
     	resloc = vorbit.val_point(l, xi_f, theta_ms, phi_ms) ;
     	if ( resloc * resloc_old < double(0) ) {
			xi_min = xi_f - 0.01 ;
			xi_max = xi_f ;
			l_ms = l ; 
			find_status = true ;  
			break ;
		} 
		
  	  }
      if (find_status) break ; 
    } 
  	
    Param par_ms ;
    par_ms.add_int(l_ms) ;
    par_ms.add_scalar(orbit) ;
    


  	if(find_status) { 
  	    
        cout << "l_ms : " << l_ms << endl ; 
        cout << "xi_min, xi_max : " << xi_min << " , " << xi_max << endl ; 
        
     	double precis_ms = 1.e-12 ;    // precision in the determination of xi_ms
									   
     	int nitermax_ms = 100 ;	       // max number of iterations

     	int niter ;
     	xi_ms = zerosec(funct_compobj_QI_isco, par_ms, xi_min, xi_max,
     					precis_ms, nitermax_ms, niter) ;     					
  		if (ost != 0x0) {
     		* ost <<
     		"    number of iterations used in zerosec to locate the ISCO : "
	  		 << niter << endl ;
     		*ost << "    zero found at xi = " << xi_ms << endl ;
        }

      	r_ms = mp.val_r(l_ms, xi_ms, theta_ms, phi_ms) ;
  	
	} else { 
		
		// ISCO not found
		r_ms = -1 ; 
	        xi_ms = -1 ; 
	        l_ms = lmin ; 
	    
	} 
  	
    p_r_isco = new double (r_ms) ;

//     p_r_isco = new double (
// 		(bbb.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms) * r_ms
//   					  ) ;

	// Determination of the frequency at the marginally stable orbit
  	// -------------------------------------------------------------				

	ucor_plus.std_spectral_base() ;
	double ucor_msplus = (ucor_plus.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms,
												  phi_ms) ;
	double nobrs = (bsn.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms) ;
	double nphirs = (nphi.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms) ;
	
 	p_f_isco = new double ( ( ucor_msplus / nobrs / r_ms + nphirs ) /
                       			(double(2) * M_PI) ) ;

	// Specific angular momentum on ms orbit
	// -------------------------------------				
	p_lspec_isco=new double (ucor_msplus/sqrt(1.-ucor_msplus*ucor_msplus)*
	((bbb.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms)) * r_ms );

	// Specific energy on ms orbit
	// ---------------------------			
	p_espec_isco=new double (( 1./nobrs / r_ms / ucor_msplus  + nphirs) *
	(ucor_msplus/sqrt(1.-ucor_msplus*ucor_msplus)*
	((bbb.get_spectral_va()).val_point(l_ms, xi_ms, theta_ms, phi_ms)) * r_ms ));
	

    }  // End of computation

    return *p_r_isco ;

}



//=============================================================================
//		f_isco()
//=============================================================================

double Compobj_QI::f_isco(int lmin) const {

    if (p_f_isco == 0x0) {    // a new computation is required

    	r_isco(lmin) ; 		// f_isco is computed in the method r_isco()

    	assert(p_f_isco != 0x0) ;
    }

    return *p_f_isco ;

}

//=============================================================================
//		lspec_isco()
//=============================================================================

double Compobj_QI::lspec_isco(int lmin) const {

    if (p_lspec_isco == 0x0) {    // a new computation is required

    	r_isco(lmin) ;    	// lspec_isco is computed in the method r_isco()

    	assert(p_lspec_isco != 0x0) ;
    }

    return *p_lspec_isco ;

}

//=============================================================================
//		espec_isco()
//=============================================================================

double Compobj_QI::espec_isco(int lmin) const {

    if (p_espec_isco == 0x0) {    // a new computation is required

    	r_isco(lmin) ;   	// espec_isco is computed in the method r_isco()

    	assert(p_espec_isco != 0x0) ;
    }

    return *p_espec_isco ;

}





//=============================================================================
//		r_mb()
//=============================================================================

double Compobj_QI::r_mb(int lmin, ostream* ost) const {

    if (p_r_mb == 0x0) {    // a new computation is required

    // Coefficients of the effective potential (A) and its derivative (B)
    // ------------------------------------------------
     
    int nzm1 = mp.get_mg()->get_nzone() - 1 ;   
    Scalar r(mp) ;
    r = mp.r ;
    Scalar r2 = r*r ;
    r2.annule_domain(nzm1) ;

    Scalar ndn = nn*nn.dsdr() ;
    ndn.annule_domain(nzm1) ;
     

    // Scalar V_eff  = A1 + A2 * E^2 + A3 * E * L + A4 * L^2 ;
    // Scalar dV_eff = B1 + B2 * E^2 + B3 * E * L + B4 * L^2 ;

    Scalar A1 =-bbb*bbb * nn*nn * r2 ;
    Scalar A2 = bbb*bbb * r2 ;
    Scalar A3 =-2. * bbb*bbb *  r2 * nphi ;
    Scalar A4 =-nn*nn + bbb*bbb * r2 * nphi*nphi ;

    Scalar B1 =-2.*r * bbb*bbb * nn*nn - 2.*r2 * bbb*bbb.dsdr() * nn*nn - 2.*r2 * bbb*bbb * nn*nn.dsdr() ;
    Scalar B2 = 2.*r * bbb*bbb + 2.*r2 * bbb*bbb.dsdr() ;
    Scalar B3 =-2.*nphi*B2 - 2.*r2 * bbb*bbb * nphi.dsdr() ;
    Scalar B4 = 2.*r * bbb*bbb * nphi*nphi + 2.*r2 * bbb*bbb.dsdr() * nphi*nphi - 2.*ndn + 2.*r2 * bbb*bbb * nphi*nphi.dsdr() ;

    
    Scalar C1 = (A1 * B3 - A3 * B1) ;
    Scalar C2 = (A2 * B3 - A3 * B2) ;
    Scalar C3 = (A4 * B3 - A3 * B4) ;


    Scalar D1 = B4 * C1 - B1 * C3 ; 
    Scalar D2 = B4 * C2 - B2 * C3 ; 
    Scalar D3 = B3 * B3 * C1 * C3 ;
    Scalar D4 = B3 * B3 * C2 * C3 ;



     


    // Constructing the orbital energy of a particle corotating with the star
    // ----------------------------------------------------------------

    /* B3 * V_eff - A3 * dV_eff = 0. ;         // solve eq. for L

     Scalar L = sqrt((C1 + C2 * E2) / C3) ;    // substitute into the eq. dVeff=0, then solve for E

     Scalar EE = (-(2.*D1*D2 + D3) + sqrt((2.*D1*D2 + D3) * (2.*D1*D2 + D3) - 4.*D1*D1 * (D2*D2 + 
            D4))) / (2.*(D2*D2 + D4)) ;        // solve eq. EE = 1 for r 
   */


    Scalar bound_orbit = -(2.*D1*D2 + D3) - sqrt((2.*D1*D2 + D3) * (2.*D1*D2 + D3) - 4.*D1*D1 * 
                             (D2*D2 + D4)) - 2.*(D2*D2 + D4) ;


    cout << "bound_orbit :" << bound_orbit << endl ;

    bound_orbit.std_spectral_base() ;
    


    // Search for the zero
    // -------------------

    const int noz(10) ;          // number of zeros    
    double zeros[2][noz] ; // define array for zeros
    int i = 0 ;            // counter
    int l ;                // number of domain

    double rmb, theta_mb, phi_mb, xi_mb;
    double xi_min = -1, xi_max = 1 ; 
  
    const Valeur& vorbit = bound_orbit.get_spectral_va() ;
        

    // Preliminary location of the zero
    // of the bound_orbit function with an error = dx

    double dx = 0.001 ; 
    
    theta_mb = M_PI / 2. ; 
    phi_mb = 0. ;


    for(l = lmin; l <= nzm1; l++) { 

        xi_min = -1. ;

        double resloc_old ;
        double xi_f = xi_min;
    
        double resloc = vorbit.val_point(l, xi_f, theta_mb, phi_mb) ;

        while(xi_f <= xi_max) { 

		    xi_f = xi_f + dx ;

     	    resloc_old = resloc ;
     	    resloc = vorbit.val_point(l, xi_f, theta_mb, phi_mb) ;

         	if ( resloc * resloc_old < double(0) ) {

                zeros[0][i] = xi_f ;   // xi_max
                zeros[1][i] = l ;      // domain number l  
                i++ ;         

            }

        }         

    }
 


    int number_of_zeros = i ;

    cout << "number of zeros: " << number_of_zeros << endl ;


    double precis_mb = 1.e-9 ;   // precision in the determination of xi_mb: 1.e-12
				                  

    int nitermax_mb = 100 ;	  // max number of iterations


    for(int i = 0; i < number_of_zeros; i++) {

	cout << i << " " << zeros[0][i] << " " << zeros[1][i] << endl ;

	int l_mb = int(zeros[1][i]) ;
	xi_max = zeros[0][i] ; 
	xi_min = xi_max - dx ;


        Param par_mb ;
        par_mb.add_scalar(bound_orbit) ;
        par_mb.add_int(l_mb) ; 

     	int niter ;
     	xi_mb = zerosec(funct_compobj_QI_rmb, par_mb, xi_min, xi_max, precis_mb, nitermax_mb, niter) ;     					
  	

	if (niter < nitermax_mb) {

		double zero_mb = mp.val_r(l_mb, xi_mb, theta_mb, phi_mb) ;
		//double r_hor = radius_hor(0); // set to 1 in the condition below
	 	double r_ms = r_isco(0) ;
		if (zero_mb < (1 + r_ms)/2){

		rmb = zero_mb ;

 		cout << "number of iterations used in zerosec to locate the marginally bound orbit : " 
                     << niter << endl ;
                cout << "    zero found at xi = " << xi_mb << endl ;
                cout << "    radial distance: " << rmb << endl ; 
        	}
	} 
 
    }	

    p_r_mb = new double (rmb) ;

    //delete [] zeros ; not used, causes "core dump" in Code kerr_qi
    
    }  // End of computation

    
    return *p_r_mb ;
   
}







//=============================================================================
//	Function used to locate the MS orbit
//=============================================================================


double funct_compobj_QI_isco(double xi, const Param& par){

    // Retrieval of the information:
    int l_ms = par.get_int() ;
    const Scalar& orbit = par.get_scalar() ;
    const Valeur& vorbit = orbit.get_spectral_va() ;

    // Value et the desired point
    double theta = M_PI / 2. ;
    double phi = 0 ;
    return vorbit.val_point(l_ms, xi, theta, phi) ;

}



//=============================================================================
//	Function used to locate the MB orbit
//=============================================================================

double funct_compobj_QI_rmb(double zeros, const Param& par){

    // Retrieval of the information:
    int l_mb = par.get_int() ;
    const Scalar& orbit = par.get_scalar() ;
    const Valeur& vorbit = orbit.get_spectral_va() ;

    // Value et the desired point
    double theta = M_PI / 2. ;
    double phi = 0 ;
    return vorbit.val_point(l_mb, zeros, theta, phi) ;

}




