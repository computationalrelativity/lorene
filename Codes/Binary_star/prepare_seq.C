/*
 *  Prepare a file describing a sequence from a list of resformat.d files
 *  
 */

/*
 *   Copyright (c) 2003  Eric Gourgoulhon
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

char prepare_seq_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2005/09/13 19:47:28  f_limousin
 * Reintroduction of the resolution of the equations in cartesian coordinates.
 *
 * Revision 1.1  2004/09/16 12:15:01  f_limousin
 * *** empty log message ***
 *
 * Revision 1.1  2003/09/16 09:15:33  e_gourgoulhon
 * First version.
 *
 *
 *
 * $Header$
 *
 */

// C++ headers
#include <iostream>
#include <fstream>
using namespace std ; 

// C headers
#include <math.h>
#include "unites.h"


#define maxval(a,b) a < b ? b : a
#define minval(a,b) a > b ? b : a

int main() {

    using namespace Unites ;

    // For the display :
    char display_bold[]="x[1m" ; display_bold[0] = 27 ;
    char display_normal[] = "x[0m" ; display_normal[0] = 27 ;

	ifstream fichlist("list_prepare_seq.d") ; 
	if ( !fichlist.good() ) {
		cout << "Problem with opening the file list_prepare_seq.d ! " << endl ;
		abort() ;
	}

	ofstream outfich("sequence.d") ; 

	//-----------------------------------------------------
	//  ***** Start of 1st modifiable part ******
	//-----------------------------------------------------

	double m0_half ;
	cout << 
	"Gravitational mass of a single star at r->oo [in sollar mass units]? " << endl ;
	cin >> m0_half ; 
	double m0 = 2 * m0_half ;
	
	outfich << "# d [km]  f_GW [Hz]   M_ADM [Msol]   Omega*M_inf  M_B [Msol]   v_m   a2/a1    a3/a1    j/m0^2"  << endl ; 
	
	//-----------------------------------------------------
	//  ***** End of 1st modifiable part ******
	//-----------------------------------------------------

	//-------------------------------------------------------
	//  Loop of the files seq.d
	//-------------------------------------------------------
	
	char filename[120] ; 
	fichlist.getline(filename, 120) ; 

	while( !fichlist.eof() ) {

		cout << filename << endl ; 
		if ( (filename[0] != ' ') && (filename[0] != '#') ) {

			ifstream infich(filename) ;
			if ( !infich.good() ) {
		    	cout << "Problem with opening the file " << filename << " ! " << endl ;
		    	abort() ;
			}
	
			double ve_m ;
			double d_km, d_g_km, d_rel, f_hz, m_adm, m_adm_vol ;
			double m_kom, m_kom_vol, jj ; 
			double entc1, enerc1, m_b1, r_eq1, a2sa11, a3sa11 ; 
			double entc2, enerc2, m_b2, r_eq2, a2sa12, a3sa12 ; 
	
			infich.ignore(1000,'\n') ; // skip first two lines
			infich.ignore(1000,'\n') ; //
			infich >> ve_m ;
			infich.ignore(1000,'\n') ;  
			infich.ignore(1000,'\n') ;  
			infich >> d_km ; 
			infich >> d_g_km ;
			infich >> d_rel ;
			infich >>  f_hz ;
			infich >>  m_adm ;
			infich >>  m_adm_vol ;
			infich >>  m_kom ;
			infich >>  m_kom_vol ;
			infich >>  jj ;
			infich.ignore(1000,'\n') ;
			infich.ignore(1000,'\n') ;
			infich >> entc1 ;
			infich >>  enerc1 ;
			infich >>  m_b1 ;
			infich >>  r_eq1 ;
			infich >>  a2sa11 ;
			infich >>  a3sa11 ; 
			infich.ignore(1000,'\n') ;
			infich.ignore(1000,'\n') ;
			infich >> entc2 ;
			infich >>  enerc2 ;
			infich >>  m_b2 ;
			infich >>  r_eq2 ;
			infich >>  a2sa12 ;
			infich >>  a3sa12 ; 
			infich.ignore(1000,'\n') ;

			infich.close() ; 
			// ---------------------------------------------------
			//  End of file reading
			// ---------------------------------------------------

			    cout.precision(8);
			cout << "d = " << d_km << display_bold 
                             << " f_gw = " << 2*f_hz
                             << " m_adm = " << m_adm
                             << " m_b1 = " << m_b1 
                             << display_normal << " j = " << jj ; 
			    cout.precision(8) ;             
			    
			//-----------------------------------------------------
			//  ***** Start of 2nd modifiable part ******
			//-----------------------------------------------------

/* 			// double m0 = 2 * 2.347359 ; 
			// double gm0 = 6.6726E-01 / pow(2.99792458E+8, 3) *
			// 			1.989E+20 * m0 ;
			// double om1 = 2.* M_PI * f_hz * gm0 ; 
			// double ebind = m_adm / m0 - 1. ;

                                                

			double om1 = omega_poly * m0 ; 
			double ebind = m_adm_poly / m0 - 1. ;
			double jj1 = jj_poly / (m0*m0) ; 
*/
			    
			    double kappa = 0.0332 ;
			    double gamma = 2. ;

			double om1 = f_hz * 2 * M_PI ; 
			double jj1 = jj / (m0*m0) ; 
//                      dM=m_adm-Mg(f=0)+4.056e-4*x^(2/3)
			outfich.setf(ios::scientific) ; 
                        outfich.precision(4) ; 
			outfich << d_km  << " ";
			outfich.precision(12) ;  
                        outfich << f_hz  << " "; 
                        outfich << (m_adm_vol - m0)  << " "; 
			outfich << om1/f_unit*m0*msol*ggrav  << " "; 
			outfich << m_b1  << " ";
                        outfich.precision(5) ;  
			outfich << ve_m << " "; 
			outfich << a2sa11 << " "; 
			outfich << a3sa11 << " "; 
			outfich << jj1  << endl ; 

			//-----------------------------------------------------
			//  ***** End of 2nd modifiable part ***** 
			//-----------------------------------------------------
		}

		fichlist.getline(filename, 120) ; 	// next file

	}
	
	outfich.close() ; 
		
	return 0 ; 
}
