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
 * Revision 1.3  2014/10/13 08:53:54  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.2  2014/10/06 15:09:43  j_novak
 * Modified #include directives to use c++ syntax.
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
#include <cmath>
#include <cstdlib>

int main() {

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
	"Value of the gravitational mass of a single star at infinity [poly. unit] ? " << endl ;
	cin >> m0_half ; 
	double m0 = 2 * m0_half ;
	
	outfich << "#   M_0 Omega           (M_ADM - M_0)/M_0       J/M_0^2" << endl ; 
	
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
	
			double ve_m, ve_gb, ve_fus ; 
			double d_km, d_g_km, d_rel, f_hz, m_adm, jj ; 
			double entc1, enerc1, m_b1, r_eq1, a2sa11, a3sa11 ; 
			double entc2, enerc2, m_b2, r_eq2, a2sa12, a3sa12 ; 
			double d_poly, d_g_poly, omega_poly, m_adm_poly, jj_poly, 
					m_b1_poly, m_b2_poly ; 

			infich.ignore(1000,'\n') ; // skip first two lines
			infich.ignore(1000,'\n') ; //
			infich >> ve_m ;
			infich >>  ve_gb ;
			infich >>  ve_fus ;
			infich.ignore(1000,'\n') ;  
			infich.ignore(1000,'\n') ;  
			infich >> d_km ; 
			infich >> d_g_km ;
			infich >> d_rel ;
			infich >>  f_hz ;
			infich >>  m_adm ;
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

			if (infich.peek() == '#') {	// There exist a line with
									    // polytropic units 
				infich.ignore(1000,'\n') ;  
				infich >> d_poly ; 
				infich >> d_g_poly ;
				infich >> omega_poly ; 
				infich >> m_adm_poly ;
				infich >> jj_poly ;
				infich >> m_b1_poly ; 
				infich >> m_b2_poly ; 			
			}
			else {
				d_poly = 0 ; 
				d_g_poly = 0 ;
				omega_poly = 0 ;
				m_adm_poly = 0 ;
				jj_poly = 0 ;
				m_b1_poly = 0 ;
				m_b2_poly = 0 ;
			}

			infich.close() ; 
			// ---------------------------------------------------
			//  End of file reading
			// ---------------------------------------------------

			cout << "d = " << d_km << " f = " << f_hz << " m_b1 = " 
				<< m_b1 << endl ; 
			cout << "d_poly = " << d_poly << " omega_poly = " 
				 << omega_poly << " m_adm_poly = " << m_adm_poly << endl ;


			//-----------------------------------------------------
			//  ***** Start of 2nd modifiable part ******
			//-----------------------------------------------------

			// double m0 = 2 * 2.347359 ; 
			// double gm0 = 6.6726E-01 / pow(2.99792458E+8, 3) *
			// 			1.989E+20 * m0 ;
			// double om1 = 2.* M_PI * f_hz * gm0 ; 
			// double ebind = m_adm / m0 - 1. ;


			double om1 = omega_poly * m0 ; 
			double ebind = m_adm_poly / m0 - 1. ;
			double jj1 = jj_poly / (m0*m0) ; 

			outfich.precision(14) ;  
			outfich.setf(ios::scientific) ; 
			outfich.width(20) ;  
			outfich << om1 ; outfich.width(23) ; 
			outfich << ebind ; outfich.width(23) ; 
			outfich << jj1 << endl ; 

			//-----------------------------------------------------
			//  ***** End of 2nd modifiable part ***** 
			//-----------------------------------------------------
		}

		fichlist.getline(filename, 120) ; 	// next file

	}
	
	outfich.close() ; 
		
	return 0 ; 
}
