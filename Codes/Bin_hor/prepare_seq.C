/*
 *  Prepare a file describing a sequence from a list of resformat.d files
 *  
 */

/*
 *   Copyright (c) 2005  Francois Limousin
 *                       Jose Luis Jaramillo
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
 * Revision 1.1  2005/04/29 14:08:46  f_limousin
 * First version
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


int main() {

    ifstream fichlist("list_prepare_seq.d") ; 
    if ( !fichlist.good() ) {
	cout << "Problem with opening the file list_prepare_seq.d ! " << endl ;
	abort() ;
    }
    
    ofstream outfich("sequence.dat") ; 

    outfich << "# Beta  omega  M_ADM  M_Komar  M_area  J_ADM  J_hor  M_ADM/M_area  J_ADM/M_area2  omega * M_area"  << endl ; 

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
		cout << "Problem with opening the file " << filename << " ! " 
		     << endl ;
		abort() ;
	    }
	    
	    double beta, omega, mass_adm, mass_area, mass_komar ; 
	    double j_adm, j_hor, madm_area, jadm_area2, omega_marea ;

	    infich.ignore(1000,'\n') ; // skip first line
	    infich >> beta ;
	    infich >> omega ;
	    infich >> mass_adm ;
	    infich >> mass_komar ;
	    infich >> mass_area ;
	    infich >> j_adm ;
	    infich >> j_hor ;
	    infich.ignore(1000,'\n') ; 
	    infich.ignore(1000,'\n') ; 
	    infich >> madm_area ;
	    infich >> jadm_area2 ;
	    infich >> omega_marea ;
	    
	    infich.close() ; 

	    // ---------------------------------------------------
	    //  End of file reading
	    // ---------------------------------------------------

	    cout.precision(8);
	    cout << "m_adm = " << mass_adm
		 << "  j_adm = " << j_adm 
		 << "  omega = " << omega
		 << "  beta = " << beta ;


	    //-----------------------------------------------------
	    //  ***** Start of modifiable part ******
	    //-----------------------------------------------------

	    outfich.setf(ios::scientific) ; 
	    outfich.precision(8) ;
	    outfich << beta << " " ;
	    outfich << omega << " " ;
	    outfich << mass_adm << " " ;
	    outfich << mass_komar << " " ;
	    outfich << mass_area << " " ;
	    outfich << j_adm << " " ;
	    outfich << j_hor << " " ;
	    outfich << madm_area << " " ;
	    outfich << jadm_area2 << " " ;
	    outfich << omega_marea << endl ;

	    //-----------------------------------------------------
	    //  ***** End of modifiable part ***** 
	    //-----------------------------------------------------
	}
	
	fichlist.getline(filename, 120) ; 	// next file
	
    }
    
    outfich.close() ; 
    
    return 0 ; 
}
