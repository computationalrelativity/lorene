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
 * Revision 1.3  2014/10/06 15:09:45  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2005/09/12 12:34:09  f_limousin
 * Compilation Warning - Change of convention for the angular velocity
 * Add Berlin boundary condition in the case of binary horizons.
 *
 * Revision 1.1  2005/04/06 14:44:45  f_limousin
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
#include <cmath>


int main() {

    ifstream fichlist("list_prepare_seq.dat") ; 
    if ( !fichlist.good() ) {
	cout << "Problem with opening the file list_prepare_seq.d ! " << endl ;
	abort() ;
    }
    
    ofstream outfich("sequence.dat") ; 

    outfich << "# M_ADM   J_ADM   J/M    J/M2   Eps_a   M_hor   J_hor  r_hor  omega_hor"  << endl ; 

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
	    
	    double ham_constr, mom_constr, r_hor, j_hor, m_hor ;
	    double kappa_hor, omega_hor, m_adm, j_adm, jsm, jsm2 ;
	    double eps_a, diff_m, diff_j ;

	    infich.ignore(1000,'\n') ; // skip first line
	    infich >> ham_constr ;
	    infich >> mom_constr ;
	    infich.ignore(1000,'\n') ;  
	    infich.ignore(1000,'\n') ; 
	    infich >> r_hor ;
	    infich >> j_hor ;
	    infich >> m_hor ;
	    infich >> kappa_hor ;
	    infich >> omega_hor ;
	    infich.ignore(1000,'\n') ;  
	    infich.ignore(1000,'\n') ; 
	    infich >> m_adm ;
	    infich >> j_adm ;
	    infich >> jsm ;
	    infich >> jsm2 ;
	    infich >> eps_a ;
	    infich.ignore(1000,'\n') ;  
	    infich.ignore(1000,'\n') ; 
	    infich >> diff_m ;
	    infich >> diff_j ;
	    
	    infich.close() ; 

	    // ---------------------------------------------------
	    //  End of file reading
	    // ---------------------------------------------------

	    cout.precision(8);
	    cout << "m_adm = " << ham_constr
		 << "  j_adm = " << j_adm 
		 << "  J/M = " << jsm
		 << "  J/M2 = " << jsm2 << endl 
		 << "eps_a = " << eps_a
		 << "  j_hor = " << j_hor
		 << "  m_hor = " << m_hor << endl ;


	    //-----------------------------------------------------
	    //  ***** Start of modifiable part ******
	    //-----------------------------------------------------

	    outfich.setf(ios::scientific) ; 
//	    outfich.width(23) ;  
	    outfich.precision(8) ;
	    outfich << m_adm << " " ;
	    outfich << j_adm << " " ;
	    outfich << jsm << " " ;
	    outfich << jsm2 << " " ;
	    outfich << eps_a << " " ;
	    outfich << m_hor << " " ;
	    outfich << j_hor << " " ;
	    outfich << r_hor << " " ;
	    outfich << omega_hor << endl ;

	    //-----------------------------------------------------
	    //  ***** End of 2nd modifiable part ***** 
	    //-----------------------------------------------------
	}
	
	fichlist.getline(filename, 120) ; 	// next file
	
    }
    
    outfich.close() ; 
    
    return 0 ; 
}
