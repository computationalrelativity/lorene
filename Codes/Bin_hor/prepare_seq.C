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

 

/*
 * $Id$
 * $Log$
 * Revision 1.11  2016/12/05 16:18:22  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.10  2014/10/13 08:53:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.9  2014/10/06 15:09:42  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.8  2007/04/13 15:30:58  f_limousin
 * Lots of improvements, generalisation to an arbitrary state of
 * rotation, implementation of the spatial metric given by Samaya.
 *
 * Revision 1.6  2006/06/29 08:54:52  f_limousin
 * Boundary conditions and grid writen in resformat.dat
 *
 * Revision 1.5  2006/06/28 13:36:52  f_limousin
 * Convergence to a given irreductible mass
 *
 * Revision 1.4  2006/05/24 16:59:08  f_limousin
 * New version
 *
 * Revision 1.2  2005/06/09 16:17:21  f_limousin
 * Many different changes.
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
    
    ofstream outfich("sequence.dat") ; 

    outfich << "# Beta  omega  M_ADM  M_Komar  M_area  J_ADM  J_hor  M_ADM/M_area  J_ADM/M_area2  omega * M_area  M_ADM/M_ih  J_ADM/M_ih2  omega * M_ih  hor-adm   adm-smarr   hor-smarr "  << endl ; 

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
	    double mass_ih1, mass_ih2, mass_ih, j1, j2, omega1, omega2 ;
	    double hor_adm, adm_smarr, hor_smarr ;

	    infich.ignore(1000,'\n') ; // skip first line
	    infich.ignore(1000,'\n') ; // skip second line
            infich.ignore(1000,'\n') ; // skip third line
	    infich >> beta ;
	    infich >> omega ;
	    infich >> mass_adm ;
	    infich >> mass_komar ;
	    infich >> mass_area ;
	    infich >> j_adm ;
	    infich >> j_hor ;
	    infich.ignore(1000,'\n') ; 
	    infich.ignore(1000,'\n') ; 
	    infich >> mass_ih1 ;
	    infich >> mass_ih2 ;
	    infich >> mass_ih ;
	    infich >> j1 ;
	    infich >> j2 ;
	    infich >> omega1 ;
	    infich >> omega2 ;
	    infich.ignore(1000,'\n') ; 
	    infich.ignore(1000,'\n') ; 
	    infich >> madm_area ;
	    infich >> jadm_area2 ;
	    infich >> omega_marea ;
	    infich.ignore(1000,'\n') ;
            infich.ignore(1000,'\n') ;
	    infich >> hor_adm ;
            infich >> adm_smarr ;
            infich >> hor_smarr ;

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
	    outfich << omega_marea << " "  ;
	    outfich << mass_adm/mass_ih << " " ;
	    outfich << j_adm/mass_ih/mass_ih << " " ;
	    outfich << omega*mass_ih << " " ;
	    outfich << hor_adm << " " ;
            outfich << adm_smarr << " " ;
            outfich << hor_smarr << endl ;

	    //-----------------------------------------------------
	    //  ***** End of modifiable part ***** 
	    //-----------------------------------------------------
	}
	
	fichlist.getline(filename, 120) ; 	// next file
	
    }
    
    outfich.close() ; 
    
    return 0 ; 
}
