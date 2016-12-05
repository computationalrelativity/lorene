/*
 * Methods of class Bin_BH (binary black hole exportation)
 *
 * (see file bin_bh.h for documentation).
 */

/*
 *   Copyright (c) 2001  Eric Gourgoulhon
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
 * Revision 1.5  2016/12/05 16:18:30  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:54:05  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2006/09/12 08:04:05  j_novak
 * Removal of the include path Export/C++/Include, updating of the relevant
 * source files in Export/C++/Source.
 *
 * Revision 1.2  2002/03/20 08:24:56  e_gourgoulhon
 * Added the derivatives of Psi.
 *
 * Revision 1.1  2001/12/18 22:27:04  e_gourgoulhon
 * Exportation of Lorene structures
 *
 *
 *
 * $Header$
 *
 */


#include "../Include/bin_bh.h"

		    //-----------------------------------------//
		    //	    Constructor from a binary file     //
		    //-----------------------------------------//

namespace Lorene {
Bin_BH::Bin_BH(FILE* fich) {
    
 	fread(&omega, sizeof(double), 1, fich) ;	
 	fread(&dist, sizeof(double), 1, fich) ;	
 	fread(&radius2, sizeof(double), 1, fich) ;	
 	fread(&np, sizeof(int), 1, fich) ;
	
	alloc_memory() ; 
	
	fread(xx, sizeof(double), np, fich) ; 
	fread(yy, sizeof(double), np, fich) ; 
	fread(zz, sizeof(double), np, fich) ; 
	fread(nnn, sizeof(double), np, fich) ; 
	fread(beta_x, sizeof(double), np, fich) ; 
	fread(beta_y, sizeof(double), np, fich) ; 
	fread(beta_z, sizeof(double), np, fich) ; 
	fread(g_xx, sizeof(double), np, fich) ; 
	fread(g_xy, sizeof(double), np, fich) ; 
	fread(g_xz, sizeof(double), np, fich) ; 
	fread(g_yy, sizeof(double), np, fich) ; 
	fread(g_yz, sizeof(double), np, fich) ; 
	fread(g_zz, sizeof(double), np, fich) ; 
	fread(k_xx, sizeof(double), np, fich) ; 
	fread(k_xy, sizeof(double), np, fich) ; 
	fread(k_xz, sizeof(double), np, fich) ; 
	fread(k_yy, sizeof(double), np, fich) ; 
	fread(k_yz, sizeof(double), np, fich) ; 
	fread(k_zz, sizeof(double), np, fich) ;
	fread(dpsi_x, sizeof(double), np, fich) ;
	fread(dpsi_y, sizeof(double), np, fich) ;
	fread(dpsi_z, sizeof(double), np, fich) ;
	fread(d2psi_xx, sizeof(double), np, fich) ;
	fread(d2psi_xy, sizeof(double), np, fich) ;
	fread(d2psi_xz, sizeof(double), np, fich) ;
	fread(d2psi_yy, sizeof(double), np, fich) ;
	fread(d2psi_yz, sizeof(double), np, fich) ;
	fread(d2psi_zz, sizeof(double), np, fich) ;

    
}

		    //--------------------------------------------//
		    //	    Constructor from a formatted file     //
		    //--------------------------------------------//

Bin_BH::Bin_BH(ifstream& fich) {
    

    char comment[120] ; 

    fich >> omega ; fich.getline(comment, 120) ; 
    fich >> dist ; fich.getline(comment, 120) ; 
    fich >> radius2 ; fich.getline(comment, 120) ; 
    fich >> np ; fich.getline(comment, 120) ; 

    alloc_memory() ; 
	
    double* field[] = {xx, yy, zz, nnn, beta_x, beta_y, beta_z, 
			     g_xx,  g_xy, g_xz, g_yy, g_yz, g_zz, 
			     k_xx,  k_xy, k_xz, k_yy, k_yz, k_zz,
			     dpsi_x, dpsi_y, dpsi_z,
			     d2psi_xx, d2psi_xy, d2psi_xz,
			     d2psi_yy, d2psi_yz, d2psi_zz} ;
			        
    for (int j=0; j<28; j++) {

	double* pdata = field[j] ;
	
	for (int i=0; i<np; i++) {
 	    fich >> *pdata ;
	    pdata++ ;  
	}	
    }
    
}




		    //--------------------------//
		    //	    Destructor		//
		    //--------------------------//


Bin_BH::~Bin_BH() {

    delete [] xx ; 
    delete [] yy ; 
    delete [] zz ; 
   
    delete [] nnn ;

    delete [] beta_x ;
    delete [] beta_y ;
    delete [] beta_z ;

    delete [] g_xx ;
    delete [] g_xy ;
    delete [] g_xz ;
    delete [] g_yy ;
    delete [] g_yz ;
    delete [] g_zz ;

    delete [] k_xx ;
    delete [] k_xy ;
    delete [] k_xz ;
    delete [] k_yy ;
    delete [] k_yz ;
    delete [] k_zz ;

    delete [] dpsi_x ;
    delete [] dpsi_y ;
    delete [] dpsi_z ;

    delete [] d2psi_xx ;
    delete [] d2psi_xy ;
    delete [] d2psi_xz ;
    delete [] d2psi_yy ;
    delete [] d2psi_yz ;
    delete [] d2psi_zz ;

}


		    //--------------------------//
		    //	    Memory allocation 	//
		    //--------------------------//

void Bin_BH::alloc_memory() {
    
    xx = new double[np] ; 
    yy = new double[np] ; 
    zz = new double[np] ; 

    nnn = new double[np] ; 

    beta_x = new double[np] ; 
    beta_y = new double[np] ; 
    beta_z = new double[np] ; 
    
    g_xx = new double[np] ; 
    g_xy = new double[np] ; 
    g_xz = new double[np] ; 
    g_yy = new double[np] ; 
    g_yz = new double[np] ; 
    g_zz = new double[np] ; 
    
    k_xx = new double[np] ; 
    k_xy = new double[np] ; 
    k_xz = new double[np] ; 
    k_yy = new double[np] ; 
    k_yz = new double[np] ; 
    k_zz = new double[np] ;

    dpsi_x = new double[np] ;
    dpsi_y = new double[np] ;
    dpsi_z = new double[np] ;

    d2psi_xx = new double[np] ;
    d2psi_xy = new double[np] ;
    d2psi_xz = new double[np] ;
    d2psi_yy = new double[np] ;
    d2psi_yz = new double[np] ;
    d2psi_zz = new double[np] ;
}

		    //--------------------------//
		    //	  	Outputs         //
		    //--------------------------//
		
ostream& operator<<(ostream& ost, const Bin_BH& bb) {

	ost << "Binary black hole system :" << endl ;
	ost << "--------------------------" << endl ;
	ost << "  Separation d/a :       " << bb.dist << endl ;
	ost << "  Omega :                " << bb.omega << " / a" << endl ;
    	ost << "  Size of black hole 2 : " << bb.radius2 << " a" << endl ;

	return ost ; 
}
		
		
		    //----------------------------------//
		    //	    Save in a binary file 	//
		    //----------------------------------//

void Bin_BH::save_bin(FILE* fresu) const {
    
 	fwrite(&omega, sizeof(double), 1, fresu) ;	
 	fwrite(&dist, sizeof(double), 1, fresu) ;	
 	fwrite(&radius2, sizeof(double), 1, fresu) ;	
 	fwrite(&np, sizeof(int), 1, fresu) ;	
	fwrite(xx, sizeof(double), np, fresu) ; 
	fwrite(yy, sizeof(double), np, fresu) ; 
	fwrite(zz, sizeof(double), np, fresu) ; 
	fwrite(nnn, sizeof(double), np, fresu) ; 
	fwrite(beta_x, sizeof(double), np, fresu) ; 
	fwrite(beta_y, sizeof(double), np, fresu) ; 
	fwrite(beta_z, sizeof(double), np, fresu) ; 
	fwrite(g_xx, sizeof(double), np, fresu) ; 
	fwrite(g_xy, sizeof(double), np, fresu) ; 
	fwrite(g_xz, sizeof(double), np, fresu) ; 
	fwrite(g_yy, sizeof(double), np, fresu) ; 
	fwrite(g_yz, sizeof(double), np, fresu) ; 
	fwrite(g_zz, sizeof(double), np, fresu) ; 
	fwrite(k_xx, sizeof(double), np, fresu) ; 
	fwrite(k_xy, sizeof(double), np, fresu) ; 
	fwrite(k_xz, sizeof(double), np, fresu) ; 
	fwrite(k_yy, sizeof(double), np, fresu) ; 
	fwrite(k_yz, sizeof(double), np, fresu) ; 
	fwrite(k_zz, sizeof(double), np, fresu) ;
	fwrite(dpsi_x, sizeof(double), np, fresu) ;
	fwrite(dpsi_y, sizeof(double), np, fresu) ;
	fwrite(dpsi_z, sizeof(double), np, fresu) ;
	fwrite(d2psi_xx, sizeof(double), np, fresu) ;
	fwrite(d2psi_xy, sizeof(double), np, fresu) ;
	fwrite(d2psi_xz, sizeof(double), np, fresu) ;
	fwrite(d2psi_yy, sizeof(double), np, fresu) ;
	fwrite(d2psi_yz, sizeof(double), np, fresu) ;
	fwrite(d2psi_zz, sizeof(double), np, fresu) ;

}

		
		    //----------------------------------//
		    //	    Save in a formatted file 	//
		    //----------------------------------//

void Bin_BH::save_form(ofstream& fich) const {

    fich.precision(13) ; 
    fich.setf(ios::scientific, ios::floatfield) ; 

    fich << omega << " Omega [a^{-1}]" << endl ; 
    fich << dist << " dist [a]" << endl ; 
    fich << radius2 << " radius2 [a]" << endl ; 
    fich << np << " Number of grid points" << endl ; 
    
    int dpl = 6 ;    // number of data per line 
    
    int nlines = np / dpl ;   // number of filled lines
    int reste = np - nlines * dpl ;	// number of remaining data
    
    const double* field[] = {xx, yy, zz, nnn, beta_x, beta_y, beta_z, 
			     g_xx,  g_xy, g_xz, g_yy, g_yz, g_zz, 
			     k_xx,  k_xy, k_xz, k_yy, k_yz, k_zz,
			     dpsi_x, dpsi_y, dpsi_z,
			     d2psi_xx, d2psi_xy, d2psi_xz,
			     d2psi_yy, d2psi_yz, d2psi_zz} ;

			        
    for (int j=0; j<28; j++) {

	const double* pdata = field[j] ;
    
	for (int line = 0; line < nlines; line++) {
	    for (int i=0; i<dpl; i++) {
		fich << *pdata << "  " ;
		pdata++ ;  
	    }
	    fich << endl ; 
	}
	for (int i=0; i<reste; i++) {
	    fich << *pdata << "  " ;
	    pdata++ ;  
	}
	fich << endl ; 
	
    }

}

}
