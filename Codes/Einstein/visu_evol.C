/*
 *  Main code for reading a series of time slices Sigma_t stored in files. 
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon & Jerome Novak
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

char visu_evol_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/06/02 21:34:39  e_gourgoulhon
 * Added creation of file anime.dxcont for OpenDX.
 *
 * Revision 1.1  2004/05/31 20:34:20  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

// Lorene headers
#include "time_slice.h"
#include "param.h"
#include "utilitaires.h"
#include "graphique.h"

int main(){
    
    ifstream fpar("par_visu_evol.d") ;
    if ( !fpar.good() ) {
        cout << "Problem with opening the file par_visu_evol.d ! " << endl ;
        abort() ;
    }
    
    char rootname[80], rootdxname[80], section_type ; 
    int jmin, jmax, jstep, nu, nv ;
    double a_section, umin, umax, vmin, vmax ;  
    fpar.ignore(1000,'\n') ;    // skip title
    fpar.ignore(1000,'\n') ;    // skip comment
    fpar.getline(rootname, 80) ;
    fpar >> jmin ; fpar.ignore(1000,'\n') ;
    fpar >> jmax ; fpar.ignore(1000,'\n') ;
    fpar >> jstep ; fpar.ignore(1000,'\n') ;
    fpar.ignore(1000,'\n') ;    // skip comment
    fpar.getline(rootdxname, 80) ;
    fpar.get(section_type) ; fpar.ignore(1000,'\n') ;
    fpar >> a_section ; fpar.ignore(1000,'\n') ;
    fpar >> umin ;  fpar >> umax ; fpar.ignore(1000,'\n') ;
    fpar >> vmin ;  fpar >> vmax ; fpar.ignore(1000,'\n') ;
    fpar >> nu ; fpar.ignore(1000,'\n') ;
    fpar >> nv ; fpar.ignore(1000,'\n') ;

    cout << "Root name of files to be read : " << rootname << endl ;         
    cout << "jmin = " << jmin << endl ; 
    cout << "jmax = " << jmax << endl ; 
    cout << "jstep = " << jstep << endl ; 
    cout << "section_type = " << section_type << endl ; 
    cout << "a_section = " << a_section << endl ; 
    cout << "umin, umax = " << umin << ", " << umax << endl ; 
    cout << "vmin, vmax = " << vmin << ", " << vmax << endl ; 
    cout << "nu x nv = " << nu << " x " << nv << endl ; 
    arrete() ; 
    
    ofstream fdx("anime.dxcont") ; 
    fdx << "fieldname     jmin   jmax  jstep  ampli " << endl ; 
    fdx << rootdxname << "     " << jmin << "   " << jmax << "   " << jstep ;
         
    for (int j=jmin; j<=jmax; j += jstep) {         
 
        char* filename = new char[ strlen(rootname)+10 ] ; 
        strcpy(filename, rootname) ; 
        char nomj[7] ; 
        sprintf(nomj, "%06d", j) ; 
        strcat(filename, nomj) ; 
        strcat(filename, ".d") ; 
        
        
        FILE* fich = fopen(filename, "r") ; 
        if (fich == 0x0) {
    	cout << "Problem in opening the file " << filename << " ! " << endl ; 
		perror(" reason") ; 
		abort() ; 
        }
    
        Mg3d mgrid(fich) ;
        Map_af map(mgrid, fich) ;
    
        Base_vect* ptriad = Base_vect::bvect_from_file(fich) ; 
        
        // Flat metric f 
        // -------------

        const Metric_flat& ff = map.flat_met_spher() ; 
    
        int depth ; 
	fread_be(&depth, sizeof(int), 1, fich) ;	
    
        Tslice_dirac_max sigma(map, *ptriad, ff, fich, false, depth) ; 
         
        assert(sigma.get_latest_j() == j) ; 
        
        double tc = sigma.get_time()[j] ;     

        cout << 
        "==============================================================\n"
        << "  step: " << j << "   time = " << tc << endl  
        << "==============================================================\n" ;
    

        cout << sigma << endl ; 
        
        bool start_dx = ( (j >= jmax - jstep) && (j != jmax) ) ;
        // bool start_dx = false ;
         
        // sigma.mu().spectral_display("mu") ; 
         
        sigma.mu().visu_section_anim(section_type, a_section, umin, umax,
		                      vmin, vmax, j, tc, 1, "mu", rootdxname, 
                                      start_dx, nu, nv) ; 
                                      
        if (j==jmin) {
            double ampli = 1. / max(maxabs(sigma.mu())) ; 
            fdx << "   " << ampli << endl ; 
            fdx.close() ; 
        }

        delete ptriad ; 

    }       
    
    return EXIT_SUCCESS ; 
}
