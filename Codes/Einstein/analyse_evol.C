/*
 *  Main code for reading a series of time slices Sigma_t stored in files
 *   and performoning some analysis of it
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

char analyse_evol_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/06/24 20:37:44  e_gourgoulhon
 * First version.
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
    
    ifstream fpar("par_ana_evol.d") ;
    if ( !fpar.good() ) {
        cout << "Problem with opening the file par_ana_evol.d ! " << endl ;
        abort() ;
    }
    
    char rootname[80] ; 
    int jmin, jmax, jstep ;
    fpar.ignore(1000,'\n') ;    // skip title
    fpar.ignore(1000,'\n') ;    // skip comment
    fpar.getline(rootname, 80) ;
    fpar >> jmin ; fpar.ignore(1000,'\n') ;
    fpar >> jmax ; fpar.ignore(1000,'\n') ;
    fpar >> jstep ; fpar.ignore(1000,'\n') ;

    cout << "Root name of files to be read : " << rootname << endl ;         
    cout << "jmin = " << jmin << endl ; 
    cout << "jmax = " << jmax << endl ; 
    cout << "jstep = " << jstep << endl << endl ; 
    
    ofstream fresu1("psidot.d") ; 
    ofstream fresu2("psidot_diff.d") ; 
    ofstream fresu3("psidot_diff_rel.d") ; 
    
    fresu1 << "#  t           max(abs(d/dt ln Psi)) in each domain" << endl ; 
    fresu2 << 
    "#  t            Absolute error on the d/dt ln Psi relation in each domain" 
    << endl ; 
    fresu3 << 
    "#  t            Relative error on the d/dt ln Psi relation in each domain" 
    << endl ; 
                 
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
    
        Base_vect::bvect_from_file(fich) ;  //##  skip the triad stored in file 
        const Base_vect* ptriad = &(map.get_bvect_spher()) ; 
        
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
        
        int nz = mgrid.get_nzone() ; 
        Tbl tlnpsi_dot(nz) ; 
        Tbl tdiff(nz) ; 
        Tbl tdiff_rel(nz) ; 
        sigma.check_psi_dot(tlnpsi_dot, tdiff, tdiff_rel) ; 
        
        fresu1 << tc ; 
        for (int l = 0; l<nz ; l++) {
            fresu1 << "  " << tlnpsi_dot(l) ; 
        }
        fresu1 << endl ; 
        
        fresu2 << tc ; 
        for (int l = 0; l<nz ; l++) {
            fresu2 << "  " << tdiff(l) ; 
        }
        fresu2 << endl ; 
        
        fresu3 << tc ; 
        for (int l = 0; l<nz ; l++) {
            fresu3 << "  " << tdiff_rel(l) ; 
        }
        fresu3 << endl ; 
        
        delete ptriad ; 

    }   
    
    fresu1.close() ; 
    fresu2.close() ; 
    fresu3.close() ; 
    
    return EXIT_SUCCESS ; 
}
