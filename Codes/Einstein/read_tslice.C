/*
 *  Main code for reading a time slice Sigma_t stored in file. 
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

char read_tslice_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/05/27 15:25:38  e_gourgoulhon
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

// Lorene headers
#include "time_slice.h"
#include "param.h"
#include "utilitaires.h"
#include "graphique.h"

int main(int argc, char** argv){
    
    if (argc < 2) {
		cout << 
		"read_tslice : the name of a file containing a configuration"
		<< endl << " must be given in argument !" << endl ; 
		abort() ; 
    }
    
    char* nomresu = argv[1] ; 
    cout << "Name of the file to be read : " << nomresu << endl ;         
    
    FILE* fich = fopen(nomresu, "r") ; 
    if (fich == 0x0) {
    	cout << "Problem in opening the file " << nomresu << " ! " << endl ; 
		perror(" reason") ; 
		abort() ; 
    }
    
    Mg3d mgrid(fich) ;
    Map_af map(mgrid, fich) ;
    
    Base_vect* ptriad = Base_vect::bvect_from_file(fich) ; 
    
    cout << "Computational grid :\n" 
         << "------------------ \n" 
         << "  " << mgrid << endl ; 

    cout << "Mapping computational grid --> physical space :\n" 
         << "---------------------------------------------\n" 
         << "  " << map << endl ;  
    
    // Flat metric f
    // -------------

    const Metric_flat& ff = map.flat_met_spher() ; 
    
    int depth ; 
	fread_be(&depth, sizeof(int), 1, fich) ;	
    
    Tslice_dirac_max sigmat(map, *ptriad, ff, fich, false, depth) ;      

    cout << sigmat << endl ; 


    return EXIT_SUCCESS ; 
}
