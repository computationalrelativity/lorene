/*
 *  save_profile function
 *
 *    (see file graphique.h for documentation).
 *
 */

/*
 *   Copyright (c) 2011  Eric Gourgoulhon 
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

char name_of_this_file_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2011/03/27 16:36:41  e_gourgoulhon
 * New function save_profile.
 *
 * Revision 1.4  2003/10/19 20:01:10  e_gourgoulhon
 * Template file
 *
 * $Header$
 *
 */

// C++ headers
#include <fstream>

// Lorene headers
#include "scalar.h"

void save_profile(const Scalar& uu, double r_min, double r_max, 
		     double theta, double phi, const char* filename) {
  
    const int npt = 400 ;   // Number of points along the axis
        
    double hr = (r_max - r_min) / double(npt-1) ; 
    
    ofstream file(filename) ;

    for (int i=0; i<npt; i++) {
    
	double r = hr * i + r_min ; 
	
	file << r << "  " << uu.val_point(r, theta, phi) << endl ; 
    }
    
    file.close() ; 
    
} 

