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

 

/*
 * $Id$
 * $Log$
 * Revision 1.4  2016/12/05 16:18:07  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.3  2016/09/19 15:26:23  j_novak
 * Correction of several bugs preventing the shared library compilation.
 *
 * Revision 1.2  2014/10/13 08:53:23  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
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

namespace Lorene {
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

}
