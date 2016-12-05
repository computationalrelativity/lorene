/*
 * Main code for testing the class Eos_Fermi 
 * 
 */

/*
 *   Copyright (c) 2012 Eric Gourgoulhon
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
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
 * Revision 1.3  2016/12/05 16:18:27  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.2  2014/10/13 08:54:00  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1  2012/10/26 14:10:04  e_gourgoulhon
 * Simple code to test the class Eos
 *
 *
 * $Header$
 *
 */

// C headers
#include <cstdlib>
#include <cmath>
#include <cstring>

// Lorene headers
#include "eos.h"
#include "utilitaires.h"


//******************************************************************************

using namespace Lorene ;

int main(){

    Eos_Fermi eos(5.11e5) ; 
    
    cout << eos << endl ; 
    
    double hmin = 0. ; 
    double hmax = 4 ; 
    int np = 100 ; 
    double dh = (hmax - hmin) / double(np-1) ; 
    for (int i = 0; i<np; i++) {
        double hh = hmin + i * dh ; 
        double n = eos.nbar_ent_p(hh) ; 
        double ener = eos.ener_ent_p(hh) ; 
        double press = eos.press_ent_p(hh) ; 
        cout << hh << "   " << n << "   " << ener << "   " << press << endl ; 
    }
}
