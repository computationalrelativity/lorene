/*
 *  Formatted file output for double's and Tbl's.
 *
 *
 */

/*
 *   Copyright (c) 2004  Eric Gourgoulhon
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

char write_formatted_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2004/05/13 21:31:06  e_gourgoulhon
 * First version.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// Lorene headers
#include "tbl.h"

// double version 
// --------------
void write_formatted(const double& x, ostream& ost) { 

    ost.width(23) ; ost << x ; 
    
}


// Tbl version
// -----------
void write_formatted(const Tbl& tb, ostream& ost) { 

    assert(tb.get_ndim() == 1) ; 
    
    for (int i=0; i<tb.get_taille(); i++) {
        ost.width(23) ; ost << tb(i) ; 
    }
    
}
