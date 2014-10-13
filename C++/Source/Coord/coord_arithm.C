/*
 *  Arithmetical operations for class Coord
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char coord_arithm_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2014/10/13 08:52:50  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/02/25  10:24:40  eric
 * Remplacement de la variable globale nom_C (!) par arithm_coord_C
 *
 * Revision 1.1  1999/10/15  13:57:58  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "coord.h"
#include "mtbl.h"

namespace Lorene {

/************************************************************************/
/*			operations sur Coord -> Mtbl			*/
/************************************************************************/

			//********************//
			// OPERATEURS UNAIRES //
			//********************//
			
Mtbl operator+(const Coord& co) {
    
    if (co.c == 0x0) co.fait() ; 
    return *(co.c) ;
    
}			

Mtbl operator-(const Coord& co) {
    
    if (co.c == 0x0) co.fait() ; 
    return -(*(co.c)) ;
    
}			

			//**********//
			// ADDITION //
			//**********//

Mtbl operator+(const Coord& c1, const Coord& c2) {

    // Sont-elles a jour ?
    if (c1.c == 0x0) c1.fait() ;
    if (c2.c == 0x0) c2.fait() ;
    
    // Termine
    return (*(c1.c)) + (*(c2.c)) ;
}

Mtbl operator+(const Coord& co, const Mtbl& mt) {

    if (co.c == 0x0) co.fait() ;
    
    return (*(co.c)) + mt ;
}

Mtbl operator+(const Mtbl& mt, const Coord& co) {

    if (co.c == 0x0) co.fait() ;
    
    return mt + (*(co.c)) ;
}

			//**************//
			// SOUSTRACTION //
			//**************//

Mtbl operator-(const Coord& c1, const Coord& c2) {

    // Sont-elles a jour ?
    if (c1.c == 0x0) c1.fait() ;
    if (c2.c == 0x0) c2.fait() ;
    
    // Termine
    return (*(c1.c)) - (*(c2.c)) ;
}

Mtbl operator-(const Coord& co, const Mtbl& mt) {

    if (co.c == 0x0) co.fait() ;
    
    return (*(co.c)) - mt ;
}

Mtbl operator-(const Mtbl& mt, const Coord& co) {

    if (co.c == 0x0) co.fait() ;
    
    return mt - (*(co.c)) ;
}

			//****************//
			// MULTIPLICATION //
			//****************//

Mtbl operator*(const Coord& c1, const Coord& c2) {

    // Sont-elles a jour ?
    if (c1.c == 0x0) c1.fait() ;
    if (c2.c == 0x0) c2.fait() ;
    
    // Termine
    return (*(c1.c)) * (*(c2.c)) ;
}

Mtbl operator*(const Mtbl& m1, const Coord& c2) {

    // A jour ?
    if (c2.c == 0x0) c2.fait() ;
    
    // Termine
    return (m1) * (*(c2.c)) ;
}

Mtbl operator*(const Coord& c2, const Mtbl& m1) {

    // A jour ?
    if (c2.c == 0x0) c2.fait() ;
    
    // Termine
    return (m1) * (*(c2.c)) ;
}


}
