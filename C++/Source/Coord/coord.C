/*
 *  Methods of class Coord
 *
 *   (see file coord.h for documentation)
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


char coord_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/06 15:13:04  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2008/02/18 13:53:39  j_novak
 * Removal of special indentation instructions.
 *
 * Revision 1.2  2002/10/16 14:36:34  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  1999/10/15  09:16:11  eric
 * Depoussierage.
 *
 * Revision 2.2  1999/03/01  15:07:25  eric
 * *** empty log message ***
 *
 * Revision 2.1  1999/02/23  14:57:03  hyc
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// Fichier includes
#include <cstdlib>
#include <cstdio>
#include "coord.h"
#include "mtbl.h"

			//---------------//
			// Constructeurs //
    			//---------------//

// Constructeur par defaut
Coord::Coord() : mp(0x0), met_fait(0x0), c(0x0) {}

// Constructeur
Coord::Coord(const Map* mpi, Mtbl* (*construit)(const Map* ) ) : mp(mpi), 
							   met_fait(construit), 
							   c(0x0)
{}

			//--------------//
			// Destructeur  //
			//--------------//

Coord::~Coord() {
    delete c ;
}

			//------------//
			// Impression //
			//------------//

// Operateurs <<
ostream& operator<<(ostream& o, const Coord & ci) {

    if (ci.c == 0x0) {
	o << "La coordonnee n'est pas a jour, je la fais." << endl ;
	ci.fait() ;
    }
    o << "Coordonnee: " << endl ;
    o << *(ci.c) << endl ;
    return o ;
}
    
			//----------//
			// Methodes //
			//----------//

void Coord::fait() const {
    delete c ;
    c = met_fait(mp) ;
}
    
			//-----------------//
			// Gestion memoire //
			//-----------------//

void Coord::del_t() const {
    delete c ;
    c = 0x0 ;
}

			//--------------------//
			// Fonctions diverses //
			//--------------------//

void Coord::set(const Map* mpi, Mtbl* (*construit)(const Map*) ) {
    mp = mpi ;
    met_fait = construit ;
}
