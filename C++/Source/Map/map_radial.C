/*
 *  Methods of class Map_radial
 *
 *   (see file map.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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


char map_radial_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/15 10:41:10  e_gourgoulhon
 * Added new Coord's drdt and stdrdp.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  1999/11/22  10:37:24  eric
 * Fonction del_coord() rebaptisee reset_coord().
 *
 * Revision 2.5  1999/10/15  14:27:20  eric
 * Suppression de l'appel a del_coord() dans le destructeur de Map_radial.
 *
 * Revision 2.4  1999/10/15  09:23:14  eric
 * *** empty log message ***
 *
 * Revision 2.3  1999/10/14  14:27:49  eric
 * Depoussierage.
 *
 * Revision 2.2  1999/03/01  16:57:35  eric
 * Suppression de l'operateur <<
 *
 * Revision 2.1  1999/03/01  14:59:25  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers Lorene
#include "map.h"


			//---------------//
			// Constructeurs //
			//---------------//

// Constructor from a grid
// -----------------------
Map_radial::Map_radial(const Mg3d& mg) : Map(mg)
{
        // The Coord's are constructed by the default Coord constructor
}

// Copy constructor 
// ----------------
Map_radial::Map_radial(const Map_radial & mp) : Map(mp)
{
        // The Coord's are constructed by the default Coord constructor
}
	
// Constructor from file
// ---------------------
Map_radial::Map_radial(const Mg3d& mgi, FILE* fd) : Map(mgi, fd) {}

			//--------------//
			// Destructeurs //
			//--------------//

// Destructeur
Map_radial::~Map_radial() {}

			//------------//
			// Sauvegarde //
			//------------//

void Map_radial::sauve(FILE* fd) const {

    Map::sauve(fd) ; 
    
}
    
			//-----------------//
			// Gestion memoire //
			//-----------------//

void Map_radial::reset_coord() {

    // Les Coord communs a toutes les classes derivees de Map : 
    
    Map::reset_coord() ;
    
    // Les Coord specifiques a Map_radial : 
    
    xsr.del_t() ; 
    dxdr.del_t() ; 
    drdt.del_t() ; 
    stdrdp.del_t() ; 
    srdrdt.del_t() ; 
    srstdrdp.del_t() ; 
    sr2drdt.del_t() ; 
    sr2stdrdp.del_t() ; 
    d2rdx2.del_t() ; 
    lapr_tp.del_t() ; 
    d2rdtdx.del_t() ; 
    sstd2rdpdx.del_t() ; 
    sr2d2rdt2.del_t() ; 
     
}
    
