/*
 *  Methods not yet implemented in class Map_log
 * 
 *   (see file map.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004 Philippe Grandclement
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


char map_log_pas_fait_C[] = "$Header $" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/11/23 12:54:45  f_limousin
 * Function poisson_frontiere(...) has two new default arguments,
 * to deal with the case of a Dirichlet + Neumann boundary condition.
 *
 * Revision 1.1  2004/06/22 08:49:58  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * 
 * $Header$
 *
 */

// headers C
#include <math.h>

// headers Lorene
#include "itbl.h"
#include "tbl.h"
#include "coord.h"
#include "grilles.h"
#include "map.h"

void pas_fait() {
  cout << "Function not implemented for Map_log..." << endl ;
  abort() ;
}

 void Map_log::homothetie (double) {
  pas_fait() ;
}
	
 void Map_log::resize (int, double) {
  pas_fait() ;
}

 void Map_log::adapt (const Cmp&, const Param&) {
  pas_fait(); 
}
	
 void Map_log::dsdr (const Cmp&, Cmp&) const {
  pas_fait() ;
}

 void Map_log::srdsdt (const Cmp&, Cmp&) const {
  pas_fait() ;
}

 void Map_log::srstdsdp (const Cmp&, Cmp&) const {
  pas_fait() ;
}

 void Map_log::srdsdt (const Scalar&, Scalar&) const {
  pas_fait() ;
}

 void Map_log::srstdsdp (const Scalar&, Scalar&) const {
  pas_fait() ; 
}

 void Map_log::dsdt (const Scalar&, Scalar&) const {
  pas_fait() ;
}

 void Map_log::stdsdp (const Scalar&, Scalar&) const {
  pas_fait() ;
}

 void Map_log::laplacien (const Cmp&, int, Cmp&) const {
  pas_fait() ;
}

 void Map_log::lapang (const Scalar&, Scalar&) const {
  pas_fait() ;
}

 Tbl* Map_log::integrale (const Cmp&) const {
  pas_fait() ;
  return 0x0 ;
}

 void Map_log::poisson (const Cmp&, Param&, Cmp&) const {
  pas_fait() ;
}

 void Map_log::poisson_regular (const Cmp&, int, int, double, Param&, Cmp&, Cmp&, Cmp&, 
				      Tenseur&, Cmp&, Cmp&) const {
  pas_fait() ;
}

 void Map_log::poisson_angu (const Scalar&, Param&, Scalar&) const {
  pas_fait() ;
}

 Param* Map_log::donne_para_poisson_vect (Param&, int) const {
  pas_fait() ;
  return 0x0 ;
}

 void Map_log::poisson_frontiere (const Cmp&, const Valeur&, int, int, Cmp&, double, double) const {
  pas_fait() ;
}

 void Map_log::poisson_frontiere_double (const Cmp&, const Valeur&, const Valeur&, int, Cmp&) const {
  pas_fait() ;
}

 void Map_log::poisson_interne (const Cmp&, const Valeur&, Param&, Cmp&) const {
  pas_fait() ;
}

 void Map_log::poisson2d (const Cmp&, const Cmp&, Param&, Cmp&) const {
  pas_fait() ;
}

 void Map_log::dalembert (Param&, Scalar&, const Scalar&, const Scalar&, const Scalar&) const {
  pas_fait() ;
}

