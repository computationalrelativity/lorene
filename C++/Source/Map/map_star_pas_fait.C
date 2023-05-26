/*
 *  Methods not yet implemented in class Map_star
 * 
 *   (see file map.h for documentation)
 *
 */

/*
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

// headers C
#include <cmath>

// headers Lorene
#include "itbl.h"
#include "tbl.h"
#include "coord.h"
#include "grilles.h"
#include "map.h"

namespace Lorene {
void pas_fait_star() {
  cout << "Function not implemented for Map_star..." << endl ;
  abort() ;
}

 void Map_star::homothetie (double) {
  pas_fait_star() ;
}
	
 void Map_star::resize (int, double) {
  pas_fait_star() ;
}

 void Map_star::adapt (const Cmp&, const Param&, int) {
  pas_fait_star(); 
}
	
 void Map_star::dsdr (const Cmp&, Cmp&) const {
  pas_fait_star() ;
}
	
 void Map_star::dsdxi (const Cmp&, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::dsdxi (const Scalar&, Scalar&) const {
  pas_fait_star() ;
}

 void Map_star::srdsdt (const Cmp&, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::srstdsdp (const Cmp&, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::laplacien (const Cmp&, int, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::laplacien (const Scalar&, int, Scalar&) const {
  pas_fait_star() ;
}

 void Map_star::lapang (const Scalar&, Scalar&) const {
  pas_fait_star() ;
}

 Tbl* Map_star::integrale (const Cmp&) const {
  pas_fait_star() ;
  return 0x0 ;
}

 void Map_star::poisson (const Cmp&, Param&, Cmp&) const {
  pas_fait_star() ;
}

void Map_star::poisson_tau (const Cmp&, Param&, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::poisson_regular (const Cmp&, int, int, double, Param&, Cmp&, Cmp&, Cmp&, 
				      Tenseur&, Cmp&, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::poisson_angu (const Scalar&, Param&, Scalar&, double) const {
  pas_fait_star() ;
}

 void Map_star::poisson_angu (const Cmp&, Param&, Cmp&, double) const {
  pas_fait_star() ;
}

 Param* Map_star::donne_para_poisson_vect (Param&, int) const {
  pas_fait_star() ;
  return 0x0 ;
}

 void Map_star::poisson_frontiere (const Cmp&, const Valeur&, int, int, Cmp&, double, double) const {
  pas_fait_star() ;
}

 void Map_star::poisson_frontiere_double (const Cmp&, const Valeur&, const Valeur&, int, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::poisson_interne (const Cmp&, const Valeur&, Param&, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::poisson2d (const Cmp&, const Cmp&, Param&, Cmp&) const {
  pas_fait_star() ;
}

 void Map_star::dalembert (Param&, Scalar&, const Scalar&, const Scalar&, const Scalar&) const {
  pas_fait_star() ;
}

const Map_af& Map_star::mp_angu(int) const {
    pas_fait_star() ;
}

void Map_star::primr(const Scalar&, Scalar&, bool) const {
  pas_fait_star() ;
}

void Map_star::poisson_falloff(const Cmp&, Param&, Cmp&, int) const {
  pas_fait_star() ;
}

void Map_star::poisson_ylm(const Cmp&, Param&, Cmp&, int, double*) const {
  pas_fait_star() ;
}

void Map_star::dsdradial (const Scalar&, Scalar&) const{
  pas_fait_star() ;
}

}