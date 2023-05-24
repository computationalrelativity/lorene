/*
 *  Methods not yet implemented in class Map_eps
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
void pas_fait_eps() {
  cout << "Function not implemented for Map_eps..." << endl ;
  abort() ;
}

 void Map_eps::homothetie (double) {
  pas_fait_eps() ;
}
	
 void Map_eps::resize (int, double) {
  pas_fait_eps() ;
}

 void Map_eps::adapt (const Cmp&, const Param&, int) {
  pas_fait_eps(); 
}
	
 void Map_eps::dsdr (const Cmp&, Cmp&) const {
  pas_fait_eps() ;
}
	
 void Map_eps::dsdxi (const Cmp&, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::dsdxi (const Scalar&, Scalar&) const {
  pas_fait_eps() ;
}

 void Map_eps::srdsdt (const Cmp&, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::srstdsdp (const Cmp&, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::laplacien (const Cmp&, int, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::laplacien (const Scalar&, int, Scalar&) const {
  pas_fait_eps() ;
}

 void Map_eps::lapang (const Scalar&, Scalar&) const {
  pas_fait_eps() ;
}

 Tbl* Map_eps::integrale (const Cmp&) const {
  pas_fait_eps() ;
  return 0x0 ;
}

 void Map_eps::poisson (const Cmp&, Param&, Cmp&) const {
  pas_fait_eps() ;
}

void Map_eps::poisson_tau (const Cmp&, Param&, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::poisson_regular (const Cmp&, int, int, double, Param&, Cmp&, Cmp&, Cmp&, 
				      Tenseur&, Cmp&, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::poisson_angu (const Scalar&, Param&, Scalar&, double) const {
  pas_fait_eps() ;
}

 void Map_eps::poisson_angu (const Cmp&, Param&, Cmp&, double) const {
  pas_fait_eps() ;
}

 Param* Map_eps::donne_para_poisson_vect (Param&, int) const {
  pas_fait_eps() ;
  return 0x0 ;
}

 void Map_eps::poisson_frontiere (const Cmp&, const Valeur&, int, int, Cmp&, double, double) const {
  pas_fait_eps() ;
}

 void Map_eps::poisson_frontiere_double (const Cmp&, const Valeur&, const Valeur&, int, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::poisson_interne (const Cmp&, const Valeur&, Param&, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::poisson2d (const Cmp&, const Cmp&, Param&, Cmp&) const {
  pas_fait_eps() ;
}

 void Map_eps::dalembert (Param&, Scalar&, const Scalar&, const Scalar&, const Scalar&) const {
  pas_fait_eps() ;
}

const Map_af& Map_eps::mp_angu(int) const {
    pas_fait_eps() ;
}

void Map_eps::primr(const Scalar&, Scalar&, bool) const {
  pas_fait_eps() ;
}

void Map_eps::poisson_falloff(const Cmp&, Param&, Cmp&, int) const {
  pas_fait_eps() ;
}

void Map_eps::poisson_ylm(const Cmp&, Param&, Cmp&, int, double*) const {
  pas_fait_eps() ;
}

void Map_eps::dsdradial (const Scalar&, Scalar&) const{
  pas_fait_eps() ;
}

}