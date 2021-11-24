/*
 *
 * Multiplication by (x+1)
 *
 * for:
 *   - Valeur
 *   - Mtbl_cf
 */

/*
 *   Copyright (c) 2021 Gaël Servignat
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


// Headers C
#include <cassert>

// Headers Lorene
#include "mtbl_cf.h"
#include "valeur.h"

// Local prototypes
namespace Lorene {
void _mult_xp1_identite(Tbl*, int&) ;
void _mult_xp1_cheb(Tbl*, int&) ;


void Valeur::mult_xp1_shell(int lz) {

    // Peut-etre ne rien faire ?
    if (etat==ETATZERO) {
	return ; 
    }

    assert(etat==ETATQCQ) ; 

    // Calcul des coef.
    coef() ;
    
    // Division par (x-1) dans la ZEC 
    c_cf->mult_xp1_shell(lz) ;
    set_etat_cf_qcq() ;

    base = c_cf->base ; // On remonte la base de sortie au niveau Valeur
    
}


/*
 * Fonction membre de la classe Mtbl_cf pour la multiplication par (x+1) 
 * dans une coquille applique a this
 *
 */

void Mtbl_cf::mult_xp1_shell(int lz)	   
{

  if (mg->get_type_r(lz) != FIN) {
    cerr << "Mtbl_cf::mult_xp1_shell() : not called on a shell!" << endl ;
    abort() ;
  }
  
  // Peut-etre ne rien faire ?
  if (etat==ETATZERO) {
    return ; 
  }

  assert(etat==ETATQCQ) ; 

  _mult_xp1_cheb(t[lz], base.b[lz]) ;
    
}

}