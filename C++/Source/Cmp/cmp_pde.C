/*
 * Methods of the class Cmp for various partial differential equations
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
 *   Copyright (c) 2000-2001 Jerome Novak
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


char cmp_pde_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.8  2001/10/16  09:59:00  novak
 * deletion of source(t=jm1) from the argument list of d'Alembert solvers
 *
 * Revision 1.7  2001/07/19 14:04:07  novak
 * new argument list for Cmp::avance_dalembert
 *
 * Revision 1.6  2000/12/04 15:07:55  novak
 * *** empty log message ***
 *
 * Revision 1.5  2000/10/19 14:16:43  novak
 * Ajout de Cmp::avance_dalembert (etat experimental)
 *
 * Revision 1.4  1999/12/21 14:47:30  eric
 * *** empty log message ***
 *
 * Revision 1.3  1999/12/21  13:04:00  eric
 * Changement de prototype de la routine poisson avec Param& : la solution est
 * desormais passee en argument (et non plus en valeur de retour)
 * pour permettre l'initialisation de methodes de resolution iteratives.
 *
 * Revision 1.2  1999/12/21  10:07:25  eric
 * Il y a desormais deux versions de poisson: une sans Param et une
 * avec Param.
 *
 * Revision 1.1  1999/12/02  14:30:13  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Header Lorene:
#include "map.h"
#include "cmp.h"
#include "param.h"

		    //-----------------------------------//
		    //      Scalar Poisson equation	 //
		    //-----------------------------------//

// Version without parameters
// --------------------------

Cmp Cmp::poisson() const {
    
    Param bidon ;
    Cmp resu(*mp) ; 
    
    mp->poisson(*this, bidon, resu) ; 

    return resu ;          
}

// Version with parameters
// -----------------------

void Cmp::poisson(Param& par, Cmp& uu) const {
    
    mp->poisson(*this, par, uu) ;     
    
}



		    //-----------------------------------//
		    //      Scalar d'Alembert equation	 //
		    //-----------------------------------//

Cmp Cmp::avance_dalembert(Param& par, Cmp& fjm1, Cmp& source) const {
  
  Cmp fjp1(*mp) ;
  
  mp->dalembert(par, fjp1, *this, fjm1, source) ;

  return fjp1 ;
}
