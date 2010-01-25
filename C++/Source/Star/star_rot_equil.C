/*
 * Method Star_rot::equilibrium
 *
 * (see file star_h.h for documentation)
 *
 */

/*
 *   Copyright (c) 2010 Eric Gourgoulhon
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


char star_rot_equil_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2010/01/25 18:15:52  e_gourgoulhon
 * First version.
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <math.h>

// Headers Lorene
#include "star_rot.h"
#include "param.h"

#include "graphique.h"
#include "utilitaires.h"
#include "unites.h"

void Star_rot::equilibrium(double ent_c, double omega0, double fact_omega, 
			     int nzadapt, const Tbl& ent_limit, const Itbl& icontrol,
			     const Tbl& control, double mbar_wanted, 
			     double aexp_mass, Tbl& diff, Param*) {
			     
	cout << "Star_rot::equilibrium : not implemented yet ! " << endl ; 
	abort() ; 

}
