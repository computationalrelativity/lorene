/*
 *  Lorene's Electro-Magnetic units
 *
 */

/*
 *   Copyright (c) 2002 Jerome Novak
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



/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/09/19 14:12:37  e_gourgoulhon
 * Modif documentation for LaTeX compliance.
 *
 * Revision 1.1  2002/05/15 09:53:59  j_novak
 * First operational version
 *
 *
 * $Header$
 *
 */
// Headers Lorene
#include "unites.h"

/** @name Fundamental constants (in SI units)
 */
//@{     
const double mu_si = 1.2566370614359173e-6 ;///Magnetic vacuum permeability
//@}

/** @name Definition of Lorene's units
 */
//@{
const double j_unit = 1e11 ; ///Lorene's current density unit [$A/m^2$]
/// Lorene's units for magnetic field [$10^9$ T]
const double mag_unit = mu_si * r_unit * j_unit / 1e9 ;
/// Lorene's unit for electric field [$10^{12}$ V/m]
const double elec_unit = mag_unit * c_si / 1e3 ;
//@}

/** @name Fundamental constants (in Lorene's units)
 */
//@{
/// $\mu_0$ in Lorene's units
const double mu0 = mu_si * pow(j_unit ,2) * pow(r_unit,2) 
     / (rho_unit*pow(c_si,2)) ; 
	//@}
    
