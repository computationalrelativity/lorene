/*
 *  Lorene's units
 *
 */

/*
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


/*
 * $Id$
 * $Log$
 * Revision 1.2  2004/03/22 13:12:43  j_novak
 * Modification of comments to use doxygen instead of doc++
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.2  2000/03/17  15:28:29  eric
 * Ajout de ggrav (G).
 *
 * Revision 1.1  1999/12/06  13:35:10  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

/** \defgroup std_unit Standard units (related to space, time and mass).
 * \ingroup (unites)
 * @{
 */

const double g_si = 6.6726E-11 ;	 ///< Newton gravitational constant [SI]
const double c_si = 2.99792458E+8 ;	 ///< Velocity of light [m/s]
const double rhonuc_si = 1.66E+17 ;	 ///< Nuclear density [kg/m3]
const double km_si = 1.E+3 ;	 ///< One kilometer [m]
const double msol_si = 1.989E+30 ;	 ///< Solar mass [kg]
const double mev_si = 1.6021892E-13 ;   ///< One MeV [J]

const double r_unit = 1.e4 ;  ///< Lorene's unit of length = 10 km
const double v_unit = c_si ; ///< Lorene's unit of velocity = c 
const double rho_unit = rhonuc_si ;	///< Lorene's unit of mass density
const double t_unit = r_unit/v_unit ; ///< Lorene's unit of time
const double m_unit = rho_unit * pow(r_unit, 3.) ;  ///< Lorene's unit of mass
const double g_unit = 1./(rho_unit*t_unit*t_unit) ; ///< Lorene's unit for G
const double f_unit = 1./t_unit ;	///< Lorene's unit of frequency
    
const double ggrav = g_si / g_unit ;  ///< G in Lorene's units
const double qpig = 4 * M_PI * ggrav ; ///< 4 Pi G in Lorene's units
const double msol = msol_si/m_unit ; ///< Solar mass in Lorene's units
const double km = km_si/r_unit ;	///< One kilometer in Lorene's units
/// 1 MeV/fm3 in Lorene's units
const double mevpfm3 = mev_si/( 1.66E-27 * v_unit *v_unit) *10 ;  

/** @} */

    
