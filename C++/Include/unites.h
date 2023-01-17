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
 * Revision 1.10  2023/01/17 15:01:25  j_novak
 * Update of the value of g_si from CODATA 2018. New definition of rhonuc_si
 * as 0.1 unified atomic mass unit per Fermi^3. Note that results with Lorene
 * may change because of these updates.
 *
 * Revision 1.9  2016/12/05 15:28:52  j_novak
 * Suppressed the use of 'pow' to avoid compilation warnings.
 *
 * Revision 1.8  2015/03/17 14:20:00  j_novak
 * New class Hot_eos to deal with temperature-dependent EOSs.
 *
 * Revision 1.7  2014/10/13 08:52:37  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/06/30 14:34:57  j_novak
 * Update of the values of some constants (G, M_sol and MeV).
 *
 * Revision 1.5  2014/05/13 10:06:12  j_novak
 * Change of magnetic units, to make the Lorene unit system coherent. Magnetic field is now expressed in Lorene units. Improvement on the comments on units.
 *
 * Revision 1.4  2004/12/01 12:28:32  p_grandclement
 * Include math.h in unite.h
 *
 * Revision 1.3  2004/03/25 10:28:56  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
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

                  //--------------------------//
                  //  Standard LORENE units   //
                  //--------------------------//     

namespace Lorene {
  /** \namespace Unites
   *  \brief Standard units of space, time and mass.
   *
   * These are the units used in LORENE for space, time and mass. They
   * are mainly designed to study compact objects like neutron stars 
   * and black holes. Last update to comply with CODATA 2018.\ingroup (unites)
   */
  namespace Unites {

    //*** G constant and density unit before January 2023 ***
    // To reproduce pre-2023 results, please change the difinitions of
    // g_si and rhonuc_si by uncommenting these two lines and commenting
    // the new corresponding ones.
    //const double g_si = 6.6738E-11 ;	 ///< Newton gravitational constant [SI]
    //const double rhonuc_si = 1.66E+17 ; ///< Nuclear density [kg/m3]
    
    const double g_si = 6.6743E-11 ;	 ///< Newton gravitational constant [SI]
    const double c_si = 2.99792458E+8 ;	 ///< Velocity of light [m/s]
    const double kB_si = 1.380649E-23 ; ///< Boltzmann constant [J/K]
    const double m_u_si = 1.6605390666E-27 ;  ///< atomic mass unit [kg]
    const double mev_si = 1.602176634E-13 ;   ///< One MeV [J]
    /// Nuclear density [kg/m3] := 0.1 atomic  mass / Fermi^3 (arbitrary)
    const double rhonuc_si = m_u_si * 1.e44 ;
    const double km_si = 1.E+3 ;	 ///< One kilometer [m]
    const double msol_si = 1.9885E+30 ;	 ///< Solar mass [kg]
    
    const double r_unit = 1.e4 ;  ///< Lorene's unit of length = 10 km
    const double v_unit = c_si ; ///< Lorene's unit of velocity = c 
    const double rho_unit = rhonuc_si ;	///< Lorene's unit of mass density
    const double t_unit = r_unit/v_unit ; ///< Lorene's unit of time
    const double m_unit = rho_unit * r_unit*r_unit*r_unit ;  ///< Lorene's unit of mass
    const double g_unit = 1./(rho_unit*t_unit*t_unit) ; ///< Lorene's unit for G
    const double f_unit = 1./t_unit ;	///< Lorene's unit of frequency
    
    const double ggrav = g_si / g_unit ;  ///< G in Lorene's units
    const double qpig = 4 * M_PI * ggrav ; ///< 4 Pi G in Lorene's units
    const double msol = msol_si/m_unit ; ///< Solar mass in Lorene's units
    const double km = km_si/r_unit ;	///< One kilometer in Lorene's units
    /// 1 MeV/fm3 in Lorene's units
    const double mevpfm3 = mev_si/( rho_unit * v_unit *v_unit) *1.e45 ;  
    /// Atomic mass conversion from Lorene's units to MeV
    const double m_u_mev = m_u_si *c_si*c_si / mev_si ;
  }
  
  
  //----------------------------------//
  //  Electro-magnetic LORENE units   //
  //----------------------------------//     
  
  
  
  /** \namespace Unites_mag
   *  \brief Standard electro-magnetic units.
   *
   * \ingroup (unites)
   */
  namespace Unites_mag {
    using namespace Unites ;
    const double mu_si = 1.2566370614359173e-6 ;///<Magnetic vacuum permeability
    
    const double j_unit = 1e11 ; ///<Lorene's current density unit [\f$A/m^2\f$]
    
    /// Lorene's units for magnetic field 
    const double mag_unit = rho_unit*v_unit*v_unit/ (r_unit * j_unit) ;
    /// Lorene's unit for electric field 
    const double elec_unit = mag_unit * v_unit ;
    /// Lorene's unit for \f$\mu_0\f$
    const double mu0_unit = rho_unit*v_unit*v_unit / (j_unit*j_unit*r_unit*r_unit);
    /// \f$\mu_0\f$ in Lorene's units
    const double mu0 = mu_si / mu0_unit ;
  }
  
}
