/*
 * Methods for magnetized axisymmetric rotating neutron stars.
 *
 * See the file et_rot_mag.h for documentation
 *
 */

/*
 *   Copyright (c) 2002 Emmanuel Marcq
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

char et_rot_mag_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.9  2002/05/20 15:44:55  e_marcq
 *
 * Dimension errors corrected, parmag.d input file created and read
 *
 * Revision 1.8  2002/05/20 08:27:59  j_novak
 * *** empty log message ***
 *
 * Revision 1.7  2002/05/17 15:08:01  e_marcq
 *
 * Rotation progressive plug-in, units corrected, Q and a_j new member data
 *
 * Revision 1.6  2002/05/16 13:27:11  j_novak
 * *** empty log message ***
 *
 * Revision 1.5  2002/05/16 11:54:11  j_novak
 * Fixed output pbs
 *
 * Revision 1.4  2002/05/15 09:53:59  j_novak
 * First operational version
 *
 * Revision 1.3  2002/05/14 13:38:36  e_marcq
 *
 *
 * Unit update, new outputs
 *
 * Revision 1.1  2002/05/10 09:26:52  j_novak
 * Added new class Et_rot_mag for magnetized rotating neutron stars (under development)
 *
 *
 * $Header$
 *
 */

// Headers C
#include "math.h"

// Headers Lorene
#include "et_rot_mag.h"
#include "utilitaires.h"
#include "unites_mag.h"

			    //--------------//
			    // Constructors //
			    //--------------//
// Standard constructor
// --------------------


Et_rot_mag::Et_rot_mag(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i)
  : Etoile_rot(mp_i, nzet_i, relat, eos_i),
    A_t(mp_i),
    A_phi(mp_i),
    j_t(mp_i),
    j_phi(mp_i),
    E_em(mp_i),
    Jp_em(mp_i),
    Srr_em(mp_i),
    Spp_em(mp_i)

{

  A_t = 0;
  A_phi = 0; 
  j_t = 0 ;
  j_phi = 0 ;

  Q = 0 ;
  a_j = 0 ;

set_der_0x0() ;  
}


// Copy constructor
// ----------------

Et_rot_mag::Et_rot_mag(const Et_rot_mag& et)
  : Etoile_rot(et),
  A_t(et.A_t),
  A_phi(et.A_phi),
  j_phi(et.j_phi),
  j_t(et.j_t),
  E_em(et.E_em),
  Jp_em(et.Jp_em),
  Srr_em(et.Srr_em),
  Spp_em(et.Spp_em)

{
  set_der_0x0() ;
}


			    //------------//
			    // Destructor //
			    //------------//

Et_rot_mag::~Et_rot_mag(){
  del_deriv() ;
}


		//----------------------------------//
		// Management of derived quantities //
		//----------------------------------//

void Et_rot_mag::del_deriv() const {

  Etoile_rot::del_deriv() ;

  // Quelles quantites derivees ?

  set_der_0x0() ;

}


void Et_rot_mag::set_der_0x0() const {
  Etoile_rot::set_der_0x0() ;

  // Quelles quantites derivees ?

}


void Et_rot_mag::del_hydro_euler() {
  Etoile_rot::del_hydro_euler() ;
  // Nouvelles qtes hydro a remettre en nondef
  del_deriv() ;
}


// Assignment to another Et_rot_mag
// --------------------------------

void Et_rot_mag::operator=(const Et_rot_mag& et) {

  // Assignement of quantities common to all the derived classes of Etoile
  Etoile_rot::operator=(et) ;
  A_t    = et.A_t    ;
  A_phi  = et.A_phi  ;
  j_t    = et.j_t    ;
  j_phi  = et.j_phi  ;
  E_em   = et.E_em   ;
  Jp_em  = et.Jp_em  ;
  Srr_em = et.Srr_em ;
  Spp_em = et.Spp_em ;
  Q      = et.Q      ;
  a_j    = et.a_j    ;

  del_deriv() ;

}

			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Printing
// --------


ostream& Et_rot_mag::operator>>(ostream& ost) const {

  //#include "unites.h"
  // Compilation warnings ?

  Etoile_rot::operator>>(ost) ;
  ost << endl ;
  ost << "Electromagnetic quantities" << endl ;
  ost << "----------------------" << endl ;
  ost << "In construction. ALL OUTPUTS IN SI UNITS HERE" << endl ;
  ost << endl ;
  ost << "Prescribed charge : " << Q*j_unit*pow(r_unit,3)/v_unit << endl;
  ost << "Prescribed current amplitude : " << a_j*j_unit << endl ;
  ost << "Magnetic Momentum : " << MagMom() << endl ;
  ost << "Radial magnetic field polar value : " << 
    Magn()(0).va.val_point(l_surf()(0,0),xi_surf()(0,0),0.,0.) << endl;

  ost << "Magnetic pressure (Srr_em) polar value : " << 
    Srr_em().va.val_point(l_surf()(0,0),xi_surf()(0,0),0.,0.) << endl ;

  ost << "Radial electric field polar value : " << 
    Elec()(0).va.val_point(l_surf()(0,0),xi_surf()(0,0),0.,0.) << endl;

  int theta_eq = mp.get_mg()->get_nt(nzet-1)-1 ;

  ost << "Central pressure : "<< 1/(2*mu_si)*(pow(Magn()(0)(0,0,0,0),2)+pow(Magn()(1)(0,0,0,0),2)+pow(Magn()(2)(0,0,0,0),2))*rho_unit*pow(v_unit,2) << endl ;
  ost << "Tangent magnetic field equatorial value : " << 
  Magn()(1).va.val_point(l_surf()(0,theta_eq),xi_surf()(0,theta_eq),M_PI_2,0.) 
      << endl;
  ost << "Magnetic pressure equatorial value : " << 
  Srr_em().va.val_point(l_surf()(0,theta_eq),xi_surf()(0,theta_eq),M_PI_2,0.) 
      << endl ;
  ost << "Radial electric field equatorial value : " << 
  Elec()(0).va.val_point(l_surf()(0,theta_eq),xi_surf()(0,theta_eq),M_PI_2,0.) 
      << endl;
  ost << "Computed charge : " << Q_comput() << endl ;
  ost << "Gyromagnetic ratio : " << GyroMag() << endl ;


  return ost ;
}
















