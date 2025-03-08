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

 

/*
 * $Id$
 * $Log$
 * Revision 1.25  2016/12/05 16:17:54  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.24  2014/10/13 08:52:58  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.23  2014/05/13 15:36:54  j_novak
 * *** empty log message ***
 *
 * Revision 1.22  2014/05/13 10:06:13  j_novak
 * Change of magnetic units, to make the Lorene unit system coherent. Magnetic field is now expressed in Lorene units. Improvement on the comments on units.
 *
 * Revision 1.21  2013/11/25 13:52:11  j_novak
 * New class Et_magnetisation to include magnetization terms in the stress energy tensor.
 *
 * Revision 1.20  2013/11/14 16:12:55  j_novak
 * Corrected a mistake in the units.
 *
 * Revision 1.19  2012/08/12 17:48:35  p_cerda
 * Magnetstar: New classes for magnetstar. Allowing for non-equatorial symmetry in Etoile et al. Adding B_phi in Et_rot_mag.
 *
 * Revision 1.18  2011/10/06 14:55:36  j_novak
 * equation_of_state() is now virtual to be able to call to the magnetized
 * Eos_mag.
 *
 * Revision 1.17  2005/06/02 11:35:30  j_novak
 * Added members for sving to a file and reading from it.
 *
 * Revision 1.16  2004/03/25 10:29:06  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.15  2002/09/30 14:21:21  j_novak
 * *** empty log message ***
 *
 * Revision 1.14  2002/09/30 08:50:37  j_novak
 * Output for central magnetic field value
 *
 * Revision 1.13  2002/06/05 15:15:59  j_novak
 * The case of non-adapted mapping is treated.
 * parmag.d and parrot.d have been merged.
 *
 * Revision 1.12  2002/06/03 13:00:45  e_marcq
 *
 * conduc parameter read in parmag.d
 *
 * Revision 1.10  2002/05/22 12:20:17  j_novak
 * *** empty log message ***
 *
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
#include "cmath"

// Headers Lorene
#include "et_rot_mag.h"
#include "utilitaires.h"
#include "unites.h"
#include "param.h"
#include "eos.h"

			    //--------------//
			    // Constructors //
			    //--------------//
// Standard constructor
// --------------------


namespace Lorene {
Et_rot_mag::Et_rot_mag(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
		       const int cond)
  : Etoile_rot(mp_i, nzet_i, relat, eos_i),
    A_t(mp_i),
    A_phi(mp_i),
    B_phi(mp_i),
    j_t(mp_i),
    j_phi(mp_i),
    E_em(mp_i),
    Jp_em(mp_i),
    Srr_em(mp_i),
    Spp_em(mp_i)

{

  A_t = 0;
  A_phi = 0; 
  B_phi = 0;
  j_t = 0 ;
  j_phi = 0 ;

  Q = 0 ;
  a_j = 0 ;
  conduc = cond ;

set_der_0x0() ;  
}



Et_rot_mag::Et_rot_mag(Map& mp_i, const Eos& eos_i, FILE* fich, int withbphi)
    : Etoile_rot(mp_i, eos_i, fich), 
      A_t(mp_i),
      A_phi(mp_i),
      B_phi(mp_i),
      j_t(mp_i),
      j_phi(mp_i),
      E_em(mp_i),
      Jp_em(mp_i),
      Srr_em(mp_i),
      Spp_em(mp_i)
{

    // Etoile parameters
    // -----------------

    fread_be(&conduc, sizeof(int), 1, fich) ;		
    fread_be(&Q, sizeof(double), 1, fich) ;		
    fread_be(&a_j, sizeof(double), 1, fich) ;		
   
    // Read of the saved fields:
    // ------------------------
    
    Cmp A_t_file(mp, *mp.get_mg(), fich) ;
    A_t = A_t_file ;

    Cmp A_phi_file(mp, *mp.get_mg(), fich) ;
    A_phi = A_phi_file ;
  
    Cmp j_t_file(mp, *mp.get_mg(), fich) ;
    j_t = j_t_file ;

    Cmp j_phi_file(mp, *mp.get_mg(), fich) ;
    j_phi = j_phi_file ;

    Tenseur E_em_file(mp, fich) ;
    E_em = E_em_file ;

    Tenseur Jp_em_file(mp, fich) ;
    Jp_em = Jp_em_file ;

    Tenseur Srr_em_file(mp, fich) ;
    Srr_em = Srr_em_file ;

    Tenseur Spp_em_file(mp, fich) ;
    Spp_em = Spp_em_file ;

    if ( withbphi == 0 ) {
      B_phi = 0;
    }else{
      Cmp B_phi_file(mp, *mp.get_mg(), fich) ;
      B_phi = B_phi_file ;
    }

    // Pointers of derived quantities initialized to zero 
    // --------------------------------------------------
    set_der_0x0() ;
    
}


// Copy constructor
// ----------------

Et_rot_mag::Et_rot_mag(const Et_rot_mag& et)
  : Etoile_rot(et),
  A_t(et.A_t),
  A_phi(et.A_phi),
  B_phi(et.B_phi),
  j_t(et.j_t),
  j_phi(et.j_phi),
  E_em(et.E_em),
  Jp_em(et.Jp_em),
  Srr_em(et.Srr_em),
  Spp_em(et.Spp_em)

{
  Q = et.Q ;
  a_j = et.a_j ;
  conduc = et.conduc ;
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

  set_der_0x0() ;

}


void Et_rot_mag::set_der_0x0() const {
  Etoile_rot::set_der_0x0() ;

}


void Et_rot_mag::del_hydro_euler() {
  Etoile_rot::del_hydro_euler() ;

  del_deriv() ;
}


// Assignment to another Et_rot_mag
// --------------------------------

void Et_rot_mag::operator=(const Et_rot_mag& et) {

  // Assignement of quantities common to all the derived classes of Etoile
  Etoile_rot::operator=(et) ;
  A_t    = et.A_t    ;
  A_phi  = et.A_phi  ;
  B_phi  = et.B_phi  ;
  j_t    = et.j_t    ;
  j_phi  = et.j_phi  ;
  E_em   = et.E_em   ;
  Jp_em  = et.Jp_em  ;
  Srr_em = et.Srr_em ;
  Spp_em = et.Spp_em ;
  Q      = et.Q      ;
  a_j    = et.a_j    ;
  conduc = et.conduc ;

  del_deriv() ;

}


			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_rot_mag::sauve(FILE* fich) const {
    
    Etoile_rot::sauve(fich) ; 
    
    fwrite_be(&conduc, sizeof(int), 1, fich) ;
    fwrite_be(&Q, sizeof(double), 1, fich) ;
    fwrite_be(&a_j, sizeof(double), 1, fich) ;

    A_t.sauve(fich) ;
    A_phi.sauve(fich) ;
    j_t.sauve(fich) ;
    j_phi.sauve(fich) ;
    E_em.sauve(fich) ;
    Jp_em.sauve(fich) ;
    Srr_em.sauve(fich) ;
    Spp_em.sauve(fich) ;
    B_phi.sauve(fich) ;
    
}


// Printing
// --------


ostream& Et_rot_mag::operator>>(ostream& ost) const {

  using namespace Unites_mag ;

  Etoile_rot::operator>>(ost) ;
  int theta_eq = mp.get_mg()->get_nt(nzet-1)-1 ;
  ost << endl ;
  ost << "Electromagnetic quantities" << endl ;
  ost << "----------------------" << endl ;
  ost << endl ;
  if (is_conduct()) {
    ost << "***************************" << endl ;
    ost << "**** perfect conductor ****" << endl ;
    ost << "***************************" << endl ;    
    ost << "Prescribed charge : " << Q*j_unit*pow(r_unit,3)/v_unit 
	<< " [C]" << endl;
  }
  else {
    ost << "***************************" << endl ;
    ost << "****     isolator      ****" << endl ;
    ost << "***************************" << endl ;  
  }  
  ost << "Prescribed current amplitude : " << a_j*j_unit 
      << " [A/m2]" << endl ;
  ost << "Magnetic Momentum : " << MagMom()/1.e32
      << " [10^32 Am^2]" << endl ;
  ost << "Radial magnetic field polar value : " << 
    Magn()(0).va.val_point(l_surf()(0,0),xi_surf()(0,0),0.,0.)*mag_unit/1.e9 
      << " [GT]" << endl;

  ost << "Tangent magnetic field equatorial value : " << 
  Magn()(1).va.val_point(l_surf()(0,theta_eq),xi_surf()(0,theta_eq),M_PI_2,0.)
    *mag_unit/1.e9
      << " [GT]" << endl;

  ost << "Central magnetic field values : " << 
    Magn()(0)(0,0,0,0)*mag_unit/1.e9  
      << " [GT]" << endl;

   ost << "Radial electric field polar value : " << 
    Elec()(0).va.val_point(l_surf()(0,0),xi_surf()(0,0),0.,0.)
     *elec_unit/1.e12
      << " [TV]" << endl;

  ost << "Radial electric field equatorial value : " << 
  Elec()(0).va.val_point(l_surf()(0,theta_eq),xi_surf()(0,theta_eq),M_PI_2,0.) 
    *elec_unit/1.e12
      << " [TV]" << endl;

  ost << "Magnetic/fluid pressure : "
      << 1/(2*mu_si)*(pow(Magn()(0)(0,0,0,0),2) + 
		      pow(Magn()(1)(0,0,0,0),2) + 
		      pow(Magn()(2)(0,0,0,0),2))*mag_unit*mag_unit
    / (press()(0,0,0,0)*rho_unit*pow(v_unit,2)) << endl ;
  ost << "Computed charge : " << Q_comput() << " [C]" << endl ;
  ost << "Interior charge : " << Q_int() << " [C]" << endl ;
  ost << "Q^2/M^2 :" << mu_si*c_si*c_si*Q_comput()*Q_comput()/(4*M_PI*g_si*mass_g()*mass_g()*m_unit*m_unit) 
      << endl ;
  ost << "Gyromagnetic ratio : " << GyroMag() << endl ;

  return ost ;
}
















}
