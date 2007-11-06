/*
 *  Methods of class Star_rot_Dirac
 *
 *    (see file star_rot_dirac.h for documentation).
 *
 */

/*
 *   Copyright (c) 2005  Lap-Ming Lin & Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char star_rot_dirac_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2007/11/06 16:23:59  j_novak
 * Added the flag spectral_filter giving the order of possible spectral filtering
 * of the hydro sources of metric equations (some members *_euler). The filtering
 * is done in strot_dirac_hydro, if this flag is non-zero.
 *
 * Revision 1.5  2007/11/06 10:15:19  j_novak
 * Change the order of updates in the constructor from a file, to avoid
 * inconsistencies.
 *
 * Revision 1.4  2007/10/30 16:55:23  j_novak
 * Completed the input/ouput to a file
 *
 * Revision 1.3  2005/02/09 13:37:37  lm_lin
 *
 * Add pointers p_tsw, p_aplat, and p_r_circ; add more screen output
 * information.
 *
 * Revision 1.2  2005/02/02 09:22:29  lm_lin
 *
 * Add the GRV3 error to screen output
 *
 * Revision 1.1  2005/01/31 08:51:48  j_novak
 * New files for rotating stars in Dirac gauge (still under developement).
 *
 *
 * $Header$
 *
 */


// C headers
#include <math.h>
#include <assert.h>

// Lorene headers
#include "star_rot_dirac.h"
#include "unites.h" 
#include "utilitaires.h"


                   //--------------//
                   // Constructors //
                   //--------------//

// Standard constructor
//-------------------------
Star_rot_Dirac::Star_rot_Dirac(Map& mpi, int nzet_i, const Eos& eos_i, int filter)
                   : Star(mpi, nzet_i, eos_i),
		     spectral_filter(filter),
		     psi4(mpi),
		     psi2(mpi),
		     qqq(mpi),
		     ln_psi(mpi),
		     j_euler(mpi, CON, mpi.get_bvect_spher()),
		     v2(mpi),
		     flat(mpi.flat_met_spher()),
		     tgamma(flat),
		     aa(mpi, CON, mpi.get_bvect_spher()),
		     taa(mpi, COV, mpi.get_bvect_spher()),
		     aa_quad(mpi),
		     hh(mpi, mpi.get_bvect_spher(), flat) 
{
    assert (spectral_filter>=0) ;
    assert (spectral_filter<1000) ;

  // Initialization to a static state
  omega = 0 ;
  v2 = 0 ;

  // All the matter quantities are initialized to zero
  j_euler.set_etat_zero() ;

  // Initialization to a flat case
  psi4 = 1 ;
  psi2 = 1 ;
  qqq = 1 ;
  ln_psi = 0 ;
  aa.set_etat_zero() ;
  taa.set_etat_zero() ;
  aa_quad = 0 ;
  hh.set_etat_zero() ;

  // Pointers of derived quantities initialized to zero : 
  set_der_0x0() ;

} 
		     

// Copy constructor
//-----------------
Star_rot_Dirac::Star_rot_Dirac(const Star_rot_Dirac& star)
                   : Star(star),
		     spectral_filter(star.spectral_filter),
		     psi4(star.psi4),
		     psi2(star.psi2),
		     qqq(star.qqq),
		     ln_psi(star.ln_psi),
		     j_euler(star.j_euler),
		     v2(star.v2),
		     flat(star.flat),
		     tgamma(star.tgamma),
		     aa(star.aa),
		     taa(star.taa),
		     aa_quad(star.aa_quad),
		     hh(star.hh)
{

  omega = star.omega ;

  // Pointers of derived quantities initialized to zero : 
  set_der_0x0() ;

}


//Constructor from a file 
//------------------------
Star_rot_Dirac::Star_rot_Dirac(Map& mpi, const Eos& eos_i, FILE* fich)
                  : Star(mpi, eos_i, fich),
		    psi4(mpi),
		    psi2(mpi),
		    qqq(mpi, *(mpi.get_mg()), fich),
		    ln_psi(mpi),
		    j_euler(mpi, CON, mpi.get_bvect_spher()),
		    v2(mpi),
		    flat(mpi.flat_met_spher()),
		    tgamma(flat),
		    aa(mpi, CON, mpi.get_bvect_spher()),
		    taa(mpi, COV, mpi.get_bvect_spher()),
		    aa_quad(mpi),
		    hh(mpi, mpi.get_bvect_spher(), flat, fich)
{

  // Pointers of derived quantities initialized to zero 
  //----------------------------------------------------
  set_der_0x0() ;

  fread_be(&spectral_filter, sizeof(int), 1, fich) ;

  // Metric fields are read in the file:
  fread_be(&omega, sizeof(double), 1, fich) ;
  Vector shift_tmp(mpi, mpi.get_bvect_spher(), fich) ;
  beta = shift_tmp ;

  update_metric() ;

  equation_of_state() ;

  hydro_euler() ;




}


                      //------------// 
                      // Destructor //
                      //------------//

Star_rot_Dirac::~Star_rot_Dirac(){

  Star_rot_Dirac::del_deriv() ;
  
}


               //----------------------------------//
               // Management of derived quantities //
               //----------------------------------//

void Star_rot_Dirac::del_deriv() const {

       if (p_angu_mom != 0x0) delete p_angu_mom ;
       if (p_grv2 != 0x0) delete p_grv2 ;
       if (p_grv3 != 0x0) delete p_grv3 ;
       if (p_tsw != 0x0) delete p_tsw ;
       if (p_aplat != 0x0) delete p_aplat ; 
       if (p_r_circ != 0x0) delete p_r_circ ;

       set_der_0x0() ;

       Star::del_deriv() ;

}


void Star_rot_Dirac::set_der_0x0() const {

       p_angu_mom = 0x0 ;
       p_grv2 = 0x0 ;
       p_grv3 = 0x0 ;
       p_tsw = 0x0 ;
       p_aplat = 0x0 ;
       p_r_circ = 0x0 ;

}


void Star_rot_Dirac::del_hydro_euler() {

    j_euler.set_etat_nondef() ;
    v2.set_etat_nondef() ;

    del_deriv() ;
    
    Star::del_hydro_euler() ;

}



                    //---------------//
                    //  Assignment   //
                    //---------------//   

// Assignment to another Star_rot_Dirac
// ------------------------------------

void Star_rot_Dirac::operator=(const Star_rot_Dirac& star) {

     // Assignment of quantities common to all the derived classes of Star
     Star::operator=(star) ;

     // Assignment of proper quantities of class Star_rot_Dirac
     spectral_filter = star.spectral_filter ;
     omega = star.omega ;
     psi4 = star.psi4 ;
     psi2 = star.psi2 ;
     qqq = star.qqq ;
     ln_psi = star.ln_psi ;
     j_euler = star.j_euler ;
     v2 = star.v2 ;
     tgamma = star.tgamma ;
     aa = star.aa ;
     aa_quad = star.aa_quad ;
     hh = star.hh ;

     assert(&flat == &star.flat) ;

     del_deriv() ;    // Deletes all derived quantities

}
     

                      //-----------//
                      //  Outputs  //
                      //-----------//

// Save in a file
// --------------

void Star_rot_Dirac::sauve(FILE* fich) const {

      Star::sauve(fich) ;

      qqq.sauve(fich) ;
      hh.sauve(fich) ;
      fwrite_be(&spectral_filter, sizeof(int), 1, fich) ;
      fwrite_be(&omega, sizeof(double), 1, fich) ;
      beta.sauve(fich) ;

}


// Printing
// ---------

ostream& Star_rot_Dirac::operator>>(ostream& ost) const {

  using namespace Unites ;

     Star::operator>>(ost) ;

     ost << "Rotating star in Dirac gauge" << endl ;

     // Only uniformly rotating star for the moment....
     ost << endl ;
     ost << "Uniformly rotating star" << endl ;
     ost << "-----------------------" << endl ;
     if (spectral_filter > 0)
	 ost << "hydro sources of equations are filtered\n"
	     << "with " << spectral_filter << "-order exponential filter" << endl ;

     double freq = omega/ (2.*M_PI) ;
     ost << "Omega : " << omega * f_unit
         << " rad/s    f : " << freq * f_unit << " Hz" << endl ;
     ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
         << endl ;

     ost << "Error on the virial identity GRV2 : " << endl ;
     ost << "GRV2 = " << grv2() << endl ;
     ost << "Error on the virial identity GRV3 : " << endl ;
     ost << "GRV3 = " << grv3() << endl ;

     ost << "Angular momentum J :    "
         << angu_mom()/( qpig / (4*M_PI) *msol*msol) << " G M_sol^2 / c"
         << endl ;
     ost << "c J / (G M^2) :         "
         << angu_mom()/( qpig / (4*M_PI) * pow(mass_g(), 2.) ) << endl ;

     if (omega != 0.) {
       double mom_iner = angu_mom() / omega ; 
       double mom_iner_38si = mom_iner * rho_unit * (pow(r_unit, double(5.)) 
         / double(1.e38) ) ; 
       ost << "Moment of inertia:       " << mom_iner_38si << " 10^38 kg m^2"
	   << endl ; 
     }

     ost << "Ratio T/W :              " << tsw() << endl ;
     ost << "Circumferential equatorial radius R_circ :     "
	 << r_circ()/km << " km" << endl ;
     ost << "Coordinate equatorial radius r_eq : " << ray_eq()/km << " km"
	 << endl ;
     ost << "Flattening r_pole/r_eq :  " << aplat() << endl ;

     double compact = qpig/(4.*M_PI) * mass_g() / r_circ() ;
     ost << "Compaction parameter M_g / R_circ : " << compact << endl ; 
     

     // More to come here.....

     return ost ;

}
