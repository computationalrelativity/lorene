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
Star_rot_Dirac::Star_rot_Dirac(Map& mpi, int nzet_i, const Eos& eos_i)
                   : Star(mpi, nzet_i, eos_i),
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


//Constructor from a file //## to be more general...
//------------------------
Star_rot_Dirac::Star_rot_Dirac(Map& mpi, const Eos& eos_i, FILE* fich)
                  : Star(mpi, eos_i, fich),
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

  // Star_rot_Dirac parameter
  //-------------------------

  // omega is read in the file:
  fread_be(&omega, sizeof(double), 1, fich) ;


  // Initialized to a static state
  //------------------------------
  
  v2 = 0 ;
  j_euler.set_etat_zero() ;

  // Initialized to a flat case
  //---------------------------

  psi4 = 1 ;
  psi2 = 1 ;
  qqq = 1 ;
  ln_psi = 0 ;
  aa.set_etat_zero() ;
  taa.set_etat_zero() ;
  aa_quad.set_etat_zero() ;
  hh.set_etat_zero() ;

  // Pointers of derived quantities initialized to zero 
  //----------------------------------------------------
  set_der_0x0() ;

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

       set_der_0x0() ;

       Star::del_deriv() ;

}


void Star_rot_Dirac::set_der_0x0() const {

       p_angu_mom = 0x0 ;
       p_grv2 = 0x0 ;

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

      fwrite_be(&omega, sizeof(double), 1, fich) ;

      // What else to save? //## to be more general ...

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

     double freq = omega/ (2.*M_PI) ;
     ost << "Omega : " << omega * f_unit
         << " rad/s    f : " << freq * f_unit << " Hz" << endl ;
     ost << "Rotation period : " << 1000. / (freq * f_unit) << " ms"
         << endl ;

     // More to come here.....

     return ost ;

}
