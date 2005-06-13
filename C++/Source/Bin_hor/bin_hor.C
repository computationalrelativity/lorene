/*
 *  Methods of class Bin_hor
 *
 *   (see file bin_hor.h for documentation)
 *
 */

/*
 *   Copyright (c) 2004-2005 Francois Limousin
 *                           Jose Luis Jaramillo
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


char bin_hor_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2005/06/13 15:47:29  jl_jaramillo
 * Add some quatities in write_global()
 *
 * Revision 1.4  2005/06/09 16:12:04  f_limousin
 * Implementation of amg_mom_adm().
 *
 * Revision 1.3  2005/04/29 14:02:44  f_limousin
 * Important changes : manage the dependances between quantities (for
 * instance psi and psi4). New function write_global(ost).
 *
 * Revision 1.2  2005/03/04 09:38:41  f_limousin
 * Implement the constructor from a file, operator>>, operator<<
 * and function sauve.
 *
 * Revision 1.1  2004/12/29 16:11:02  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// Lorene
#include "nbr_spx.h"
#include "tenseur.h"
#include "tensor.h"
#include "isol_hor.h"
#include "proto.h"
#include "utilitaires.h"
#include "graphique.h"

// Standard constructor
// --------------------

Bin_hor::Bin_hor (Map_af& mp1, Map_af& mp2, int depth_in) :
	hole1(mp1, depth_in), hole2(mp2, depth_in), omega(0){

    holes[0] = &hole1 ;
    holes[1] = &hole2 ;
}

// Copy constructor
// ----------------

Bin_hor::Bin_hor (const Bin_hor& source) :
	    hole1(source.hole1), hole2(source.hole2), omega(source.omega) {
    
    holes[0] = &hole1 ;
    holes[1] = &hole2 ;
    }

// Constructor from a file
// -----------------------
    
Bin_hor::Bin_hor(Map_af& mp1, Map_af& mp2, FILE* fich, 
		   bool partial_read, int depth_in)
    : hole1(mp1, fich, partial_read, depth_in),
      hole2(mp2, fich, partial_read, depth_in),
      omega(0) {

    fread_be(&omega, sizeof(double), 1, fich) ;
    holes[0] = &hole1 ;
    holes[1] = &hole2 ;

}

			    //--------------//
			    //  Destructor  //
			    //--------------//

Bin_hor::~Bin_hor () {
}

                    //-----------------------//
                    // Mutators / assignment //
                    //-----------------------//

void Bin_hor::operator= (const Bin_hor& source) {    
    hole1 = source.hole1 ;
    hole2 = source.hole2 ;
    
    omega = source.omega ;
}

                //------------------//
                //      output      //
                //------------------//

// Printing
// --------
ostream& operator<<(ostream& flux, const Bin_hor& bibi)  {
    bibi >> flux ;
    return flux ;
}
    
ostream& Bin_hor::operator>>(ostream& flux) const {

    flux << "black hole 1" << '\n' ;
    flux << "----------------------------" << '\n' ;
    flux << hole1 << '\n' << '\n' ;
    flux << "black hole 2" << '\n' ;
    flux << "----------------------------" << '\n' ;
    flux << hole2 << '\n' << '\n' ;
    
    cout << "orbital angular velocity  : " << omega << '\n' ;

    return flux ;

}

                //--------------------------//
                //      Save in a file      //
                //--------------------------//

void Bin_hor::sauve(FILE* fich, bool partial_save) const {

    hole1.sauve(fich, partial_save) ;
    hole2.sauve(fich, partial_save) ;
    fwrite_be (&omega, sizeof(double), 1, fich) ;
   
}


//Initialisation : Sum of two static BH
void Bin_hor::init_bin_hor() {
    set_omega (0) ;
    hole1.init_bhole() ;
    hole2.init_bhole() ;
    
    hole1.psi_comp(hole2) ;
    hole2.psi_comp(hole1) ;
    
    hole1.n_comp(hole2) ;
    hole2.n_comp(hole1) ;
    
    decouple() ;
}


void Bin_hor::write_global(ostream& ost) const {

  double beta = hole1.get_mp().get_ori_x() - hole2.get_mp().get_ori_x() ;
  beta /= hole1.get_radius() ;
  double mass_adm = adm_mass() ;
  double mass_komar = komar_mass() ;
  double mass_area = sqrt(hole1.area_hor()/16/M_PI) + 
      sqrt(hole2.area_hor()/16/M_PI) ;
  double J_adm = ang_mom_adm() ;
  double J_hor = ang_mom_hor() ; //hole1.ang_mom_hor() + hole2.ang_mom_hor() ;
  double j1 = hole1.ang_mom_hor() ;
  double j2 = hole2.ang_mom_hor() ;
  double mass_ih1 = hole1.mass_hor() ;
  double mass_ih2 = hole2.mass_hor() ;
  double mass_ih = mass_ih1 + mass_ih2 ;
  double omega1 = hole1.omega_hor() ;
  double omega2 = hole2.omega_hor() ;


  ost.precision(8) ;
  ost << "# beta  omega  Mass_ADM  Mass_K  M_area  J_ADM  J_hor" << endl ;
  ost << beta << " " ;
  ost << omega << " " ;
  ost << mass_adm << " " ;
  ost << mass_komar << " " ;
  ost << mass_area << " " ;
  ost << J_adm << " " ;
  ost << J_hor << endl ;
  ost << "# mass_ih1  mass_ih2  mass_ih  j1  J2  omega1  omega2" << endl ;
  ost << mass_ih1 << " " ;
  ost << mass_ih2 << " " ;
  ost << mass_ih << " " ;
  ost << j1 << " " ;
  ost << j2 << " " ;
  ost << omega1 << " " ;
  ost << omega2 << endl ;
  ost << "# ADM_mass/M_area  J/M_area2  omega*M_area" << endl ;
  ost << mass_adm / mass_area << " " ;
  ost << J_adm /mass_area / mass_area << " " ;
  ost << omega * mass_area << endl ;

}
      
