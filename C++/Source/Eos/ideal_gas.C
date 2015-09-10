/*
 * Methods of the class Ideal_gas.
 *
 * (see file hoteos.h for documentation).
 */

/*
 *   Copyright (c) 2015 Jerome Novak
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


char ideal_gas_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2015/09/10 13:54:04  j_novak
 * Allows for negative entropy in the temperature function.
 *
 * Revision 1.1  2015/03/17 14:20:00  j_novak
 * New class Hot_eos to deal with temperature-dependent EOSs.
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cmath>

// Headers Lorene
#include "hoteos.h"
#include "eos.h"
#include "utilitaires.h"
#include "unites.h"

namespace Lorene {

			//--------------//
			// Constructors //
			//--------------//

  // Standard constructor
  // --------------------
  Ideal_gas::Ideal_gas(double gam0, double kappa, double mass) :
    Hot_eos("Ideal gas EOS"), gam(gam0), kap(kappa), m_0(mass) {
    
    set_auxiliary() ;
    
  }
  
  // Copy constructor
  // ----------------
  Ideal_gas::Ideal_gas(const Ideal_gas& eosi) : 
    Hot_eos(eosi), gam(eosi.gam), kap(eosi.kap), m_0(eosi.m_0) 
  {
    set_auxiliary() ; 
  }  
  

  // Constructor from binary file
  // ----------------------------
  Ideal_gas::Ideal_gas(FILE* fich) : 
    Hot_eos(fich) {
    
    fread_be(&gam, sizeof(double), 1, fich) ;		
    fread_be(&kap, sizeof(double), 1, fich) ;		
    fread_be(&m_0, sizeof(double), 1, fich) ;
    
    set_auxiliary() ; 
    
  }


  // Constructor from a formatted file
  // ---------------------------------
  Ideal_gas::Ideal_gas(ifstream& fich) : 
    Hot_eos(fich) {
    
    fich >> gam ; fich.ignore(80, '\n') ;
    fich >> kap ; fich.ignore(80, '\n') ;
    fich >> m_0 ; fich.ignore(80, '\n') ;

    set_auxiliary() ; 

  }
			//--------------//
			//  Destructor  //
			//--------------//

  Ideal_gas::~Ideal_gas(){}

			//--------------//
			//  Assignment  //
			//--------------//

  void Ideal_gas::operator=(const Ideal_gas& eosi) {
    
    name = eosi.name ; 
    
    gam = eosi.gam ; 
    kap = eosi.kap ; 
    m_0 = eosi.m_0 ;
    
    set_auxiliary() ;
    
  }
  
  
		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

  void Ideal_gas::set_auxiliary() {
    
    gam1 = gam - double(1) ;
    
    unsgam1 = double(1) / gam1 ;
    
    gam1sgamkap = m_0 * gam1 / (gam * kap) ;
    
  }
  
  double Ideal_gas::get_gam() const {
    return gam ;
  }
  
  double Ideal_gas::get_kap() const {
    return kap ; 
  }
  
  double Ideal_gas::get_m_0() const {
    return m_0 ;
  }


			//-------------------------------//
			//  The corresponding cold Eos   //
			//-------------------------------//

  const Eos& Ideal_gas::new_cold_Eos() const {
    
    if (p_cold_eos == 0x0) {
      p_cold_eos = new Eos_poly(gam, kap, m_0) ;
    }
    
    return *p_cold_eos ;
  }


			//------------------------//
			//  Comparison operators  //
			//------------------------//


  bool Ideal_gas::operator==(const Hot_eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
      cout << "The second EOS is not of type Ideal_gas !" << endl ; 
      resu = false ; 
    }
    else {
      
      const Ideal_gas& eos = dynamic_cast<const Ideal_gas&>( eos_i ) ; 
      
      if (eos.gam != gam) {
	cout 
	  << "The two Ideal_gas have different gamma : " << gam << " <-> " 
	  << eos.gam << endl ; 
	resu = false ; 
      }
      
      if (eos.kap != kap) {
	cout 
	  << "The two Ideal_gas have different kappa : " << kap << " <-> " 
	  << eos.kap << endl ; 
	resu = false ; 
      }
      
      if (eos.m_0 != m_0) {
	cout
	  << "The two Ideal_gas have different m_0 : " << m_0 << " <-> "
	  << eos.m_0 << endl ;
	resu = false ;
      }
    }    
    return resu ; 
  }


  bool Ideal_gas::operator!=(const Hot_eos& eos_i) const {
    return !(operator==(eos_i)) ;  
  }
  
  
			//------------//
			//  Outputs   //
			//------------//

  void Ideal_gas::sauve(FILE* fich) const {
    
    Hot_eos::sauve(fich) ; 
    
    fwrite_be(&gam, sizeof(double), 1, fich) ;	
    fwrite_be(&kap, sizeof(double), 1, fich) ;
    fwrite_be(&m_0, sizeof(double), 1, fich) ;

  }

  ostream& Ideal_gas::operator>>(ostream & ost) const {
    
    ost << "Hot EOS of class Ideal_gas (relativistic ideal gas) : " << endl ; 
    ost << "   Adiabatic index gamma :      " << gam << endl ; 
    ost << "   Pressure coefficient kappa : " << kap << 
      " rho_nuc c^2 / n_nuc^gamma" << endl ; 
    ost << "   Mean particle mass : " << m_0 << " m_B" << endl ;

    return ost ;
}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

  double Ideal_gas::nbar_Hs_p(double ent, double sb) const {

    if ( ent > 0. ) {
      return pow( gam1sgamkap * ( exp(ent) - 1. ), unsgam1 ) * exp(-sb) ;
    }
    else {
      return 0 ;
    }
  }

// Energy density from enthalpy
//------------------------------

  double Ideal_gas::ener_Hs_p(double ent, double sb) const {

    if ( ent > 0. ) {
      double nn = pow( gam1sgamkap * ( exp(ent) - 1. ), unsgam1 ) * exp(-sb) ;
      double pp = kap * pow( nn, gam ) * exp(gam1*sb) ;

      return  unsgam1 * pp + m_0 * nn ;
    }
    else {
      return 0. ;
    }
  }

// Pressure from enthalpy
//------------------------

  double Ideal_gas::press_Hs_p(double ent, double sb) const {

    if ( ent > 0. ) {
      double nn = pow( gam1sgamkap * ( exp(ent) - 1. ), unsgam1 ) * exp(-sb) ;

      return kap * pow( nn, gam ) * exp(gam1*sb) ;
    }
    else {
      return 0. ;
    }
  }

// Temperature from enthalpy
//---------------------------

  double Ideal_gas::temp_Hs_p(double ent, double) const {

    using namespace Unites ;

    //    if ( ent > 0. ) {
      return kap * gam1sgamkap * ( exp(ent) - 1. ) ;
      //      return m_u_mev * kap * gam1sgamkap * ( exp(ent) - 1. ) ;
      //}
      //else {
      //return 0 ;
      //}
  }

}
