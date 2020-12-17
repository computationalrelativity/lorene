/*
 * Methods of the class Dyn_eos_poly.
 *
 * (see file dyneos.h for documentation).
 */

/*
 *   Copyright (c) 2019 Jerome Novak
 *             (c) 2000 Eric Gourgoulhon for Eos classes
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
 * Revision 1.3  2020/12/17 17:00:27  j_novak
 * Output of sound speed squared, instead of sound speed.
 *
 * Revision 1.2  2019/12/19 13:38:37  e_declerck
 * Correction pour la formaule de la vitesse du son dans dyneos_poly.C
 *
 * Revision 1.1  2019/12/06 14:30:50  j_novak
 * New classes Dyn_eos... for cold Eos's with baryon density as input.
 *
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cstring>
#include <cmath>

// Headers Lorene
#include "dyneos.h"
#include "utilitaires.h"

namespace Lorene {

			//--------------//
			// Constructors //
			//--------------//

  // Standard constructor with m_0 = 1 and mu_0 = 1
  // -----------------------------------------------
  Dyn_eos_poly::Dyn_eos_poly(double gam0, double kappa) : 
    Dyn_eos("Relativistic polytropic EOS"), 
    gam(gam0), kap(kappa), m_0(double(1)), mu_0(double(1))
  {  
    set_auxiliary() ;  
  }  

  // Standard constructor with mu_0 = 1
  // ----------------------------------
  Dyn_eos_poly::Dyn_eos_poly(double gam0, double kappa, double mass) :
    Dyn_eos("Relativistic polytropic EOS"),
    gam(gam0), kap(kappa), m_0(mass), mu_0(double(1))
  {
    set_auxiliary() ;
  }

  // Standard constructor with mu_0 = 1
  // ----------------------------------
  Dyn_eos_poly::Dyn_eos_poly(double gam0, double kappa, double mass, double mu_zero) :
    Dyn_eos("Relativistic polytropic EOS"),
    gam(gam0), kap(kappa), m_0(mass), mu_0(mu_zero)
  {
    set_auxiliary() ;
  }

  // Copy constructor
  // ----------------
  Dyn_eos_poly::Dyn_eos_poly(const Dyn_eos_poly& eosi) : 
    Dyn_eos(eosi), 
    gam(eosi.gam), kap(eosi.kap), m_0(eosi.m_0), mu_0(eosi.mu_0)
  {
    set_auxiliary() ; 
  }  
  

  // Constructor from binary file
  // ----------------------------
  Dyn_eos_poly::Dyn_eos_poly(FILE* fich) : Dyn_eos(fich)
  {      
    fread_be(&gam, sizeof(double), 1, fich) ;		
    fread_be(&kap, sizeof(double), 1, fich) ;		
    fread_be(&m_0, sizeof(double), 1, fich) ;
    
    if (m_0 < 0) // to ensure compatibility with previous version (revision <= 1.2)
      {       	 //  of Eos_poly
        m_0 = fabs( m_0 ) ;
        fread_be(&mu_0, sizeof(double), 1, fich) ;
      }
    else
      {
        mu_0 = double(1) ;
      }
    set_auxiliary() ; 
  }


  // Constructor from a formatted file
  // ---------------------------------
  Dyn_eos_poly::Dyn_eos_poly(ifstream& fich) : Dyn_eos(fich) {

    fich >> gam ; fich.ignore(1000, '\n') ; 
    fich >> kap ; fich.ignore(1000, '\n') ; 
    fich >> m_0 ; fich.ignore(1000, '\n') ;

    if (m_0 < 0)  // to ensure compatibility with previous version (revision <= 1.2)
      {           //  of Eos_poly
        m_0 = fabs( m_0 ) ;
        fich >> mu_0 ; fich.ignore(1000, '\n') ; 
      }
    else
      {
        mu_0 = double(1) ;
      }
    set_auxiliary() ; 
  }
			//--------------//
			//  Destructor  //
			//--------------//

  Dyn_eos_poly::~Dyn_eos_poly(){}  // does nothing
    
			//--------------//
			//  Assignment  //
			//--------------//

void Dyn_eos_poly::operator=(const Dyn_eos_poly& eosi) {
    
    set_name(eosi.name) ; 
    
    gam = eosi.gam ; 
    kap = eosi.kap ; 
    m_0 = eosi.m_0 ;
    mu_0 = eosi.mu_0 ;

    set_auxiliary() ;

}


		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

  void Dyn_eos_poly::set_auxiliary()
  {  
    gam1 = gam - double(1) ; 
    kapsgam1 = kap / gam1 ;
    gamkapsgam1 = gam * kap / (gam1*m_0) ;
    rel_mu_0 = mu_0 / m_0 ;
  }

  double Dyn_eos_poly::get_gam() const
  {
    return gam ;
  }

  double Dyn_eos_poly::get_kap() const
  {
    return kap ; 
  }

  double Dyn_eos_poly::get_m_0() const
  {
    return m_0 ;
  }

  double Dyn_eos_poly::get_mu_0() const
  {
    return mu_0 ;
  }


			//------------------------//
			//  Comparison operators  //
			//------------------------//


  bool Dyn_eos_poly::operator==(const Dyn_eos& eos_i) const
  {  
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() )
      {
	cout << "The second EOS is not of type Dyn_eos_poly !" << endl ; 
	resu = false ; 
      }
    else
      {
	const Dyn_eos_poly& eos = dynamic_cast<const Dyn_eos_poly&>( eos_i ) ; 
	if (eos.gam != gam)
	  {
	    cout 
	      << "The two Dyn_eos_poly have different gamma : " << gam << " <-> " 
	      << eos.gam << endl ; 
	    resu = false ; 
	  }
	if (eos.kap != kap)
	  {
	    cout 
	      << "The two Dyn_eos_poly have different kappa : " << kap << " <-> " 
	      << eos.kap << endl ; 
	    resu = false ; 
	  }
	if (eos.m_0 != m_0)
	  {
	    cout
	      << "The two Dyn_eos_poly have different m_0 : " << m_0 << " <-> "
	      << eos.m_0 << endl ;
	    resu = false ;
	  }
	if (eos.mu_0 != mu_0)
	  {
	    cout
	      << "The two Dyn_eos_poly have different mu_0 : " << mu_0 << " <-> "
	      << eos.mu_0 << endl ;
	    resu = false ;
	  }
      }
    return resu ; 
  }

  bool Dyn_eos_poly::operator!=(const Dyn_eos& eos_i) const
  {
    return !(operator==(eos_i)) ;    
  }


			//------------//
			//  Outputs   //
			//------------//

  void Dyn_eos_poly::sauve(FILE* fich) const
  {
    Dyn_eos::sauve(fich) ; 
    
    fwrite_be(&gam, sizeof(double), 1, fich) ;	
    fwrite_be(&kap, sizeof(double), 1, fich) ;
    double tempo = - m_0 ;
    fwrite_be(&tempo, sizeof(double), 1, fich) ;
    fwrite_be(&mu_0, sizeof(double), 1, fich) ;
  }

  ostream& Dyn_eos_poly::operator>>(ostream & ost) const
  {  
    ost << "EOS of class Dyn_eos_poly (relativistic polytrope) : " << endl ; 
    ost << "   Adiabatic index gamma :      " << gam << endl ; 
    ost << "   Pressure coefficient kappa : " << kap << 
	   " rho_nuc c^2 / n_nuc^gamma" << endl ; 
    ost << "   Mean particle mass : " << m_0 << " m_B" << endl ;
    ost << "   Relativistic chemical potential at zero pressure : "
	<< mu_0 << " m_B c^2" << endl ;
    return ost ;
  }


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

  // Baryon log-enthalpy from baryon density
  //-----------------------------------------
  double Dyn_eos_poly::ent_nbar_p(double nbar, const Param* ) const
  {
    if ( nbar > 0. )
      return log( gamkapsgam1*pow(nbar, gam1) + rel_mu_0 ) ;
    else
      return 0. ;
  }

  // Energy density from baryon density
  //------------------------------------
  double Dyn_eos_poly::ener_nbar_p(double nbar, const Param* ) const
  {
    if ( nbar > 0. )
      return kapsgam1*pow(nbar, gam) + mu_0*nbar ;
    else
      return 0. ;
  }

  // Pressure from baryon density
  //------------------------------
  double Dyn_eos_poly::press_nbar_p(double nbar, const Param* ) const
  {
    if ( nbar > 0. )
      return kap * pow( nbar, gam ) ;
    else
      return 0 ;
  }

  // Sound speed from baryon density
  //---------------------------------
  double Dyn_eos_poly::csound_square_nbar_p(double nbar, const Param* ) const
  {
    if ( nbar > 0. )
      {
	double ngam = pow(nbar, gam1) ;
	return kap*gam*ngam / ( gam*kapsgam1*ngam + mu_0 ) ;
      }
    else
      return 0. ;
  }


}
