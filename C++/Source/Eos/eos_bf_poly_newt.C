/*
 * Methods of the class Eos_bf_poly_newt.
 *
 * (see file eos_bifluid.h for documentation).
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


char eos_bf_poly_newt_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/01/11 14:09:34  j_novak
 * Added newtonian version for 2-fluid stars
 *
 *
 * $Header$
 *
 */

// Headers C++
#include <iostream.h>
#include <fstream.h>

// Headers C
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Headers Lorene
#include "eos_bifluid.h"
#include "eos.h"
#include "cmp.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor with gam1 = gam2 = 2, 
// gam3 = gam4 = gam5 = gam6 = 1, m_1 = 1 and m_2 =1
// -------------------------------------------------
Eos_bf_poly_newt::Eos_bf_poly_newt(double kappa1, double kappa2, double kappa3,
				   double bet) :
  Eos_bf_poly(kappa1, kappa2, kappa3, bet) {
  set_name("bi-fluid polytropic non-relativistic EOS") ;
}  

// Standard constructor with everything specified
// -----------------------------------------------
Eos_bf_poly_newt::Eos_bf_poly_newt(double gamma1, double gamma2, double gamma3,
			 double gamma4, double gamma5, double gamma6,
			 double kappa1, double kappa2, double kappa3,
			 double bet, double mass1, double mass2) : 
  Eos_bf_poly(gamma1, gamma2, gamma3, gamma4, gamma5, gamma6,
	      kappa1, kappa2, kappa3, bet, mass1, mass2) {
  set_name("bi-fluid polytropic non-relativistic EOS") ;
} 
  
// Copy constructor
// ----------------
Eos_bf_poly_newt::Eos_bf_poly_newt(const Eos_bf_poly_newt& eosi) : 
  Eos_bf_poly(eosi) {} 
  

// Constructor from binary file
// ----------------------------
Eos_bf_poly_newt::Eos_bf_poly_newt(FILE* fich) : 
  Eos_bf_poly(fich) {} 

// Constructor from a formatted file
// ---------------------------------
Eos_bf_poly_newt::Eos_bf_poly_newt(ifstream& fich) : 
  Eos_bf_poly(fich) {} 

			//--------------//
			//  Destructor  //
			//--------------//

Eos_bf_poly_newt::~Eos_bf_poly_newt(){
    
    // does nothing
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Eos_bf_poly_newt::operator=(const Eos_bf_poly_newt& eosi) {
    
    set_name(eosi.name) ; 
    
    gam1 = eosi.gam1 ; 
    gam2 = eosi.gam2 ; 
    gam3 = eosi.gam3 ; 
    kap1 = eosi.kap1 ; 
    kap2 = eosi.kap2 ; 
    kap3 = eosi.kap3 ;
    beta = eosi.beta ;
    m_1 = eosi.m_1 ; 
    m_2 = eosi.m_2 ; 
    
    set_auxiliary() ; 
    
}

			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_bf_poly_newt::operator==(const Eos_bifluid& eos_i) const {
    
  bool resu = true ; 
  
  if ( eos_i.identify() != identify() ) {
    cout << "The second EOS is not of type Eos_bf_poly_newt !" << endl ; 
    resu = false ; 
  }
  else{
    
    const Eos_bf_poly_newt& eos = 
      dynamic_cast<const Eos_bf_poly_newt&>( eos_i ) ; 
    
    if ((eos.gam1 != gam1)||(eos.gam2 != gam2)||(eos.gam3 != gam3)
	||(eos.gam4 != gam4)||(eos.gam5 != gam5)||(eos.gam6 != gam6)) {
      cout 
	<< "The two Eos_bf_poly_newt have different gammas : " << gam1 << " <-> " 
	<< eos.gam1 << ", " << gam2 << " <-> " 
	<< eos.gam2 << ", " << gam3 << " <-> " 
	<< eos.gam3 << ", " << gam4 << " <-> " 
	<< eos.gam4 << ", " << gam5 << " <-> " 
	<< eos.gam5 << ", " << gam6 << " <-> " 
	<< eos.gam6 << endl ; 
      resu = false ; 
    }
	
    if ((eos.kap1 != kap1)||(eos.kap2 != kap2)|| (eos.kap3 != kap3)){
      cout 
	<< "The two Eos_bf_poly_newt have different kappas : " << kap1 << " <-> " 
	<< eos.kap1 << ", " << kap2 << " <-> " 
	<< eos.kap2 << ", " << kap3 << " <-> " 
	<< eos.kap3 << endl ; 
      resu = false ; 
    }
    
    if (eos.beta != beta) {
      cout 
	<< "The two Eos_bf_poly_newt have different betas : " << beta << " <-> " 
	<< eos.beta << endl ; 
      resu = false ; 
    }

    if ((eos.m_1 != m_1)||(eos.m_2 != m_2)) {
      cout 
	<< "The two Eos_bf_poly_newt have different masses : " << m_1 << " <-> " 
	<< eos.m_1 << ", " << m_2 << " <-> " 
	<< eos.m_2 << endl ; 
      resu = false ; 
    }
    
  }
  
  return resu ; 
  
}

bool Eos_bf_poly_newt::operator!=(const Eos_bifluid& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}


			//------------//
			//  Outputs   //
			//------------//

void Eos_bf_poly_newt::sauve(FILE* fich) const {

    Eos_bf_poly::sauve(fich) ; 
}

ostream& Eos_bf_poly_newt::operator>>(ostream & ost) const {
    
    ost << "EOS of class Eos_bf_poly_newt (non-relativistic polytrope) : " << endl ; 
    ost << "   Adiabatic index gamma1 :      " << gam1 << endl ; 
    ost << "   Adiabatic index gamma2 :      " << gam2 << endl ; 
    ost << "   Adiabatic index gamma3 :      " << gam3 << endl ; 
    ost << "   Adiabatic index gamma4 :      " << gam4 << endl ; 
    ost << "   Adiabatic index gamma5 :      " << gam5 << endl ; 
    ost << "   Adiabatic index gamma6 :      " << gam6 << endl ; 
    ost << "   Pressure coefficient kappa1 : " << kap1 << 
	   " rho_nuc c^2 / n_nuc^gamma" << endl ; 
    ost << "   Pressure coefficient kappa2 : " << kap2 << 
	   " rho_nuc c^2 / n_nuc^gamma" << endl ; 
    ost << "   Pressure coefficient kappa3 : " << kap3 << 
	   " rho_nuc c^2 / n_nuc^gamma" << endl ; 
    ost << "   Coefficient beta : " << beta << 
	   " rho_nuc c^2 / n_nuc^gamma" << endl ; 
    ost << "   Mean particle 1 mass : " << m_1 << " m_B" << endl ;
    ost << "   Mean particle 2 mass : " << m_2 << " m_B" << endl ;
    
    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon densities from enthalpies 
//---------------------------------

void  Eos_bf_poly_newt::nbar_ent_p(const double ent1, const double ent2, 
			      const double delta2, double& nbar1, 
			      double& nbar2, bool tronc) const {  

  if ((gam1 == double(2)) && (gam2 == double(2)) && (gam3 == double(1))
      && (gam4 == double(1))) {
    
    assert ((gam5 == double(1)) || (gam5 == double(0))) ;
    assert ((gam6 == double(1)) || (gam6 == double(0))) ;

    if ((gam5 == double(1))&&(gam6 == double(1))) {
      double kpd = kap3+beta*delta2 ;
      double determ = kap1*kap2 - kpd*kpd ;
    
      nbar1 = (kap2*m_1*ent1 - kpd*m_2*ent2) / determ ;
      nbar2 = (kap1*m_2*ent2 - kpd*m_1*ent1) / determ ;
      
    }
    else {
      double determ = kap1*kap2 - kap3*kap3 ;
      double mu_1 = ent1*m_1 - gam5*beta*delta2 ;
      double mu_2 = ent2*m_2 - gam6*beta*delta2 ;

      nbar1 = (kap2*mu_1 - kap3*mu_2) / determ ;
      nbar2 = (kap1*mu_2 - kap3*mu_1) / determ ;
    }
      
    if (tronc) {
      if (nbar1 < 0.) {
	nbar1 = 0. ;
	nbar2 = ent2*m_2 / kap2 ;
	nbar2 = nbar2 < 0. ? 0. : nbar2 ;
	return ;
      }
      if (nbar2 < 0.) {
	nbar2 = 0. ;
	nbar1 = ent1*m_1 / kap1 ;
	nbar1 = nbar1 < 0. ? 0. : nbar1 ;
	return ;
      }
    }
    return ;
  }
  else {
    cout << "Eos_bf_poly_newt::nbar_ent_p: the case gamma_i != 2" << endl;
    cout << " is not implemented yet. Sorry!" << endl ;
    abort() ;
  }
}

// Energy density from baryonic densities
//---------------------------------------

double Eos_bf_poly_newt::ener_nbar_p(const double nbar1, const double nbar2, 
				const double delta2) const {
    
    if (( nbar1 > double(0) ) || ( nbar2 > double(0))) {
      
      double n1 = (nbar1>double(0) ? nbar1 : double(0)) ;
      double n2 = (nbar2>double(0) ? nbar2 : double(0)) ;
      double x2 = ((nbar1>double(0))&&(nbar2>double(0))) ? delta2 : 0 ;

      double resu = 0.5*kap1*pow(n1, gam1) + 0.5*kap2*pow(n2,gam2)
	+ kap3*pow(n1,gam3)*pow(n2,gam4) 
	+ x2*beta*pow(n1,gam5)*pow(n2,gam6) ;
      return resu ;
    }
    else return 0 ;
}

// Pressure from baryonic densities
//---------------------------------

double Eos_bf_poly_newt::press_nbar_p(const double nbar1, const double nbar2,
				const double delta2) const {
    
  if (( nbar1 > double(0) ) || ( nbar2 > double(0))) {
    
    double n1 = (nbar1>double(0) ? nbar1 : double(0)) ;
    double n2 = (nbar2>double(0) ? nbar2 : double(0)) ;
    double x2 = ((nbar1>double(0))&&(nbar2>double(0))) ? delta2 : 0 ;
    
    double resu = 0.5*gam1m1*kap1*pow(n1,gam1) + 0.5*gam2m1*kap2*pow(n2,gam2)
      + gam34m1*kap3*pow(n1,gam3)*pow(n2,gam4) + 
      x2*gam56m1*beta*pow(n1,gam5)*pow(n2,gam6) ;
    
    return resu ;
  }
  else return 0 ;
}

// Derivatives of energy
//----------------------

double Eos_bf_poly_newt::get_K11(const double n1, const double n2, const
			       double)  const 
{
  double xx ;
  if (n1 <= 0.) xx = 0. ;
  else xx = m_1/n1 -2*beta*pow(n1,gam5-2)*pow(n2,gam6) ;

  return xx ;
}

double Eos_bf_poly_newt::get_K22(const double n1, const double n2, const
			       double ) const  
{
  double xx ;
  if (n2 <= 0.) xx = 0. ;
  else xx = m_2/n2 - 2*beta*pow(n1,gam5)*pow(n2,gam6-2) ;

  return xx ;
}

double Eos_bf_poly_newt::get_K12(const double n1, const double n2, const
			       double) const  
{
  double xx ;
  if ((n1 <= 0.) || (n2 <= 0.)) xx = 0.; 
  else xx = 2*beta*pow(n1,gam5-1)*pow(n2,gam6-1) ;

  return xx ;
}

// Conversion functions
// ---------------------

Eos* Eos_bf_poly_newt::trans2Eos() const {

  Eos_poly_newt* eos_simple = new Eos_poly_newt(gam1, kap1) ;
  return eos_simple ;

}
       
