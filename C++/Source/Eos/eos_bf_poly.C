/*
 * Methods of the class Eos_bf_poly.
 *
 * (see file eos_bifluid.h for documentation).
 */

/*
 *   Copyright (c) 2001 Jerome Novak
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


char eos_bf_poly_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.2  2001/11/29 15:05:26  j_novak
 * The entrainment term in 2-fluid eos is modified
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.4  2001/08/31  15:48:50  novak
 * The flag tronc has been added to nbar_ent_p
 *
 * Revision 1.3  2001/08/27 09:52:21  novak
 * New version of formulas
 *
 * Revision 1.2  2001/06/22 15:36:46  novak
 * Modification de trans2Eos
 *
 * Revision 1.1  2001/06/21 15:24:46  novak
 * Initial revision
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
Eos_bf_poly::Eos_bf_poly(double kappa1, double kappa2, double kappa3,
			 double bet) :
  Eos_bifluid("bi-fluid polytropic EOS"), 
  gam1(2), gam2(2),gam3(1),gam4(1),gam5(1),
  gam6(1),kap1(kappa1), kap2(kappa2), kap3(kappa3),beta(bet),
  m_1(1),m_2(1) {

  set_auxiliary() ; 

}  

// Standard constructor with everything specified
// -----------------------------------------------
Eos_bf_poly::Eos_bf_poly(double gamma1, double gamma2, double gamma3,
			 double gamma4, double gamma5, double gamma6,
			 double kappa1, double kappa2, double kappa3,
			 double bet, double mass1, double mass2) : 
  Eos_bifluid("bi-fluid polytropic EOS"), 
  gam1(gamma1),gam2(gamma2),gam3(gamma3),gam4(gamma4),gam5(gamma5),
  gam6(gamma6),kap1(kappa1),kap2(kappa2),kap3(kappa3),beta(bet), 
  m_1(mass1),m_2(mass2) {

    set_auxiliary() ; 

}  
  
// Copy constructor
// ----------------
Eos_bf_poly::Eos_bf_poly(const Eos_bf_poly& eosi) : 
	Eos_bifluid(eosi), 
	gam1(eosi.gam1), gam2(eosi.gam2), gam3(eosi.gam3),
	gam4(eosi.gam4), gam5(eosi.gam5), gam6(eosi.gam6),
	kap1(eosi.kap1), kap2(eosi.kap2), kap3(eosi.kap3),
	beta(eosi.beta),m_1(eosi.m_1), m_2(eosi.m_2) {

    set_auxiliary() ; 

}  
  

// Constructor from binary file
// ----------------------------
Eos_bf_poly::Eos_bf_poly(FILE* fich) : 
	Eos_bifluid(fich) {
        
    fread_be(&gam1, sizeof(double), 1, fich) ;		
    fread_be(&gam2, sizeof(double), 1, fich) ;		
    fread_be(&gam3, sizeof(double), 1, fich) ;		
    fread_be(&gam4, sizeof(double), 1, fich) ;		
    fread_be(&gam5, sizeof(double), 1, fich) ;		
    fread_be(&gam6, sizeof(double), 1, fich) ;		
    fread_be(&kap1, sizeof(double), 1, fich) ;		
    fread_be(&kap2, sizeof(double), 1, fich) ;		
    fread_be(&kap3, sizeof(double), 1, fich) ;		
    fread_be(&beta, sizeof(double), 1, fich) ;		
    fread_be(&m_1, sizeof(double), 1, fich) ;		
    fread_be(&m_2, sizeof(double), 1, fich) ;		
    
    set_auxiliary() ; 

}


// Constructor from a formatted file
// ---------------------------------
Eos_bf_poly::Eos_bf_poly(ifstream& fich) : 
	Eos_bifluid(fich) {

    char blabla[80] ;
        
    fich >> gam1 ; fich.getline(blabla, 80) ;
    fich >> gam2 ; fich.getline(blabla, 80) ;
    fich >> gam3 ; fich.getline(blabla, 80) ;
    fich >> gam4 ; fich.getline(blabla, 80) ;
    fich >> gam5 ; fich.getline(blabla, 80) ;
    fich >> gam6 ; fich.getline(blabla, 80) ;
    fich >> kap1 ; fich.getline(blabla, 80) ;
    fich >> kap2 ; fich.getline(blabla, 80) ;
    fich >> kap3 ; fich.getline(blabla, 80) ;
    fich >> beta ; fich.getline(blabla, 80) ;
    fich >> m_1 ; fich.getline(blabla, 80) ;
    fich >> m_2 ; fich.getline(blabla, 80) ;
    cout << m_2 << endl ;
    
    set_auxiliary() ; 

}
			//--------------//
			//  Destructor  //
			//--------------//

Eos_bf_poly::~Eos_bf_poly(){
    
    // does nothing
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Eos_bf_poly::operator=(const Eos_bf_poly& eosi) {
    
    set_name(eosi.name) ; 
    
    gam1 = eosi.gam1 ; 
    gam2 = eosi.gam2 ; 
    gam3 = eosi.gam3 ; 
    kap1 = eosi.kap1 ; 
    kap2 = eosi.kap2 ; 
    kap3 = eosi.kap3 ; 
    m_1 = eosi.m_1 ; 
    m_2 = eosi.m_2 ; 
    
    set_auxiliary() ; 
    
}


		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

void Eos_bf_poly::set_auxiliary() {
    
    gam1m1 = gam1 - double(1) ; 
    gam2m1 = gam2 - double(1) ; 
    gam34m1 = gam3 + gam4 - double(1) ; 
    gam56m1 = gam5 + gam6 - double(1) ;

    if (fabs(kap3*kap3-kap2*kap1) < 1.e-15) {
      cout << "WARNING!: Eos_bf_poly: the parameters are degenerate!" << endl ;
      abort() ;
    }

}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_bf_poly::operator==(const Eos_bifluid& eos_i) const {
    
  bool resu = true ; 
    
  if ( eos_i.identify() != identify() ) {
    cout << "The second EOS is not of type Eos_bf_poly !" << endl ; 
    resu = false ; 
  }
  else{
    
    const Eos_bf_poly& eos = dynamic_cast<const Eos_bf_poly&>( eos_i ) ; 
    
    if ((eos.gam1 != gam1)||(eos.gam2 != gam2)||(eos.gam3 != gam3)
	||(eos.gam4 != gam4)||(eos.gam5 != gam5)||(eos.gam6 != gam6)) {
      cout 
	<< "The two Eos_bf_poly have different gammas : " << gam1 << " <-> " 
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
	<< "The two Eos_bf_poly have different kappas : " << kap1 << " <-> " 
	<< eos.kap1 << ", " << kap2 << " <-> " 
	<< eos.kap2 << ", " << kap3 << " <-> " 
	<< eos.kap3 << endl ; 
      resu = false ; 
    }
    
    if (eos.beta != beta) {
      cout 
	<< "The two Eos_bf_poly have different betas : " << beta << " <-> " 
	<< eos.beta << endl ; 
      resu = false ; 
    }

    if ((eos.m_1 != m_1)||(eos.m_2 != m_2)) {
      cout 
	<< "The two Eos_bf_poly have different masses : " << m_1 << " <-> " 
	<< eos.m_1 << ", " << m_2 << " <-> " 
	<< eos.m_2 << endl ; 
      resu = false ; 
    }
    
  }
  
  return resu ; 
  
}

bool Eos_bf_poly::operator!=(const Eos_bifluid& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}


			//------------//
			//  Outputs   //
			//------------//

void Eos_bf_poly::sauve(FILE* fich) const {

    Eos_bifluid::sauve(fich) ; 
    
    fwrite_be(&gam1, sizeof(double), 1, fich) ;	
    fwrite_be(&gam2, sizeof(double), 1, fich) ;	
    fwrite_be(&gam3, sizeof(double), 1, fich) ;	
    fwrite_be(&gam4, sizeof(double), 1, fich) ;	
    fwrite_be(&gam5, sizeof(double), 1, fich) ;	
    fwrite_be(&gam6, sizeof(double), 1, fich) ;	
    fwrite_be(&kap1, sizeof(double), 1, fich) ;	
    fwrite_be(&kap2, sizeof(double), 1, fich) ;	
    fwrite_be(&kap3, sizeof(double), 1, fich) ;	
    fwrite_be(&beta, sizeof(double), 1, fich) ;	
    fwrite_be(&m_1, sizeof(double), 1, fich) ;	
    fwrite_be(&m_2, sizeof(double), 1, fich) ;	
   
}

ostream& Eos_bf_poly::operator>>(ostream & ost) const {
    
    ost << "EOS of class Eos_bf_poly (relativistic polytrope) : " << endl ; 
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

void  Eos_bf_poly::nbar_ent_p(const double ent1, const double ent2, 
			      const double delta2, double& nbar1, 
			      double& nbar2, bool tronc) const {  

  if ((gam1 == double(2)) && (gam2 == double(2)) && (gam3 == double(1))
      && (gam4 == double(1)) && (gam5 == double(1)) 
      && (gam6 == double(1))) {

    double kpd = kap3+beta*delta2 ;
    double determ = kap1*kap2 - kpd*kpd ;
    
    nbar1 = (kap2*(exp(ent1) - m_1) - kpd*(exp(ent2) - m_2)) / determ ;
    nbar2 = (kap1*(exp(ent2) - m_2) - kpd*(exp(ent1) - m_1)) / determ ;

    if (tronc) {
      if (nbar1 < 0.) {
	nbar1 = 0. ;
	nbar2 = (exp(ent2) - m_2)/kap2 ;
	nbar2 = nbar2 < 0. ? 0. : nbar2 ;
	return ;
      }
      if (nbar2 < 0.) {
	nbar2 = 0. ;
	nbar1 = (exp(ent1) - m_1)/kap1 ;
	nbar1 = nbar1 < 0. ? 0. : nbar1 ;
	return ;
      }
    }
    return ;
  }
  else {
    cout << "Eos_bf_poly::nbar_ent_p: the case gamma_i != 2" << endl;
    cout << " is not implemented yet. Sorry!" << endl ;
    abort() ;
  }
}

// Energy density from baryonic densities
//---------------------------------------

double Eos_bf_poly::ener_nbar_p(const double nbar1, const double nbar2, 
				const double delta2) const {
    
    if (( nbar1 > double(0) ) || ( nbar2 > double(0))) {
      
      double n1 = (nbar1>double(0) ? nbar1 : double(0)) ;
      double n2 = (nbar2>double(0) ? nbar2 : double(0)) ;
      double x2 = ((nbar1>double(0))&&(nbar2>double(0))) ? delta2 : 0 ;

      double resu = 0.5*kap1*pow(n1, gam1) + 0.5*kap2*pow(n2,gam2)
	+ kap3*pow(n1,gam3)*pow(n2,gam4) + m_1*n1 + m_2*n2
	+ x2*beta*pow(n1,gam5)*pow(n2,gam6) ;
      return resu ;
    }
    else return 0 ;
}

// Pressure from baryonic densities
//---------------------------------

double Eos_bf_poly::press_nbar_p(const double nbar1, const double nbar2,
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

double Eos_bf_poly::get_K11(const double n1, const double n2, const
			       double delta2)  const 
{
  double xx ;
  if (n1 <= 0.) {
    xx = 0. ;
  }
  else {
    xx = 0.5*gam1*kap1 * pow(n1,gam1 - 2) + m_1/n1 + 
      gam3*kap3 * pow(n1,gam3 - 2) * pow(n2,gam4) + 
      delta2*gam5*beta * pow(n1,gam5 - 2)*pow(n2, gam6) ;
  }
  return xx ;
}

double Eos_bf_poly::get_K22(const double n1, const double n2, const
			       double delta2) const  
{
  double xx ;
  if (n2 <= 0.) {
    xx = 0. ;
  }
  else {
    xx = 0.5*gam2*kap2 * pow(n2,gam2 - 2) + m_2/n2 + 
      gam3*kap3 * pow(n2,gam4 - 2) * pow(n1,gam3) + 
      delta2*gam6*beta * pow(n1, gam5) * pow(n2,gam6 - 2) ;
  }
  return xx ;
}

double Eos_bf_poly::get_K12(const double n1, const double n2, const
			       double delta2) const  
{
  double xx ;
  if ((n1 <= 0.) || (n2 <= 0.)) { xx = 0.; }
  else { 
    double gamma_delta3 = pow(1-delta2,-1.5) ;
    xx = 2*beta / gamma_delta3 ;
  }
  return xx ;
}

// Conversion functions
// ---------------------

Eos* Eos_bf_poly::trans2Eos(bool relat) const {

  if (relat) {
    Eos_poly* eos_simple = new Eos_poly(gam1, kap1, m_1) ;
    return eos_simple ;
  }
  else {
    Eos_poly_newt* eos_simple = new Eos_poly_newt(gam1, kap1) ;
    return eos_simple ;
  }
}
       
