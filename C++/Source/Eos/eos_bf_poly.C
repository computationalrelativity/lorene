/*
 * Methods of the class Eos_bf_poly.
 *
 * (see file eos_bifluid.h for documentation).
 */

/*
 *   Copyright (c) 2001-2002 Jerome Novak
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
 * Revision 1.13  2003/11/18 18:28:38  r_prix
 * moved particle-masses m_1, m_2 of the two fluids into class eos_bifluid (from eos_bf_poly)
 *
 * Revision 1.12  2003/03/17 10:28:04  j_novak
 * Removed too strong asserts on beta and kappa
 *
 * Revision 1.11  2003/02/07 10:47:43  j_novak
 * The possibility of having gamma5 xor gamma6 =0 has been introduced for
 * tests. The corresponding parameter files have been added.
 *
 * Revision 1.10  2002/10/16 14:36:34  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.9  2002/05/31 16:13:36  j_novak
 * better inversion for eos_bifluid
 *
 * Revision 1.8  2002/05/10 09:26:52  j_novak
 * Added new class Et_rot_mag for magnetized rotating neutron stars (under development)
 *
 * Revision 1.7  2002/05/02 15:16:22  j_novak
 * Added functions for more general bi-fluid EOS
 *
 * Revision 1.6  2002/04/05 09:09:36  j_novak
 * The inversion of the EOS for 2-fluids polytrope has been modified.
 * Some errors in the determination of the surface were corrected.
 *
 * Revision 1.5  2002/01/16 15:03:28  j_novak
 * *** empty log message ***
 *
 * Revision 1.4  2002/01/11 14:09:34  j_novak
 * Added newtonian version for 2-fluid stars
 *
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

// Headers C
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Headers Lorene
#include "eos_bifluid.h"
#include "param.h"
#include "eos.h"
#include "cmp.h"
#include "utilitaires.h"

//****************************************************************************
// local prototypes
double puis(double, double) ;
double enthal1(const double x, const Param& parent) ;
double enthal23(const double x, const Param& parent) ;
double enthal(const double x, const Param& parent) ;      
//****************************************************************************

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor with gam1 = gam2 = 2, 
// gam3 = gam4 = gam5 = gam6 = 1, m_1 = 1 and m_2 =1
// -------------------------------------------------
Eos_bf_poly::Eos_bf_poly(double kappa1, double kappa2, double kappa3, double bet) :
  Eos_bifluid("bi-fluid polytropic EOS", 1, 1), 
  gam1(2), gam2(2),gam3(1),gam4(1),gam5(1),
  gam6(1),kap1(kappa1), kap2(kappa2), kap3(kappa3),beta(bet),
  relax(0.5), precis(1.e-9), ecart(1.e-8) 
{

  set_auxiliary() ; 
  determine_type() ;

}  

// Standard constructor with everything specified
// -----------------------------------------------
Eos_bf_poly::Eos_bf_poly(double gamma1, double gamma2, double gamma3,
			 double gamma4, double gamma5, double gamma6,
			 double kappa1, double kappa2, double kappa3,
			 double bet, double mass1, double mass2, 
			 double rel, double prec, double ec) : 
  Eos_bifluid("bi-fluid polytropic EOS", mass1, mass2), 
  gam1(gamma1),gam2(gamma2),gam3(gamma3),gam4(gamma4),gam5(gamma5),
  gam6(gamma6),kap1(kappa1),kap2(kappa2),kap3(kappa3),beta(bet), 
  relax(rel), precis(prec), ecart(ec) 
{

  set_auxiliary() ; 
  determine_type() ;
}  
  
// Copy constructor
// ----------------
Eos_bf_poly::Eos_bf_poly(const Eos_bf_poly& eosi) : 
	Eos_bifluid(eosi), 
	gam1(eosi.gam1), gam2(eosi.gam2), gam3(eosi.gam3),
	gam4(eosi.gam4), gam5(eosi.gam5), gam6(eosi.gam6),
	kap1(eosi.kap1), kap2(eosi.kap2), kap3(eosi.kap3),
	beta(eosi.beta),
        typeos(eosi.typeos), relax(eosi.relax), precis(eosi.precis),
        ecart(eosi.ecart) {
  
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
    fread_be(&relax, sizeof(double), 1, fich) ;	
    fread_be(&precis, sizeof(double), 1, fich) ;	
    fread_be(&ecart, sizeof(double), 1, fich) ;	
    
    set_auxiliary() ; 
    determine_type() ;

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
    fich >> relax ; fich.getline(blabla, 80) ;
    fich >> precis ; fich.getline(blabla, 80) ;
    fich >> ecart ; fich.getline(blabla, 80) ;

    set_auxiliary() ; 
    determine_type() ;

}
			//--------------//
			//  Destructor  //
			//--------------//

Eos_bf_poly::~Eos_bf_poly(){
  
  //Does nothing ...
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Eos_bf_poly::operator=(const Eos_bf_poly& eosi) {
    
  // Assignment of quantities common to all the derived classes of Eos_bifluid
  Eos_bifluid::operator=(eosi) ;	    
    
    gam1 = eosi.gam1 ; 
    gam2 = eosi.gam2 ; 
    gam3 = eosi.gam3 ; 
    kap1 = eosi.kap1 ; 
    kap2 = eosi.kap2 ; 
    kap3 = eosi.kap3 ;
    beta = eosi.beta ;
    typeos = eosi.typeos ;
    relax = eosi.relax ;
    precis = eosi.precis ;
    ecart = eosi.ecart ;
    
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

void Eos_bf_poly::determine_type() {
    
  if ((gam1 == double(2)) && (gam2 == double(2)) && (gam3 == double(1))
      && (gam4 == double(1)) && (gam5 == double(1)) 
      && (gam6 == double(0))) {
    typeos = -1 ;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl ;
    cout << "WARNING!! The entrainment factor does not depend on" <<endl ;
    cout << " density 2! This may be incorrect and should only be used"<<endl ;
    cout << " for testing purposes..." << endl ;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl ;
    return ;
  }
  if ((gam1 == double(2)) && (gam2 == double(2)) && (gam3 == double(1))
      && (gam4 == double(1)) && (gam5 == double(0)) 
      && (gam6 == double(1))) {
    typeos = -2 ;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl ;
    cout << "WARNING!! The entrainment factor does not depend on" << endl ;
    cout << " density 1! This may be incorrect and should only be used"<<endl ;
    cout << " for testing purposes..." << endl ;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl ;
    return ;
  }
  if ((gam1 == double(2)) && (gam2 == double(2)) && (gam3 == double(1))
      && (gam4 == double(1)) && (gam5 == double(1)) 
      && (gam6 == double(1))) {
    typeos = 0 ;
    return ;
  }
  if ((gam3 == double(1)) && (gam4 == double(1)) && (gam5 == double(1)) 
      && (gam6 == double(1))) {
    typeos = 1 ;
    return ;
  }
  if ((gam3 == double(1)) && (gam5 == double(1))) {
    typeos = 2 ;
    return ;
  }
  if ((gam4 == double(1)) && (gam6 == double(1))) {
    typeos = 3 ;
    return ;
  }
  typeos = 4 ;
  return ;  
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
    fwrite_be(&relax, sizeof(double), 1, fich) ;
    fwrite_be(&precis, sizeof(double), 1, fich) ;
    fwrite_be(&ecart, sizeof(double), 1, fich) ;
   
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
    ost << "   Type for inversion : " << typeos << endl ;
    ost << "   Parameters for inversion (used if typeos = 4) : " << endl ;
    ost << "   relaxation : " << relax << endl ;
    ost << "   precision for zerosec_b : " << precis << endl ;
    ost << "   final discrepancy in the result : " << ecart << endl ;
    
    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon densities from enthalpies 
//---------------------------------

bool Eos_bf_poly::nbar_ent_p(const double ent1, const double ent2, 
			      const double delta2, double& nbar1, 
			      double& nbar2) const {  

  bool one_fluid = ((ent1<=0.)||(ent2<=0.)) ;
  if (!one_fluid) {
    switch (typeos) {
    case 0: {
      double kpd = kap3+beta*delta2 ;
      double determ = kap1*kap2 - kpd*kpd ;
    
      nbar1 = (kap2*(exp(ent1) - m_1) - kpd*(exp(ent2) - m_2)) / determ ;
      nbar2 = (kap1*(exp(ent2) - m_2) - kpd*(exp(ent1) - m_1)) / determ ;
      one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)) ;
      break ;
    }
    case 1: {
      double b1 = exp(ent1) - m_1 ;
      double b2 = exp(ent2) - m_2 ;
      double denom = kap3 + beta*delta2 ;
      if (fabs(denom) < 1.e-15) {
	nbar1 = puis(2*b1/(gam1*kap1), 1./gam1m1) ;
	nbar2 = puis(2*b2/(gam2*kap2), 1./gam2m1) ;
	one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)) ;
      }
      else {
	double coef0 = 0.5*gam2*kap2/pow(denom, gam2m1) ;
	double coef1 = 0.5*kap1*gam1 ;
	Param parent ;
	parent.add_double(coef0, 0) ;
	parent.add_double(b1, 1) ;
	parent.add_double(coef1, 2) ;
	parent.add_double(gam1m1,3) ;
	parent.add_double(gam2m1,4) ;
	parent.add_double(denom, 5) ;
	parent.add_double(b2, 6) ;

	double xmin, xmax ;
	double f0 = enthal1(0.,parent) ;
	if (fabs(f0)<1.e-15) {
	  nbar1 = 0. ;}
	else {
	  double cas = (gam1m1*gam2m1 - 1.)*f0;
	  assert (fabs(cas) > 1.e-15) ;
	  xmin = 0. ;
	  xmax = cas/fabs(cas) ;
	  do {
	    xmax *= 1.1 ;
	    if (fabs(xmax) > 1.e10) {
	      cout << "Problem in Eos_bf_poly::nbar_ent_p!" << endl ;
	      cout << f0 << ", " << cas << endl ; // to be removed!!
	      cout << "typeos = 1" << endl ;
	      abort() ;
	    }
	  } while (enthal1(xmax,parent)*f0 > 0.) ;
	  double l_precis = 1.e-12 ;
	  int nitermax = 400 ;
	  int niter = 0 ;
	  nbar1 = zerosec_b(enthal1, parent, xmin, xmax, l_precis, 
			    nitermax, niter) ;
	}
	nbar2 = (b1 - coef1*puis(nbar1, gam1m1))/denom ;
	double resu1 = coef1*puis(nbar1,gam1m1) + denom*nbar2 ;
	double resu2 = 0.5*gam2*kap2*puis(nbar2,gam2m1) + denom*nbar1 ;
	double erreur = fabs((log(fabs(1+resu1))-ent1)/ent1) + 
	  fabs((log(fabs(1+resu2))-ent2)/ent2) ;
	one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)||(erreur > 1.e-4)) ;
      }
      break ;
    }
    case 2: {
      double b1 = exp(ent1) - m_1 ;
      double b2 = exp(ent2) - m_2 ;
      double denom = kap3 + beta*delta2 ;
      if (fabs(denom) < 1.e-15) {
	nbar1 = puis(2*b1/(gam1*kap1), 1./gam1m1) ;
	nbar2 = puis(2*b2/(gam2*kap2), 1./gam2m1) ;
	one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)) ;
      }
      else {
	double coef0 = beta*delta2 ;
	double coef1 = 0.5*kap1*gam1 ;
	double coef2 = 0.5*kap2*gam2 ;
	double coef3 = 1./gam1m1 ;
	Param parent ;
	parent.add_double(b1, 0) ;
	parent.add_double(kap3, 1) ;
	parent.add_double(gam4, 2) ;
	parent.add_double(coef0, 3) ;
	parent.add_double(gam6,4) ;
	parent.add_double(coef1, 5) ;
	parent.add_double(coef3, 6) ;
	parent.add_double(coef2, 7) ;
	parent.add_double(gam2m1, 8) ;
	parent.add_double(b2, 9) ;

	double xmin, xmax ;
	double f0 = enthal23(0.,parent) ;
	if (fabs(f0)<1.e-15) {
	  nbar2 = 0. ;}
	else {
	  double pmax = (fabs(kap3) < 1.e-15 ? 0. : gam4*(gam4-1) ) ;
	  double ptemp = (fabs(kap3*coef0) < 1.e-15  ? 0. 
			  : gam6*(gam4-1) ) ;
	  pmax = (pmax>ptemp ? pmax : ptemp) ;
	  ptemp = (fabs(kap3*coef0) < 1.e-15 ? 0. : gam4*(gam6-1) ) ;
	  pmax = (pmax>ptemp ? pmax : ptemp) ;
	  ptemp = (fabs(coef0) < 1.e-15 ? 0. : gam6*(gam6-1) ) ;
	  pmax = (pmax>ptemp ? pmax : ptemp) ;
	  double cas = (pmax - gam1m1*gam2m1)*f0;
	  //	  cout << "Pmax, cas: " << pmax << ", " << cas << endl ;
	  assert (fabs(cas) > 1.e-15) ;
	  xmin = 0. ;
	  xmax = cas/fabs(cas) ;
	  do {
	    xmax *= 1.1 ;
	    if (fabs(xmax) > 1.e10) {
	      cout << "Problem in Eos_bf_poly::nbar_ent_p!" << endl ;
	      cout << "typeos = 2" << endl ;
	      abort() ;
	    }
	  } while (enthal23(xmax,parent)*f0 > 0.) ;
	  double l_precis = 1.e-12 ;
	  int nitermax = 400 ;
	  int niter = 0 ;
	  nbar2 = zerosec_b(enthal23, parent, xmin, xmax, l_precis, 
			    nitermax, niter) ;
	}
	nbar1 = (b1 - kap3*puis(nbar2,gam4) - coef0*puis(nbar2,gam6))
	  / coef1 ;
	nbar1 = puis(nbar1,coef3) ;
	double resu1 = coef1*puis(nbar1,gam1m1) + kap3*puis(nbar2,gam4)
	  + coef0*puis(nbar2, gam6) ;
	double resu2 = coef2*puis(nbar2,gam2m1) 
	  + gam4*kap3*puis(nbar2, gam4-1)*nbar1
	  + gam6*coef0*puis(nbar2, gam6-1)*nbar1 ;
	double erreur = fabs((log(fabs(1+resu1))-ent1)/ent1) + 
	  fabs((log(fabs(1+resu2))-ent2)/ent2) ;
	//cout << "Erreur d'inversion: " << erreur << endl ;
	one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)||(erreur > 1.e-4)) ;
      }
      break ;
    }
    case 3: {
      double b1 = exp(ent1) - m_1 ;
      double b2 = exp(ent2) - m_2 ;
      double denom = kap3 + beta*delta2 ;
      if (fabs(denom) < 1.e-15) {
	nbar1 = puis(2*b1/(gam1*kap1), 1./gam1m1) ;
	nbar2 = puis(2*b2/(gam2*kap2), 1./gam2m1) ;
	one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)) ;
      }
      else {
	double coef0 = beta*delta2 ;
	double coef1 = 0.5*kap1*gam1 ;
	double coef2 = 0.5*kap2*gam2 ;
	double coef3 = 1./gam2m1 ;
	Param parent ;
	parent.add_double(b2, 0) ;
	parent.add_double(kap3, 1) ;
	parent.add_double(gam3, 2) ;
	parent.add_double(coef0, 3) ;
	parent.add_double(gam5, 4) ;
	parent.add_double(coef2, 5) ;
	parent.add_double(coef3, 6) ;
	parent.add_double(coef1, 7) ;
	parent.add_double(gam1m1, 8) ;
	parent.add_double(b1, 9) ;
	
	double xmin, xmax ;
	double f0 = enthal23(0.,parent) ;
	if (fabs(f0)<1.e-15) {
	  nbar1 = 0. ;}
	else {
	  double pmax = (fabs(kap3) < 1.e-15 ? 0. : gam3*(gam3-1) ) ;
	  double ptemp = (fabs(kap3*coef0) < 1.e-15  ? 0. 
			  : gam5*(gam3-1) ) ;
	  pmax = (pmax>ptemp ? pmax : ptemp) ;
	  ptemp = (fabs(kap3*coef0) < 1.e-15 ? 0. : gam3*(gam5-1) ) ;
	  pmax = (pmax>ptemp ? pmax : ptemp) ;
	  ptemp = (fabs(coef0) < 1.e-15 ? 0. : gam5*(gam5-1) ) ;
	  pmax = (pmax>ptemp ? pmax : ptemp) ;
	  double cas = (pmax - gam1m1*gam2m1)*f0;
	  //	  cout << "Pmax, cas: " << pmax << ", " << cas << endl ;
	  assert (fabs(cas) > 1.e-15) ;
	  xmin = 0. ;
	  xmax = cas/fabs(cas) ;
	  do {
	    xmax *= 1.1 ;
	    if (fabs(xmax) > 1.e10) {
	      cout << "Problem in Eos_bf_poly::nbar_ent_p!" << endl ;
	      cout << "typeos = 3" << endl ;
	      abort() ;
	    }
	  } while (enthal23(xmax,parent)*f0 > 0.) ;
	  double l_precis = 1.e-12 ;
	  int nitermax = 400 ;
	  int niter = 0 ;
	  nbar1 = zerosec_b(enthal23, parent, xmin, xmax, l_precis, 
			    nitermax, niter) ;
	}
	nbar2 = (b2 - kap3*puis(nbar1,gam3) - coef0*puis(nbar1,gam5))
	  / coef2 ;
	nbar2 = puis(nbar2,coef3) ;
	double resu1 = coef1*puis(nbar1,gam1m1) 
	  + gam3*kap3*puis(nbar1,gam3-1)*nbar2
	  + coef0*gam5*puis(nbar1, gam5-1)*nbar2 ;
	double resu2 = coef2*puis(nbar2,gam2m1) 
	  + kap3*puis(nbar1, gam3) + coef0*puis(nbar1, gam5);
	double erreur = fabs((log(fabs(1+resu1))-ent1)/ent1) + 
	  fabs((log(fabs(1+resu2))-ent2)/ent2) ;
	one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)||(erreur > 1.e-4)) ;
      }
      break ;
    }    
    case 4:{
      double b1 = exp(ent1) - m_1 ; 
      double b2 = exp(ent2) - m_2 ; 
      double denom = kap3 + beta*delta2 ;
      if (fabs(denom) < 1.e-15) {
	nbar1 = puis(2*b1/(gam1*kap1), 1./gam1m1) ;
	nbar2 = puis(2*b2/(gam2*kap2), 1./gam2m1) ;
	one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)) ;
      }
      else {
	int nitermax = 200 ;
	int niter = 0 ;
	int nmermax = 800 ;
	
	double a11 = 0.5*gam1*kap1 ;
	double a13 = gam3*kap3 ;
	double a14 = beta*gam5*delta2 ;
	
	double a22 = 0.5*gam2*kap2 ;
	double a23 = gam4*kap3 ;
	double a24 = beta*gam6*delta2 ;
	
	double n1l, n2l, n1s, n2s ;
	
	double delta = a11*a22 - (a13+a14)*(a23+a24) ;
	n1l = (a22*b1 - (a13+a14)*b2)/delta ;
	n2l = (a11*b2 - (a23+a24)*b1)/delta ;
	n1s = puis(b1/a11, 1./(gam1m1)) ;
	n2s = puis(b2/a22, 1./(gam2m1)) ;
	
	double n1m = (n1l<0. ? n1s : sqrt(n1l*n1s)) ;
	double n2m = (n2l<0. ? n2s : sqrt(n2l*n2s)) ;
	
	Param parf1 ;
	Param parf2 ;
	double c1, c2, c3, c4, c5, c6, c7 ;
	c1 = gam1m1 ;
	c2 = gam3 - 1. ;
	c3 = gam5 - 1. ;
	c4 = a11 ;
	c5 = a13*puis(n2m,gam4) ;
	c6 = a14*puis(n2m,gam6) ;
	c7 = b1 ; 
	parf1.add_double_mod(c1,0) ;
	parf1.add_double_mod(c2,1) ;
	parf1.add_double_mod(c3,2) ;
	parf1.add_double_mod(c4,3) ;
	parf1.add_double_mod(c5,4) ;
	parf1.add_double_mod(c6,5) ;
	parf1.add_double_mod(c7,6) ;
	
	double d1, d2, d3, d4, d5, d6, d7 ;
	d1 = gam2m1 ;
	d2 = gam4 - 1. ;
	d3 = gam6 - 1. ;
	d4 = a22 ;
	d5 = a23*puis(n1m,gam3) ;
	d6 = a24*puis(n1m,gam5) ;
	d7 = b2 ; 
	parf2.add_double_mod(d1,0) ;
	parf2.add_double_mod(d2,1) ;
	parf2.add_double_mod(d3,2) ;
	parf2.add_double_mod(d4,3) ;
	parf2.add_double_mod(d5,4) ;
	parf2.add_double_mod(d6,5) ;
	parf2.add_double_mod(d7,6) ;
	
	double xmin = -3*(n1s>n2s ? n1s : n2s) ;
	double xmax = 3*(n1s>n2s ? n1s : n2s) ;
	
	double n1 = n1m ;
	double n2 = n2m ;
	bool sortie = true ;
	int mer = 0 ;
	
	//cout << "Initial guess: " << n1 << ", " << n2 << endl ;
	n1s *= 0.1 ;
	n2s *= 0.1 ;
	do {
	  //cout << "n1, n2: " << n1 << ", " << n2 << endl ;
	  n1 = zerosec_b(&enthal, parf1, xmin, xmax, precis, nitermax, niter) ;
	  n2 = zerosec_b(&enthal, parf2, xmin, xmax, precis, nitermax, niter) ;
	  
	  sortie = (fabs((n1m-n1)/(n1m+n1)) + fabs((n2m-n2)/(n2m+n2)) > ecart) ;
	  n1m = relax*n1 + (1.-relax)*n1m ;
	  n2m = relax*n2 + (1.-relax)*n2m ;
	  if (n2m>-n2s) {
	    parf1.get_double_mod(4) = a13*puis(n2m,gam4) ;
	    parf1.get_double_mod(5) = a14*puis(n2m,gam6) ;
	  }
	  else {
	    parf1.get_double_mod(4) = a13*puis(-n2s,gam4) ;
	    parf1.get_double_mod(5) = a14*puis(-n2s,gam6) ;
	  }
	  
	  if (n1m>-n1s) {
	    parf2.get_double_mod(4) = a23*puis(n1m,gam3) ;
	    parf2.get_double_mod(5) = a24*puis(n1m,gam5) ;
	  }
	  else {
	    parf2.get_double_mod(4) = a23*puis(-n1s,gam3) ;
	    parf2.get_double_mod(5) = a24*puis(-n1s,gam5) ;
	  }
	  
	  mer++ ;
	} while (sortie&&(mer<nmermax)) ;
	nbar1 = n1m ;
	nbar2 = n2m ;
	
//  	double resu1 = a11*puis(n1,gam1m1) + a13*puis(n1,gam3-1.)*puis(n2,gam4)
//  	  +a14*puis(n1,gam5-1.)*puis(n2,gam6) ;
//  	double resu2 = a22*puis(n2,gam2m1) + a23*puis(n1,gam3)*puis(n2,gam4-1.)
//  	  +a24*puis(n1,gam5)*puis(n2,gam6-1.) ;
	//cout << "Nbre d'iterations: " << mer << endl ;
	//cout << "Resus: " << n1m << ", " << n2m << endl ;
	//cout << "Verification: " << log(fabs(1+resu1)) << ", " 
	//	   << log(fabs(1+resu2)) << endl ; 
	//    cout << "Erreur: " << fabs(enthal(n1, parf1)/b1) + 
	//      fabs(enthal(n2, parf2)/b2) << endl ;
	//cout << "Erreur: " << fabs((log(fabs(1+resu1))-ent1)/ent1) + 
	//fabs((log(fabs(1+resu2))-ent2)/ent2) << endl ;
      }
      break ;
    }
    case -1: {
      double determ = kap1*kap2 - kap3*kap3 ;
    
      nbar1 = (kap2*(exp(ent1) - m_1 - beta*delta2) 
	       - kap3*(exp(ent2) - m_2)) / determ ;
      nbar2 = (kap1*(exp(ent2) - m_2) - kap3*(exp(ent1) - m_1 - beta*delta2))
	/ determ ;
      one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)) ;
      break ;
    }
    case -2: {
      double determ = kap1*kap2 - kap3*kap3 ;
    
      nbar1 = (kap2*(exp(ent1) - m_1 ) 
	       - kap3*(exp(ent2) - m_2 - beta*delta2)) / determ ;
      nbar2 = (kap1*(exp(ent2) - m_2 - beta*delta2) 
	       - kap3*(exp(ent1) - m_1)) / determ ;
      one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)) ;
      break ;
    }
    }
  }
  return one_fluid ;
}
// One fluid sub-EOSs
//-------------------

double Eos_bf_poly::nbar_ent_p1(const double ent1) const {
  return puis(2*(exp(ent1) - m_1)/(gam1*kap1), 1./gam1m1) ;
}

double Eos_bf_poly::nbar_ent_p2(const double ent2) const {
  return puis(2*(exp(ent2) - m_2)/(gam2*kap2), 1./gam2m1) ;
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
      (delta2*(gam5 + 2) - 2)*beta * pow(n1,gam5 - 2)*pow(n2, gam6) ;
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
      gam4*kap3 * pow(n2,gam4 - 2) * pow(n1,gam3) + 
      (delta2*(gam6 + 2) - 2)*beta * pow(n1, gam5) * pow(n2,gam6 - 2) ;
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
    xx = 2*beta*pow(n1,gam5-1)*pow(n2,gam6-1) / gamma_delta3 ;
  }
  return xx ;
}

// Conversion functions
// ---------------------

Eos* Eos_bf_poly::trans2Eos() const {

  Eos_poly* eos_simple = new Eos_poly(gam1, kap1, m_1) ;
  return eos_simple ;
}
/***************************************************************************
 *
 *                         Non-class functions
 *
 ***************************************************************************/

// New "pow"
//-----------

 double puis(double x, double p) {
  assert(p>=0.) ;
  if (p==0.) return (x>=0 ? 1 : -1) ;
  //if (p==0.) return 1 ;
  if (x<0.) return (-pow(-x,p)) ;
  else return pow(x,p) ;
}

// Auxilliary functions for nbar_ent_p
//------------------------------------

double enthal1(const double x, const Param& parent) {
  assert(parent.get_n_double() == 7) ;

  return parent.get_double(0)*puis(parent.get_double(1) - parent.get_double(2)
	   *puis(x,parent.get_double(3)), parent.get_double(4)) 
    + parent.get_double(5)*x - parent.get_double(6) ;

}

double enthal23(const double x, const Param& parent) {
  assert(parent.get_n_double() == 10) ;

  double nx = (parent.get_double(0) - parent.get_double(1)*
	       puis(x,parent.get_double(2)) - parent.get_double(3)*
	       puis(x,parent.get_double(4)) )/parent.get_double(5) ;
  nx = puis(nx,parent.get_double(6)) ;
  return parent.get_double(7)*puis(x,parent.get_double(8)) 
    + parent.get_double(1)*parent.get_double(2)*nx*
    puis(x,parent.get_double(2) - 1) 
    + parent.get_double(3)*parent.get_double(4)*nx*
    puis(x,parent.get_double(4) - 1) 
    - parent.get_double(9) ;

}

double enthal(const double x, const Param& parent) {
  assert(parent.get_n_double_mod() == 7) ;

  double alp1 = parent.get_double_mod(0) ;
  double alp2 = parent.get_double_mod(1) ;
  double alp3 = parent.get_double_mod(2) ;
  double cc1 = parent.get_double_mod(3) ;
  double cc2 = parent.get_double_mod(4) ;
  double cc3 = parent.get_double_mod(5) ;
  double cc4 = parent.get_double_mod(6) ;

  return (cc1*puis(x,alp1) + cc2*puis(x,alp2) + cc3*puis(x,alp3) - cc4) ;

}

