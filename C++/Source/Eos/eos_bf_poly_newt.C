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
 * Revision 1.4  2002/05/02 15:16:22  j_novak
 * Added functions for more general bi-fluid EOS
 *
 * Revision 1.3  2002/04/05 09:09:36  j_novak
 * The inversion of the EOS for 2-fluids polytrope has been modified.
 * Some errors in the determination of the surface were corrected.
 *
 * Revision 1.2  2002/01/16 15:03:28  j_novak
 * *** empty log message ***
 *
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
			 double bet, double mass1, double mass2, double
			 relax, double precis, double ecart) : 
  Eos_bf_poly(gamma1, gamma2, gamma3, gamma4, gamma5, gamma6,
	      kappa1, kappa2, kappa3, bet, mass1, mass2, relax, precis, 
	      ecart) {
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
    relax = eosi.relax ;
    precis = eosi.precis ;
    ecart = eosi.ecart ;
    
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
    ost << "   Parameters for inversion (used if typeos = 4 : " << endl ;
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

void  Eos_bf_poly_newt::nbar_ent_p(const double ent1, const double ent2, 
			      const double delta2, double& nbar1, 
			      double& nbar2, bool tronc) const {  

  bool one_fluid = ((ent1<=0.)||(ent2<=0.)) ;
  if (!one_fluid) {
    switch (typeos) {
    case 0: {
      double kpd = kap3+beta*delta2 ;
      double determ = kap1*kap2 - kpd*kpd ;
    
      nbar1 = (kap2*ent1*m_1 - kpd*ent2*m_2) / determ ;
      nbar2 = (kap1*ent2*m_2 - kpd*ent1*m_1) / determ ;
      one_fluid = ((nbar1 < 0.)||(nbar2 < 0.)) ;
      break ;
    }
    case 1: {
      double b1 = ent1*m_1 ;
      double b2 = ent2*m_2 ;
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
	  double precis = 1.e-12 ;
	  int nitermax = 400 ;
	  int niter = 0 ;
	  nbar1 = zerosec_b(enthal1, parent, xmin, xmax, precis, 
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
      double b1 = ent1*m_1 ;
      double b2 = ent2*m_2 ;
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
	  double precis = 1.e-12 ;
	  int nitermax = 400 ;
	  int niter = 0 ;
	  nbar2 = zerosec_b(enthal23, parent, xmin, xmax, precis, 
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
      double b1 = ent1*m_1 ;
      double b2 = ent2*m_2 ;
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
	  double precis = 1.e-12 ;
	  int nitermax = 400 ;
	  int niter = 0 ;
	  nbar1 = zerosec_b(enthal23, parent, xmin, xmax, precis, 
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
      double b1 = ent1*m_1 ; 
      double b2 = ent2*m_2 ; 
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
    }
  }
  if (one_fluid) {
    nbar1 = puis(2*ent1*m_1/(gam1*kap1), 1./gam1m1) ;
    nbar2 = puis(2*ent2*m_2/(gam2*kap2), 1./gam2m1) ;
    if (tronc) { 
      nbar1 = nbar1 < 0. ? 0. : nbar1 ;
      nbar2 = nbar2 < 0. ? 0. : nbar2 ;
    }
    if ((nbar1 > 0.)&&(nbar2 > 0.)) {
      if (nbar1>nbar2) nbar2 = 0. ;
      if (nbar2>nbar1) nbar1 = 0. ;
    }
    return ;
  }
  return ;
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
       
