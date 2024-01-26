/*
 * Methods of the class Pseudo_polytrope_1D.
 *
 * (see file eos.h for documentation).
 */

/*
 *   Copyright (c) 2023 Gaël Servignat
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
 
 */

// Headers C
#include <cstdlib>
#include <cstring>
#include <cmath>

// Headers Lorene
#include "eos.h"
#include "cmp.h"
#include "unites.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor 
// --------------------
namespace Lorene {
Pseudo_polytrope_1D::Pseudo_polytrope_1D(const Tbl& coefs0, double n_lim0, double m_n0): 
	      n_lim1(n_lim0), n_lim2(n_lim0), m_n(m_n0) {
	
	cerr << "Deprecated constructor, please do not use ! Aborting..." << endl; abort() ;
    set_name("Pseudo-polytropic fit of cold, beta-equilibrated EoS") ; 
    
	assert(coefs0.get_ndim() == 1) ;
    coefs   = new Tbl(coefs0) ;
	n_coefs = coefs->get_taille() ;
	
	double x_lim = log(n_lim1*0.1) ;
	double alpha = (*coefs)(n_coefs-1) ;
	double sum_poly = 0., sum_der_poly = 0. ;
	for (int i=0; i<n_coefs-1; i++){
		sum_poly     += (*coefs)(i)*pow(x_lim, double(i)) ;
		sum_der_poly += (i == 0) ? 0. : (*coefs)(i)*double(i)*pow(x_lim, double(i-1)) ;
	} 
  ent_lim1 = log1p(exp(alpha*x_lim)*((alpha+1)*sum_poly + sum_der_poly)) ;

}  

  
// Copy constructor
// ----------------
Pseudo_polytrope_1D::Pseudo_polytrope_1D(const Pseudo_polytrope_1D& eosi): Eos(eosi), coefs(eosi.coefs), gamma_low(eosi.gamma_low), kappa_low(eosi.kappa_low),
																			 n_lim1(eosi.n_lim1), n_lim2(eosi.n_lim2), ent_lim1(eosi.ent_lim1), ent_lim2(eosi.ent_lim2), 
                                       m_n(eosi.m_n), Lambda(eosi.Lambda), ddd(eosi.ddd), Kappa_GPP(eosi.Kappa_GPP),
                                       Gamma_GPP(eosi.Gamma_GPP), Lambda_GPP(eosi.Lambda_GPP), ddd_GPP(eosi.ddd_GPP), n_coefs(eosi.n_coefs) {}
  

// Constructor from a binary file
// ------------------------------
Pseudo_polytrope_1D::Pseudo_polytrope_1D(FILE* fich): Eos(fich) {
    
    fread_be(&n_coefs, sizeof(int), 1, fich) ;
	coefs = new Tbl(fich) ;
    fread_be(&kappa_low, sizeof(double), 1, fich) ;
    fread_be(&gamma_low, sizeof(double), 1, fich) ;
	fread_be(&n_lim1, sizeof(double), 1, fich) ;
	fread_be(&ent_lim1, sizeof(double), 1, fich) ;
	fread_be(&n_lim2, sizeof(double), 1, fich) ;
	fread_be(&ent_lim2, sizeof(double), 1, fich) ;
	fread_be(&m_n, sizeof(double), 1, fich) ;
	fread_be(&Lambda, sizeof(double), 1, fich) ;
	fread_be(&ddd, sizeof(double), 1, fich) ;
	fread_be(&Kappa_GPP, sizeof(double), 1, fich) ;
	fread_be(&Gamma_GPP, sizeof(double), 1, fich) ;
	fread_be(&Lambda_GPP, sizeof(double), 1, fich) ;
	fread_be(&ddd_GPP, sizeof(double), 1, fich) ;
	
	eos_low = new Eos_poly(gamma_low, kappa_low) ;
}
	       
// Constructor from a formatted file
// ---------------------------------
Pseudo_polytrope_1D::Pseudo_polytrope_1D(ifstream& fich): Eos(fich) {

	fich >> n_coefs ; fich.ignore(1000, '\n') ;
	coefs = new Tbl(n_coefs) ; coefs->set_etat_qcq() ;

	for (int i=0; i < n_coefs; i++)
		fich >> coefs->set(i) ;
	fich.ignore(1000, '\n') ;
	
	fich >> n_lim1 ; fich.ignore(1000, '\n') ;
	fich >> n_lim2 ; fich.ignore(1000, '\n') ;
	
	fich >> kappa_low ; fich.ignore(1000, '\n') ;
	
	fich >> gamma_low ; fich.ignore(1000, '\n') ;
	
	fich >> m_n ; fich.ignore(1000, '\n') ;
	
	fich >> Lambda ; fich.ignore(1000, '\n') ;

	fich >> ddd ; fich.ignore(1000, '\n') ;
  
  fich >> Kappa_GPP ; fich.ignore(1000, '\n');
  
  fich >> Gamma_GPP ; fich.ignore(1000, '\n');

	fich >> Lambda_GPP ; fich.ignore(1000, '\n') ;

	fich >> ddd_GPP ; fich.ignore(1000, '\n') ;
	
	double x_lim2 = log(n_lim2*0.1) ;
	double alpha = (*coefs)(n_coefs-1) ;
	double sum_poly = 0., sum_der_poly = 0. ;
	for (int i=0; i<n_coefs-1; i++){
		sum_poly     += (*coefs)(i)*pow(x_lim2, double(i)) ;
		sum_der_poly += (i == 0) ? 0. : (*coefs)(i)*double(i)*pow(x_lim2, double(i-1)) ;
	}
  cout << setprecision(16) << "nlim1: " << n_lim1*0.1 << " n_lim2: " << n_lim2*0.1 << endl << " ddd: " << ddd << " et: " ;
  ent_lim2 = log1p(ddd + exp(alpha*x_lim2)*((alpha+1)*sum_poly + sum_der_poly)) ;
  ent_lim1 = log1p(ddd_GPP + Gamma_GPP*Kappa_GPP/(Gamma_GPP-1.) * pow(n_lim1, Gamma_GPP-1.)) ;
	cout << ent_lim1 << endl;
	cout << ent_lim2 << endl;
	
	eos_low = new Eos_poly(gamma_low, kappa_low) ;
	
}
	       

			//--------------//
			//  Destructor  //
			//--------------//

Pseudo_polytrope_1D::~Pseudo_polytrope_1D(){
    
    if (coefs != 0x0)   delete coefs ;
	if (eos_low != 0x0) delete eos_low ;
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Pseudo_polytrope_1D::operator=(const Pseudo_polytrope_1D& eosi) {
    
    set_name(eosi.name) ; 
    
    coefs     = eosi.coefs ; 
    n_lim1     = eosi.n_lim1 ; 
    n_lim2     = eosi.n_lim2 ; 
	n_coefs   = eosi.n_coefs ;
	gamma_low = eosi.gamma_low ;
	kappa_low = eosi.kappa_low ;
	m_n       = eosi.m_n ;
	ent_lim1   = eosi.ent_lim1 ;
	ent_lim2   = eosi.ent_lim2 ;
	Lambda    = eosi.Lambda ;
	ddd    = eosi.ddd ;
	Kappa_GPP    = eosi.Kappa_GPP ;
	Gamma_GPP    = eosi.Gamma_GPP ;
	Lambda_GPP    = eosi.Lambda_GPP ;
	ddd_GPP    = eosi.ddd_GPP ;
	eos_low   = eosi.eos_low ;
    
}
			//------------------------//
			//  Comparison operators  //
			//------------------------//

bool Pseudo_polytrope_1D::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Pseudo_polytrope_1D !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Pseudo_polytrope_1D& eos = dynamic_cast<const Pseudo_polytrope_1D&>( eos_i ) ; 
	
	for (int i=0; i < n_coefs; i++)
		if ((*eos.coefs)(i) != (*coefs)(i)) {
			cout 
			<< "The two Pseudo_polytrope_1D have different coefficients: " << (*coefs)(i) << " <-> " 
			<< (*eos.coefs)(i) << endl ; 
			resu = false ; 
		}

	if (eos.n_lim1 != n_lim1) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different limiting densities n_lim1: " << n_lim1 << " <-> " 
		<< eos.n_lim1 << endl ; 
	    resu = false ; 
	}
	
	if (eos.ent_lim1 != ent_lim1) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different limiting enthalpies ent_lim1: " << ent_lim1 << " <-> " 
		<< eos.ent_lim1 << endl ; 
	    resu = false ; 
	}

	if (eos.n_lim2 != n_lim2) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different limiting densities n_lim2: " << n_lim2 << " <-> " 
		<< eos.n_lim2 << endl ; 
	    resu = false ; 
	}
	
	if (eos.ent_lim2 != ent_lim2) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different limiting enthalpies ent_lim2: " << ent_lim2 << " <-> " 
		<< eos.ent_lim2 << endl ; 
	    resu = false ; 
	}
	
	if (eos.n_coefs != n_coefs) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different number of coefficients: " << n_coefs << " <-> " 
		<< eos.n_coefs << endl ; 
	    resu = false ; 
	}
	
	if (eos.kappa_low != kappa_low) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different kappa in low polytrope: " << kappa_low << " <-> " 
		<< eos.kappa_low << endl ; 
	    resu = false ; 
	}
	
	if (eos.gamma_low != gamma_low) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different kappa in low polytrope: " << gamma_low << " <-> " 
		<< eos.gamma_low << endl ; 
	    resu = false ; 
	}
	
	if (eos.ddd!= ddd) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different ddd: " << ddd << " <-> " 
		<< eos.ddd << endl ; 
	    resu = false ; 
	}

	if (eos.Lambda!= Lambda) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different Lambda: " << Lambda << " <-> " 
		<< eos.Lambda << endl ; 
	    resu = false ; 
	}

	if (eos.Kappa_GPP!= Kappa_GPP) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different Kappa_GPP: " << Kappa_GPP << " <-> " 
		<< eos.Kappa_GPP << endl ; 
	    resu = false ; 
	}

	if (eos.Gamma_GPP!= Gamma_GPP) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different Gamma_GPP: " << Gamma_GPP << " <-> " 
		<< eos.Gamma_GPP << endl ; 
	    resu = false ; 
	}

	if (eos.Lambda_GPP!= Lambda_GPP) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different Lambda_GPP: " << Lambda_GPP << " <-> " 
		<< eos.Lambda_GPP << endl ; 
	    resu = false ; 
	}
	
  	if (eos.ddd_GPP!= ddd_GPP) {
	    cout 
	    << "The two Pseudo_polytrope_1D have different ddd_GPP: " << ddd_GPP << " <-> " 
		<< eos.ddd_GPP << endl ; 
	    resu = false ; 
	}

    }
    
    return resu ; 
    
}

bool Pseudo_polytrope_1D::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}


			//------------//
			//  Outputs   //
			//------------//

void Pseudo_polytrope_1D::sauve(FILE* fich) const {

    Eos::sauve(fich) ; 
    
    fwrite_be(&n_coefs, sizeof(int), 1, fich) ;	
    coefs->sauve(fich) ;
    fwrite_be(&kappa_low, sizeof(double), 1, fich) ;
    fwrite_be(&gamma_low, sizeof(double), 1, fich) ;
	fwrite_be(&n_lim1, sizeof(double), 1, fich) ;
	fwrite_be(&ent_lim1, sizeof(double), 1, fich) ;
	fwrite_be(&n_lim2, sizeof(double), 1, fich) ;
	fwrite_be(&ent_lim2, sizeof(double), 1, fich) ;
	fwrite_be(&m_n, sizeof(double), 1, fich) ;
	fwrite_be(&Lambda, sizeof(double), 1, fich) ;
	fwrite_be(&ddd, sizeof(double), 1, fich) ;
	fwrite_be(&Kappa_GPP, sizeof(double), 1, fich) ;
	fwrite_be(&Gamma_GPP, sizeof(double), 1, fich) ;
	fwrite_be(&Lambda_GPP, sizeof(double), 1, fich) ;
	fwrite_be(&ddd_GPP, sizeof(double), 1, fich) ;
	
       
}

ostream& Pseudo_polytrope_1D::operator>>(ostream & ost) const {
    
    ost << setprecision(16) << "EOS of class Pseudo_polytrope_1D (analytical fit of cold EoS): " << '\n' ; 
    ost << "   Coefficients:            " << *coefs << '\n' ; 
	ost << setprecision(16) << "   Baryon mass:             " << m_n << " [MeV]" << '\n' ;
	ost << setprecision(16) << "   ddd: " << ddd << " and Lambda: " << Lambda << '\n' ;
	ost << setprecision(16) << "   Kappa_GPP: " << Kappa_GPP << " and Gamma_GPP: " << Gamma_GPP << '\n' ;
	ost << setprecision(16) << "   Lambda_GPP: " << ddd << " and Lambda_GPP: " << Lambda << '\n' ;
	ost << "Low densities polytrope :" << '\n' ;
	cout << *eos_low << '\n' ;
    ost << setprecision(16) << "   Limiting density n_lim1:  " << 0.1*n_lim1 << " [fm^-3]" << '\n' ; 
	ost << setprecision(16) << "   Limiting enthalpy h_lim1: " << ent_lim1 << " [c^2]" << endl;
    ost << setprecision(16) << "   Limiting density n_lim2:  " << 0.1*n_lim2 << " [fm^-3]" << '\n' ; 
	ost << setprecision(16) << "   Limiting enthalpy h_lim2: " << ent_lim2 << " [c^2]" << endl;
    
    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy 
//------------------------------

double Pseudo_polytrope_1D::nbar_ent_p(double ent, const Param* ) const {

    if (ent >= ent_lim2 ) {
		
		auto ent_nb_p = [&](double nbar) {
			double xx = log(nbar*0.1) ;
			double alpha = (*coefs)(n_coefs-1) ;
			double sum_poly = 0.;
			for (int i=0; i < n_coefs-1; i++){
				double cc1 = (alpha+1.)*(*coefs)(i) ;
				double cc2 = (i == n_coefs - 2) ? 0. : double(i+1)*(*coefs)(i+1) ;
				sum_poly  += (cc1 + cc2)*pow(xx, double(i));
			}
			double arg = exp(alpha*xx)*sum_poly + ddd ;
			return log1p(arg) - ent ;
		} ;
		
		double a = n_lim2/2., b = max(100.*ent, 50.*n_lim2) ;
		double f0 = ent_nb_p(a), f1 = ent_nb_p(b) ;
		double c=1., c_old=2., f2 ;
		double eps = 5e-16 ;
		assert( f0 * f1 < 0.) ;
		while(fabs((c-c_old)/c)>eps) {
			c_old = c ;
			c     = (b*f0 - a*f1)/(f0 - f1) ;
			f2    = ent_nb_p(c) ;
			
			if (f2*f0 < 0.){
				b  = c  ;
				f1 = f2 ;
			}
			else{
				a  = c  ;
				f0 = f2 ;
			}
		}
	    return c ;
    }
    else if ( (ent_lim1 <= ent) && (ent < ent_lim2)){

      return pow((expm1(ent)-ddd_GPP)*(Gamma_GPP-1.)/(Kappa_GPP*Gamma_GPP), 1./(Gamma_GPP-1.));
    }
    else if ( (0 < ent) && (ent < ent_lim1)){

		return eos_low->nbar_ent_p(ent) ;
    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Pseudo_polytrope_1D::ener_ent_p(double ent, const Param* ) const {
    if (ent >= ent_lim2 ) {
		
		double nbar = nbar_ent_p(ent) ;
		double xx = log(nbar*0.1) ;
		
		double alpha = (*coefs)(n_coefs-1) ;
		double sum_poly = 0.;
		for (int i=0; i < n_coefs-1; i++)
			sum_poly += (*coefs)(i)*pow(xx, double(i));
	    return m_n*Unites::mevpfm3*(exp(xx)*(1.+ddd) + exp((alpha+1) * xx)*sum_poly) - Lambda ;
    }
    else if ( (ent_lim1 <= ent) && (ent < ent_lim2)){
      double nbar = pow((expm1(ent)-ddd_GPP)*(Gamma_GPP-1.)/(Kappa_GPP*Gamma_GPP), 1./(Gamma_GPP-1.));
      return nbar*(1 + (Gamma_GPP-1.)/Gamma_GPP * (1./(Gamma_GPP-1.) * expm1(ent) + ddd_GPP)) - Lambda_GPP;
    }
    else if ( (0 < ent) && (ent < ent_lim1)){
		
		return eos_low->ener_ent_p(ent) ;
    }
    else{
	return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Pseudo_polytrope_1D::press_ent_p(double ent, const Param* ) const {
    if (ent >= ent_lim2 ) {
		
		double nbar = nbar_ent_p(ent) ;

		double xx = log(nbar*0.1) ;
		
		double alpha = (*coefs)(n_coefs-1) ;
		double sum_poly = 0.;
		for (int i=0; i < n_coefs-1; i++){
			double cc1 = alpha*(*coefs)(i) ;
			double cc2 = (i == n_coefs - 2) ? 0. : double(i+1)*(*coefs)(i+1) ;
			sum_poly  += (cc1 + cc2)*pow(xx, double(i));
		}
		
	    return m_n*Unites::mevpfm3*(exp((alpha+1.)*xx) * sum_poly) + Lambda;
    }
    else if ( (ent_lim1 <= ent) && (ent < ent_lim2)){
      return Kappa_GPP*pow((expm1(ent)-ddd_GPP)*(Gamma_GPP-1.)/(Kappa_GPP*Gamma_GPP), Gamma_GPP/(Gamma_GPP-1.)) + Lambda_GPP;
    }
    else if ( (0 < ent) && (ent < ent_lim1)){
		  return eos_low->press_ent_p(ent) ;
    }
    else{
	return 0 ;
    }
}

// dln(n)/ln(h) from enthalpy
//---------------------------

double Pseudo_polytrope_1D::der_nbar_ent_p(double , const Param* ) const {
	c_est_pas_fait("Pseudo_polytrope_1D::der_nbar_ent_p") ;
	return 0.;
}

// dln(e)/ln(h) from enthalpy
//---------------------------

double Pseudo_polytrope_1D::der_ener_ent_p(double, const Param* ) const {
    c_est_pas_fait("Pseudo_polytrope_1D::der_ener_ent_p") ;
    return 0.;
}

// dln(p)/ln(h) from enthalpy 
//---------------------------

double Pseudo_polytrope_1D::der_press_ent_p(double, const Param* ) const {
	c_est_pas_fait("Pseudo_polytrope_1D::der_press_ent_p") ;
	return 0. ;
}

double Pseudo_polytrope_1D::csound_square_ent_p(double ent, const Param*) const {
	if (ent >= ent_lim2 ) {
		double nbar = nbar_ent_p(ent) ;
		double xx = log(nbar*0.1) ;
		
		double alpha = (*coefs)(n_coefs-1) ;
		double sum_poly = 0., sum_der_poly = 0., sum_der2_poly = 0.;
		for (int i=0; i < n_coefs-1; i++){
			sum_poly      += (*coefs)(i)*pow(xx, double(i));
			sum_der_poly  += (i == 0) ? 0. : double(i) * (*coefs)(i) * pow(xx, double(i-1)) ;
			sum_der2_poly += (i < 2)  ? 0. : double(i) * double(i-1) * (*coefs)(i) * pow(xx, double(i-2)) ;
		}
		double num   = (alpha*(alpha+1.) * sum_poly + (2.*alpha+1.) * sum_der_poly + sum_der2_poly)*exp(alpha*xx) ;
		double denom = 1. + ddd + ((alpha+1.) * sum_poly + sum_der_poly)*exp(alpha*xx) ;
		
		return num/denom ;
	}
  else if ( (ent_lim1 <= ent) && (ent < ent_lim2)){
    return (Gamma_GPP-1.)*(1.-exp(-ent)*(1.+ddd_GPP));
  }
	else if ( (0 < ent) && (ent < ent_lim1)){
		
		return eos_low->csound_square_ent_p(ent) ;
	}
	else{
	return 0 ;
	}   
}

}
