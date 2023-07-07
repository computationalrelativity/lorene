/*
 * Methods of the class Piecewise_polytrope_1D.
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
void huntm(const Tbl& xx, double& x, int& i_low) ;

Piecewise_polytrope_1D::Piecewise_polytrope_1D(const Tbl& Gamma, const Tbl& Kappa, const Tbl& Lambda, 
                                               const Tbl& AA, const Tbl& N_lim, const Tbl& Ent_lim, double gamma, double kappa, double n_lim0): 
	      gamma_high(0x0), kappa_high(0x0), Lambda_high(0x0), a_high(0x0), n_lim_high(0x0), ent_lim_high(0x0), gamma_low(gamma), kappa_low(kappa), n_lim(n_lim0) {

    set_name("Pseudo-polytropic fit of cold, beta-equilibrated EoS") ; 
    
	assert(Gamma.get_ndim() == 1) ;
	assert(Kappa.get_ndim() == 1) ;
	assert(Lambda.get_ndim() == 1) ;
	assert(AA.get_ndim() == 1) ;
	assert(N_lim.get_ndim() == 1) ;
	assert(Ent_lim.get_ndim() == 1) ;
  gamma_high   = new Tbl(Gamma) ;
	n_param_high = gamma_high->get_taille() ;
  assert(Kappa.get_taille() == n_param_high) ;
  kappa_high = new Tbl(Kappa) ;
  assert(Lambda.get_taille() == n_param_high) ;
  Lambda_high = new Tbl(Lambda) ;
  assert(AA.get_taille() == n_param_high) ;
  a_high = new Tbl(AA) ;
  assert(N_lim.get_taille() == n_param_high) ;
  n_lim_high = new Tbl(N_lim) ;
  assert(Ent_lim.get_taille() == n_param_high) ;
  ent_lim_high = new Tbl(Ent_lim) ; // maybe ent_lim should be computed here ? Anyway this constructor is very likely not to be used a lot so maybe it is not a problem.
	
  eos_low = new Eos_poly(gamma_low, kappa_low) ;
}  

  
// Copy constructor
// ----------------
Piecewise_polytrope_1D::Piecewise_polytrope_1D(const Piecewise_polytrope_1D& eosi): Eos(eosi), gamma_high(eosi.gamma_high), kappa_high(eosi.kappa_high),
                                      Lambda_high(eosi.Lambda_high), a_high(eosi.a_high), n_lim_high(eosi.n_lim_high), ent_lim_high(eosi.ent_lim_high), gamma_low(eosi.gamma_low),
																			kappa_low(eosi.kappa_low), n_lim(eosi.n_lim), ent_lim(eosi.ent_lim), n_param_high(eosi.n_param_high) {}
  

// Constructor from a binary file
// ------------------------------
Piecewise_polytrope_1D::Piecewise_polytrope_1D(FILE* fich): Eos(fich) {
    
    fread_be(&n_param_high, sizeof(int), 1, fich) ;
	gamma_high = new Tbl(fich) ;
	kappa_high = new Tbl(fich) ;
	Lambda_high = new Tbl(fich) ;
	a_high = new Tbl(fich) ;
	n_lim_high = new Tbl(fich) ;
	ent_lim_high = new Tbl(fich) ;
    fread_be(&kappa_low, sizeof(double), 1, fich) ;
    fread_be(&gamma_low, sizeof(double), 1, fich) ;
	fread_be(&n_lim, sizeof(double), 1, fich) ;
	fread_be(&ent_lim, sizeof(double), 1, fich) ;
		
	eos_low = new Eos_poly(gamma_low, kappa_low) ;
}
	       
// Constructor from a formatted file
// ---------------------------------
Piecewise_polytrope_1D::Piecewise_polytrope_1D(ifstream& fich): Eos(fich) {
	
	fich >> n_param_high ; fich.ignore(1000, '\n') ;
  
	gamma_high = new Tbl(n_param_high) ; gamma_high->set_etat_qcq() ;
	for (int i=0; i < n_param_high; i++)
		fich >> gamma_high->set(i) ;
	fich.ignore(1000, '\n') ;

	kappa_high = new Tbl(n_param_high) ; kappa_high->set_etat_qcq() ;
	for (int i=0; i < n_param_high; i++)
		fich >> kappa_high->set(i) ;
	fich.ignore(1000, '\n') ;

	Lambda_high = new Tbl(n_param_high) ; Lambda_high->set_etat_qcq() ;
	for (int i=0; i < n_param_high; i++)
		fich >> Lambda_high->set(i) ;
	fich.ignore(1000, '\n') ;

	a_high = new Tbl(n_param_high) ; a_high->set_etat_qcq() ;
	for (int i=0; i < n_param_high; i++)
		fich >> a_high->set(i) ;
	fich.ignore(1000, '\n') ;

	n_lim_high = new Tbl(n_param_high) ; n_lim_high->set_etat_qcq() ;
	for (int i=0; i < n_param_high; i++)
		fich >> n_lim_high->set(i) ;
	fich.ignore(1000, '\n') ;

	ent_lim_high = new Tbl(n_param_high) ; ent_lim_high->set_etat_qcq() ;
	for (int i=0; i < n_param_high; i++)
		fich >> ent_lim_high->set(i) ;
	fich.ignore(1000, '\n') ;
	
  n_lim = (*n_lim_high)(0) ;
  ent_lim = (*ent_lim_high)(0) ;
  
	fich >> kappa_low ; fich.ignore(1000, '\n') ;
	
	fich >> gamma_low ; fich.ignore(1000, '\n') ;
	
  // double Gamma1 = (*gamma_high)(0) ;
  // double Kappa1 = (*kappa_high)(0) ;
  // double a1     = (*a_high)(0) ;
  // ent_lim = log(1. + a1 + Kappa1*Gamma1/(Gamma1-1.) * pow(0.1*n_lim,Gamma1-1.)) ;
		
	eos_low = new Eos_poly(gamma_low, kappa_low) ;
	
}
	       

			//--------------//
			//  Destructor  //
			//--------------//

Piecewise_polytrope_1D::~Piecewise_polytrope_1D(){
    
    if (gamma_high != 0x0)   delete gamma_high ;
    if (kappa_high != 0x0)   delete kappa_high ;
    if (Lambda_high != 0x0)   delete Lambda_high ;
    if (a_high != 0x0)   delete a_high ;
    if (n_lim_high != 0x0)   delete n_lim_high ;
    if (ent_lim_high != 0x0)   delete ent_lim_high ;
	if (eos_low != 0x0) delete eos_low ;
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Piecewise_polytrope_1D::operator=(const Piecewise_polytrope_1D& eosi) {
    
    set_name(eosi.name) ; 
    
    gamma_high     = eosi.gamma_high ; 
    kappa_high     = eosi.kappa_high ; 
    Lambda_high    = eosi.Lambda_high ; 
    a_high         = eosi.a_high ; 
    n_lim_high     = eosi.n_lim_high ; 
    ent_lim_high   = eosi.ent_lim_high ; 
    n_lim          = eosi.n_lim ; 
	  n_param_high   = eosi.n_param_high ;
	  gamma_low      = eosi.gamma_low ;
	  kappa_low      = eosi.kappa_low ;
  	ent_lim        = eosi.ent_lim ;
	  eos_low        = eosi.eos_low ;
    
}
			//------------------------//
			//  Comparison operators  //
			//------------------------//

bool Piecewise_polytrope_1D::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Piecewise_polytrope_1D !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Piecewise_polytrope_1D& eos = dynamic_cast<const Piecewise_polytrope_1D&>( eos_i ) ; 
	
	for (int i=0; i < n_param_high; i++){
		if ((*eos.gamma_high)(i) != (*gamma_high)(i)) {
			cout 
			<< "The two Piecewise_polytrope_1D have different coefficients: " << (*gamma_high)(i) << " <-> " 
			<< (*eos.gamma_high)(i) << endl ; 
			resu = false ; 
		}
    
    if ((*eos.kappa_high)(i) != (*kappa_high)(i)) {
			cout 
			<< "The two Piecewise_polytrope_1D have different coefficients: " << (*kappa_high)(i) << " <-> " 
			<< (*eos.kappa_high)(i) << endl ; 
			resu = false ; 
		}
    
    if ((*eos.Lambda_high)(i) != (*Lambda_high)(i)) {
			cout 
			<< "The two Piecewise_polytrope_1D have different coefficients: " << (*Lambda_high)(i) << " <-> " 
			<< (*eos.Lambda_high)(i) << endl ; 
			resu = false ; 
		}
    
    if ((*eos.a_high)(i) != (*a_high)(i)) {
			cout 
			<< "The two Piecewise_polytrope_1D have different coefficients: " << (*a_high)(i) << " <-> " 
			<< (*eos.a_high)(i) << endl ; 
			resu = false ; 
		}
    
    if ((*eos.n_lim_high)(i) != (*n_lim_high)(i)) {
			cout 
			<< "The two Piecewise_polytrope_1D have different coefficients: " << (*n_lim_high)(i) << " <-> " 
			<< (*eos.n_lim_high)(i) << endl ; 
			resu = false ; 
		}
    
    if ((*eos.ent_lim_high)(i) != (*ent_lim_high)(i)) {
			cout 
			<< "The two Piecewise_polytrope_1D have different coefficients: " << (*ent_lim_high)(i) << " <-> " 
			<< (*eos.ent_lim_high)(i) << endl ; 
			resu = false ; 
		}
  }

	if (eos.n_lim != n_lim) {
	    cout 
	    << "The two Piecewise_polytrope_1D have different limiting densities n_lim: " << n_lim << " <-> " 
		<< eos.n_lim << endl ; 
	    resu = false ; 
	}
	
	if (eos.ent_lim != ent_lim) {
	    cout 
	    << "The two Piecewise_polytrope_1D have different limiting enthalpies ent_lim: " << ent_lim << " <-> " 
		<< eos.ent_lim << endl ; 
	    resu = false ; 
	}
	
	if (eos.n_param_high != n_param_high) {
	    cout 
	    << "The two Piecewise_polytrope_1D have different number of coefficients: " << n_param_high << " <-> " 
		<< eos.n_param_high << endl ; 
	    resu = false ; 
	}
	
	if (eos.kappa_low != kappa_low) {
	    cout 
	    << "The two Piecewise_polytrope_1D have different kappa in low polytrope: " << kappa_low << " <-> " 
		<< eos.kappa_low << endl ; 
	    resu = false ; 
	}
	
	if (eos.gamma_low != gamma_low) {
	    cout 
	    << "The two Piecewise_polytrope_1D have different kappa in low polytrope: " << gamma_low << " <-> " 
		<< eos.gamma_low << endl ; 
	    resu = false ; 
	}
	
    }
    
    return resu ; 
    
}

bool Piecewise_polytrope_1D::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}


			//------------//
			//  Outputs   //
			//------------//

void Piecewise_polytrope_1D::sauve(FILE* fich) const {

    Eos::sauve(fich) ; 
    
    fwrite_be(&n_param_high, sizeof(int), 1, fich) ;	
    gamma_high->sauve(fich) ;
    kappa_high->sauve(fich) ;
    Lambda_high->sauve(fich) ;
    a_high->sauve(fich) ;
    n_lim_high->sauve(fich) ;
    ent_lim_high->sauve(fich) ;
    fwrite_be(&kappa_low, sizeof(double), 1, fich) ;
    fwrite_be(&gamma_low, sizeof(double), 1, fich) ;
	fwrite_be(&n_lim, sizeof(double), 1, fich) ;
	fwrite_be(&ent_lim, sizeof(double), 1, fich) ;	
       
}

ostream& Piecewise_polytrope_1D::operator>>(ostream & ost) const {
    
    ost << setprecision(16) << "EOS of class Piecewise_polytrope_1D (analytical fit of cold EoS): " << '\n' ; 
    ost << "   Gamma_high coefficients:            " << *gamma_high << '\n' ; 
    ost << "   Kappa_high coefficients:            " << *kappa_high << '\n' ; 
    ost << "   Lambda_high coefficients:            " << *Lambda_high << '\n' ; 
    ost << "   a_high coefficients:            " << *a_high << '\n' ; 
    ost << "   n_lim_high coefficients:            " << *n_lim_high << '\n' ; 
    ost << "   ent_lim_high coefficients:            " << *ent_lim_high << '\n' ; 
	ost << "Low densities polytrope :" << '\n' ;
	cout << *eos_low << '\n' ;
    ost << setprecision(16) << "   Limiting density n_lim:  " << 0.1*n_lim << " [fm^-3]" << '\n' ; 
	ost << setprecision(16) << "   Limiting enthalpy h_lim: " << ent_lim << " [c^2]" << endl;
    
    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy 
//------------------------------

double Piecewise_polytrope_1D::nbar_ent_p(double ent, const Param* ) const {

    if (ent > ent_lim ) {
      
      int i_low = 0;
      huntm(*ent_lim_high, ent, i_low) ;
      Tbl& Gamma  = *gamma_high ;
      Tbl& Kappa  = *kappa_high ;
      Tbl& AA     = *a_high ;
      double nn = 1./(Gamma(i_low)-1.) ;

	    return pow((exp(ent)-1.-AA(i_low))/Kappa(i_low)/(nn+1.), nn) ;
      
    }
    else if ( (0 < ent) && (ent < ent_lim)){
		
		return eos_low->nbar_ent_p(ent) ;
    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Piecewise_polytrope_1D::ener_ent_p(double ent, const Param* ) const {
    if (ent > ent_lim ) {
      
      int i_low = 0;
      huntm(*ent_lim_high, ent, i_low) ;
      Tbl& Gamma  = *gamma_high ;
      Tbl& Kappa  = *kappa_high ;
      Tbl& AA     = *a_high ;
      Tbl& Lambda = *Lambda_high ;
      double nn = 1./(Gamma(i_low)-1.) ;
      double nbar = pow((exp(ent)-1.-AA(i_low))/Kappa(i_low)/(nn+1.), nn) ;
      
		
	    return Kappa(i_low)*nn * pow(nbar, Gamma(i_low)) + (1.+AA(i_low))*nbar - Lambda(i_low) ;
    }
    else if ( (0 < ent) && (ent < ent_lim)){
		
		return eos_low->ener_ent_p(ent) ;
    }
    else{
	return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Piecewise_polytrope_1D::press_ent_p(double ent, const Param* ) const {
	if (ent > ent_lim ) {
		
      int i_low = 0;
      huntm(*ent_lim_high, ent, i_low) ;
      Tbl& Gamma  = *gamma_high ;
      Tbl& Kappa  = *kappa_high ;
      Tbl& AA     = *a_high ;
      Tbl& Lambda = *Lambda_high ;
      double nn = 1./(Gamma(i_low)-1.) ;
      double nbar = pow((exp(ent)-1.-AA(i_low))/Kappa(i_low)/(nn+1.), nn) ;

	    return Kappa(i_low) * pow(nbar, Gamma(i_low)) + Lambda(i_low) ;
    }
    else if ( (0 < ent) && (ent < ent_lim)){
		
		return eos_low->press_ent_p(ent) ;
    }
    else{
	return 0 ;
    }
}

// dln(n)/ln(h) from enthalpy
//---------------------------

double Piecewise_polytrope_1D::der_nbar_ent_p(double , const Param* ) const {
	c_est_pas_fait("Piecewise_polytrope_1D::der_nbar_ent_p") ;
	return 0. ;
}

// dln(e)/ln(h) from enthalpy
//---------------------------

double Piecewise_polytrope_1D::der_ener_ent_p(double, const Param* ) const {
    c_est_pas_fait("Piecewise_polytrope_1D::der_ener_ent_p") ;
    return 0.;
}

// dln(p)/ln(h) from enthalpy 
//---------------------------

double Piecewise_polytrope_1D::der_press_ent_p(double, const Param* ) const {
	c_est_pas_fait("Piecewise_polytrope_1D::der_press_ent_p") ;
	return 0. ;
}

double Piecewise_polytrope_1D::csound_square_ent_p(double ent, const Param*) const {
	if (ent > ent_lim ) {

    int i_low = 0;
    huntm(*ent_lim_high, ent, i_low) ;
    Tbl& AA    = *a_high ;
    Tbl& Gamma = *gamma_high ;
      
		return (Gamma(i_low)-1.)*(1.-exp(-ent)*(1.+AA(i_low))) ;
	}
	else if ( (0 < ent) && (ent < ent_lim)){
		
		return eos_low->csound_square_ent_p(ent) ;
	}
	else{
	return 0 ;
	}   
}

}
