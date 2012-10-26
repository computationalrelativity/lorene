/*
 * Methods of the class Eos_Fermi.
 *
 * (see file eos.h for documentation).
 */

/*
 *   Copyright (c) 2012 Eric Gourgoulhon
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


char eos_fermi_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2012/10/26 14:09:33  e_gourgoulhon
 * Added new class Eos_Fermi
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
#include "eos.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor with g_s = 2
// ----------------------------------

Eos_Fermi::Eos_Fermi(double mass) : 
	Eos("Degenerate ideal Fermi gas"), 
	m_0(mass), g_s(2) {

    set_auxiliary() ; 

}  

// Standard constructor
// --------------------

Eos_Fermi::Eos_Fermi(double mass, int g_degen) : 
	Eos("Degenerate ideal Fermi gas"), 
	m_0(mass), g_s(g_degen) {

    set_auxiliary() ; 

}  


// Copy constructor
// ----------------
Eos_Fermi::Eos_Fermi(const Eos_Fermi& eosi) : 
	Eos(eosi), 
	m_0(eosi.m_0), g_s(eosi.g_s) {

    set_auxiliary() ; 

}  
  

// Constructor from binary file
// ----------------------------
Eos_Fermi::Eos_Fermi(FILE* fich) : 
	Eos(fich) {
        
    fread_be(&m_0, sizeof(double), 1, fich) ;		
    fread_be(&g_s, sizeof(int), 1, fich) ;		

    set_auxiliary() ; 

}


// Constructor from a formatted file
// ---------------------------------
Eos_Fermi::Eos_Fermi(ifstream& fich) : 
	Eos(fich) {

    char blabla[80] ;
        
    fich >> m_0 ; fich.getline(blabla, 80) ;
    fich >> g_s ; fich.getline(blabla, 80) ;
 
    set_auxiliary() ; 

}
			//--------------//
			//  Destructor  //
			//--------------//

Eos_Fermi::~Eos_Fermi(){
    
    // does nothing
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Eos_Fermi::operator=(const Eos_Fermi& eosi) {
    
    set_name(eosi.name) ; 
    
    m_0 = eosi.m_0 ;
    g_s = eosi.g_s ;

    set_auxiliary() ;

}


		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

void Eos_Fermi::set_auxiliary() {

    n_0 =  2.19780778967127e-26 * double(g_s) * m_0 * m_0 * m_0 ;

    p_0 = 1.34236578162651e-10 * n_0 * m_0 ;
    
    ener_0 = double(3) * p_0 ; 
    
    cout << "n_0 = " << n_0 << endl ; 
    cout << "p_0 = " << p_0 << endl ; 
    cout << "ener_0 = " << ener_0 << endl ; 

}

double Eos_Fermi::get_m() const {
    return m_0 ;
}

int Eos_Fermi::get_g_degen() const {
    return g_s ; 
}



			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_Fermi::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_Fermi !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Eos_Fermi& eos = dynamic_cast<const Eos_Fermi&>( eos_i ) ; 

	if (eos.m_0 != m_0) {
	    cout 
	    << "The two Eos_Fermi have different m_0 : " << m_0 << " <-> " 
		<< eos.m_0 << endl ; 
	    resu = false ; 
	}

	if (eos.g_s != g_s) {
	    cout 
	    << "The two Eos_Fermi have different g_s : " << g_s << " <-> " 
		<< eos.g_s << endl ; 
	    resu = false ; 
	}

	
    }
    
    return resu ; 
    
}

bool Eos_Fermi::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}


			//------------//
			//  Outputs   //
			//------------//

void Eos_Fermi::sauve(FILE* fich) const {

    Eos::sauve(fich) ; 
    
    fwrite_be(&m_0, sizeof(double), 1, fich) ;	
    fwrite_be(&g_s, sizeof(int), 1, fich) ;

}

ostream& Eos_Fermi::operator>>(ostream & ost) const {
    
    ost << "EOS of class Eos_Fermi (degenerate ideal Fermi gas) : " << endl ; 
    ost << "   Fermion mass : " << m_0 << " eV/c2" << endl ;
    ost << "   Degeneracy factor : " << g_s << endl ;

    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double Eos_Fermi::nbar_ent_p(double ent, const Param* ) const {

    if ( ent > 0 ) {
        double x2 =  exp(double(2)*ent) - double(1) ; 
        double x = sqrt( x2 ) ; 
        // cout << "nbar: ent,  x = " << ent << " " << x << endl ; 
        return n_0 * x2 * x ;
    }
    else{
        return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Eos_Fermi::ener_ent_p(double ent, const Param* ) const {

    if ( ent > 0 ) {

        double hh = exp(ent) ; 
        double x2 = hh*hh - double(1)  ; 
        double x = sqrt(x2) ; 
        // cout << "ener: ent,  x = " << ent << " " << x << endl ; 
        return ener_0 * (x*(double(2)*x2 + double(1))*hh - log(x+hh)) ;
    }
    else{
        return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Eos_Fermi::press_ent_p(double ent, const Param* ) const {

    if ( ent > 0 ) {

        double hh = exp(ent) ; 
        double x2 = hh*hh - double(1)  ; 
        double x = sqrt(x2) ; 
        // cout << "press: ent,  x = " << ent << " " << x << endl ; 
        return p_0 * (x*(double(2)*x2 - double(3))*hh + double(3)*log(x+hh)) ;

    }
    else{
        return 0 ;
    }
}

// dln(n)/ln(H) from enthalpy
//---------------------------

double Eos_Fermi::der_nbar_ent_p(double , const Param* ) const {

    cerr << "Eos_Fermi::der_nbar_ent_p : not implemented yet ! " << endl ; 
    abort() ; 
 
}

// dln(e)/ln(H) from enthalpy
//---------------------------

double Eos_Fermi::der_ener_ent_p(double , const Param* ) const {

    cerr << "Eos_Fermi::der_ener_ent_p : not implemented yet ! " << endl ; 
    abort() ; 

 }

// dln(p)/ln(H) from enthalpy
//---------------------------

double Eos_Fermi::der_press_ent_p(double , const Param* ) const {

    cerr << "Eos_Fermi::der_press_ent_p : not implemented yet ! " << endl ; 
    abort() ; 
    
}

