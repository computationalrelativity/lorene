/*
 * Methods for the class Eos_strange
 *
 * (see file eos.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000 J. Leszek Zdunik
 *   Copyright (c) 2000-2001 Eric Gourgoulhon
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


char eos_strange_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2004/03/25 10:29:02  j_novak
 * All LORENE's units are now defined in the namespace Unites (in file unites.h).
 *
 * Revision 1.4  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/04/09 14:32:15  e_gourgoulhon
 * 1/ Added extra parameters in EOS computational functions (argument par)
 * 2/ New class MEos for multi-domain EOS
 *
 * Revision 1.2  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.2  2001/02/07  09:49:47  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *       der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 *
 * Revision 2.1  2000/10/25  10:54:41  eric
 * Correction erreur dans la densite d'energie (conversion d'unite).
 *
 * Revision 2.0  2000/10/24  15:29:11  eric
 * *** empty log message ***
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
#include "eos.h"
#include "cmp.h"
#include "utilitaires.h"
#include "unites.h"

		    //------------------------------------//
		    //		Constructors		  //
		    //------------------------------------//

// Standard constructor
// --------------------
Eos_strange::Eos_strange(double n0_b60_i, double b60_i, double ent0_i, 
		    double eps_fit_i, double rho0_b60_i) :
	     Eos("Strange matter EOS from Zdunik (2000)"), 
	     n0_b60(n0_b60_i), 
	     b60(b60_i), 
	     ent0(ent0_i), 
	     eps_fit(eps_fit_i), 
	     rho0_b60(rho0_b60_i) {

    set_auxiliary() ; 
		    
}		    

// Copy constructor
// -----------------

Eos_strange::Eos_strange(const Eos_strange& eos_i) : 
	    Eos(eos_i), 
	    n0_b60(eos_i.n0_b60), 
	    b60(eos_i.b60), 
	    ent0(eos_i.ent0), 
	    eps_fit(eos_i.eps_fit), 
	    rho0_b60(eos_i.rho0_b60)	{

    set_auxiliary() ; 

}


// Constructor from binary file
// ----------------------------
Eos_strange::Eos_strange(FILE* fich) : 
	Eos(fich) {
        
    fread_be(&n0_b60, sizeof(double), 1, fich) ;		
    fread_be(&b60, sizeof(double), 1, fich) ;		
    fread_be(&ent0, sizeof(double), 1, fich) ;		
    fread_be(&eps_fit, sizeof(double), 1, fich) ;		
    fread_be(&rho0_b60, sizeof(double), 1, fich) ;		
    
    set_auxiliary() ; 

}

// Constructor from a formatted file
// ---------------------------------
Eos_strange::Eos_strange(ifstream& fich) : 
	Eos(fich) {

    char blabla[80] ;
        
    fich >> n0_b60 ; fich.getline(blabla, 80) ;
    fich >> b60 ; fich.getline(blabla, 80) ;
    fich >> ent0 ; fich.getline(blabla, 80) ;
    fich >> eps_fit ; fich.getline(blabla, 80) ;
    fich >> rho0_b60 ; fich.getline(blabla, 80) ;
    
    set_auxiliary() ; 

}
			//--------------//
			//  Destructor  //
			//--------------//

Eos_strange::~Eos_strange(){
    
    // does nothing
        
}

			//--------------//
			//  Assignment  //
			//--------------//

void Eos_strange::operator=(const Eos_strange& eosi) {
    
    set_name(eosi.name) ; 
    
    n0_b60 = eosi.n0_b60 ; 
    b60 = eosi.b60 ; 
    ent0 = eosi.ent0 ; 
    eps_fit = eosi.eps_fit ; 
    rho0_b60 = eosi.rho0_b60 ; 
    
    set_auxiliary() ; 
    
}


		  //-----------------------//
		  //	Miscellaneous	   //
		  //-----------------------//

void Eos_strange::set_auxiliary() {
    
  using namespace Unites ;
    
    rho0 = b60 * rho0_b60 * mevpfm3 ;
    
    b34 = pow(b60, double(0.75)) ;

    n0 = b34 * n0_b60 * double(10) ; 	// 10 : fm^{-3} --> 0.1 fm^{-3}
    
    fach = (double(4) + eps_fit) / (double(1) + eps_fit) ;
   
}


			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool Eos_strange::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_strange !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Eos_strange& eos = dynamic_cast<const Eos_strange&>( eos_i ) ; 

	if (eos.n0_b60 != n0_b60) {
	    cout 
	    << "The two Eos_strange have different n0_b60 : " << n0_b60 << " <-> " 
		<< eos.n0_b60 << endl ; 
	    resu = false ; 
	}

	if (eos.b60 != b60) {
	    cout 
	    << "The two Eos_strange have different b60 : " << b60 << " <-> " 
		<< eos.b60 << endl ; 
	    resu = false ; 
	}

	if (eos.ent0 != ent0) {
	    cout 
	    << "The two Eos_strange have different ent0 : " << ent0 << " <-> " 
		<< eos.ent0 << endl ; 
	    resu = false ; 
	}

	if (eos.eps_fit != eps_fit) {
	    cout 
	    << "The two Eos_strange have different eps_fit : " << eps_fit 
		<< " <-> " << eos.eps_fit << endl ; 
	    resu = false ; 
	}

	if (eos.rho0_b60 != rho0_b60) {
	    cout 
	    << "The two Eos_strange have different rho0_b60 : " << rho0_b60 
		<< " <-> " << eos.rho0_b60 << endl ; 
	    resu = false ; 
	}

	
    }
    
    return resu ; 
    
}

bool Eos_strange::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}

			//------------//
			//  Outputs   //
			//------------//

void Eos_strange::sauve(FILE* fich) const {

    Eos::sauve(fich) ; 
    
    fwrite_be(&n0_b60, sizeof(double), 1, fich) ;	
    fwrite_be(&b60, sizeof(double), 1, fich) ;	
    fwrite_be(&ent0, sizeof(double), 1, fich) ;	
    fwrite_be(&eps_fit, sizeof(double), 1, fich) ;	
    fwrite_be(&rho0_b60, sizeof(double), 1, fich) ;	

}		    

ostream& Eos_strange::operator>>(ostream & ost) const {
    
    ost << 
    "EOS of class Eos_strange (Strange matter EOS from Zdunik (2000)) : " 
    << endl ; 
    ost << "   Baryon density at zero pressure : " << n0_b60 
	<< " * B_{60}^{3/4}" << endl ; 
    ost << "   Bag constant B :  " << b60 << " * 60 MeV/fm^3"<< endl ;
    ost << 
    "   Log-enthalpy threshold for setting the energy density to non-zero: " 
	<< endl << "              " << ent0 << endl ; 
    ost << "   Fitting parameter eps_fit : " << eps_fit << endl ; 
    ost << "   Energy density at zero pressure : " << rho0_b60 
	<< " * B_{60}  MeV/fm^3" << endl ; 
    
    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy 
//------------------------------

double Eos_strange::nbar_ent_p(double ent, const Param* ) const {

    if ( ent > ent0 ) {

	return n0 * exp( double(3) * ent / (double(1) + eps_fit))  ;

    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Eos_strange::ener_ent_p(double ent, const Param* ) const {


    if ( ent > ent0 ) {

	double pp = ( exp(fach * ent) - 1) / fach  * rho0 ;

	return rho0 + double(3) * pp / (double(1) + eps_fit) ;

    }
    else{
	return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Eos_strange::press_ent_p(double ent, const Param* ) const {

    if ( ent > ent0 ) {

	return ( exp(fach * ent) - 1) / fach  * rho0 ;

    }
    else{
	return 0 ;
    }
}



// dln(n)/ln(H) from enthalpy
//---------------------------

double Eos_strange::der_nbar_ent_p(double ent, const Param* ) const {

    if ( ent > ent0 ) {

	return double(3) * ent / ( double(1) +  eps_fit ) ;

    }
    else{
	return 0 ;
    }
}

// dln(e)/ln(H) from enthalpy
//---------------------------

double Eos_strange::der_ener_ent_p(double ent, const Param* ) const {

    if ( ent > ent0 ) {

	double xx = fach * ent ;

	return xx / ( double(1) +
		( double(1) + eps_fit ) / double(3) * exp(-xx) ) ;

    }
    else{
	return 0 ;
    }
}

// dln(p)/ln(H) from enthalpy
//---------------------------

double Eos_strange::der_press_ent_p(double ent, const Param* ) const {
    
    if ( ent > ent0 ) {

	double xx = fach * ent ; 

	return xx / ( double(1) - exp(-xx) ) ;

    }
    else{
	return 0 ;
    }
}

