/*
 * Methods of the class Eos_incomp_newt.
 *
 * (see file eos.h for documentation).
 */

/*
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2016/12/05 16:17:51  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.5  2014/10/13 08:52:53  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.4  2014/10/06 15:13:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.3  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/04/09 14:32:15  e_gourgoulhon
 * 1/ Added extra parameters in EOS computational functions (argument par)
 * 2/ New class MEos for multi-domain EOS
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.5  2001/02/07  09:48:48  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *      der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 *
 * Revision 2.4  2000/02/14  14:49:55  eric
 * Modif affichage.
 *
 * Revision 2.3  2000/02/14  14:34:02  eric
 * Ajout du constructeur par lecture de fichier formate.
 *
 * Revision 2.2  2000/01/21  15:18:30  eric
 * Ajout des operateurs de comparaison == et !=
 *
 * Revision 2.1  2000/01/19  08:53:39  eric
 * Ajout du set_name dans les constructeurs standards.
 *
 * Revision 2.0  2000/01/18  16:11:38  eric
 * *** empty log message ***
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
#include "cmp.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor with ent0 = 1
// ---------------------------------
namespace Lorene {
Eos_incomp_newt::Eos_incomp_newt(double rho_c) : Eos_incomp(rho_c) {

    set_name("Newtonian EOS for incompressible matter") ;
     
}

// Standard constructor with ent0 specified
// ---------------------------------------
Eos_incomp_newt::Eos_incomp_newt(double rho_c, double ent_c) : 
		    Eos_incomp(rho_c, ent_c) {

    set_name("Newtonian EOS for incompressible matter") ;
     
}
  
// Copy constructor
// ----------------
Eos_incomp_newt::Eos_incomp_newt(const Eos_incomp_newt& eosi) : 
		    Eos_incomp(eosi) {} 
  

// Constructor from a binary file
// ------------------------------
Eos_incomp_newt::Eos_incomp_newt(FILE* fich) : Eos_incomp(fich) {}

// Constructor from a formatted file
// ---------------------------------
Eos_incomp_newt::Eos_incomp_newt(ifstream& fich) : Eos_incomp(fich) {}

			//--------------//
			//  Destructor  //
			//--------------//

Eos_incomp_newt::~Eos_incomp_newt(){
    
    // does nothing
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Eos_incomp_newt::operator=(const Eos_incomp_newt& eosi) {
    
    set_name(eosi.name) ; 
    
    rho0 = eosi.rho0 ; 
    ent0 = eosi.ent0 ; 

}

			//------------------------//
			//  Comparison operators  //
			//------------------------//

bool Eos_incomp_newt::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_incomp_newt !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Eos_incomp_newt& eos = 
			    dynamic_cast<const Eos_incomp_newt&>( eos_i ) ; 

	if (eos.rho0 != rho0) {
	    cout 
	    << "The two Eos_incomp_newt have different rho0 : " << rho0 << " <-> " 
		<< eos.rho0 << endl ; 
	    resu = false ; 
	}

	if (eos.ent0 != ent0) {
	    cout 
	    << "The two Eos_incomp_newt have different ent0 : " << ent0 << " <-> " 
		<< eos.ent0 << endl ; 
	    resu = false ; 
	}

    }
    
    return resu ; 
    
}

bool Eos_incomp_newt::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}



			//------------//
			//  Outputs   //
			//------------//

void Eos_incomp_newt::sauve(FILE* fich) const {

    Eos_incomp::sauve(fich) ; 
    
}

ostream& Eos_incomp_newt::operator>>(ostream & ost) const {
    
    ost << "EOS of class Eos_incomp_newt (Newtonian incompressible matter) : " 
	<< endl ; 
    ost << "   Constant density : " << rho0 << " rho_nuc" << endl ; 
    ost << "   Log-enthalpy threshold for non-zero density : " << ent0 
	<< " c^2" <<  endl ; 
    
    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy 
//------------------------------

double Eos_incomp_newt::nbar_ent_p(double ent, const Param* ) const {

    if ( ent >= ent0 ) {

	return rho0 ;
    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Eos_incomp_newt::ener_ent_p(double ent, const Param* ) const {

    if ( ent >= ent0 ) {

	return rho0 ;
    }
    else{
	return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Eos_incomp_newt::press_ent_p(double ent, const Param* ) const {

    if ( ent >= ent0 ) {

	return rho0 * ent ;
    }
    else{
	return 0 ;
    }
}

// dln(n)/ln(h) from enthalpy
//---------------------------

double Eos_incomp_newt::der_nbar_ent_p(double ent, const Param* ) const {

    if ( ent >= ent0 ) {

	return 0 ;
    }
    else{
	return 0 ;
    }
}

// dln(e)/ln(h) from enthalpy
//---------------------------

double Eos_incomp_newt::der_ener_ent_p(double ent, const Param* ) const {

    if ( ent >= ent0 ) {

	return 0 ;
    }
    else{
	return 0 ;
    }
}

// dln(p)/ln(h) from enthalpy
//---------------------------

double Eos_incomp_newt::der_press_ent_p(double ent, const Param* ) const {
    
    if ( ent >= ent0 ) {

	return double(1) ;

    }
    else{
	return 0 ;
    }
}
}
