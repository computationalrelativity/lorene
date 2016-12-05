/*
 * Methods of the class Eos_poly_newt.
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
 * Revision 2.5  2001/02/23  15:16:51  eric
 * Continuite en ent=0 des quantites derivees.
 *
 * Revision 2.4  2001/02/07  09:50:11  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *      der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 *
 * Revision 2.3  2000/02/14  14:49:41  eric
 * Modif affichage.
 *
 * Revision 2.2  2000/02/14  14:33:33  eric
 * Ajout du constructeur par lecture de fichier formate.
 *
 * Revision 2.1  2000/01/21  15:18:56  eric
 * Ajout des operateurs de comparaison == et !=
 *
 * Revision 2.0  2000/01/18  15:14:12  eric
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

// Standard constructor 
// --------------------
namespace Lorene {
Eos_poly_newt::Eos_poly_newt(double gamma, double kappa) : 
	       Eos_poly(gamma, kappa) {

    set_name("Newtonian polytropic EOS") ; 

}  

  
// Copy constructor
// ----------------
Eos_poly_newt::Eos_poly_newt(const Eos_poly_newt& eosi) : Eos_poly(eosi) {}
  

// Constructor from a binary file
// ------------------------------
Eos_poly_newt::Eos_poly_newt(FILE* fich) : Eos_poly(fich) {}
	       
// Constructor from a formatted file
// ---------------------------------
Eos_poly_newt::Eos_poly_newt(ifstream& fich) : Eos_poly(fich) {}
	       

			//--------------//
			//  Destructor  //
			//--------------//

Eos_poly_newt::~Eos_poly_newt(){
    
    // does nothing
        
}
			//--------------//
			//  Assignment  //
			//--------------//

void Eos_poly_newt::operator=(const Eos_poly_newt& eosi) {
    
    set_name(eosi.name) ; 
    
    gam = eosi.gam ; 
    kap = eosi.kap ; 
    m_0 = eosi.m_0 ; 
    
    set_auxiliary() ; 
    
}
			//------------------------//
			//  Comparison operators  //
			//------------------------//

bool Eos_poly_newt::operator==(const Eos& eos_i) const {
    
    bool resu = true ; 
    
    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type Eos_poly_newt !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Eos_poly_newt& eos = dynamic_cast<const Eos_poly_newt&>( eos_i ) ; 

	if (eos.gam != gam) {
	    cout 
	    << "The two Eos_poly_newt have different gamma : " << gam << " <-> " 
		<< eos.gam << endl ; 
	    resu = false ; 
	}

	if (eos.kap != kap) {
	    cout 
	    << "The two Eos_poly_newt have different kappa : " << kap << " <-> " 
		<< eos.kap << endl ; 
	    resu = false ; 
	}

	if (eos.m_0 != m_0) {
	    cout 
	    << "The two Eos_poly_newt have different m_0 : " << m_0 << " <-> " 
		<< eos.m_0 << endl ; 
	    resu = false ; 
	}
	
    }
    
    return resu ; 
    
}

bool Eos_poly_newt::operator!=(const Eos& eos_i) const {
 
    return !(operator==(eos_i)) ; 
       
}


			//------------//
			//  Outputs   //
			//------------//

void Eos_poly_newt::sauve(FILE* fich) const {

    Eos_poly::sauve(fich) ; 
       
}

ostream& Eos_poly_newt::operator>>(ostream & ost) const {
    
    ost << "EOS of class Eos_poly_newt (Newtonian polytrope) : " << endl ; 
    ost << "   Adiabatic index gamma :      " << gam << endl ; 
    ost << "   Pressure coefficient kappa : " << kap << 
	   " rho_nuc c^2 / n_nuc^gamma" << endl ; 
    
    return ost ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy 
//------------------------------

double Eos_poly_newt::nbar_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

	return pow( gam1sgamkap * ent, unsgam1 ) ;
    }
    else{
	return 0 ;
    }
}

// Energy density from enthalpy
//------------------------------

double Eos_poly_newt::ener_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

	double nn = pow( gam1sgamkap * ent, unsgam1 ) ;

	double pp = kap * pow( nn, gam ) ;

	return  unsgam1 * pp + m_0 * nn ;
    }
    else{
	return 0 ;
    }
}

// Pressure from enthalpy
//------------------------

double Eos_poly_newt::press_ent_p(double ent, const Param* ) const {

    if ( ent > double(0) ) {

	double nn = pow( gam1sgamkap * ent, unsgam1 ) ;

	return kap * pow( nn, gam ) ;

    }
    else{
	return 0 ;
    }
}

// dln(n)/ln(h) from enthalpy
//---------------------------

double Eos_poly_newt::der_nbar_ent_p(double , const Param* ) const {

    return double(1) / gam1 ;

}

// dln(e)/ln(h) from enthalpy
//---------------------------

double Eos_poly_newt::der_ener_ent_p(double ent, const Param* ) const {
    
    if ( ent > double(0) ) {


	double nn = pow( gam1sgamkap * ( exp(ent) - double(1) ),  
				     unsgam1 ) ;

	double pp = kap * pow( nn, gam ) ;

	double ee =  unsgam1 * pp + m_0 * nn ; 
	

	return ( double(1) + pp / ee) / gam1 ;

    }
    else{
	return double(1) / gam1 ;   //  to ensure continuity at ent=0
    }
}

// dln(p)/ln(h) from enthalpy 
//---------------------------

double Eos_poly_newt::der_press_ent_p(double, const Param* ) const {
    
    return gam / gam1 ;

}

}
