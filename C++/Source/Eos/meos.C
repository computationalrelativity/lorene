/*
 *  Methods of class MEos
 *
 */

/*
 *   Copyright (c) 2002 Michal Bejger, Eric Gourgoulhon & Leszek Zdunik
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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

char meos_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/16 14:36:35  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/04/09 14:42:29  e_gourgoulhon
 * Dummy argument in assignment.
 *
 * Revision 1.1  2002/04/09 14:40:36  e_gourgoulhon
 * Methods for new class MEos
 *
 *
 *
 * $Header$
 *
 */

// C headers
#include <stdlib.h>

// Lorene headers
#include "headcpp.h"
#include "eos.h"
#include "utilitaires.h"
#include "param.h"


                        //--------------------------------------//
                        //              Constructors            //
                        //--------------------------------------//

MEos::MEos(int ndom_i, const Eos** mono_eos_i) : ndom(ndom_i) ,
                               constructed_from_file(false) {


        mono_eos = new const Eos* [ndom] ;

        for (int l=0; l<ndom; l++) {
                mono_eos[l] =  mono_eos_i[l] ;
        }

}


MEos::MEos(const Eos& eos1, const Eos& eos2) : ndom(2) ,
                               constructed_from_file(false) {

        mono_eos = new const Eos* [ndom] ;

        mono_eos[0] = &eos1 ;
        mono_eos[1] = &eos2 ;

}

MEos::MEos(const Eos& eos1, const Eos& eos2, const Eos& eos3) : ndom(3) ,
                               constructed_from_file(false) {

        mono_eos = new const Eos* [ndom] ;

        mono_eos[0] = &eos1 ;
        mono_eos[1] = &eos2 ;
        mono_eos[2] = &eos3 ;

}

MEos::MEos(const Eos& eos1, const Eos& eos2, const Eos& eos3, const Eos& eos4) : ndom(4) ,
                               constructed_from_file(false) {

        mono_eos = new const Eos* [ndom] ;

        mono_eos[0] = &eos1 ;
        mono_eos[1] = &eos2 ;
        mono_eos[2] = &eos3 ;
        mono_eos[3] = &eos4 ;

}

// Copy constructor
MEos::MEos(const MEos& meos) : ndom(meos.ndom),
                               constructed_from_file(false) {

        mono_eos = new const Eos* [ndom] ;

        for (int l=0; l<ndom; l++) {
                mono_eos[l] =  meos.mono_eos[l] ;
        }

}


//  Constructor from a binary file
MEos::MEos(FILE* fich) : Eos(fich), constructed_from_file(true) {

    fread_be(&ndom, sizeof(int), 1, fich) ;

    mono_eos = new const Eos* [ndom] ;

    for (int l=0; l<ndom; l++) {
        mono_eos[l] = Eos::eos_from_file(fich) ;
    }
}

//  Constructor from a formatted  file
MEos::MEos(ifstream& fich) : Eos(fich),
                             constructed_from_file(true) {

    char blabla[80] ;

    fich >> ndom ; fich.getline(blabla, 80) ;

    mono_eos = new const Eos* [ndom] ;

    for (int l=0; l<ndom; l++) {
        mono_eos[l] = Eos::eos_from_file(fich) ;
    }
}



// Destructor
MEos::~MEos() {

        if (constructed_from_file) {
                for (int l=0; l<ndom; l++) {
                        delete mono_eos[l] ;
                }
        }

        delete [] mono_eos ;

}

			//--------------//
			//  Assignment  //
			//--------------//

void MEos::operator=(const MEos& ) {

        cout << "MEos::operator=  : not implemented yet !" << endl ;
                abort() ;

}


                     //---------------------------------------//
                        //              Outputs                  //
                        //---------------------------------------//

void MEos::sauve(FILE* fich) const {

    Eos::sauve(fich) ;

    fwrite_be(&ndom, sizeof(int), 1, fich) ;

    for (int l=0; l<ndom; l++) {
        mono_eos[l]->sauve(fich) ;
    }

}

ostream& MEos::operator>>(ostream & ost) const {

    ost << "EOS of class MEos (multi-domain equation of state) : " << endl ;
    ost << "   Number of domains :      " << ndom << endl ;

    for (int l=0; l<ndom; l++) {
        ost << "Equation of state in domain " << l << " : " << endl ;
        ost << "-------------------------------" << endl ; 
        ost << *(mono_eos[l]) ;
    }

    return ost ;

}

			//------------------------//
			//  Comparison operators  //
			//------------------------//


bool MEos::operator==(const Eos& eos_i) const {

    bool resu = true ;

    if ( eos_i.identify() != identify() ) {
	cout << "The second EOS is not of type MEos !" << endl ;
	resu = false ;
    }
    else{

	const MEos& eos = dynamic_cast<const MEos&>( eos_i ) ;

        if (eos.ndom != ndom) {
                cout <<  "The two MEos have different number of domains" << endl ;
                resu = false ;

        }
        else {
                for (int l=0; l<ndom; l++) {
                        resu = resu && ( *(mono_eos[l]) == *(eos.mono_eos[l]) )  ;
                }
        }
    }

    return resu ;

}

bool MEos::operator!=(const Eos& eos_i) const {

    return !(operator==(eos_i)) ;

}


			//------------------------------//
			//    Computational routines    //
			//------------------------------//

// Baryon density from enthalpy
//------------------------------

double MEos::nbar_ent_p(double ent, const Param* par) const {

        int l0 = par->get_int() ;        // index of the domain

        return mono_eos[l0]->nbar_ent_p(ent) ;

}

// Energy density from enthalpy
//------------------------------

double MEos::ener_ent_p(double ent, const Param* par) const {

        int l0 = par->get_int() ;        // index of the domain

        return mono_eos[l0]->ener_ent_p(ent) ;
}

// Pressure from enthalpy
//------------------------

double MEos::press_ent_p(double ent, const Param* par) const {

        int l0 = par->get_int() ;        // index of the domain

        return mono_eos[l0]->press_ent_p(ent) ;
}

// dln(n)/ln(H) from enthalpy
//---------------------------

double MEos::der_nbar_ent_p(double ent, const Param* par) const {

        int l0 = par->get_int() ;        // index of the domain

        return mono_eos[l0]->der_nbar_ent_p(ent) ;
}

// dln(e)/ln(H) from enthalpy
//---------------------------

double MEos::der_ener_ent_p(double ent, const Param* par) const {

        int l0 = par->get_int() ;        // index of the domain

        return mono_eos[l0]->der_ener_ent_p(ent) ;
}

// dln(p)/ln(H) from enthalpy
//---------------------------

double MEos::der_press_ent_p(double ent, const Param* par) const {

        int l0 = par->get_int() ;        // index of the domain

        return mono_eos[l0]->der_press_ent_p(ent) ;
}



