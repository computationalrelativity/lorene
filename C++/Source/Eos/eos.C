/*
 * Methods of class Eos.
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
 * Revision 1.10  2021/05/14 15:39:22  g_servignat
 * Added sound speed computation from enthalpy to Eos class and tabulated+polytropic derived classes
 *
 * Revision 1.9  2017/12/15 15:36:38  j_novak
 * Improvement of the MEos class. Implementation of automatic offset computation accross different EoSs/domains.
 *
 * Revision 1.8  2016/12/05 16:17:50  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.7  2014/10/13 08:52:51  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.6  2014/10/06 15:13:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.5  2004/01/14 15:59:42  f_limousin
 * Add methos calcule, nbar_ent, der_bar_ent, der_press_ent and press_ent
 * for Scalar's.
 *
 * Revision 1.4  2002/10/16 14:36:34  j_novak
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
 * Revision 2.5  2001/10/10  13:45:06  eric
 * Modif Joachim : &(Eos::der_press_ent_p) -> &Eos::der_press_ent_p, etc...
 *  pour conformite au compilateur HP.
 *
 * Revision 2.4  2001/02/07  09:50:49  eric
 * Suppression de la fonction derent_ent_p.
 * Ajout des fonctions donnant les derivees de l'EOS:
 *      der_nbar_ent_p
 *      der_ener_ent_p
 *      der_press_ent_p
 * ainsi que des fonctions Cmp associees.
 *
 * Revision 2.3  2000/02/14  14:33:05  eric
 *  Ajout des constructeurs par lecture de fichier formate.
 *
 * Revision 2.2  2000/01/21  15:17:28  eric
 * fonction sauve: on ecrit en premier l'identificateur.
 *
 * Revision 2.1  2000/01/18  13:47:07  eric
 * Premiere version operationnelle
 *
 * Revision 2.0  2000/01/18  10:46:15  eric
 * /
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cstring>

// Headers Lorene
#include "eos.h"
#include "cmp.h"
#include "scalar.h"
#include "utilitaires.h"
#include "param.h"

			//--------------//
			// Constructors //
			//--------------//


// Standard constructor without name
// ---------------------------------
namespace Lorene {
Eos::Eos(){
        
    set_name("") ; 
    
}

// Standard constructor with name
// ---------------------------------
Eos::Eos(const char* name_i){
        
    set_name(name_i) ; 
    
}

// Copy constructor
// ----------------
Eos::Eos(const Eos& eos_i){
        
    set_name(eos_i.name) ; 
    
}

// Constructor from a binary file
// ------------------------------
Eos::Eos(FILE* fich){
        
    fread(name, sizeof(char), 100, fich) ;		
    
}

// Constructor from a formatted file
// ---------------------------------
Eos::Eos(ifstream& fich){
        
    fich.getline(name, 100) ;
    
}



			//--------------//
			//  Destructor  //
			//--------------//

Eos::~Eos(){
    
    // does nothing
        
}

			//-------------------------//
			//  Manipulation of name   //
			//-------------------------//
			
			
void Eos::set_name(const char* name_i) {

    strncpy(name, name_i,  100) ; 
    
}

const char* Eos::get_name() const {
    
    return name ; 
    
}

			//------------//
			//  Outputs   //
			//------------//

void Eos::sauve(FILE* fich) const {

    int ident = identify() ; 
    fwrite_be(&ident, sizeof(int), 1, fich) ;	
    	
    fwrite(name, sizeof(char), 100, fich) ;		
   
}
    



ostream& operator<<(ostream& ost, const Eos& eqetat)  {
    ost << eqetat.get_name() << endl ; 
    eqetat >> ost ;
    return ost ;
}


			//-------------------------------//
			// Generic computational routine //
			//-------------------------------//


void Eos::calcule(const Cmp& ent, int nzet, int l_min,  
		       double (Eos::*fait)(double, const Param*) const, Param* par, Cmp& resu) const {
    
    assert(ent.get_etat() != ETATNONDEF) ; 
    
    const Map* mp = ent.get_mp() ;	// Mapping
    
    
    if (ent.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
	return ; 
    }

    assert(ent.get_etat() == ETATQCQ) ;
    const MEos* tryMEos = dynamic_cast<const MEos*>(this) ;
    bool isMEos = (tryMEos != 0x0) ;
    int l_index = 0 ; // Index of the current domain in case of MEos
    if (isMEos) par->add_int_mod(l_index) ;
    
    const Valeur& vent = ent.va ;
    vent.coef_i() ;	// the values in the configuration space are required
   
    const Mg3d* mg = mp->get_mg() ;	// Multi-grid
    
    int nz = mg->get_nzone() ;		// total number of domains
    
    // Preparations for a point by point computation:
    resu.set_etat_qcq() ;
    Valeur& vresu = resu.va ; 
    vresu.set_etat_c_qcq() ;
    vresu.c->set_etat_qcq() ;

    // Loop on domains where the computation has to be done :
    for (int l = l_min; l< l_min + nzet; l++) {
	
	assert(l>=0) ; 
	assert(l<nz) ;

	l_index = l ; // The domain index is passed to the 'fait' function
	
	Tbl* tent = vent.c->t[l] ; 
	Tbl* tresu = vresu.c->t[l] ; 
	
        if (tent->get_etat() == ETATZERO) {
	    tresu->set_etat_zero() ; 
	}
	else {
	    assert( tent->get_etat() == ETATQCQ ) ; 
	    tresu->set_etat_qcq() ;

	    for (int i=0; i<tent->get_taille(); i++) {
		    
		tresu->t[i] = (this->*fait)( tent->t[i], par ) ;
	    }  
	    
	}  // End of the case where ent != 0 in the considered domain  
	
    }  // End of the loop on domains where the computation had to be done

    // resu is set to zero in the other domains :
    
    if (l_min > 0) {
	resu.annule(0, l_min-1) ; 
    }
    
    if (l_min + nzet < nz) {
	resu.annule(l_min + nzet, nz - 1) ; 
    }
}



void Eos::calcule(const Scalar& ent, int nzet, int l_min,  
		       double (Eos::*fait)(double, const Param*) const, Param* par, Scalar& resu) const {
    
    assert(ent.get_etat() != ETATNONDEF) ; 
    
    const Map* mp = &(ent.get_mp()) ;	// Mapping
    
    
    if (ent.get_etat() == ETATZERO) {
	resu.set_etat_zero() ; 
	return ; 
    }

    assert(ent.get_etat() == ETATQCQ) ;
    const MEos* tryMEos = dynamic_cast<const MEos*>(this) ;
    bool isMEos = (tryMEos != 0x0) ;
    int l_index = 0 ; // Index of the current domain in case of MEos
    if (isMEos) par->add_int_mod(l_index) ;

    const Valeur& vent = ent.get_spectral_va() ;
    vent.coef_i() ;	// the values in the configuration space are required
   
    const Mg3d* mg = mp->get_mg() ;	// Multi-grid
    
    int nz = mg->get_nzone() ;		// total number of domains
    
    // Preparations for a point by point computation:
    resu.set_etat_qcq() ;
    Valeur& vresu = resu.set_spectral_va() ; 
    vresu.set_etat_c_qcq() ;
    vresu.c->set_etat_qcq() ;

    // Loop on domains where the computation has to be done :
    for (int l = l_min; l< l_min + nzet; l++) {
	
	assert(l>=0) ; 
	assert(l<nz) ;

	l_index = l ;  // The domain index is passed to the 'fait' function
	
	Tbl* tent = vent.c->t[l] ; 
	Tbl* tresu = vresu.c->t[l] ; 
	
        if (tent->get_etat() == ETATZERO) {
	    tresu->set_etat_zero() ; 
	}
	else {
	    assert( tent->get_etat() == ETATQCQ ) ; 
	    tresu->set_etat_qcq() ;

	    for (int i=0; i<tent->get_taille(); i++) {
		    
		tresu->t[i] = (this->*fait)( tent->t[i], par ) ;
	    }  
	    
	}  // End of the case where ent != 0 in the considered domain  
	
    }  // End of the loop on domains where the computation had to be done

    // resu is set to zero in the other domains :
    
    if (l_min > 0) {
	resu.annule(0, l_min-1) ; 
    }
    
    if (l_min + nzet < nz) {
	resu.annule(l_min + nzet, nz - 1) ; 
    }
}
			//---------------------------------//
			//	Public  functions	   //
			//---------------------------------//
			

// Baryon density from enthalpy 
//------------------------------

Cmp Eos::nbar_ent(const Cmp& ent, int nzet, int l_min, Param* par) const {
    
    Cmp resu(ent.get_mp()) ; 
    
    calcule(ent, nzet, l_min, &Eos::nbar_ent_p, par, resu) ;
    
    return resu ; 
    
}

Scalar Eos::nbar_ent(const Scalar& ent, int nzet, int l_min, Param* par) const {
    
    Scalar resu(ent.get_mp()) ; 
    
    calcule(ent, nzet, l_min, &Eos::nbar_ent_p, par, resu) ;
    
    return resu ; 
    
}



// Energy density from enthalpy 
//------------------------------

Cmp Eos::ener_ent(const Cmp& ent, int nzet, int l_min, Param* par) const {
    
    Cmp resu(ent.get_mp()) ; 
    
    calcule(ent, nzet, l_min, &Eos::ener_ent_p, par, resu) ;
    
    return resu ; 
    
}

Scalar Eos::ener_ent(const Scalar& ent, int nzet, int l_min, Param* par) const {
    
    Scalar resu(ent.get_mp()) ; 
    
    calcule(ent, nzet, l_min, &Eos::ener_ent_p, par, resu) ;
    
    return resu ; 
    
}
// Pressure from enthalpy 
//-----------------------

Cmp Eos::press_ent(const Cmp& ent, int nzet, int l_min, Param* par) const {
    
    Cmp resu(ent.get_mp()) ; 
    
    calcule(ent, nzet, l_min, &Eos::press_ent_p, par, resu) ;

    return resu ;

}

Scalar Eos::press_ent(const Scalar& ent, int nzet, int l_min, Param* par) const {
    
    Scalar resu(ent.get_mp()) ; 
    
    calcule(ent, nzet, l_min, &Eos::press_ent_p, par, resu) ;

    return resu ;

}
// Derivative of baryon density from enthalpy
//-------------------------------------------

Cmp Eos::der_nbar_ent(const Cmp& ent, int nzet, int l_min, Param* par) const {

    Cmp resu(ent.get_mp()) ;

    calcule(ent, nzet, l_min, &Eos::der_nbar_ent_p, par, resu) ;

    return resu ;

}

Scalar Eos::der_nbar_ent(const Scalar& ent, int nzet, int l_min, Param* par) const {

    Scalar resu(ent.get_mp()) ;

    calcule(ent, nzet, l_min, &Eos::der_nbar_ent_p, par, resu) ;

    return resu ;

}

// Derivative of energy density from enthalpy
//-------------------------------------------

Cmp Eos::der_ener_ent(const Cmp& ent, int nzet, int l_min, Param* par) const {

    Cmp resu(ent.get_mp()) ;

    calcule(ent, nzet, l_min, &Eos::der_ener_ent_p, par, resu) ;

    return resu ;

}

Scalar Eos::der_ener_ent(const Scalar& ent, int nzet, int l_min, Param* par) const {

    Scalar resu(ent.get_mp()) ;

    calcule(ent, nzet, l_min, &Eos::der_ener_ent_p, par, resu) ;

    return resu ;

}
// Derivative of pressure from enthalpy
//-------------------------------------------

Cmp Eos::der_press_ent(const Cmp& ent, int nzet, int l_min, Param* par) const {

    Cmp resu(ent.get_mp()) ;

    calcule(ent, nzet, l_min, &Eos::der_press_ent_p, par, resu) ;

    return resu ;

}

Scalar Eos::der_press_ent(const Scalar& ent, int nzet, int l_min, Param* par) const {

    Scalar resu(ent.get_mp()) ;

    calcule(ent, nzet, l_min, &Eos::der_press_ent_p, par, resu) ;

    return resu ;

}

// Sound speed from enthalpy
//---------------------------------

Scalar Eos::csound_square_ent(const Scalar& ent, int nzet,
				      int l_min, Param* par) const {
    Scalar resu(ent.get_mp()) ;
    
    calcule(ent, nzet, l_min, &Eos::csound_square_ent_p, par, resu) ;
    
    return resu ;
  }

}


