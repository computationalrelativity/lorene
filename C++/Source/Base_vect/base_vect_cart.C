/*
 *  Methods of class Base_vect_cart
 *
 *   (see file base_vect.h for documentation)
 *
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


char base_vect_cart_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2001/12/04 21:27:52  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.6  2000/09/12  12:15:15  eric
 * Fonction change_basis : ajout du cas spherique --> cartesiennes
 *  (appel des fonctions Map::comp_*_from_spherical).
 *
 * Revision 2.5  2000/02/10  17:04:13  eric
 * *** empty log message ***
 *
 * Revision 2.4  2000/02/10  15:38:51  eric
 * Premiere version de change_basis.
 *
 * Revision 2.3  2000/02/09  13:24:29  eric
 * REFONTE COMPLETE DE LA CLASSE
 * L'identification n'est plus base sur un membre statique (numero
 * d'instance) mais sur les caracteres physiques (rot_phi, etc...)
 * Ajout des constructeurs par copie et lecture de fichier.
 *
 * Revision 2.2  2000/01/10  15:43:29  eric
 * Methode change_basis.
 *
 * Revision 2.1  2000/01/10  13:27:32  eric
 * Ajout de la fonction set_rot_phi.
 *
 * Revision 2.0  2000/01/10  12:43:29  eric
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Headers C++
#include <iostream.h>

// Headers C
#include <math.h>
#include <stdlib.h>

// Headers Lorene
#include "base_vect.h"
#include "tenseur.h"
#include "utilitaires.h"

			//--------------//
			// Constructors //
			//--------------//

// Standard constructor without name
// ---------------------------------
Base_vect_cart::Base_vect_cart(double rot_phi_i) 
			      : rot_phi(rot_phi_i) {
    
    set_align() ;     
}



// Standard constructor with name
// ------------------------------
Base_vect_cart::Base_vect_cart(double rot_phi_i, const char* name_i)  
		    : Base_vect(name_i), 
		      rot_phi(rot_phi_i) {
    
    set_align() ;     
}


// Copy constructor
// ----------------
Base_vect_cart::Base_vect_cart(const Base_vect_cart& bi) 
			      : Base_vect(bi), 
			        rot_phi(bi.rot_phi) {
    
    set_align() ;     
}

// Constructor from file
// ---------------------
Base_vect_cart::Base_vect_cart(FILE* fich) 
			       : Base_vect(fich) {
        
    fread_be(&rot_phi, sizeof(double), 1, fich) ;		

    set_align() ;     
}


			//--------------//
			//  Destructor  //
			//--------------//

Base_vect_cart::~Base_vect_cart(){
    
    // does nothing
        
}

 
			//--------------//
			//   Mutators   //
			//--------------//

// Assignment 
//-----------
void Base_vect_cart::operator=(const Base_vect_cart& bi) {
    
    set_name(bi.name) ; 
    
    rot_phi = bi.rot_phi ;
    
    set_align() ;     
    
}



// Sets the value of the indicator align
//--------------------------------------
void Base_vect_cart::set_align() {

    if (rot_phi == double(0)) {
	    align = 1 ; 
    }
    else{
	if (rot_phi == M_PI) {
	    align = - 1 ; 
	}
	else{

	    if (fabs(rot_phi) < 1.e-14) {
		cout << 
		"WARNING : Base_vect_cart::set_align : " << endl ; 
		cout << "   rot_phi is close to zero !" << endl ;
		arrete() ; 
	    }
	    
	    if (fabs(rot_phi - M_PI) < 1.e-14) {
		cout << 
		"WARNING : Base_vect_cart::set_align : " << endl ; 
		cout << "   rot_phi is close to Pi !" << endl ;
		arrete() ; 
	    }
	    
	    align = 0 ; 
	    
	}

    }

}

void Base_vect_cart::set_rot_phi(double rot_phi_i) {
    
    rot_phi = rot_phi_i ;
     
    set_align() ;   
}

			//-----------------------//
			//  Comparison operator  //
			//-----------------------//


bool Base_vect_cart::operator==(const Base_vect& bi) const {
    
    bool resu = true ; 
    
    if ( bi.identify() != identify() ) {
	cout << "The second Base_vect is not of type Base_vect_cart !" << endl ; 
	resu = false ; 
    }
    else{
	
	const Base_vect_cart& bic = dynamic_cast<const Base_vect_cart&>( bi ) ; 

	if (bic.rot_phi != rot_phi) {
	    cout 
	    << "The two Base_vect_cart have different rot_phi : " << rot_phi 
		<< " <-> " << bic.rot_phi << endl ; 
	    resu = false ; 
	}

	
    }
    
    return resu ; 
    
}

			//------------//
			//  Outputs   //
			//------------//

void Base_vect_cart::sauve(FILE* fich) const {

    Base_vect::sauve(fich) ; 
    
    fwrite_be(&rot_phi, sizeof(double), 1, fich) ;	
   
}

ostream& Base_vect_cart::operator>>(ostream & ost) const {
    
    ost << "Azimuthal angle with respect to the Absolute frame : " 
	<< rot_phi << endl ; 
    
    return ost ;

}






		    //--------------------------------------//
		    //		Change of basis		    //
		    //--------------------------------------//

void Base_vect_cart::change_basis(Tenseur& ti) const {


    assert(ti.get_etat() != ETATNONDEF) ; 

    const Base_vect* triad_i = ti.get_triad() ; 

    assert(triad_i != 0x0) ; 
    
    if (ti.get_etat() == ETATZERO) {
	ti.set_triad(*this) ; 
	return ; 
    }
    
    assert(ti.get_etat() == ETATQCQ) ; 
    
    const Base_vect_cart* bvc = dynamic_cast<const Base_vect_cart*>(triad_i) ; 
    const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad_i) ; 
    
    // ---------------------------------------------
    // Case where the input triad is a Cartesian one
    // ---------------------------------------------
    if (bvc != 0x0) {	
	assert(bvs == 0x0) ; 
	
	int ind = align * (bvc->align) ; 
	
	switch (ind) {
	    
	    case 1 : {	// the two basis are aligned : nothing to do
			// -----------------------------------------
	    
		break ; 		
	    }
	    
	    case - 1 : {    // the two basis are anti-aligned 
			    // ------------------------------
	    
		switch (ti.get_valence()) {
		    
		    case 1 : {	// vector
			ti.set(0) = - ti(0) ;	 // V^x --> - V^x
			ti.set(1) = - ti(1) ; 	 // V^y --> - V^y
						 // V^z unchanged
			break ; 	
		    }
		    
		    case 2 : { // rank 2 tensor
			ti.set(0, 2) = - ti(0, 2) ; // {xz} --> - {xz}
			ti.set(1, 2) = - ti(1, 2) ; // {yz} --> - {yz}
			ti.set(2, 0) = - ti(2, 0) ; // {zx} --> - {zx}
			ti.set(2, 1) = - ti(2, 1) ; // {zy} --> - {zy}
					// all other components are unchanged
			break ; 
		    }
		    
		    default : {
			cout << 
			"Base_vect_cart::change_basis : the case of valence "
			<< ti.get_valence() << " is not treated !" << endl ;
			abort() ; 
			break ; 
		    }
		}
		break ; 		
	    }		    // end of case ind = -1
	    
	    case 0 : {	// the two basis have not a special relative orientation
			// -----------------------------------------------------
		cout << 
		"Base_vect_cart::change_basis : general value of rot_phi "
		<< " not contemplated yet, sorry !" << endl ;
		abort() ; 
		break ; 		
	    }
	    
	    default : {	    // error
		cout << 
		"Base_vect_cart::change_basis : unexpected value of ind !" << endl ;
		cout << "  ind = " << ind << endl ; 
		abort() ; 
		break ; 		
	    }
	}

    }	// end of the Cartesian basis case


    // ---------------------------------------------
    // Case where the input triad is a spherical one
    // ---------------------------------------------
    if (bvs != 0x0) {	

	assert(bvc == 0x0) ; 
	
	switch (ti.get_valence()) {
		    
	    case 1 : {	// vector
	    
		// The triads should be the same as that associated 
		// with the mapping :
		const Map* mp = ti.get_mp() ; 
		assert( *this == mp->get_bvect_cart() ) ; 
		assert( *bvs == mp->get_bvect_spher() ) ; 

		Cmp vr = ti(0) ; 
		Cmp vt = ti(1) ; 
		Cmp vp = ti(2) ; 
	    
		mp->comp_x_from_spherical(vr, vt, vp, ti.set(0)) ; 
		mp->comp_y_from_spherical(vr, vt, vp, ti.set(1)) ; 
		mp->comp_z_from_spherical(vr, vt, ti.set(2)) ; 
		
		break ; 	
	    }
		    
	    default : {
		cout << 
		"Base_vect_cart::change_basis : the case of valence "
		<< ti.get_valence() << " is not treated !" << endl ;
		abort() ; 
		break ; 
	    }
	}
	
	
    }	// end of the spherical basis case

    ti.set_triad(*this) ; 

} 
