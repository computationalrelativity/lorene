/*
 *  Methods of class Base_vect_spher
 *
 *   (see file base_vect.h for documentation)
 *
 */

/*
 *   Copyright (c) 2000-2002 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Philippe Grandclement
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

char base_vect_spher_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/01/15 09:09:49  e_gourgoulhon
 * Suppression of cout printing in the comparison operator
 *
 * Revision 1.2  2001/12/04 21:27:52  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.5  2000/09/19  16:09:40  phil
 * *** empty log message ***
 *
 * Revision 2.4  2000/09/19  15:33:23  phil
 * ajout passage de cartesinne en spherique
 *
 * Revision 2.3  2000/02/09  13:24:47  eric
 * REFONTE COMPLETE DE LA CLASSE
 * L'identification n'est plus base sur un membre statique (numero
 * d'instance) mais sur les caracteres physiques (rot_phi, etc...)
 * Ajout des membres ori_x, ori_y, ori_z
 * Ajout des constructeurs par copie et lecture de fichier.
 *
 * Revision 2.2  2000/01/10  15:43:57  eric
 * Methode change_basis.
 *
 * Revision 2.1  2000/01/10  13:27:09  eric
 * Ajout de la fonction set_rot_phi.
 *
 * Revision 2.0  2000/01/10  12:43:33  eric
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
Base_vect_spher::Base_vect_spher(double xa0, double ya0, double za0,
				 double rot_phi_i) 
				: ori_x(xa0), 
				  ori_y(ya0), 
				  ori_z(za0), 
				  rot_phi(rot_phi_i) {} 





// Standard constructor with name
// ------------------------------
Base_vect_spher::Base_vect_spher(double xa0, double ya0, double za0,
				 double rot_phi_i, const char* name_i)  
				: Base_vect(name_i), 
				  ori_x(xa0), 
				  ori_y(ya0), 
				  ori_z(za0), 
				  rot_phi(rot_phi_i) {} 
				
// Copy constructor
// ----------------
Base_vect_spher::Base_vect_spher(const Base_vect_spher& bi)  
				: Base_vect(bi), 
				  ori_x(bi.ori_x), 
				  ori_y(bi.ori_y), 
				  ori_z(bi.ori_z), 
				  rot_phi(bi.rot_phi) {} 

// Constructor from file
// ---------------------
Base_vect_spher::Base_vect_spher(FILE* fich) 
			       : Base_vect(fich) {
        
    fread_be(&ori_x, sizeof(double), 1, fich) ;		
    fread_be(&ori_y, sizeof(double), 1, fich) ;		
    fread_be(&ori_z, sizeof(double), 1, fich) ;		
    fread_be(&rot_phi, sizeof(double), 1, fich) ;		

}


			//--------------//
			//  Destructor  //
			//--------------//

Base_vect_spher::~Base_vect_spher(){
    
    // does nothing
        
}

			//--------------//
			//   Mutators   //
			//--------------//

// Assignment 
//-----------
void Base_vect_spher::operator=(const Base_vect_spher& bi) {
    
    set_name(bi.name) ; 
    
    ori_x = bi.ori_x ; 
    ori_y = bi.ori_y ; 
    ori_z = bi.ori_z ; 
    rot_phi = bi.rot_phi ;
    
}

// Change of the elements
// ----------------------
void Base_vect_spher::set_ori(double xa0, double ya0, double za0) {

    ori_x = xa0 ; 
    ori_y = ya0 ; 
    ori_z = za0 ; 

}


void Base_vect_spher::set_rot_phi(double rot_phi_i) {
    
    rot_phi = rot_phi_i ;
     
}

			//-----------------------//
			//  Comparison operator  //
			//-----------------------//


bool Base_vect_spher::operator==(const Base_vect& bi) const {
    
    bool resu = true ; 
    
    if ( bi.identify() != identify() ) {
	// cout << "The second Base_vect is not of type Base_vect_spher !" << endl ;
	resu = false ; 
    }
    else{
	
	const Base_vect_spher& bis = dynamic_cast<const Base_vect_spher&>( bi ) ; 

	if (bis.ori_x != ori_x) {
	    cout 
	    << "The two Base_vect_spher have different X origin : " << ori_x 
		<< " <-> " << bis.ori_x << endl ; 
	    resu = false ; 
	}

	if (bis.ori_y != ori_y) {
	    cout 
	    << "The two Base_vect_spher have different Y origin : " << ori_y 
		<< " <-> " << bis.ori_y << endl ; 
	    resu = false ; 
	}

	if (bis.ori_z != ori_z) {
	    cout 
	    << "The two Base_vect_spher have different Z origin : " << ori_z 
		<< " <-> " << bis.ori_z << endl ; 
	    resu = false ; 
	}

	if (bis.rot_phi != rot_phi) {
	    cout 
	    << "The two Base_vect_spher have different rot_phi : " << rot_phi 
		<< " <-> " << bis.rot_phi << endl ; 
	    resu = false ; 
	}

	
    }
    
    return resu ; 
    
}

			//------------//
			//  Outputs   //
			//------------//

void Base_vect_spher::sauve(FILE* fich) const {

    Base_vect::sauve(fich) ; 
    
    fwrite_be(&ori_x, sizeof(double), 1, fich) ;	
    fwrite_be(&ori_y, sizeof(double), 1, fich) ;	
    fwrite_be(&ori_z, sizeof(double), 1, fich) ;	
    fwrite_be(&rot_phi, sizeof(double), 1, fich) ;	
   
}

ostream& Base_vect_spher::operator>>(ostream & ost) const {
    
    ost << "Absolute coordinates (X,Y,Z) of the origin : " 
	<< ori_x << "  " << ori_y << "  " << ori_z << endl ; 
    ost << "Azimuthal angle with respect to the Absolute frame : " 
	<< rot_phi << endl ; 
    
    return ost ;

}





		    //--------------------------------------//
		    //		Change of basis		    //
		    //--------------------------------------//

void Base_vect_spher::change_basis(Tenseur& ti) const {
    

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
		
	switch (ti.get_valence()) {
		    
	    case 1 : {	// vector
	    
		// The triads should be the same as that associated 
		// with the mapping :
		const Map* mp = ti.get_mp() ; 
		assert( *this == mp->get_bvect_spher() ) ; 
		assert( *bvc == mp->get_bvect_cart() ) ; 

		Cmp vx = ti(0) ; 
		Cmp vy = ti(1) ; 
		Cmp vz = ti(2) ; 
	    
		mp->comp_r_from_cartesian(vx, vy, vz, ti.set(0)) ; 
		mp->comp_t_from_cartesian(vx, vy, vz, ti.set(1)) ; 
		mp->comp_p_from_cartesian(vx, vy, ti.set(2)) ; 
		
		break ; 	
	    }
		    
	    default : {
		cout << 
		"Base_vect_sphere::change_basis : the case of valence "
		<< ti.get_valence() << " is not treated !" << endl ;
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
	
	cout << "Base_vect_spher::change_basis : case not treated yet !" << endl ;
	abort() ; 
    }	// end of the spherical basis case

    ti.set_triad(*this) ; 
}
