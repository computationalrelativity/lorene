/*
 *  Methods of class Param
 *
 *   (see file param.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 2000-2001 Jerome Novak
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


char param_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2002/10/16 14:37:12  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2002/09/19 09:52:42  j_novak
 * Added objects Qtenseur and Qmetrique for 4D tensor and metric handling.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.8  2001/10/11  07:44:27  eric
 * Ajout du stokage des Etoile's
 *
 * Revision 1.7  2000/10/24  14:55:20  novak
 * Added the function clean_all()
 *
 * Revision 1.6  2000/05/25 12:40:47  eric
 *  MODIFICATION MAJEURE: pour les int et les double, ce sont desormais les
 * dresses qui sont stokees, et non plus les nombres eux-memes
 * (le traitement des int et des double est donc desormais completement
 * aligne sur celui des Tbl, Cmp, etc...)
 *
 * Revision 1.5  1999/12/29  13:10:54  eric
 *  Ajout du stokage des Mtbl_cf.
 *
 * Revision 1.4  1999/12/27  12:17:02  eric
 * Ajout du stokage des mappings (class Map).
 *
 * Revision 1.3  1999/12/16  10:28:25  eric
 * Ajout des membres modifiables.
 * Par defaut, les objets listes sont const.
 *
 * Revision 1.2  1999/12/15  16:23:22  eric
 * Changement de l'ordre des arguments dans add_*
 * Argument par defaut: position = 0
 * Ajout du stokage des int et des double.
 *
 * Revision 1.1  1999/12/13  14:36:00  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <assert.h>

// Headers Lorene
#include "param.h"
#include "tbl.h"
#include "itbl.h"
#include "cmp.h"
#include "qtenseur.h"

			//------------------------//
			//	Constructor	  //
			//------------------------//

Param::Param() : n_int(0), 
                 n_int_mod(0),
		 n_double(0), 
		 n_double_mod(0), 
		 n_tbl(0), 
		 n_tbl_mod(0), 
		 n_itbl(0), 
		 n_itbl_mod(0), 
		 n_cmp(0), 
		 n_cmp_mod(0), 
		 n_tenseur(0),
		 n_tenseur_mod(0), 
		 n_qtenseur(0),
		 n_qtenseur_mod(0), 
		 n_map(0), 
		 n_mtbl_cf(0), 
		 n_etoile(0)
		 {}


			//----------------------//
			//	Destructor	//
			//----------------------//

Param::~Param(){

    if (n_int > 0)  delete [] p_int ; 
    if (n_int_mod > 0)  delete [] p_int_mod ; 
    if (n_double > 0)  delete [] p_double ; 
    if (n_double_mod > 0)  delete [] p_double_mod ; 
    if (n_tbl > 0)  delete [] p_tbl ; 
    if (n_tbl_mod > 0)  delete [] p_tbl_mod ; 
    if (n_itbl > 0) delete [] p_itbl ; 
    if (n_itbl_mod > 0) delete [] p_itbl_mod ; 
    if (n_cmp > 0)  delete [] p_cmp ; 
    if (n_cmp_mod > 0)  delete [] p_cmp_mod ; 
    if (n_tenseur > 0) delete [] p_tenseur ; 
    if (n_tenseur_mod > 0) delete [] p_tenseur_mod ; 
    if (n_qtenseur > 0) delete [] p_qtenseur ; 
    if (n_qtenseur_mod > 0) delete [] p_qtenseur_mod ; 
    if (n_map > 0)  delete [] p_map ; 
    if (n_mtbl_cf > 0)  delete [] p_mtbl_cf ; 
    if (n_etoile > 0)  delete [] p_etoile ; 

}
 
		    //------------------------------------//
		    //	      cleaning the memory 	  //
		    //------------------------------------//

void Param::clean_all() {

  for (int i=0; i<n_int_mod; i++) 
    if (p_int_mod[i] != 0x0) {
      delete p_int_mod[i] ;
      p_int_mod[i] = 0x0 ;
    }

  for (int i=0; i<n_double_mod; i++) 
    if (p_double_mod[i] != 0x0) {
      delete p_double_mod[i] ;
      p_double_mod[i] = 0x0 ;
    }

  for (int i=0; i<n_tbl_mod; i++) 
    if (p_tbl_mod[i] != 0x0) { 
      delete p_tbl_mod[i] ;
      p_tbl_mod[i] = 0x0 ;
    }

  for (int i=0; i<n_itbl_mod; i++) 
    if (p_itbl_mod[i] != 0x0) {
      delete p_itbl_mod[i] ;
      p_itbl_mod[i] = 0x0 ;
    }

  for (int i=0; i<n_cmp_mod; i++) 
    if (p_cmp_mod[i] != 0x0) {
      delete p_cmp_mod[i] ;
      p_cmp_mod[i] = 0x0 ;
    }

  for (int i=0; i<n_tenseur_mod; i++) 
    if (p_tenseur_mod[i] != 0x0) {
      delete p_tenseur_mod[i] ;
      p_tenseur_mod[i] = 0x0 ;
    }

  for (int i=0; i<n_qtenseur_mod; i++) 
    if (p_qtenseur_mod[i] != 0x0) {
      delete p_qtenseur_mod[i] ;
      p_qtenseur_mod[i] = 0x0 ;
    }
}


		    //------------------------------------//
		    //		int storage		  //
		    //------------------------------------//

// Total number of stored int
// --------------------------

int Param::get_n_int() const {
    return n_int ; 
}

// Addition  
// --------
		    
void Param::add_int(const int& ti, int index){
    
	if (index >= n_int) {    // p_int must be rescaled
	    	    
	    int n_int_nouveau = index + 1 ; 
	    const int** p_int_nouveau = new const int*[n_int_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_int; i++) {
		p_int_nouveau[i] = p_int[i] ; 
	    }
	    
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_int; i<index; i++) {
		p_int_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_int_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_int > 0) delete [] p_int ; 
	    p_int = p_int_nouveau ; 
	    n_int = n_int_nouveau ; 
	    
	}
	else {
	
	    if (p_int[index] != 0x0) {
		cout << "Param::add_int : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_int[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const int& Param::get_int(int index) const {

    assert(index >= 0) ;
    assert(index < n_int) ; 
    
    return *(p_int[index]) ; 

} 
		    
		    //------------------------------------//
		    //		double storage		  //
		    //------------------------------------//

// Total number of stored doubles
// ------------------------------

int Param::get_n_double() const {
    return n_double ; 
}

// Addition  
// --------
		    
void Param::add_double(const double& ti, int index){
    
	if (index >= n_double) {    // p_double must be rescaled
	    	    
	    int n_double_nouveau = index + 1 ; 
	    const double** p_double_nouveau = 
				    new const double*[n_double_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_double; i++) {
		p_double_nouveau[i] = p_double[i] ; 
	    }
	    	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_double; i<index; i++) {
		p_double_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_double_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_double > 0) delete [] p_double ; 
	    p_double = p_double_nouveau ; 
	    n_double = n_double_nouveau ; 
	    
	}
	else {
	
	    if (p_double[index] != 0x0) {
		cout << "Param::add_double : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_double[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const double& Param::get_double(int index) const {

    assert(index >= 0) ;
    assert(index < n_double) ; 
    
    return *(p_double[index]) ; 

} 
		    

		    //------------------------------------//
		    //	     Modifiable int storage	  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_int_mod() const {
    return n_int_mod ; 
}

// Addition  
// --------
		    
void Param::add_int_mod(int& ti, int index){
    
	if (index >= n_int_mod) {    // p_int_mod must be rescaled
	    	    
	    int n_int_nouveau = index + 1 ; 
	    int** p_int_nouveau = new int*[n_int_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_int_mod; i++) {
		p_int_nouveau[i] = p_int_mod[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_int_mod; i<index; i++) {
		p_int_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_int_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_int_mod > 0) delete [] p_int_mod ; 
	    p_int_mod = p_int_nouveau ; 
	    n_int_mod = n_int_nouveau ; 
	    
	}
	else {
	
	    if (p_int_mod[index] != 0x0) {
		cout << "Param::add_int_mod : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_int_mod[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
int& Param::get_int_mod(int index) const {

    assert(index >= 0) ;
    assert(index < n_int_mod) ; 
    
    return *(p_int_mod[index]) ; 

} 
		    
		    //------------------------------------//
		    //	 Modifiable double storage	  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_double_mod() const {
    return n_double_mod ; 
}

// Addition  
// --------
		    
void Param::add_double_mod(double& ti, int index){
    
	if (index >= n_double_mod) {    // p_double_mod must be rescaled
	    	    
	    int n_double_nouveau = index + 1 ; 
	    double** p_double_nouveau = new double*[n_double_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_double_mod; i++) {
		p_double_nouveau[i] = p_double_mod[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_double_mod; i<index; i++) {
		p_double_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_double_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_double_mod > 0) delete [] p_double_mod ; 
	    p_double_mod = p_double_nouveau ; 
	    n_double_mod = n_double_nouveau ; 
	    
	}
	else {
	
	    if (p_double_mod[index] != 0x0) {
		cout << "Param::add_double_mod : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_double_mod[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
double& Param::get_double_mod(int index) const {

    assert(index >= 0) ;
    assert(index < n_double_mod) ; 
    
    return *(p_double_mod[index]) ; 

} 
		    

		    //------------------------------------//
		    //		Tbl storage		  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_tbl() const {
    return n_tbl ; 
}

// Addition  
// --------
		    
void Param::add_tbl(const Tbl& ti, int index){
    
	if (index >= n_tbl) {    // p_tbl must be rescaled
	    	    
	    int n_tbl_nouveau = index + 1 ; 
	    const Tbl** p_tbl_nouveau = new const Tbl*[n_tbl_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_tbl; i++) {
		p_tbl_nouveau[i] = p_tbl[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_tbl; i<index; i++) {
		p_tbl_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_tbl_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_tbl > 0) delete [] p_tbl ; 
	    p_tbl = p_tbl_nouveau ; 
	    n_tbl = n_tbl_nouveau ; 
	    
	}
	else {
	
	    if (p_tbl[index] != 0x0) {
		cout << "Param::add_tbl : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_tbl[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const Tbl& Param::get_tbl(int index) const {

    assert(index >= 0) ;
    assert(index < n_tbl) ; 
    
    return *(p_tbl[index]) ; 

} 
		    

		    //------------------------------------//
		    //	    Modifiable Tbl storage	  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_tbl_mod() const {
    return n_tbl_mod ; 
}

// Addition  
// --------
		    
void Param::add_tbl_mod(Tbl& ti, int index){
    
	if (index >= n_tbl_mod) {    // p_tbl_mod must be rescaled
	    	    
	    int n_tbl_nouveau = index + 1 ; 
	    Tbl** p_tbl_nouveau = new Tbl*[n_tbl_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_tbl_mod; i++) {
		p_tbl_nouveau[i] = p_tbl_mod[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_tbl_mod; i<index; i++) {
		p_tbl_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_tbl_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_tbl_mod > 0) delete [] p_tbl_mod ; 
	    p_tbl_mod = p_tbl_nouveau ; 
	    n_tbl_mod = n_tbl_nouveau ; 
	    
	}
	else {
	
	    if (p_tbl_mod[index] != 0x0) {
		cout << "Param::add_tbl_mod : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_tbl_mod[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
Tbl& Param::get_tbl_mod(int index) const {

    assert(index >= 0) ;
    assert(index < n_tbl_mod) ; 
    
    return *(p_tbl_mod[index]) ; 

} 

		    
		    //------------------------------------//
		    //		Itbl storage		  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_itbl() const {
    return n_itbl ; 
}

// Addition  
// --------
		    
void Param::add_itbl(const Itbl& ti, int index){
    
	if (index >= n_itbl) {    // p_itbl must be rescaled
	    	    
	    int n_itbl_nouveau = index + 1 ; 
	    const Itbl** p_itbl_nouveau = new const Itbl*[n_itbl_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_itbl; i++) {
		p_itbl_nouveau[i] = p_itbl[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_itbl; i<index; i++) {
		p_itbl_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_itbl_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_itbl > 0) delete [] p_itbl ; 
	    p_itbl = p_itbl_nouveau ; 
	    n_itbl = n_itbl_nouveau ; 
	    
	}
	else {
	
	    if (p_itbl[index] != 0x0) {
		cout << "Param::add_itbl : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_itbl[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const Itbl& Param::get_itbl(int index) const {

    assert(index >= 0) ;
    assert(index < n_itbl) ; 
    
    return *(p_itbl[index]) ; 

} 
		    

		    //------------------------------------//
		    //	    Modifiable Itbl storage	  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_itbl_mod() const {
    return n_itbl_mod ; 
}

// Addition  
// --------
		    
void Param::add_itbl_mod(Itbl& ti, int index){
    
	if (index >= n_itbl_mod) {    // p_itbl_mod must be rescaled
	    	    
	    int n_itbl_nouveau = index + 1 ; 
	    Itbl** p_itbl_nouveau = new Itbl*[n_itbl_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_itbl_mod; i++) {
		p_itbl_nouveau[i] = p_itbl_mod[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_itbl_mod; i<index; i++) {
		p_itbl_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_itbl_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_itbl_mod > 0) delete [] p_itbl_mod ; 
	    p_itbl_mod = p_itbl_nouveau ; 
	    n_itbl_mod = n_itbl_nouveau ; 
	    
	}
	else {
	
	    if (p_itbl_mod[index] != 0x0) {
		cout << "Param::add_itbl_mod : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_itbl_mod[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
Itbl& Param::get_itbl_mod(int index) const {

    assert(index >= 0) ;
    assert(index < n_itbl_mod) ; 
    
    return *(p_itbl_mod[index]) ; 

} 
		    
		    
		    //------------------------------------//
		    //		Cmp storage		  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_cmp() const {
    return n_cmp ; 
}

// Addition  
// --------
		    
void Param::add_cmp(const Cmp& ti, int index){
    
	if (index >= n_cmp) {    // p_cmp must be rescaled
	    	    
	    int n_cmp_nouveau = index + 1 ; 
	    const Cmp** p_cmp_nouveau = new const Cmp*[n_cmp_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_cmp; i++) {
		p_cmp_nouveau[i] = p_cmp[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_cmp; i<index; i++) {
		p_cmp_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_cmp_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_cmp > 0) delete [] p_cmp ; 
	    p_cmp = p_cmp_nouveau ; 
	    n_cmp = n_cmp_nouveau ; 
	    
	}
	else {
	
	    if (p_cmp[index] != 0x0) {
		cout << "Param::add_cmp : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_cmp[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const Cmp& Param::get_cmp(int index) const {

    assert(index >= 0) ;
    assert(index < n_cmp) ; 
    
    return *(p_cmp[index]) ; 

} 
		    

		    //------------------------------------//
		    //	    Modifiable Cmp storage	  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_cmp_mod() const {
    return n_cmp_mod ; 
}

// Addition  
// --------
		    
void Param::add_cmp_mod(Cmp& ti, int index){
    
	if (index >= n_cmp_mod) {    // p_cmp_mod must be rescaled
	    	    
	    int n_cmp_nouveau = index + 1 ; 
	    Cmp** p_cmp_nouveau = new Cmp*[n_cmp_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_cmp_mod; i++) {
		p_cmp_nouveau[i] = p_cmp_mod[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_cmp_mod; i<index; i++) {
		p_cmp_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_cmp_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_cmp_mod > 0) delete [] p_cmp_mod ; 
	    p_cmp_mod = p_cmp_nouveau ; 
	    n_cmp_mod = n_cmp_nouveau ; 
	    
	}
	else {
	
	    if (p_cmp_mod[index] != 0x0) {
		cout << "Param::add_cmp_mod : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_cmp_mod[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
Cmp& Param::get_cmp_mod(int index) const {

    assert(index >= 0) ;
    assert(index < n_cmp_mod) ; 
    
    return *(p_cmp_mod[index]) ; 

} 

		    
		    //------------------------------------//
		    //		Tenseur storage		  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_tenseur() const {
    return n_tenseur ; 
}

// Addition  
// --------
		    
void Param::add_tenseur(const Tenseur& ti, int index){
    
	if (index >= n_tenseur) {    // p_tenseur must be rescaled
	    	    
	    int n_tenseur_nouveau = index + 1 ; 
	    const Tenseur** p_tenseur_nouveau = new const Tenseur*[n_tenseur_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_tenseur; i++) {
		p_tenseur_nouveau[i] = p_tenseur[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_tenseur; i<index; i++) {
		p_tenseur_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_tenseur_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_tenseur > 0) delete [] p_tenseur ; 
	    p_tenseur = p_tenseur_nouveau ; 
	    n_tenseur = n_tenseur_nouveau ; 
	    
	}
	else {
	
	    if (p_tenseur[index] != 0x0) {
		cout << "Param::add_tenseur : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_tenseur[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const Tenseur& Param::get_tenseur(int index) const {

    assert(index >= 0) ;
    assert(index < n_tenseur) ; 
    
    return *(p_tenseur[index]) ; 

} 
		    

		    //------------------------------------//
		    //	    Modifiable Tenseur storage	  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_tenseur_mod() const {
    return n_tenseur_mod ; 
}

// Addition  
// --------
		    
void Param::add_tenseur_mod(Tenseur& ti, int index){
    
	if (index >= n_tenseur_mod) {    // p_tenseur_mod must be rescaled
	    	    
	    int n_tenseur_nouveau = index + 1 ; 
	    Tenseur** p_tenseur_nouveau = new Tenseur*[n_tenseur_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_tenseur_mod; i++) {
		p_tenseur_nouveau[i] = p_tenseur_mod[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_tenseur_mod; i<index; i++) {
		p_tenseur_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_tenseur_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_tenseur_mod > 0) delete [] p_tenseur_mod ; 
	    p_tenseur_mod = p_tenseur_nouveau ; 
	    n_tenseur_mod = n_tenseur_nouveau ; 
	    
	}
	else {
	
	    if (p_tenseur_mod[index] != 0x0) {
		cout << "Param::add_tenseur_mod : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_tenseur_mod[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
Tenseur& Param::get_tenseur_mod(int index) const {

    assert(index >= 0) ;
    assert(index < n_tenseur_mod) ; 
    
    return *(p_tenseur_mod[index]) ; 

} 

		     
		    //------------------------------------//
		    //		Qtenseur storage		  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_qtenseur() const {
    return n_qtenseur ; 
}

// Addition  
// --------
		    
void Param::add_qtenseur(const Qtenseur& ti, int index){
    
	if (index >= n_qtenseur) {    // p_qtenseur must be rescaled
	    	    
	    int n_qtenseur_nouveau = index + 1 ; 
	    const Qtenseur** p_qtenseur_nouveau = new const Qtenseur*[n_qtenseur_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_qtenseur; i++) {
		p_qtenseur_nouveau[i] = p_qtenseur[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_qtenseur; i<index; i++) {
		p_qtenseur_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_qtenseur_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_qtenseur > 0) delete [] p_qtenseur ; 
	    p_qtenseur = p_qtenseur_nouveau ; 
	    n_qtenseur = n_qtenseur_nouveau ; 
	    
	}
	else {
	
	    if (p_qtenseur[index] != 0x0) {
		cout << "Param::add_qtenseur : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_qtenseur[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const Qtenseur& Param::get_qtenseur(int index) const {

    assert(index >= 0) ;
    assert(index < n_qtenseur) ; 
    
    return *(p_qtenseur[index]) ; 

} 
		    

		    //------------------------------------//
		    //	    Modifiable Qtenseur storage	  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_qtenseur_mod() const {
    return n_qtenseur_mod ; 
}

// Addition  
// --------
		    
void Param::add_qtenseur_mod(Qtenseur& ti, int index){
    
	if (index >= n_qtenseur_mod) {    // p_qtenseur_mod must be rescaled
	    	    
	    int n_qtenseur_nouveau = index + 1 ; 
	    Qtenseur** p_qtenseur_nouveau = new Qtenseur*[n_qtenseur_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_qtenseur_mod; i++) {
		p_qtenseur_nouveau[i] = p_qtenseur_mod[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_qtenseur_mod; i<index; i++) {
		p_qtenseur_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_qtenseur_nouveau[index] = &ti ; 
	    
	    // Update 
	    if (n_qtenseur_mod > 0) delete [] p_qtenseur_mod ; 
	    p_qtenseur_mod = p_qtenseur_nouveau ; 
	    n_qtenseur_mod = n_qtenseur_nouveau ; 
	    
	}
	else {
	
	    if (p_qtenseur_mod[index] != 0x0) {
		cout << "Param::add_qtenseur_mod : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_qtenseur_mod[index] = &ti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
Qtenseur& Param::get_qtenseur_mod(int index) const {

    assert(index >= 0) ;
    assert(index < n_qtenseur_mod) ; 
    
    return *(p_qtenseur_mod[index]) ; 

} 


		    //------------------------------------//
		    //		Map storage		  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_map() const {
    return n_map ; 
}

// Addition  
// --------
		    
void Param::add_map(const Map& mi, int index){
    
	if (index >= n_map) {    // p_map must be rescaled
	    	    
	    int n_map_nouveau = index + 1 ; 
	    const Map** p_map_nouveau = new const Map*[n_map_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_map; i++) {
		p_map_nouveau[i] = p_map[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_map; i<index; i++) {
		p_map_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_map_nouveau[index] = &mi ; 
	    
	    // Update 
	    if (n_map > 0) delete [] p_map ; 
	    p_map = p_map_nouveau ; 
	    n_map = n_map_nouveau ; 
	    
	}
	else {
	
	    if (p_map[index] != 0x0) {
		cout << "Param::add_map : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_map[index] = &mi ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const Map& Param::get_map(int index) const {

    assert(index >= 0) ;
    assert(index < n_map) ; 
    
    return *(p_map[index]) ; 

} 
		    
		    //------------------------------------//
		    //		Mtbl_cf storage		  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_mtbl_cf() const {
    return n_mtbl_cf ; 
}

// Addition  
// --------
		    
void Param::add_mtbl_cf(const Mtbl_cf& mi, int index){
    
	if (index >= n_mtbl_cf) {    // p_mtbl_cf must be rescaled
	    	    
	    int n_mtbl_cf_nouveau = index + 1 ; 
	    const Mtbl_cf** p_mtbl_cf_nouveau = 
				    new const Mtbl_cf*[n_mtbl_cf_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_mtbl_cf; i++) {
		p_mtbl_cf_nouveau[i] = p_mtbl_cf[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_mtbl_cf; i<index; i++) {
		p_mtbl_cf_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_mtbl_cf_nouveau[index] = &mi ; 
	    
	    // Update 
	    if (n_mtbl_cf > 0) delete [] p_mtbl_cf ; 
	    p_mtbl_cf = p_mtbl_cf_nouveau ; 
	    n_mtbl_cf = n_mtbl_cf_nouveau ; 
	    
	}
	else {
	
	    if (p_mtbl_cf[index] != 0x0) {
		cout << "Param::add_mtbl_cf : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_mtbl_cf[index] = &mi ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const Mtbl_cf& Param::get_mtbl_cf(int index) const {

    assert(index >= 0) ;
    assert(index < n_mtbl_cf) ; 
    
    return *(p_mtbl_cf[index]) ; 

} 
		    

		    //------------------------------------//
		    //		Etoile storage		  //
		    //------------------------------------//

// Total number of stored addresses
// --------------------------------

int Param::get_n_etoile() const {
    return n_etoile ; 
}

// Addition  
// --------
		    
void Param::add_etoile(const Etoile& eti, int index){
    
	if (index >= n_etoile) {    // p_etoile must be rescaled
	    	    
	    int n_etoile_nouveau = index + 1 ; 
	    const Etoile** p_etoile_nouveau = new const Etoile*[n_etoile_nouveau] ; 
	    
	   
	    // Copy of the previous addresses  
	    for (int i=0; i<n_etoile; i++) {
		p_etoile_nouveau[i] = p_etoile[i] ; 
	    }
	    
	    // The intermediate addresses are set to 0x0
	    for (int i=n_etoile; i<index; i++) {
		p_etoile_nouveau[i] = 0x0 ; 
	    }
	    
	    // The new address 
	    p_etoile_nouveau[index] = &eti ; 
	    
	    // Update 
	    if (n_etoile > 0) delete [] p_etoile ; 
	    p_etoile = p_etoile_nouveau ; 
	    n_etoile = n_etoile_nouveau ; 
	    
	}
	else {
	
	    if (p_etoile[index] != 0x0) {
		cout << "Param::add_etoile : the position " << index 
		     << " is already occupied !" << endl ; 
		abort() ; 
	    }
	    else{
		p_etoile[index] = &eti ; 
	    }
	    
	}   
    
}

// Extraction 
// ----------
		    
const Etoile& Param::get_etoile(int index) const {

    assert(index >= 0) ;
    assert(index < n_etoile) ; 
    
    return *(p_etoile[index]) ; 

} 
		    

		    
