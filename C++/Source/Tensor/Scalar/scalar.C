/*
 *  Methods of class Scalar
 *
 *   (see file scalar.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char scalar_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/09/23 08:52:17  e_gourgoulhon
 * new version
 *
 *
 * $Header$
 *
 */

// headers C
#include <assert.h>
#include <stdlib.h>
#include <math.h>

// headers Lorene
#include "tensor.h"
#include "type_parite.h"
#include "utilitaires.h"

			//---------------//
			// Constructors  //
			//---------------//


Scalar::Scalar(const Map& mpi) : Tensor(mpi), dzpuis(0), va(mpi.get_mg()) {

	cmp[0] = this ; 
    set_der_0x0() ;

}

Scalar::Scalar(const Map* mpi) : Tensor(*mpi), dzpuis(0), va(mpi->get_mg()) {

	cmp[0] = this ; 
    set_der_0x0() ;

}



// Copy constructor
// ----------------
Scalar::Scalar(const Scalar& sci)  : Tensor(*(sci.mp)), dzpuis(sci.dzpuis), 
			   va(sci.va) {
    
	cmp[0] = this ; 
	etat = sci.etat ; 
    set_der_0x0() ;	// On ne recopie pas les derivees

}
	
// From file
// ---------
Scalar::Scalar(const Map& mpi, const Mg3d& mgi, FILE* fd) : Tensor(mpi), 
		va(mgi, fd) {

    assert( mpi.get_mg() == &mgi ) ; 

    fread_be(&etat, sizeof(int), 1, fd) ;		    // L'etat
    fread_be(&dzpuis, sizeof(int), 1, fd) ;	    // dzpuis

	cmp[0] = this ; 

    set_der_0x0() ;	// Les derivees sont initialisees a zero

}

			//--------------//
			// Destructor  //
			//--------------//

// Destructeur
Scalar::~Scalar() {
    del_t() ;
	cmp[0] = 0x0 ;
}

			//-----------------------//
			// Memory management     //
			//-----------------------//

// Destructeur logique
void Scalar::del_t() {

    va.del_t() ;
    del_deriv() ;

}

void Scalar::del_deriv() {
    delete p_dsdr ; p_dsdr = 0x0 ;
    delete p_srdsdt ; p_srdsdt = 0x0 ;
    delete p_srstdsdp ; p_srstdsdp = 0x0 ;
    delete p_dsdx ; p_dsdx = 0x0 ;
    delete p_dsdy ; p_dsdy = 0x0 ;
    delete p_dsdz ; p_dsdz = 0x0 ;
    delete p_lap ; p_lap = 0x0 ;
    delete p_integ ; p_integ = 0x0 ;
}

void Scalar::set_der_0x0() {
    p_dsdr = 0x0 ;
    p_srdsdt = 0x0 ;
    p_srstdsdp = 0x0 ;
    p_dsdx = 0x0 ;
    p_dsdy = 0x0 ;
    p_dsdz = 0x0 ;
    p_lap = 0x0 ; 
    ind_lap = - 1 ; 
    p_integ = 0x0 ; 
}

// ETATZERO
void Scalar::set_etat_zero() {
    if (etat == ETATZERO) return ;
    del_deriv() ;
    va.set_etat_zero() ;
    etat = ETATZERO ;
}

// ETATNONDEF
void Scalar::set_etat_nondef() {
    if (etat == ETATNONDEF) return ;
    del_t() ;
    etat = ETATNONDEF ;
}

// ETATQCQ
void Scalar::set_etat_qcq() {

    if (etat == ETATQCQ) {
		del_deriv() ; 
		return ;
    }
    
    // Protection
    assert( (etat == ETATZERO) || (etat == ETATNONDEF) ) ; // sinon...
    
    del_t() ;
    
    etat = ETATQCQ ;
}


// Allocates everything
// --------------------
void Scalar::allocate_all() {
    
	set_etat_qcq() ; 
	va.set_etat_c_qcq() ;	    // allocation in configuration space
	Mtbl* mt = va.c ; 
	mt->set_etat_qcq() ;
	for (int l=0; l<mt->get_nzone(); l++) {
	    mt->t[l]->set_etat_qcq() ; 
	}
	
} 



// ZERO hard
void Scalar::annule_hard() {

    va.annule_hard() ;
    del_deriv() ; 
    etat = ETATQCQ ;
}


// Sets the Scalar to zero in several domains
// ---------------------------------------

void Scalar::annule(int l_min, int l_max) {
    
    // Cas particulier: annulation globale : 
    if ( (l_min == 0) && (l_max == va.mg->get_nzone()-1) ) {
		set_etat_zero() ;
		return ; 
    }
    
    assert( etat != ETATNONDEF ) ; 
    
    if ( etat == ETATZERO ) {
		return ;		// rien n'a faire si c'est deja zero
    }
    else {
		assert( etat == ETATQCQ ) ;	// sinon...

		va.annule(l_min, l_max) ;	// Annule la Valeur
	
		// Annulation des membres derives
		if (p_dsdr != 0x0) p_dsdr->annule(l_min, l_max) ;
		if (p_srdsdt != 0x0) p_srdsdt->annule(l_min, l_max) ;
		if (p_srstdsdp != 0x0) p_srstdsdp->annule(l_min, l_max) ;
		if (p_dsdx != 0x0) p_dsdx->annule(l_min, l_max) ;
		if (p_dsdy != 0x0) p_dsdy->annule(l_min, l_max) ;
		if (p_dsdz != 0x0) p_dsdz->annule(l_min, l_max) ;
		if (p_lap != 0x0) p_lap->annule(l_min, l_max) ;
		if (p_integ != 0x0) delete p_integ ;
    }
    
}


			//------------//
			// Assignment //
			//------------//


// From tensor
// -----------

void Scalar::operator=(const Tensor& uu) {

	assert(uu.valence == 0) ; 
	
	operator=(*(uu.cmp[0])) ; 

}

// From Scalar
// ----------
void Scalar::operator=(const Scalar& ci) {
    

    assert(&ci != this) ;    // pour eviter l'auto-affectation

    // Menage general de la Valeur, mais pas des quantites derivees !
    va.del_t() ;
    
    // Les elements fixes
    assert( mp == ci.mp ) ;
	
    dzpuis = ci.dzpuis ; 
    
    // La valeur eventuelle
    switch(ci.etat) {
	case ETATNONDEF: {
	    set_etat_nondef() ; 
	    break ;		    // valeur par defaut
	}
	
	case ETATZERO: {
	    set_etat_zero() ;
	    break ;
	}
	
	case ETATQCQ: {
	    set_etat_qcq() ;
	    va = ci.va ;

	    // On detruit les quantites derivees (seulement lorsque tout est fini !)
	    del_deriv() ; 

	    break ;
	}
	
	default: {
	    cout << "Unkwown state in Scalar::operator=(const Scalar&) !" 
		 << endl ;
	    abort() ;
	    break ;
	}
    }

}
    
// From Valeur
// -----------
void Scalar::operator=(const Valeur& vi) {

    // Traitement de l'auto-affectation :
    if (&vi == &va) {
		return ; 
    }

    // Protection
    assert(vi.get_etat() != ETATNONDEF) ;
    
    // Menage general de la Valeur, mais pas des quantites derivees !
    va.del_t() ;

    
    // La valeure eventuelle
    switch(vi.get_etat()) {

	case ETATZERO: {
	    set_etat_zero() ;
	    break ;
	}

	case ETATQCQ: {
	    set_etat_qcq() ;
	    va = vi ;
	    
	    // On detruit les quantites derivees (seulement lorsque tout est fini !)
	    del_deriv() ; 

	    break ;
	}
	
	default: {
	    cout << "Unkwown state in Scalar::operator=(const Valeur&) !" << endl ;
	    abort() ;
	    break ;
	}
    }

}

// From Mtbl
// ---------
void Scalar::operator=(const Mtbl& mi) {
    
    // Protection
    assert(mi.get_etat() != ETATNONDEF) ;

    assert(&mi != va.c) ;  // pour eviter l'auto-affectation

   
    // Menage general de la Valeur, mais pas des quantites derivees !
    va.del_t() ;

    // La valeure eventuelle
    switch(mi.get_etat()) {
	case ETATZERO: {
	    set_etat_zero() ;
	    break ;
	}
	
	case ETATQCQ: {
	    set_etat_qcq() ;
	    va = mi ;

	    // On detruit les quantites derivees (seulement lorsque tout est fini !)
	    del_deriv() ; 

	    break ;
	 }
	
	default: {
	    cout << "Unkwown state in Scalar::operator=(const Mtbl&) !" << endl ;
	    abort() ;
	    break ;
	}
    }


}

// From double
// -----------
void Scalar::operator=(double x) {
    
    if (x == double(0)) {
		set_etat_zero() ;
    }
    else {
		set_etat_qcq() ;
		del_deriv() ;
		va = x ;
    }

    dzpuis = 0 ; 
}

// From int
// --------
void Scalar::operator=(int n) {
    
    if (n == 0) {
		set_etat_zero() ;
    }
    else {
		set_etat_qcq() ;
		del_deriv() ;
		va = n ;
    }

    dzpuis = 0 ; 

}


			//------------//
			// Sauvegarde //
			//------------//

void Scalar::sauve(FILE* fd) const {

    va.sauve(fd) ;	    // la valeur (en premier pour la construction
			    //   lors de la lecture du fichier)

    fwrite_be(&etat, sizeof(int), 1, fd) ;		    // l'etat
    fwrite_be(&dzpuis, sizeof(int), 1, fd) ;	    // dzpuis

}
    

