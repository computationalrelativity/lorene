/*
 *  Methods of class Scalar
 *
 *   (see file scalar.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Philippe Grandclement (for preceding class Cmp)
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
 * Revision 1.16  2003/10/19 19:54:37  e_gourgoulhon
 * -- Modified method spectral_display: now calling Valeur::display_coef.
 *
 * Revision 1.15  2003/10/15 16:03:38  j_novak
 * Added the angular Laplace operator for Scalar.
 *
 * Revision 1.14  2003/10/15 10:43:06  e_gourgoulhon
 * Added new members p_dsdt and p_stdsdp.
 *
 * Revision 1.13  2003/10/13 13:52:40  j_novak
 * Better managment of derived quantities.
 *
 * Revision 1.12  2003/10/11 14:42:16  e_gourgoulhon
 * Suppressed unusued argument new_triad in method change_triad.
 *
 * Revision 1.11  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
 * Revision 1.10  2003/10/07 08:05:03  j_novak
 * Added an assert for the constructor from a Tensor.
 *
 * Revision 1.9  2003/10/06 16:16:03  j_novak
 * New constructor from a Tensor.
 *
 * Revision 1.8  2003/10/06 13:58:48  j_novak
 * The memory management has been improved.
 * Implementation of the covariant derivative with respect to the exact Tensor
 * type.
 *
 * Revision 1.7  2003/10/05 21:15:42  e_gourgoulhon
 * - Suppressed method std_spectral_base_scal().
 * - Added method std_spectral_base().
 *
 * Revision 1.6  2003/10/01 13:04:44  e_gourgoulhon
 * The method Tensor::get_mp() returns now a reference (and not
 * a pointer) onto a mapping.
 *
 * Revision 1.5  2003/09/29 12:52:58  j_novak
 * Methods for changing the triad are implemented.
 *
 * Revision 1.4  2003/09/24 20:55:27  e_gourgoulhon
 * Added -- constructor by conversion of a Cmp
 *       -- operator=(const Cmp&)
 *
 * Revision 1.3  2003/09/24 15:10:54  j_novak
 * Suppression of the etat flag in class Tensor (still present in Scalar)
 *
 * Revision 1.2  2003/09/24 10:21:07  e_gourgoulhon
 * added more methods
 *
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
#include "proto.h"
#include "cmp.h"

			//---------------//
			// Constructors  //
			//---------------//


Scalar::Scalar(const Map& mpi) : Tensor(mpi), etat(ETATNONDEF), dzpuis(0), 
				 va(mpi.get_mg()) {

	cmp[0] = this ; 
    set_der_0x0() ;

}


// Constructor from a Tensor
// -------------------------
Scalar::Scalar(const Tensor& ti) : Tensor(*(ti.mp)), etat(ETATNONDEF), 
				   dzpuis(0), va(ti.cmp[0]->va) {

  assert(valence == 0) ;

  cmp[0] = this ; 
  set_der_0x0() ;

}


// Copy constructor
// ----------------
Scalar::Scalar(const Scalar& sci)  : Tensor(*(sci.mp)), etat(sci.etat), 
				     dzpuis(sci.dzpuis), va(sci.va) {
    
  cmp[0] = this ; 
  set_der_0x0() ;	// On ne recopie pas les derivees

}

// Conversion of a Cmp
//--------------------
Scalar::Scalar(const Cmp& ci) : Tensor(*(ci.get_mp())),
								etat(ci.get_etat()),
								dzpuis(ci.get_dzpuis()),
								va(ci.va) {
	cmp[0] = this ;
	set_der_0x0() ; 
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
    cmp[0] = 0x0 ; //cmp[0] was set to 'this' before (now 0x0 to break a 
                   // in loop in ~Tensor!)
           
}

			//-----------------------//
			// Memory management     //
			//-----------------------//

// Destructeur logique
void Scalar::del_t() {

    va.del_t() ;
    va.set_etat_nondef() ;
    Scalar::del_deriv() ;

}

void Scalar::del_deriv() const{
    if (p_dsdr != 0x0) delete p_dsdr ;
    if (p_srdsdt != 0x0) delete p_srdsdt ;
    if (p_srstdsdp != 0x0) delete p_srstdsdp ; 
    if (p_dsdt != 0x0) delete p_dsdt ;
    if (p_stdsdp != 0x0) delete p_stdsdp ;
    if (p_dsdx != 0x0) delete p_dsdx ; 
    if (p_dsdy != 0x0) delete p_dsdy ;
    if (p_dsdz != 0x0) delete p_dsdz ;
    if (p_lap != 0x0) delete p_lap ; 
    if (p_lapang != 0x0) delete p_lapang ; 
    if (p_integ != 0x0) delete p_integ ; 
    set_der_0x0() ;

    Tensor::del_deriv() ;
}

void Scalar::set_der_0x0() const {
    p_dsdr = 0x0 ;
    p_srdsdt = 0x0 ;
    p_srstdsdp = 0x0 ;
    p_dsdt = 0x0 ;
    p_stdsdp = 0x0 ;
    p_dsdx = 0x0 ;
    p_dsdy = 0x0 ;
    p_dsdz = 0x0 ;
    p_lap = 0x0 ; 
    p_lapang = 0x0 ; 
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

// ETATUN
void Scalar::set_etat_one() {
    if (etat == ETATUN) return ;
    else {
      del_deriv() ;
      va = 1 ;
      etat = ETATUN ;
    }
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
		if (p_dsdt != 0x0) p_dsdt->annule(l_min, l_max) ;
		if (p_stdsdp != 0x0) p_stdsdp->annule(l_min, l_max) ;
		if (p_dsdx != 0x0) p_dsdx->annule(l_min, l_max) ;
		if (p_dsdy != 0x0) p_dsdy->annule(l_min, l_max) ;
		if (p_dsdz != 0x0) p_dsdz->annule(l_min, l_max) ;
		if (p_lap != 0x0) p_lap->annule(l_min, l_max) ;
		if (p_lapang != 0x0) p_lapang->annule(l_min, l_max) ;
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
	
	case ETATUN: {
	    set_etat_one() ;
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
    
// From Cmp
// --------
void Scalar::operator=(const Cmp& ci) {
    

    // Menage general de la Valeur, mais pas des quantites derivees !
    va.del_t() ;
    
    // Les elements fixes
    assert( mp == ci.get_mp() ) ;
	
    dzpuis = ci.get_dzpuis() ; 
    
    // La valeur eventuelle
    switch(ci.get_etat()) {
	
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
	    cout << "Unkwown state in Scalar::operator=(const Cmp&) !" 
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
    if (x == double(1)) {
      set_etat_one() ;
    }
    else {
      set_etat_qcq() ;
      del_deriv() ;
	va = x ;
    }
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
    if (n == 1) {
      set_etat_one() ;
    }
    else {
      set_etat_qcq() ;
      del_deriv() ;
      va = n ;
    }
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
    
			//------------//
			// Impression //
			//------------//

// Operator <<
// -----------
ostream& operator<<(ostream& ost, const Scalar& ci) {

    switch(ci.etat) {
	case ETATNONDEF: {
	    ost << "*** UNDEFINED STATE *** " << endl ;
	    break ;
	}
	
	case ETATZERO: {
		ost << "*** identically ZERO ***" << endl ; 
	    break ; 
	}
	
	case ETATUN: {
		ost << "*** identically ONE ***" << endl ; 
	    break ; 
	}
	
	case ETATQCQ: {
	    ost << "*** dzpuis = " << ci.get_dzpuis() << endl ; 
	    ost << ci.va << endl ; 
	    break ;
	}
	
	default: {
	    cout << "operator<<(ostream&, const Scalar&) : unknown state !" 
		 << endl ;
	    abort() ;
	    break ; 
	}
    }
    
    // Termine
    return ost ;
}

// spectral_display
//-----------------

void Scalar::spectral_display(double thres, int precis, ostream& ost) const {

    // Cas particuliers
    //-----------------

    if (etat == ETATNONDEF) {
		ost << "*** UNDEFINED ***" << endl << endl ; 
		return ;
    }

    if (etat == ETATZERO) {
		ost << "*** identically ZERO ***" << endl ; 
		return ;
    }

    if (etat == ETATUN) {
		ost << "*** identically ONE ***" << endl ; 
		return ;
    }

    // Cas general : on affiche la Valeur
    //------------
	   
    if (dzpuis != 0) {
		ost << "*** dzpuis = " << dzpuis << endl ; 
 	}
    va.display_coef(thres, precis, ost) ;

}



    
		    //------------------------------------//
		    //	Spectral bases of the Valeur va   //
		    //------------------------------------//
		    
void Scalar::std_spectral_base() {
      
    va.std_base_scal() ;  
                   
}    


void Scalar::set_spectral_base(const Base_val& bi) {

	va.set_base(bi) ; 

}	 


		    //--------------------------//
		    //	dzpuis manipulations    //
		    //--------------------------//
		    
void Scalar::set_dzpuis(int dzi) {
    
    dzpuis = dzi ;
    
}

bool Scalar::dz_nonzero() const {
    
    assert(etat != ETATNONDEF) ; 
    
    const Mg3d* mg = mp->get_mg() ;
    
    int nzm1 = mg->get_nzone() - 1; 
    if (mg->get_type_r(nzm1) != UNSURR) {
		return false ; 
    } 
    
    if (etat == ETATZERO) {
		return false ; 
    }
    
    assert( (etat == ETATQCQ) || (etat == ETATUN)) ;//## to be checked!!
    
    if (va.etat == ETATZERO) {
		return false ; 
    }

    assert(va.etat == ETATQCQ) ; 
    
    if (va.c != 0x0) {
	if ( (va.c)->get_etat() == ETATZERO ) {
	    return false ; 
	}
	
	assert( (va.c)->get_etat() == ETATQCQ ) ; 
	if ( (va.c)->t[nzm1]->get_etat() == ETATZERO ) {
	    return false ; 
	}
	else {
	    assert( (va.c)->t[nzm1]->get_etat() == ETATQCQ ) ; 
	    return true ; 
	}
    }
    else{
	assert(va.c_cf != 0x0) ; 
	if ( (va.c_cf)->get_etat() == ETATZERO ) {
	    return false ; 
	}
	assert( (va.c_cf)->get_etat() == ETATQCQ ) ; 
	if ( (va.c_cf)->t[nzm1]->get_etat() == ETATZERO ) {
	    return false ; 
	}
	else {
	    assert( (va.c_cf)->t[nzm1]->get_etat() == ETATQCQ ) ; 
	    return true ; 
	}
    
    } 
    
}

bool Scalar::check_dzpuis(int dzi) const {
    
    if (dz_nonzero()) {	    // the check must be done
		return (dzpuis == dzi) ; 
    }
    else{
		return true ; 
    }
    
}



		//---------------------------------------//
		//	    Value at an arbitrary point		 //
		//---------------------------------------//

double Scalar::val_point(double r, double theta, double phi) const {

    assert(etat != ETATNONDEF) ; 
    
    if (etat == ETATZERO) {
		return double(0) ; 
    }
    
    if (etat == ETATUN) {
		return double(1) ; 
    }
    
    assert(etat == ETATQCQ) ; 
    
    // 1/ Search for the domain and the grid coordinates (xi,theta',phi')
    //    which corresponds to the point (r,theta,phi)
    
    int l ; 
    double xi ; 
    
    mp->val_lx(r, theta, phi, l,  xi) ;	    // call of val_lx with default 
					    // accuracy parameters
    
    // 2/ Call to the Valeur version
    
    return va.val_point(l, xi, theta, phi) ; 

}
 
		//---------------------------------//
        //	    Multipolar spectrum	       //
		//---------------------------------//

Tbl Scalar::multipole_spectrum() {

  assert (etat != ETATNONDEF) ;

  const Mg3d* mg = mp->get_mg() ;
  int nzone = mg->get_nzone() ;
  int lmax = 0 ;
  
  for (int lz=0; lz<nzone; lz++) 
    lmax = (lmax < 2*mg->get_nt(lz) - 1 ? 2*mg->get_nt(lz) - 1 : lmax) ;

  Tbl resu(nzone, lmax) ;
  if (etat == ETATZERO) {
    resu.set_etat_zero() ;
    return resu ;
  }

  assert((etat == ETATQCQ) || (etat == ETATUN));

  va.coef() ;
  va.ylm() ;
  resu.annule_hard() ;
  const Base_val& base = va.c_cf->base ;
  int m_quant, l_quant, base_r ;
  for (int lz=0; lz<nzone; lz++) 
    for (int k=0 ; k<mg->get_np(lz) ; k++) 
      for (int j=0 ; j<mg->get_nt(lz) ; j++) {
	if (nullite_plm(j, mg->get_nt(lz), k, mg->get_np(lz), base) == 1) 
	  {
	    // quantic numbers and spectral bases
	    donne_lm(nzone, lz, j, k, base, m_quant, l_quant, base_r) ;
	    for (int i=0; i<mg->get_nr(lz); i++) resu.set(lz, l_quant) 
				     += fabs((*va.c_cf)(0, k, j, i)) ; 
	  }
      }

  return resu ;
}

void Scalar::change_triad(const Base_vect& ) {

  cout << "WARNING: Scalar::change_triad : "<< endl ;
  cout << "This method does nothing ... " << endl ;

}
