/*
 *  Methods of class Cmp
 *
 *   (see file cmp.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char cmp_C[] = "$Header$" ;


/*
 * $Id$
 * $Log$
 * Revision 1.6  2003/10/01 15:49:33  e_gourgoulhon
 * Method Scalar::get_mp() now returns a reference onto a mapping.
 *
 * Revision 1.5  2003/09/24 20:54:24  e_gourgoulhon
 * Added constructor by conversion of a Scalar.
 *
 * Revision 1.4  2003/08/26 09:46:10  j_novak
 * Added the method multipole_spectrum
 *
 * Revision 1.3  2002/10/16 14:36:33  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
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
 * Revision 2.36  2000/09/13  12:11:56  eric
 * Ajout de la fonction allocate_all().
 *
 * Revision 2.35  2000/09/06  09:45:23  keisuke
 * Ajout de l'appel a del_deriv() dans set_etat_qcq() pour le cas
 *   ou etat est deja ETATQCQ.
 *
 * Revision 2.34  2000/08/16  10:43:30  eric
 * Suppression de Mtbl_cf::dzpuis.
 *
 * Revision 2.33  2000/08/16  10:30:41  eric
 * Suppression de Mtbl::dzpuis.
 *
 * Revision 2.32  2000/03/07  16:52:58  eric
 * Modif dz_nonzero() : cas etat=ETATQCQ et va.etat=ETATZERO.
 *
 * Revision 2.31  2000/01/28  16:09:21  eric
 * Ajout des fonctions dz_nonzero et check_dzpuis.
 *
 * Revision 2.30  1999/12/22  16:44:28  eric
 * set_etat_zero() : remplacement de l'appel a del_t() par
 * 			1/ va.set_etat_zero() ;
 * 			2/ del_deriv() ;
 *
 * Revision 2.29  1999/12/10  15:59:14  eric
 * Changement de la place de del_deriv() dans l'affectation
 *   (pour permettre l'affectation a des membres derives).
 * Annulation des membres derives.
 *  dans la fonction annule(int,int).
 *
 * Revision 2.28  1999/12/09  10:50:01  eric
 * Ajout du membre p_integ.
 *
 * Revision 2.27  1999/12/07  14:53:32  eric
 * Changement ordre des arguments (phi,theta,r) --> (r,theta,phi)
 *   dans la routine val_point.
 *
 * Revision 2.26  1999/12/06  16:47:11  eric
 * Ajout de la fonction val_point.
 *
 * Revision 2.25  1999/11/30  16:30:23  eric
 * *** empty log message ***
 *
 * Revision 2.24  1999/11/30  16:26:39  eric
 * Ajout (provisoire) de l'affectation des dzpuis des Mtbl et Mtbl_cf
 *  dans set_dzpuis.
 *
 * Revision 2.23  1999/11/29  15:14:48  phil
 * *** empty log message ***
 *
 * Revision 2.22  1999/11/29  12:56:49  eric
 * Introduction des membres p_lap, ind_lap.
 *
 * Revision 2.21  1999/11/26  14:23:17  eric
 * Ajout du membre dzpuis et des fonctions de manipulation associees.
 *
 * Revision 2.20  1999/11/25  16:27:45  eric
 * Reorganisation complete du calcul et stokage des derivees partielles.
 *
 * Revision 2.19  1999/11/23  16:21:02  eric
 * Suppression du membre statique Cmp_Zero.
 * Suppression du constructeur par defaut.
 *
 * Revision 2.18  1999/11/22  16:34:56  eric
 * Ajout du constructeur prive sans argument pour Cmp_Zero.
 *
 * Revision 2.17  1999/11/22  15:41:57  eric
 * Ajout de la fonction annule(int l).
 *
 * Revision 2.16  1999/10/29  08:14:59  eric
 * Ajout de assert( mpi.get_mg() == &mgi ) dans le constructeur
 * par lecture de fichier.
 *
 * Revision 2.15  1999/10/28  09:39:00  eric
 * Ajout de la fonction affiche_seuil.
 *
 * Revision 2.14  1999/10/28  09:01:56  eric
 * Constructeur par lecture de fichier.
 * Ajout de la fonction annule(int, int).
 *
 * Revision 2.13  1999/10/27  15:38:52  eric
 * Suppression du membre c.
 *
 * Revision 2.12  1999/10/27  09:51:46  eric
 * *** empty log message ***
 *
 * Revision 2.11  1999/10/27  08:45:31  eric
 * Introduction du membre Valeur va.
 * Le pointeur Valeur* c est desormais un membre prive constant qui pointe
 *
 * sur va.
 *
 * Revision 2.10  1999/10/22  08:14:32  eric
 * Depoussierage.
 *
 * Revision 2.9  1999/10/18  16:08:15  phil
 * Correction de set_etat_qcq
 * Evite les memory leak
 *
 * Revision 2.8  1999/10/18  15:07:58  eric
 * La fonction membre Valeur::annule() est rebaptisee Valeur::annule_hard().
 *
 * Revision 2.7  1999/04/09  13:38:58  phil
 * *** empty log message ***
 *
 * Revision 2.6  1999/04/09  13:10:09  phil
 * ajout de cmp = valeur
 *
 * Revision 2.5  1999/03/03  11:16:24  hyc
 * *** empty log message ***
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
#include "cmp.h"
#include "type_parite.h"
#include "utilitaires.h"
#include "proto.h"
#include "tensor.h"

			//---------------//
			// Constructeurs //
			//---------------//


Cmp::Cmp(const Map& mpi) : mp(&mpi), etat(ETATNONDEF), dzpuis(0), 
			   va(mpi.get_mg()) {

    set_der_0x0() ;

}

Cmp::Cmp(const Map* mpi) : mp(mpi), etat(ETATNONDEF),  dzpuis(0),
			   va(mpi->get_mg()) {

    set_der_0x0() ;

}


// Copy constructor
// ----------------
Cmp::Cmp(const Cmp& ci)  : mp(ci.mp), etat(ci.etat), dzpuis(ci.dzpuis), 
			   va(ci.va) {
    
    set_der_0x0() ;	// On ne recopie pas les derivees

}

// Conversion of a Scalar
//-----------------------
Cmp::Cmp(const Scalar& uu) : mp(&(uu.get_mp())), 
							 etat(uu.get_etat()),
							 dzpuis(uu.get_dzpuis()),
							 va(uu.get_spectral_va()) {							 
	set_der_0x0() ; 
}
	
// From file
// ---------
Cmp::Cmp(const Map& mpi, const Mg3d& mgi, FILE* fd) : mp(&mpi), va(mgi, fd) {

    assert( mpi.get_mg() == &mgi ) ; 

    fread_be(&etat, sizeof(int), 1, fd) ;		    // L'etat
    fread_be(&dzpuis, sizeof(int), 1, fd) ;	    // dzpuis

    set_der_0x0() ;	// Les derivees sont initialisees a zero

}

			//--------------//
			// Destructeur  //
			//--------------//

// Destructeur
Cmp::~Cmp() {
    del_t() ;
}

			//-----------------------//
			// Gestion de la memoire //
			//-----------------------//

// Destructeur logique
void Cmp::del_t() {
    va.del_t() ;
    del_deriv() ;
    etat = ETATNONDEF ;
}

void Cmp::del_deriv() {
    delete p_dsdr ; p_dsdr = 0x0 ;
    delete p_srdsdt ; p_srdsdt = 0x0 ;
    delete p_srstdsdp ; p_srstdsdp = 0x0 ;
    delete p_dsdx ; p_dsdx = 0x0 ;
    delete p_dsdy ; p_dsdy = 0x0 ;
    delete p_dsdz ; p_dsdz = 0x0 ;
    delete p_lap ; p_lap = 0x0 ;
    delete p_integ ; p_integ = 0x0 ;
}

void Cmp::set_der_0x0() {
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
void Cmp::set_etat_zero() {
    if (etat == ETATZERO) return ;
    del_deriv() ;
    va.set_etat_zero() ;
    etat = ETATZERO ;
}

// ETATNONDEF
void Cmp::set_etat_nondef() {
    if (etat == ETATNONDEF) return ;
    del_t() ;
    etat = ETATNONDEF ;
}

// ETATQCQ
void Cmp::set_etat_qcq() {

    if (etat == ETATQCQ) {
	del_deriv() ; 
	return ;
    }
    
    // Protection
    assert( (etat == ETATZERO) || (etat == ETATNONDEF) ) ; // sinon...
    
    del_t() ;
    
    // Termine
    etat = ETATQCQ ;
}


// Allocates everything
// --------------------
void Cmp::allocate_all() {
    
	set_etat_qcq() ; 
	va.set_etat_c_qcq() ;	    // allocation in configuration space
	Mtbl* mt = va.c ; 
	mt->set_etat_qcq() ;
	for (int l=0; l<mt->get_nzone(); l++) {
	    mt->t[l]->set_etat_qcq() ; 
	}
	
} 



// ZERO hard
void Cmp::annule_hard() {

    va.annule_hard() ;
    del_deriv() ; 
    etat = ETATQCQ ;
}

// Sets the Cmp to zero in a given domain
// --------------------------------------

void Cmp::annule(int l) {
    
    annule(l, l) ;     
}


// Sets the Cmp to zero in several domains
// ---------------------------------------

void Cmp::annule(int l_min, int l_max) {
    
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

// From Cmp
// --------
void Cmp::operator=(const Cmp& ci) {
    
    assert(&ci != this) ;    // pour eviter l'auto-affectation

    // Menage general de la Valeur, mais pas des quantites derivees !
    va.del_t() ;
    
    // Les elements fixes
    mp = ci.mp ;
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
	    cout << "Unkwown state in Cmp::operator=(const Cmp&) !" 
		 << endl ;
	    abort() ;
	    break ;
	}
    }

}
    
// From Valeur
// -----------
void Cmp::operator=(const Valeur& vi) {

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
	    cout << "Unkwown state in Cmp::operator=(const Valeur&) !" << endl ;
	    abort() ;
	    break ;
	}
    }

}

// From Mtbl
// ---------
void Cmp::operator=(const Mtbl& mi) {
    
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
	    cout << "Unkwown state in Cmp::operator=(const Mtbl&) !" << endl ;
	    abort() ;
	    break ;
	}
    }


}

// From double
// -----------
void Cmp::operator=(double x) {
    
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
void Cmp::operator=(int n) {
    
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

void Cmp::sauve(FILE* fd) const {

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
ostream& operator<<(ostream& o, const Cmp& ci) {

    switch(ci.etat) {
	case ETATNONDEF: {
	    o << "*** Cmp in UNDEFINED STATE" ;
	    break ;
	}
	
	case ETATZERO: {
	    o << "*** Cmp IDENTICALLY ZERO" ;
	    break ; 
	}
	
	case ETATQCQ: {
	    o << "*** Cmp : " << endl ; 
	    o << "                        dzpuis = " << ci.get_dzpuis() << endl ; 
	    o << ci.va << endl ; 
	    break ;
	}
	
	default: {
	    cout << "operator<<(ostream&, const Cmp&) : unknown state !" 
		 << endl ;
	    abort() ;
	    break ; 
	}
    }
    
    // Termine
    return o ;
}

// affiche_seuil
//---------------

void Cmp::affiche_seuil(ostream& ost, int type, int precis, 
			double seuil) const {
    ost << "*** Cmp " << endl ;

    // Cas particuliers
    //-----------------

    if (etat == ETATNONDEF) {
	ost << "    state: UNDEFINED" << endl ;
	return ;
    }

    if (etat == ETATZERO) {
	ost << "    state: ZERO" << endl ;
	return ;
    }

    // Cas general : on affiche la Valeur
    //------------
	   
    ost << "                        dzpuis = " << dzpuis << endl ; 
    va.affiche_seuil(ost, type, precis, seuil) ;

}



    
		    //------------------------------------//
		    //	Spectral bases of the Valeur va   //
		    //------------------------------------//
		    
void Cmp::std_base_scal() {
      
    va.std_base_scal() ;  
                   
}    

		    //--------------------------//
		    //	dzpuis manipulations    //
		    //--------------------------//
		    
void Cmp::set_dzpuis(int dzi) {
    
    dzpuis = dzi ;
    
}

bool Cmp::dz_nonzero() const {
    
    assert(etat != ETATNONDEF) ; 
    
    const Mg3d* mg = mp->get_mg() ;
    
    int nzm1 = mg->get_nzone() - 1; 
    if (mg->get_type_r(nzm1) != UNSURR) {
	return false ; 
    } 
    
    if (etat == ETATZERO) {
	return false ; 
    }
    
    assert(etat == ETATQCQ) ;
    
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

bool Cmp::check_dzpuis(int dzi) const {
    
    if (dz_nonzero()) {	    // the check must be done
	return (dzpuis == dzi) ; 
    }
    else{
	return true ; 
    }
    
}



		//-----------------------------------------------//
		//	    Value at an arbitrary point		 //
		//-----------------------------------------------//

double Cmp::val_point(double r, double theta, double phi) const {

    assert(etat != ETATNONDEF) ; 
    
    if (etat == ETATZERO) {
	return double(0) ; 
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
 

		//-------------------------------------//
                //	    Multipolar spectrum	       //
		//-------------------------------------//

Tbl Cmp::multipole_spectrum() {
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

  assert(etat == ETATQCQ) ;

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

