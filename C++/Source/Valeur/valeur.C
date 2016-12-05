/*
 *  Methods of class Valeur
 *
 *   (see file valeur.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2003 Eric Gourgoulhon
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.9  2016/12/05 16:18:20  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.8  2014/10/13 08:53:49  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:13:22  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2005/10/25 08:56:40  p_grandclement
 * addition of std_spectral_base in the case of odd functions near the origin
 *
 * Revision 1.5  2004/11/23 12:45:00  f_limousin
 * Add function filtre_tp(int nn, int nz1, int nz2).
 *
 * Revision 1.4  2003/10/19 19:52:56  e_gourgoulhon
 * Added new method display_coef.
 *
 * Revision 1.3  2002/10/16 14:37:15  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.2  2001/12/04 21:27:54  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 2.37  2000/09/11  13:53:07  eric
 * Ajout des membres p_mult_cp et p_mult_sp.
 *
 * Revision 2.36  2000/09/08  10:07:13  eric
 * Ajout des methodes set_base_r, etc...
 *
 * Revision 2.35  1999/12/29  13:11:23  eric
 * Ajout de la fonction val_point_jk.
 *
 * Revision 2.34  1999/12/20  16:35:23  eric
 * Ajout de la fonction set_base.
 *
 * Revision 2.33  1999/12/10  16:09:47  eric
 * Annulation des membres derives dans la fonction annule(int,int).
 *
 * Revision 2.32  1999/12/07  14:53:00  eric
 * Changement ordre des arguments (phi,theta,xi) --> (xi,theta,phi)
 *  dans la routine val_point.
 *
 * Revision 2.31  1999/12/06  16:47:48  eric
 * Ajout de la fonction val_point.
 *
 * Revision 2.30  1999/11/30  12:41:53  eric
 * Le membre base est desormais un objet de type Base_val et non plus
 *  un pointeur vers une Base_val.
 *
 * Revision 2.29  1999/11/23  16:19:33  eric
 * Suppression du constructeur par defaut.
 *
 * Revision 2.28  1999/11/23  14:32:16  novak
 * Ajout des membres mult_ct et mult_st
 *
 * Revision 2.27  1999/11/22  15:41:33  eric
 * Ajout de la fonction annule(int l).
 *
 * Revision 2.26  1999/11/19  11:51:01  eric
 * Ajout du membre p_stdsdp.
 *
 * Revision 2.25  1999/11/19  10:53:17  eric
 * *** empty log message ***
 *
 * Revision 2.24  1999/11/19  10:34:52  eric
 * Corrections de diverses erreurs dans l'affectation (operator=):
 *  1/ l'auto-affectation est desormais interdite
 *  2/ l'affectation a des elements derives est desormais possible
 *
 * Revision 2.23  1999/11/19  09:32:15  eric
 * Ajout du membre p_lapang.
 *
 * Revision 2.22  1999/11/16  13:28:32  novak
 * Ajout de la gestion des pointeurs sur mult_x et scost
 *
 * Revision 2.21  1999/10/29  15:14:59  eric
 * *** empty log message ***
 *
 * Revision 2.20  1999/10/28  13:24:15  phil
 * copie passe avec ETATNONDEF
 *
 * Revision 2.19  1999/10/27  09:54:34  eric
 * *** empty log message ***
 *
 * Revision 2.18  1999/10/22  08:42:53  eric
 * *** empty log message ***
 *
 * Revision 2.17  1999/10/22  08:32:55  eric
 * Anglicisation de operator<<
 *
 * Revision 2.16  1999/10/21  14:33:05  eric
 * *** empty log message ***
 *
 * Revision 2.15  1999/10/21  14:21:23  eric
 * Constructeur par lecture de fichier.
 * Fonction sauve().
 *
 * Revision 2.14  1999/10/20  15:32:20  eric
 * Ajout de la routine Valeur::std_base_scal().
 * Ajout de operator=(const Mtbl_cf&).
 *
 * Revision 2.13  1999/10/19  15:30:33  eric
 * Ajout de la fonction affiche_seuil.
 *
 * Revision 2.12  1999/10/18  15:16:17  eric
 * *** empty log message ***
 *
 * Revision 2.11  1999/10/18  15:09:01  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 * Introduction de la fonction membre annule(int, int).
 *
 * Revision 2.10  1999/10/18  13:43:29  eric
 * Suppression de sxdsdx() et p_sxdsdx (car non-implementes).
 *
 * Revision 2.9  1999/10/13  15:52:35  eric
 * Depoussierage.
 *
 * Revision 2.8  1999/04/12  13:05:40  phil
 * *** empty log message ***
 *
 * Revision 2.7  1999/04/12  12:47:56  phil
 * Correction du constructeur par recopie.
 *
 * Revision 2.6  1999/04/12  12:13:26  phil
 * Correction de l'operateur Valeur = Valeur
 *
 * Revision 2.5  1999/03/01  15:10:42  eric
 * *** empty log message ***
 *
 * Revision 2.4  1999/02/26  11:43:07  hyc
 * *** empty log message ***
 *
 * Revision 2.3  1999/02/24  15:26:31  hyc
 * *** empty log message ***
 *
 * Revision 2.2  1999/02/23  11:46:40  hyc
 * *** empty log message ***
 *
 * $Header$
 *
 */

// headers C
#include <cassert>
#include <cstdlib>

// headers Lorene
#include "valeur.h"
#include "utilitaires.h"
#include "proto.h"



			//---------------//
			// Constructeurs //
			//---------------//

namespace Lorene {
Valeur::Valeur(const Mg3d& mgi) : mg(&mgi), base(mgi.get_nzone()) {
    
    // C'est nouveau
    nouveau() ;
}

Valeur::Valeur(const Mg3d* mgi) : mg(mgi), base(mgi->get_nzone()) {
    
    // C'est nouveau
    nouveau() ;
}

// Copie
Valeur::Valeur(const Valeur& vc) : base(vc.base) {
     
    // Tout par default
    mg = vc.get_mg() ;
    nouveau() ;
    
    // Que faire ?
    switch(vc.get_etat()) {
    	case ETATZERO:
	set_etat_zero() ;
	break ;
	
	case ETATNONDEF:
	set_etat_nondef() ;
	break ;
	
	case ETATQCQ:
	assert((vc.c != 0x0) || (vc.c_cf != 0x0) ) ;
	del_deriv() ;
	
	if (vc.c != 0x0) {
	    if (c == 0x0) {
		c = new Mtbl( *(vc.c) ) ;
	    }
	    else{
		*c = *(vc.c) ;
	    }
	}
	
	if (vc.c_cf != 0x0) {
	    if (c_cf == 0x0) {
		c_cf = new Mtbl_cf( *(vc.c_cf) ) ;
	    }
	    else{
		*c_cf = *(vc.c_cf) ;
	    }
	}
	
	etat = ETATQCQ ;
	break ;
	
	default:
	cout << "Etat pas possible" << endl ;
	abort() ;
	break ;
    }
}
	
// Constructeur a partir d'une grille et d'un fichier
// --------------------------------------------------
Valeur::Valeur(const Mg3d& g, FILE* fd) : mg(&g), base(fd) {

    // La multi-grille
    Mg3d* mg_tmp = new Mg3d(fd) ;	// la multi-grille d'origine
    if (*mg != *mg_tmp) {
	cout << 
	 "Valeur::Valeur(const Mg3d& g, FILE* fd): grid not consistent !" 
	 << endl ;
	abort() ;
    }
    delete mg_tmp ;

    fread_be(&etat, sizeof(int), 1, fd) ;	    // L'etat

    if (etat == ETATQCQ) {
	c = new Mtbl(g, fd) ;		    // Les valeurs ds l'espace des conf.
	
	// Tous les autres pointeurs sont mis a zero : 
	c_cf = 0x0 ; 
	set_der_0x0() ; 	
    }    
    else {
	c = 0x0 ; 
	c_cf = 0x0 ; 
	set_der_0x0() ; 	
    }

}

			//--------------//
			// Destructeurs //
			//--------------//

// Destructeur
Valeur::~Valeur() {
    del_t() ;
}

			//-------------//
			// Affectation //
			//-------------//

// From double
// ----------- 
void Valeur::operator=(double x) {

    // Cas x = 0
    if (x == 0) {
	set_etat_zero() ;
	return ;
    }

    // Cas general
    // les dependances
    set_etat_c_qcq() ;

    // Les bases
    base.set_base_nondef() ;
    
    *c = x ;
}
    
// From Valeur
// ----------- 
void Valeur::operator=(const Valeur& v) {

    // Protections:
    // ----------- 
    assert(&v != this) ;    // pour eviter l'auto-affectation
    assert(mg == v.mg) ;    // meme grille
    
    // Copie de v : 
    // ------------
    etat = v.etat ;  // l'etat 	    

    base = v.base ;	    // Les bases

    delete c ; 
    c = 0x0 ;		    // eventuellement provisoire...
	
    delete c_cf ; 
    c_cf = 0x0 ;	    // eventuellement provisoire...

    // La valeur eventuelle
    switch(v.etat) {
	case ETATNONDEF: {
	    set_etat_nondef() ; 
	    break ;		    // valeur par defaut
	}
	
	case ETATZERO: {
	    set_etat_zero() ;
	    break ;
	}
		
	case ETATQCQ: {
	    if (v.c != 0x0) {
		c = new Mtbl( *(v.c) ) ;
	    }
	
	    if (v.c_cf != 0x0) {
		c_cf = new Mtbl_cf( *(v.c_cf) ) ;
	    }
	
	    etat = ETATQCQ ;

	    // On detruit les dependances (seulement lorsque tout est fini !)
	    del_deriv() ;

	    break ;
	}
	
	default: {
	    cout << "Unkwon state in Valeur::operator=(const Valeur&) !" 
		 << endl ;
	    abort() ;
	    break ;
	}
    }
    
}
    
// From Mtbl 
// --------- 
void Valeur::operator=(const Mtbl& mt) {
    
    // Protections
    // -----------
    assert(mt.get_etat() != ETATNONDEF) ;
    assert(mg == mt.get_mg()) ;	    // meme grille
    assert(&mt != c) ;		    // pour eviter l'autoaffectation 
    
    // Les bases
    base.set_base_nondef() ;
    	
    delete c_cf ; 
    c_cf = 0x0 ;	    

    // Suivant le cas...
    switch(mt.get_etat()) {
    	case ETATZERO: {
	    set_etat_zero() ;
	    break ; 
	}
	
	case ETATQCQ: {
	    delete c ; 
	    c = new Mtbl(mt) ;

	    del_deriv() ;   // On detruit les dependances...

	    etat = ETATQCQ ; 
	    break ;
	}

	default: {
	    cout << "Unkwon state in Valeur::operator=(const Mtbl&) !" 
		 << endl ;
	    abort() ;
	    break ;
	}
    }

}

// From Mtbl_cf 
// ------------
void Valeur::operator=(const Mtbl_cf& mt) {

    // Protections
    // -----------
    assert(mt.get_etat() != ETATNONDEF) ;
    assert(mg == mt.get_mg()) ;	    // meme grille
    assert(&mt != c_cf) ;		    // pour eviter l'autoaffectation 
    
    // Les bases
    base = mt.base ;  
    	
    delete c ; 
    c = 0x0 ;	    

    // Suivant le cas...
    switch(mt.get_etat()) {
    	case ETATZERO: {
	    set_etat_zero() ;
	    break ; 
	}
	
	case ETATQCQ: {
	    delete c_cf ; 
	    c_cf = new Mtbl_cf(mt) ;

	    del_deriv() ;   // On detruit les dependances...

	    etat = ETATQCQ ; 
	    break ;
	}

	default: {
	    cout << "Unkwon state in Valeur::operator=(const Mtbl_cf&) !" 
		 << endl ;
	    abort() ;
	    break ;
	}
    }

}

			//------------//
			// Sauvegarde //
			//------------//

void Valeur::sauve(FILE* fd) const {

    base.sauve(fd) ;			    // la base 
    mg->sauve(fd) ;			    // la multi-grille
    fwrite_be(&etat, sizeof(int), 1, fd) ;	    // l'etat

    if (etat == ETATQCQ) {	
	if (c == 0x0) {		// La sauvegarde s'effectue dans l'espace
	    coef_i() ;		// des configurations
	} 
	c->sauve(fd) ;    
    }

}
    
			//------------//
			// Impression //
			//------------//

// Operator <<
// -----------
ostream& operator<<(ostream& o, const Valeur & vi) {
    
    switch(vi.etat) {
	case ETATNONDEF: {
	    o << "*** Valeur in UNDEFINED STATE" ;
	    break ;
	}
	
	case ETATZERO: {
	    o << "*** Valeur IDENTICALLY ZERO" ;
	    break ; 
	}
	
	case ETATQCQ: {
	    if (vi.c != 0x0) {
		o << "*** Valeur (configuration space) :" << endl ;
		o << *(vi.c) << endl ;
	    }
	    if (vi.c_cf != 0x0) {
		o << "*** Valeur (coefficients) :" << endl ;
		o << *(vi.c_cf) << endl ;
	    }
	    break ;
	}
	
	default: {
	    cout << "operator<<(ostream&, const Valeur&) : unknown state !" 
		 << endl ;
	    abort() ;
	    break ; 
	}

    }
    
    // Termine
    return o ;
}

// display_coef
// ------------
void Valeur::display_coef(double thres, int precis, ostream& ost) const {

    if (etat == ETATNONDEF) {
		ost << "    state: UNDEFINED" << endl ;
		return ;
    }

    if (etat == ETATZERO) {
		ost << "    state: ZERO" << endl ;
		return ;
    }

	coef() ; 	// the coefficients are required
	
	c_cf->display(thres, precis, ost) ; 
		
}


// affiche_seuil
//---------------

void Valeur::affiche_seuil(ostream& ost, int type, int precis,  
			   double seuil) const {
    ost << "*** Valeur " << endl ;


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

    switch(type) {
	case 0 : {
	    coef() ;	   // Mise a jour eventuelle des coefficients
	    ost << " Coefficients of the Valeur : " << endl ; 
	    c_cf->affiche_seuil(ost, precis, seuil) ;
	    break ; 
	}
	
	case 1 : {
	    coef_i() ;	   // Mise a jour eventuelle dans l'espace des conf.
	    ost << " Values of the Valeur at the collocation points: " << endl ; 
	    c->affiche_seuil(ost, precis, seuil) ;
	    break ; 
	}
	
	case 2 : {
	    coef() ;	   // Mise a jour eventuelle des coefficients
	    coef_i() ;	   // Mise a jour eventuelle dans l'espace des conf.
	    ost << " Coefficients of the Valeur : " << endl ; 
	    c_cf->affiche_seuil(ost, precis, seuil) ;
	    ost << " Values of the Valeur at the collocation points: " << endl ; 
	    c->affiche_seuil(ost, precis, seuil) ;
	    break ; 
	}
	
	default : {
	    cout << "Unknown type in Valeur::affiche_seuil !" << endl ; 
	    abort() ;
	    break ; 
	}
    }
 

}


    
			//-----------------------//
			// Gestion de la memoire //
			//-----------------------//

void Valeur::nouveau() {
    // Les donnees
    c = 0x0 ;
    c_cf = 0x0 ;
    set_der_0x0() ;
    // L'etat
    etat = ETATNONDEF ;
}

void Valeur::del_t() {

    delete c ; 
    c = 0x0 ;
	
    delete c_cf ; 
    c_cf = 0x0 ;

    del_deriv() ;

}

void Valeur::del_deriv() {
    // Destruction
    delete p_dsdx ;
    delete p_d2sdx2 ;
    delete p_sx ;
    delete p_sx2 ;
    delete p_mult_x ;

    delete p_dsdt ;
    delete p_d2sdt2 ;
    delete p_ssint ;
    delete p_scost ;
    delete p_mult_ct ;
    delete p_mult_st ;

    delete p_dsdp ;
    delete p_stdsdp ;
    delete p_d2sdp2 ;
    delete p_mult_cp ;
    delete p_mult_sp ;
    
    delete p_lapang ;

    // Pointeurs a 0x0
    set_der_0x0() ;
}

void Valeur::set_der_0x0() {
    p_dsdx = 0x0 ;
    p_d2sdx2 = 0x0 ;
    p_sx = 0x0 ;
    p_sx2 = 0x0 ;
    p_mult_x = 0x0 ;

    p_dsdt = 0x0 ;
    p_d2sdt2 = 0x0 ;
    p_ssint = 0x0 ;
    p_scost = 0x0 ;
    p_mult_ct = 0x0 ;
    p_mult_st = 0x0 ;

    p_dsdp = 0x0 ;
    p_stdsdp = 0x0 ;
    p_d2sdp2 = 0x0 ;
    p_mult_cp = 0x0 ;
    p_mult_sp = 0x0 ;

    p_lapang = 0x0 ;
}

// ETATZERO
void Valeur::set_etat_zero() {
    if (etat == ETATZERO) return ;
    del_t() ;
    etat = ETATZERO ;
}
// ETATNONDEF
void Valeur::set_etat_nondef() {
    if (etat == ETATNONDEF) return ;
    del_t() ;
    etat = ETATNONDEF ;
}
// ETATQCQ physique
void Valeur::set_etat_c_qcq() {
    // Detruit les dependances
    del_deriv() ;
    delete c_cf ; c_cf = 0x0 ;
        
    if (c == 0x0) {
	c = new Mtbl(mg) ;
    }
    etat = ETATQCQ ;
}
// ETATQCQ coefficients
void Valeur::set_etat_cf_qcq() {
    // Detruit les dependances
    del_deriv() ;
    delete c ; c = 0x0 ;
        
    if (c_cf == 0x0) {
	c_cf = new Mtbl_cf(mg, base) ;
    }
    etat = ETATQCQ ;
}
// ZERO hard
void Valeur::annule_hard() {
    // Detruit les dependances
    del_deriv() ;
    
    if (c == 0x0) {
	c = new Mtbl(mg) ;
    }
    if (c_cf == 0x0) {
	c_cf = new Mtbl_cf(mg, base) ;
    }
    
    c->annule_hard() ;
    c_cf->annule_hard() ;
    
    etat = ETATQCQ ;
}


// Sets the Valeur to zero in a given domain
// -----------------------------------------

void Valeur::annule(int l) {
    
    annule(l, l) ;     
}



// Sets the Valeur to zero in several domains
// ------------------------------------------

void Valeur::annule(int l_min, int l_max) {
    
    // Cas particulier: annulation globale : 
    if ( (l_min == 0) && (l_max == mg->get_nzone()-1) ) {
	set_etat_zero() ;
	return ; 
    }
    
    assert( etat != ETATNONDEF ) ; 
    
    if ( etat == ETATZERO ) {
	return ;		// rien n'a faire si c'est deja zero
    }
    else {
	assert( etat == ETATQCQ ) ;	// sinon...

	if (c != 0x0) {
	    c->annule(l_min, l_max) ;	    // Annule le Mtbl
	}

	if (c_cf != 0x0) {
	    c_cf->annule(l_min, l_max) ;    // Annule le Mtbl_cf
	}

	// Annulation des membres derives

	if (p_dsdx != 0x0) p_dsdx->annule(l_min, l_max) ;
	if (p_d2sdx2 != 0x0) p_d2sdx2->annule(l_min, l_max) ;
	if (p_sx != 0x0) p_sx->annule(l_min, l_max) ;
	if (p_sx2 != 0x0) p_sx2->annule(l_min, l_max) ;
	if (p_mult_x != 0x0) p_mult_x->annule(l_min, l_max) ;

	if (p_dsdt != 0x0) p_dsdt->annule(l_min, l_max) ;
	if (p_d2sdt2 != 0x0) p_d2sdt2->annule(l_min, l_max) ;
	if (p_ssint != 0x0) p_ssint->annule(l_min, l_max) ;
	if (p_scost != 0x0) p_scost->annule(l_min, l_max) ;
	if (p_mult_ct != 0x0) p_mult_ct->annule(l_min, l_max) ;
	if (p_mult_st != 0x0) p_mult_st->annule(l_min, l_max) ;

	if (p_dsdp != 0x0) p_dsdp->annule(l_min, l_max) ;
	if (p_stdsdp != 0x0) p_stdsdp->annule(l_min, l_max) ;
	if (p_d2sdp2 != 0x0) p_d2sdp2->annule(l_min, l_max) ;
	if (p_mult_cp != 0x0) p_mult_cp->annule(l_min, l_max) ;
	if (p_mult_sp != 0x0) p_mult_sp->annule(l_min, l_max) ;

	if (p_lapang != 0x0) p_lapang->annule(l_min, l_max) ;	
	 
    }
    
}


		    //--------------------------------------//
		    //	 Spectral bases manipulation        //
		    //--------------------------------------//
		    
void Valeur::set_base(const Base_val& base_i) {

    base = base_i ; 
    
    if (c_cf != 0x0) {
	if ( !(c_cf->base == base_i) ) {
	    delete c_cf ; 
	    c_cf = 0x0 ; 
	}
    }
    
}


void Valeur::std_base_scal() {
      
    set_base( mg->std_base_scal() ) ;  
                   
}    

void Valeur::std_base_scal_odd() {
      
    set_base( mg->std_base_scal_odd() ) ;  
                   
}    

void Valeur::set_base_r(int l, int base_r) {
    
    base.set_base_r(l, base_r) ; 
    
    if (c_cf != 0x0) {
	if ( !(c_cf->base == base) ) {
	    delete c_cf ; 
	    c_cf = 0x0 ; 
	}
    }
        
}

void Valeur::set_base_t(int base_t) {
    
    base.set_base_t(base_t) ; 
    
    if (c_cf != 0x0) {
	if ( !(c_cf->base == base) ) {
	    delete c_cf ; 
	    c_cf = 0x0 ; 
	}
    }
        
}

void Valeur::set_base_p(int base_p) {
    
    base.set_base_p(base_p) ; 
    
    if (c_cf != 0x0) {
	if ( !(c_cf->base == base) ) {
	    delete c_cf ; 
	    c_cf = 0x0 ; 
	}
    }
        
}




		//-----------------------------------------------//
		//	    Value at an arbitrary point		 //
		//-----------------------------------------------//

double Valeur::val_point(int l, double x, double theta, double phi) const {

    assert(etat != ETATNONDEF) ; 
    
    double resu ; 
    
    if (etat == ETATZERO) {
	resu = 0 ; 
    }
    else{
	assert(etat == ETATQCQ) ; 
	coef() ;			    // The coefficients are required 
	resu = c_cf->val_point(l, x, theta, phi) ;  // Call to the Mtbl_cf version
    }

    return resu ; 
}

double Valeur::val_point_jk(int l, double x, int j, int k) const {

    assert(etat != ETATNONDEF) ; 
    
    double resu ; 
    
    if (etat == ETATZERO) {
	resu = 0 ; 
    }
    else{
	assert(etat == ETATQCQ) ; 
	coef() ;			    // The coefficients are required 
	resu = c_cf->val_point_jk(l, x, j, k) ;  // Call to the Mtbl_cf version
    }

    return resu ; 
}


		//-----------------------------------------------//
		//	              Filtres	                 //
		//-----------------------------------------------//


void Valeur::filtre_tp(int nn, int nz1, int nz2) {

    int nz = mg->get_nzone() ;
    int nr = mg->get_nr(0) ;
    int nt = mg->get_nt(0) ;
    int np = mg->get_np(0) ;

     if (etat != ETATZERO) {
	assert (etat == ETATQCQ) ;
	ylm() ;
	for (int lz=nz1; lz<=nz2; lz++) {
	    if (c_cf->operator()(lz).get_etat() != ETATZERO)
		for (int k=0; k<np+1; k++)
		    for (int j=0; j<nt; j++) {
			    int l_q, m_q, base_r ;
			    if (nullite_plm(j, nt, k, np, base) == 1) {
				donne_lm(nz, lz, j, k, base, m_q, l_q,base_r) ;
				if (l_q > nn)
				    for (int i=0; i<nr; i++)
					c_cf->set(lz, k, j, i) = 0. ;
			    }
		    }
	}
	if (c != 0x0) {
	    delete c ;
	    c = 0x0 ;
	}
    }

}
}
