/*
 * Methods of class Tbl
 *
 *  (see file tbl.h for documentation)
 *
 */

/*
 *   Copyright (c) 1999-2000 Jean-Alain Marck
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char tbl_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.8  2006/09/25 10:01:50  p_grandclement
 * Addition of N-dimensional Tbl
 *
 * Revision 1.7  2003/11/03 13:53:20  j_novak
 * Yet another efficiency improvement.
 *
 * Revision 1.6  2003/10/19 19:58:56  e_gourgoulhon
 * Slightly improved operator<<.
 *
 * Revision 1.5  2003/02/26 10:47:11  j_novak
 * The copy of a Tbl to another has been improved in speed.
 *
 * Revision 1.4  2002/10/16 14:37:13  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/09/24 08:32:07  e_gourgoulhon
 *
 * Added constructor from Matrice.
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
 * Revision 2.13  1999/11/24  16:00:43  eric
 * Modif affichage dimensions dans affiche_seuil.
 *
 * Revision 2.12  1999/11/23  13:32:44  eric
 * Le constructeur Tbl::Tbl(const Dim_tbl ) devient Tbl::Tbl(const Dim_tbl& ).
 * Le constructeur Tbl::Tbl(const Grille3d* ) devient
 *   Tbl(const Grille3d& ).
 * Modif affichage.
 *
 * Revision 2.11  1999/11/15  16:36:36  eric
 * Le membre dim est desormais un Dim_tbl et non plus un pointeur sur un
 * Dim_tbl.
 *
 * Revision 2.10  1999/10/29  15:45:44  eric
 * Ajout des cas 1-D et 2-D dans operator<<
 *
 * Revision 2.9  1999/10/29  15:05:16  eric
 * Suppression des fonctions membres min() et max():
 * elles deviennent des fonctions externes.
 *
 * Revision 2.8  1999/10/21  14:37:34  eric
 * *** empty log message ***
 *
 * Revision 2.7  1999/10/18  15:06:13  eric
 * La fonction membre annule() est rebaptisee annule_hard().
 *
 * Revision 2.6  1999/10/15  13:58:39  eric
 * Modification de l'affichage (operator<<).
 *
 * Revision 2.5  1999/10/01  12:36:08  eric
 * Ajout de la fonction affiche_seuil.
 *
 * Revision 2.4  1999/10/01  10:09:52  eric
 * 0 -> double(0)
 *
 * Revision 2.3  1999/09/30  12:51:46  eric
 * Enleve l'include "grilles.h" non necessaire.
 *
 * Revision 2.2  1999/09/24  14:24:23  eric
 * Depoussierage, changement de prototypes, etc...
 *
 * Revision 2.1  1999/03/02  16:49:42  eric
 * *** empty log message ***
 *
 * Revision 2.0  1999/02/15  10:42:45  hyc
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// headers C
#include <math.h>

// headers Lorene
#include "itbl.h"
#include "tbl.h"
#include "grilles.h"
#include "matrice.h"
#include "utilitaires.h"


			//---------------//
			// Constructeurs //
			//---------------//


// Constructeur 1D
Tbl::Tbl(int n1) : etat(ETATNONDEF), dim(n1), t(0x0) {}

// Constructeur 2D
Tbl::Tbl(int n1, int n0) : etat(ETATNONDEF), dim(n1, n0), t(0x0) {}

// Constructeur 3D
Tbl::Tbl(int n2, int n1, int n0) : etat(ETATNONDEF), dim(n2, n1, n0), t(0x0) {}

// Constructeur a partir d'une grille 3D
Tbl::Tbl(const Grille3d& g) : etat(ETATNONDEF), 
			      dim(g.get_np(), g.get_nt(), g.get_nr()), 
			      t(0x0) {}
			      
// Constructeur a partir d'un Itbl
Tbl::Tbl(Itbl sizes) : etat(ETATNONDEF), 
			      dim(1), 
			      t(0x0) {
			      
	int n = sizes.get_dim(0) ;
	int* dims = new int[n] ;
	for (int i=0 ; i<n ; i++)
	    dims[i] = sizes(i) ;
	Dim_tbl new_dim (n, dims) ;
	dim = new_dim ;	
	delete [] dims ;
}
			      
// Constructeur a partir d'un Dim_tbl
Tbl::Tbl(const Dim_tbl& dt) : etat(ETATNONDEF), dim(dt), t(0x0) {}

// Copie
Tbl::Tbl(const Tbl& tc) : etat(tc.etat), dim(tc.dim) {

    // La valeur eventuelle
    if (tc.etat == ETATQCQ) {
      int n = dim.taille ;
	t = new double[n] ;
	double* tin = tc.t ; 
	double* tout = t ;
	for (int i=0 ; i<n ; i++) {
	  *tout = *tin ;
	  tout++ ; tin++ ; 
	}
    }
    else{
	t = 0x0 ; 
    }
    
}

// From file
Tbl::Tbl(FILE* fd) : dim(fd) {

    fread_be(&etat, sizeof(int), 1, fd) ;		// etat
    
    // Le tableau
    if (etat == ETATQCQ) {
	t = new double[get_taille()] ;
	fread_be(t, sizeof(double), get_taille(), fd) ;	    // le tableau
    }
    else{
	t = 0x0 ; 
    }
}

// From a matrix
Tbl::Tbl(const Matrice& aa) : etat(aa.get_etat()),
							  dim( (aa.get_array()).dim ) {

	int nbl = dim.dim[1] ;
	int nbc = dim.dim[0] ;
	
	// Special case of one row :
	if (nbl == 1)  {
		dim.ndim = 1 ;
	}
	
	// Special case of one column :
	if (nbc == 1)  {
		dim.ndim = 1 ;
		dim.dim[0] = dim.dim[1] ;
		dim.dim[1] = 1 ;
	}
	
    if (etat == ETATQCQ) {

		t = new double[get_taille()] ;
		
		Tbl taa = aa.get_array() ;
		double* ta = taa.t ;
		
		for (int i=0 ; i<get_taille() ; i++) {
			t[i] = ta[i] ;
		}
		
	}
    else{
		t = 0x0 ;
    }

}

			//-------------//
			// Destructeur //
			//-------------//

Tbl::~Tbl() {
    delete [] t ;
}

			//-------------//
			// Affectation //
			//-------------//

// From Tbl
void Tbl::operator=(const Tbl& tx)
{
    // Protection
    assert( dim == tx.dim ) ;
    assert(tx.get_etat() != ETATNONDEF) ;

    switch (tx.etat) {
    case ETATZERO:
	set_etat_zero() ;
	break ;
	
    case ETATQCQ: {
      set_etat_qcq() ;
      int n = get_taille() ;
      double* tin = tx.t ;
      double* tout = t ;
      for (int i=0 ; i<n ; i++) {
	*tout = *tin ;
	tout++ ;
	tin ++ ;
      }
      break ;
    }
    default:
      cout << "Erreur bizarre !" << endl ;
      abort() ;
      break ;
    }
}

// From double
void Tbl::operator=(double a)
{
    if ( a == double(0) ) {
	set_etat_zero() ;
    }
    else {
	int n = get_taille() ;
	set_etat_qcq() ;
	for (int i=0 ; i<n ; i++) {
	    t[i] = a ;
	}
    }
}

// From int
void Tbl::operator=(int m)
{
    if (m == 0) {
	set_etat_zero() ;
    }
    else {
	int n = get_taille() ;
	set_etat_qcq() ;
	for (int i=0 ; i<n ; i++) {
	    t[i] = m ;
	}
    }
}

    
			//------------//
			// Sauvegarde //
			//------------//

// save in a file

void Tbl::sauve(FILE* fd) const {

    dim.sauve(fd) ;	    	    	    	    // dim
    fwrite_be(&etat, sizeof(int), 1, fd) ;		    // etat
    if (etat == ETATQCQ) {
	fwrite_be(t, sizeof(double), get_taille(), fd) ;   // le tableau
    }
}
    
		    //-----------------//
    	    	    // Gestion memoire //
		    //-----------------//

// Destructeur logique
void Tbl::del_t() {
    delete [] t ;
    t = 0x0 ;
    etat = ETATNONDEF ;
}

// ETATZERO
void Tbl::set_etat_zero() {
    if (etat == ETATZERO) return ;
    del_t() ;
    etat = ETATZERO ;
}

// ETATNONDEF
void Tbl::set_etat_nondef() {
    if (etat == ETATNONDEF) return ;
    del_t() ;
    etat = ETATNONDEF ;
}

// ETATQCQ
void Tbl::set_etat_qcq() {
    if (etat == ETATQCQ) return ;

    // Protection
    assert( (etat == ETATZERO) || (etat == ETATNONDEF) ) ; // sinon...

    t = new double[get_taille()] ;
    etat = ETATQCQ ;
}

// ZERO hard
void Tbl::annule_hard() {
    if (t == 0x0) {
	t = new double[get_taille()] ;
    }
    for (int i=0 ; i<get_taille() ; i++) {
	t[i] = 0. ;
    }
    etat = ETATQCQ ;
}


			//------------------------//
			//	Display		  //
			//------------------------//
			
//-----------			
// Operator<<
//-----------			

ostream& operator<<(ostream& o, const Tbl& t) {
    
    int ndim = t.get_ndim() ;
    o.precision(4);
    o.setf(ios::showpoint);
    o << "*** Tbl " << ndim << "D" << "   size: " ; 
    for (int i = 0; i<ndim-1; i++) {
	o << t.get_dim(i) << " x " ;
    } 
    o << t.get_dim(ndim-1) << "  =  " << t.get_taille() 
      << incindent << iendl ;


    if (t.get_etat() == ETATZERO) {
	o << "Identically ZERO" << decindent << iendl ;
	return o ;
    }

    if (t.get_etat() == ETATNONDEF) {
	o << "UNDEFINED STATE" << decindent << iendl ;
	return o ;
    }

    assert(t.etat == ETATQCQ) ;
    switch (ndim) {

	case 1 : {
	    for (int i=0 ; i<t.get_dim(0) ; i++) {
		o << " " << t(i)  ;
	    }
	    o << decindent << endl ;
	    break ;
	}


	case 2 : {
	    for (int j=0 ; j<t.get_dim(1) ; j++) {
		o << " j = " << j << " : " << incindent << iendl ;
		for (int i=0 ; i<t.get_dim(0) ; i++) {
		    o << " " << t(j, i)  ;
		}
		o << decindent << iendl ;
	    }
	    o << decindent << endl ;
	    break ;
	}
		
	case 3 : {
	    for (int k=0 ; k<t.get_dim(2) ; k++) {
		o << "k = " << k << " : " << incindent << iendl ;
		for (int j=0 ; j<t.get_dim(1) ; j++) {
		    o << "j = " << j << " : " << incindent ;
		    for (int i=0 ; i<t.get_dim(0) ; i++) {
			o << " " << t(k, j, i)  ;
		    }
		    o << decindent << iendl ;
		}
		o << decindent << iendl ;
	    }
	    o << decindent << endl ;
	    break ;
	}
		
	default : {
	    cout << "operator<< Tbl : unexpected dimension !" << endl ;
	    cout << " ndim = " << ndim << endl ; 	
	    abort() ;
	    break ;
	}
    }
    return o ;
}

//---------------
// Affiche_seuil
//---------------

void Tbl::affiche_seuil(ostream& ost, int precis,  double seuil) const {

    int ndim = get_ndim() ;
    ost << "*** Tbl " << ndim << "D" << "   size: " ; 
    for (int i = 0; i<ndim-1; i++) {
	ost << get_dim(i) << " x " ;
    } 
    ost << get_dim(ndim-1) << "  =  " << get_taille() << endl ; 

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
    
    // Affichage des elements du tableau 
    //----------------------------------
    
    ost << "    threshold for display : " << seuil << endl ; 
    ost.precision(precis);
    ost.setf(ios::showpoint);

    switch (get_ndim()) {
	case 1 : {			    	// cas 1-D

	    for (int i=0; i<get_dim(0); i++) {
		ost <<  " " << setw(precis) << (*this)(i)  ;
	    }
	    ost << endl ;
	    break ;
	}

	case 2 : {				// cas 2-D

	    for (int j=0; j<get_dim(1); j++) {
		ost <<  " #j=" << j << " : "  ;
		for (int i=0; i<get_dim(0); i++){
		    ost <<  " " << setw(precis) << (*this)(j, i)  ;
		}
		ost << endl;
	    }
	    ost << endl;
	    break;
	}

	case 3 : {				// cas 3-D
	    for (int k=0; k<get_dim(2); k++) {
		for (int j=0; j<get_dim(1); j++){
		    int test_imp = 0 ;
		    for (int i=0; i<get_dim(0); i++){
			if ( fabs( (*this)(k, j, i) ) >= seuil ) 
				    test_imp = 1 ; 
		    }
		    if (test_imp == 1 ) {
			ost <<  " #k=" << k <<",j=" << j << " : "  ;
			for (int i=0; i<get_dim(0); i++){
			    ost <<  " " << setw(precis) << (*this)(k, j, i) ;
			}
			ost << endl ;
		    }
		}
	    }
	    ost << endl;
	    break;
	}
	
	default : {
	    cout << "Tbl:affiche_seuil : unexpected dimension !" << endl ;
	    cout << " get_ndim() = " << get_ndim() << endl ; 	
	    abort() ; 
	    break;
	}               
	
    }     // fin du switch sur le nombre de dimensions

    // On restaure l'etat du flot ost a ses valeurs standards:
    ost.precision(6);
    ost.unsetf(ios::showpoint);
}



