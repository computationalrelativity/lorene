/*
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


char bhole_coal_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2002/10/16 14:36:32  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.15  2001/05/07  09:12:17  phil
 * *** empty log message ***
 *
 * Revision 2.14  2001/04/26  12:23:06  phil
 * *** empty log message ***
 *
 * Revision 2.13  2001/04/26  12:04:17  phil
 * *** empty log message ***
 *
 * Revision 2.12  2001/03/22  10:49:42  phil
 * *** empty log message ***
 *
 * Revision 2.11  2001/02/28  13:23:54  phil
 * vire kk_auto
 *
 * Revision 2.10  2001/01/29  14:31:04  phil
 * ajout tuype rotation
 *
 * Revision 2.9  2001/01/22  09:29:34  phil
 * vire convergence vers bare masse
 *
 * Revision 2.8  2001/01/10  09:31:52  phil
 * ajoute fait_kk_auto
 *
 * Revision 2.7  2000/12/20  15:02:57  phil
 * *** empty log message ***
 *
 * Revision 2.6  2000/12/20  09:09:48  phil
 * ajout set_statiques
 *
 * Revision 2.5  2000/12/18  17:43:06  phil
 * ajout sortie pour le rayon
 *
 * Revision 2.4  2000/12/18  16:38:39  phil
 * ajout convergence vers une masse donneee
 *
 * Revision 2.3  2000/12/14  10:45:38  phil
 * ATTENTION : PASSAGE DE PHI A PSI
 *
 * Revision 2.2  2000/12/04  14:29:17  phil
 * changement nom omega pour eviter interference avec membre prive
 *
 * Revision 2.1  2000/11/17  10:07:14  phil
 * corrections diverses
 *
 * Revision 2.0  2000/11/17  10:04:08  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>

// Lorene
#include "tenseur.h"
#include "bhole.h"


void Bhole_binaire::set_statiques (double precis, double relax) {
    
    int nz = hole1.mp.get_mg()->get_nzone() ;
    
    set_omega(0) ;
    init_bhole_binaire() ;
      
    int indic = 1 ;
    int conte = 0 ;
 
    cout << "TROUS STATIQUES : " << endl ;
    while (indic == 1) {
	Cmp lapse_un_old (hole1.get_n_auto()()) ;
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;
	
	double erreur = 0 ;
	Tbl diff (diffrelmax (lapse_un_old, hole1.get_n_auto()())) ;
	for (int i=1 ; i<nz ; i++)
	    if (diff(i) > erreur)
		erreur = diff(i) ;
	
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	
	if (erreur < precis)
	    indic = -1 ;
	conte ++ ;
    }
}

double Bhole_binaire::coal (double angulaire, double precis, double relax, 
			double nbre_ome, const int sortie) {
    
    assert (omega == 0) ;
    int nz = hole1.mp.get_mg()->get_nzone() ;
    
    int indic = 1 ;
    int conte = 0 ;
    
    char name_iteration[20] ;
    char name_correction[20] ;
    char name_viriel[20] ;
     
    sprintf(name_iteration, "ite_%e.dat", angulaire) ;
    sprintf(name_correction, "cor_%e.dat", angulaire) ;
    sprintf(name_viriel, "vir_%e.dat", angulaire) ;
    
    ofstream fiche_iteration(name_iteration) ;
    fiche_iteration.precision(8) ; 

    ofstream fiche_correction(name_correction) ;
    fiche_correction.precision(8) ; 
    
    ofstream fiche_viriel(name_viriel) ;
    fiche_viriel.precision(8) ; 
    
    // LA BOUCLE EN AUGMENTANT OMEGA A LA MAIN PROGRESSIVEMENT : 
    cout << "OMEGA AUGMENTE A LA MAIN." << endl ;
    double homme = 0 ;
    for (int pas = 0 ; pas <nbre_ome ; pas ++) {
	
	homme += angulaire/nbre_ome ;
	set_omega (homme) ;
	Cmp shift_un_old (hole1.get_shift_auto()(0)) ;
	
	solve_shift (precis, relax) ;
	fait_tkij() ;
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;
	
	double erreur = 0 ;
	Tbl diff (diffrelmax (shift_un_old, hole1.get_shift_auto()(0))) ;
	for (int i=1 ; i<nz ; i++)
	    if (diff(i) > erreur)
		erreur = diff(i) ;
	if (sortie != 0) {
	    fiche_iteration << conte << " " << erreur << endl ;
	    fiche_correction << conte << " " << hole1.get_regul() << endl ;
	    fiche_viriel << conte << " " << viriel() << endl ;
	    }
	    
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	conte ++ ;
    }
    
    // BOUCLE AVEC BLOQUE :
    cout << "OMEGA BLOQUE" << endl ;
    indic = 1 ; 
    double erreur ;
    while (indic == 1) {
	
	Cmp shift_un_old (hole1.get_shift_auto()(0)) ;
	solve_shift (precis, relax) ;
	fait_tkij() ;
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;

	erreur = 0 ;
	Tbl diff (diffrelmax (shift_un_old, hole1.get_shift_auto()(0))) ;
	for (int i=1 ; i<nz ; i++)
	    if (diff(i) > erreur)
		erreur = diff(i) ;
	
	if (sortie != 0) {
	    fiche_iteration << conte << " " << erreur << endl ;
	    fiche_correction << conte << " " << hole1.regul << endl ;
	    fiche_viriel << conte << " " << viriel() << endl ;
	    }
	    
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	if (erreur < precis)
	    indic = -1 ;
	conte ++ ;
    }
    
    fiche_iteration.close() ;
    fiche_correction.close() ;
    fiche_viriel.close() ;
      
    return viriel() ;
}
