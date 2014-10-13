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
 * Revision 1.8  2014/10/13 08:52:40  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.7  2014/10/06 15:12:58  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.6  2005/08/31 09:48:00  m_saijo
 * Delete one <math.h>
 *
 * Revision 1.5  2005/08/31 09:06:18  p_grandclement
 * add math.h in bhole_coal.C
 *
 * Revision 1.4  2005/08/29 15:10:14  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 * Revision 1.3  2003/11/13 13:43:54  p_grandclement
 * Addition of things needed for Bhole::update_metric (const Etoile_bin&, double, double)
 *
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
#include <cmath>
#include <cstdlib>

// Lorene
#include "tenseur.h"
#include "bhole.h"

namespace Lorene {
void Bhole_binaire::set_statiques (double precis, double relax) {
    
    int nz_un = hole1.mp.get_mg()->get_nzone() ;
    int nz_deux = hole2.mp.get_mg()->get_nzone() ;
    
    set_omega(0) ;
    init_bhole_binaire() ;
    
    int indic = 1 ;
    int conte = 0 ;
 
    cout << "TROUS STATIQUES : " << endl ;
    while (indic == 1) {
	Cmp lapse_un_old (hole1.get_n_auto()()) ;
	Cmp lapse_deux_old (hole2.get_n_auto()()) ;
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;
	
	double erreur = 0 ;
	Tbl diff_un (diffrelmax (lapse_un_old, hole1.get_n_auto()())) ;
	for (int i=1 ; i<nz_un ; i++)
	    if (diff_un(i) > erreur)
		erreur = diff_un(i) ;
	
	Tbl diff_deux (diffrelmax (lapse_deux_old, hole2.get_n_auto()())) ;
	for (int i=1 ; i<nz_deux ; i++)
	    if (diff_deux(i) > erreur)
		erreur = diff_deux(i) ;
	
		
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	
	if (erreur < precis)
	    indic = -1 ;
	conte ++ ;
    }
}

void Bhole_binaire::coal (double precis, double relax, int nbre_ome, double seuil_search, double m1, double m2, const int sortie) {

    assert (omega == 0) ;
    int nz1 = hole1.mp.get_mg()->get_nzone() ;
    int nz2 = hole1.mp.get_mg()->get_nzone() ;
    
    // Distance initiale
    double distance = hole1.mp.get_ori_x()-hole2.mp.get_ori_x() ;
    set_pos_axe (distance*(hole2.rayon/(hole1.rayon+hole2.rayon))) ;
    double scale_linear = (hole1.rayon + hole2.rayon)/2.*distance ;
    
    // Omega initial :
    double angulaire = sqrt((hole1.rayon+hole2.rayon)/distance/distance/distance) ;
    
    int indic = 1 ;
    int conte = 0 ;
    
    char name_iteration[40] ;
    char name_correction[40] ;
    char name_viriel[40] ;
    char name_ome [40] ;
    char name_linear[40] ;
    char name_axe[40] ;
    char name_error_m1[40] ;
    char name_error_m2[40] ;
    char name_r1[40] ;
    char name_r2[40] ;
    
    sprintf(name_iteration, "ite.dat") ;
    sprintf(name_correction, "cor.dat") ;
    sprintf(name_viriel, "vir.dat") ;
    sprintf(name_ome, "ome.dat") ;
    sprintf(name_linear, "linear.dat") ;
    sprintf(name_axe, "axe.dat") ;
    sprintf(name_error_m1, "error_m1.dat") ;
    sprintf(name_error_m2, "error_m2.dat") ; 
    sprintf(name_r1, "r1.dat") ;
    sprintf(name_r2, "r2.dat") ; 
    
    ofstream fiche_iteration(name_iteration) ;
    fiche_iteration.precision(8) ; 

    ofstream fiche_correction(name_correction) ;
    fiche_correction.precision(8) ; 
    
    ofstream fiche_viriel(name_viriel) ;
    fiche_viriel.precision(8) ; 
    
    ofstream fiche_ome(name_ome) ;
    fiche_ome.precision(8) ; 
   
    ofstream fiche_linear(name_linear) ;
    fiche_linear.precision(8) ; 
     
    ofstream fiche_axe(name_axe) ;
    fiche_axe.precision(8) ; 
    
    ofstream fiche_error_m1 (name_error_m1) ;
    fiche_error_m1.precision(8) ;
      
    ofstream fiche_error_m2 (name_error_m2) ;
    fiche_error_m2.precision(8) ;
  
    ofstream fiche_r1 (name_r1) ;
    fiche_r1.precision(8) ;
      
    ofstream fiche_r2 (name_r2) ;
    fiche_r2.precision(8) ;
    
    // LA BOUCLE EN AUGMENTANT OMEGA A LA MAIN PROGRESSIVEMENT : 
    cout << "OMEGA AUGMENTE A LA MAIN." << endl ;
    double homme = 0 ;
    for (int pas = 0 ; pas <nbre_ome ; pas ++) {
	
	homme += angulaire/nbre_ome ;
	set_omega (homme) ;
	
	Cmp shift_un_old (hole1.get_shift_auto()(0)) ;
	Cmp shift_deux_old (hole2.get_shift_auto()(0)) ;
	
	solve_shift (precis, relax) ;
	fait_tkij() ;
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;
	
	double erreur = 0 ;
	Tbl diff_un (diffrelmax (shift_un_old, hole1.get_shift_auto()(0))) ;
	for (int i=1 ; i<nz1 ; i++)
	    if (diff_un(i) > erreur)
		erreur = diff_un(i) ;
	
	Tbl diff_deux (diffrelmax (shift_deux_old, hole2.get_shift_auto()(0))) ;
	for (int i=1 ; i<nz2 ; i++)
	    if (diff_deux(i) > erreur)
		erreur = diff_deux(i) ;
	
	double error_viriel = viriel() ;
	double error_linear = linear_momentum_systeme_inf()(1)/scale_linear ;
	double error_m1 = 1.-sqrt(hole1.area()/16./M_PI)/m1 ;
	double error_m2 = 1.-sqrt(hole2.area()/16./M_PI)/m2 ;
	double r1 = hole1.mp.val_r(0, 1, 0, 0) ;
	double r2 = hole2.mp.val_r(0, 1, 0, 0) ;
		
	if (sortie != 0) {
	    fiche_iteration << conte << " " << erreur << endl ;
	    fiche_correction << conte << " " << hole1.get_regul() << " " << hole2.get_regul() << endl ;
	    fiche_viriel << conte << " " << error_viriel << endl ;
	    fiche_ome << conte << " " << homme << endl ;
	    fiche_linear << conte << " " << error_linear << endl ;
	    fiche_axe << conte << " " << pos_axe << endl ;
	    fiche_error_m1 << conte << " " << error_m1 << endl ;
	    fiche_error_m2 << conte << " " << error_m2 << endl ;
	    fiche_r1 << conte << " " << r1 << endl ;
	    fiche_r2 << conte << " " << r2 << endl ;
	    }
	    
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	conte ++ ;
    }
    
    // BOUCLE AVEC BLOQUE :
    cout << "OMEGA VARIABLE" << endl ;
    indic = 1 ;     
    bool scale = false ;
    
    while (indic == 1) {
	
	Cmp shift_un_old (hole1.get_shift_auto()(0)) ;
	Cmp shift_deux_old (hole2.get_shift_auto()(0)) ;
	
	solve_shift (precis, relax) ;
	fait_tkij() ;
	
	solve_psi (precis, relax) ;
	solve_lapse (precis, relax) ;

	double erreur = 0 ;
	Tbl diff_un (diffrelmax (shift_un_old, hole1.get_shift_auto()(0))) ;
	for (int i=1 ; i<nz1 ; i++)
	    if (diff_un(i) > erreur)
		erreur = diff_un(i) ;
	
	Tbl diff_deux (diffrelmax (shift_deux_old, hole2.get_shift_auto()(0))) ;
	for (int i=1 ; i<nz2 ; i++)
	    if (diff_deux(i) > erreur)
		erreur = diff_deux(i) ;
	
        double error_viriel = viriel() ;
	double error_linear = linear_momentum_systeme_inf()(1)/scale_linear ;
	double error_m1 = 1.-sqrt(hole1.area()/16./M_PI)/m1 ;
	double error_m2 = 1.-sqrt(hole2.area()/16./M_PI)/m2 ;
	double r1 = hole1.mp.val_r(0, 1, 0, 0) ;
	double r2 = hole2.mp.val_r(0, 1, 0, 0) ;
	
	if (sortie != 0) {
	    fiche_iteration << conte << " " << erreur << endl ;
	    fiche_correction << conte << " " << hole1.regul << " " << hole2.regul << endl ;
	    fiche_viriel << conte << " " << error_viriel << endl ;
	    fiche_ome << conte << " " << omega << endl ;
	    fiche_linear << conte << " " << error_linear << endl ;
	    fiche_axe << conte << " " << pos_axe << endl ; 
	    fiche_error_m1 << conte << " " << error_m1 << endl ;
	    fiche_error_m2 << conte << " " << error_m2 << endl ; 
	    fiche_r1 << conte << " " << r1 << endl ;
	    fiche_r2 << conte << " " << r2 << endl ;
	    }
	    
	// On modifie omega, position de l'axe et les masses !
	if (erreur <= seuil_search)
	    scale = true ;
	if (scale) {
	    double scaling_ome = pow((2-error_viriel)/(2-2*error_viriel), 1.) ;
	    set_omega (omega*scaling_ome) ; 
	    
	    double scaling_axe = pow((2-error_linear)/(2-2*error_linear), 0.1) ;
	    set_pos_axe (pos_axe*scaling_axe) ;
	    
	    double scaling_r1 = pow((2-error_m1)/(2-2*error_m1), 0.1) ;
	    hole1.mp.homothetie_interne(scaling_r1) ;
	  
	    double scaling_r2 = pow((2-error_m2)/(2-2*error_m2), 0.1) ;
	    hole2.mp.homothetie_interne(scaling_r2) ;
	}
	
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	if (erreur < precis)
	    indic = -1 ;
	conte ++ ;
    }
    
    fiche_iteration.close() ;
    fiche_correction.close() ;
    fiche_viriel.close() ;
    fiche_ome.close() ;
    fiche_linear.close() ;
    fiche_axe.close() ;
    fiche_error_m1.close() ;
    fiche_error_m2.close() ;
    fiche_r1.close() ;
    fiche_r2.close() ;
}
}
