/*
 *   Copyright (c) 2004 Philippe Grandclement
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


char bin_ns_bh_coal_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2005/12/01 12:59:10  p_grandclement
 * Files for bin_ns_bh project
 *
 * Revision 1.2  2005/10/18 13:12:32  p_grandclement
 * update of the mixted binary codes
 *
 * Revision 1.1  2005/08/29 15:10:15  p_grandclement
 * Addition of things needed :
 *   1) For BBH with different masses
 *   2) Provisory files for the mixted binaries (Bh and NS) : THIS IS NOT
 *   WORKING YET !!!
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>

// Lorene
#include "tenseur.h"
#include "bin_ns_bh.h"
#include "unites.h"
#include "graphique.h"

void Bin_ns_bh::coal (double precis, double relax, int itemax_equil, int itemax_mp_et, double ent_c_init, double seuil_dist, double dist, double m1, double m2, const int sortie) {
    
    using namespace Unites ;
     
    int nz_bh = hole.mp.get_mg()->get_nzone() ;
    int nz_ns = star.mp.get_mg()->get_nzone() ;
    
    // Distance initiale 
    double distance = fabs(hole.mp.get_ori_x()-star.mp.get_ori_x()) ;
    double mass_ns =  star.mass_g() * ggrav;
    double mass_bh =  hole.masse_adm_seul() ;
    //double axe_pos = distance*(mass_bh/(mass_bh+mass_ns)) ;
    double axe_pos = star.mp.get_ori_x() ;
    //double angulaire = sqrt((mass_bh+mass_ns)/distance/distance/distance) ;
    double angulaire = omega ;
    double scale_linear = (mass_ns+mass_bh)/2.*distance*angulaire ;
    //star.set_mp().set_ori (axe_pos, 0., 0.) ;
    //hole.set_mp().set_ori (-distance + axe_pos, 0., 0.) ;
    
    int conte = 0 ;
  
    char name_iteration[40] ;
    char name_correction[40] ;
    char name_viriel[40] ;
    char name_ome [40] ;
    char name_linear[40] ;
    char name_axe[40] ;
    char name_error_m1[40] ;
    char name_error_m2[40] ;
    
    sprintf(name_iteration, "ite.dat") ;
    sprintf(name_correction, "cor.dat") ;
    sprintf(name_viriel, "vir.dat") ;
    sprintf(name_ome, "ome.dat") ;
    sprintf(name_linear, "linear.dat") ;
    sprintf(name_axe, "axe.dat") ;
    sprintf(name_error_m1, "error_m_bh.dat") ;
    sprintf(name_error_m2, "error_m_ns.dat") ;
    
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
    
    
    // BOUCLE AVEC BLOQUE :
    bool loop = true ;     
    bool search_dist = false ;
    double ent_c = ent_c_init ;
    
    Cmp shift_bh_old (hole.mp) ;
    Cmp shift_ns_old (star.mp) ;
	
    double erreur ;
    
    while (loop) {
    
	if (hole.get_shift_auto().get_etat() != ETATZERO)
	    shift_bh_old = hole.get_shift_auto()(0) ;
	else
	    shift_bh_old = 0 ;
	    
	if (star.get_shift_auto().get_etat() != ETATZERO)
	    shift_ns_old = star.get_shift_auto()(0) ;
	else
	    shift_ns_old = 0 ;
		
	star.kinematics(omega, x_axe) ;
        star.fait_d_psi() ;
        star.hydro_euler() ;
	
        Tbl diff (7) ;
	diff.set_etat_qcq() ;
	int ite ;
	  
	star.equilibrium_nsbh (true, ent_c, ite, itemax_equil, itemax_mp_et, relax, itemax_mp_et, relax, diff) ;
	
	hole.update_metric(star) ;    
	
	hole.equilibrium (star, precis, relax) ;
        cout << "Apres equilibrium" << endl ;   
	
        star.update_metric(hole) ;
	cout << "Apres star::update_metric" << endl ;
	
	star.update_metric_der_comp(hole) ;
	cout << "Apres star::update_metric_der_comp" << endl ;	
	fait_tkij() ;
	cout << "Apres Bin_ns_bh::fait_tkij" << endl ;     
	
	erreur = 0 ;
	Tbl diff_bh (diffrelmax (shift_bh_old, hole.get_shift_auto()(0))) ;
	for (int i=1 ; i<nz_bh ; i++)
	    if (diff_bh(i) > erreur)
		erreur = diff_bh(i) ;
	
	Tbl diff_ns (diffrelmax (shift_ns_old, star.get_shift_auto()(0))) ;
	for (int i=0 ; i<nz_ns ; i++)
	    if (diff_ns(i) > erreur)
		erreur = diff_ns(i) ;
	
	if (erreur<seuil_dist)
	    search_dist = true ;
		
	cout << "Avant viriel" << endl ;
        double error_viriel = viriel() ;
	cout << "Apres viriel" << endl ;
	double error_linear = linear_momentum_systeme_inf()(1)/scale_linear ;
	cout << "Apres linear" << endl ;
	double error_m1 = 1.-sqrt(hole.area()/16./M_PI)/m1 ;
	cout << "Apres Mbh" << endl ;
	double error_m2 = 1 - star.mass_b()/m2 ;
	cout << "Apres Mns" << endl ;
	
	if (sortie != 0) {
	    fiche_iteration << conte << " " << erreur << endl ;
	    fiche_correction << conte << " " << hole.regul << endl ;
	    fiche_viriel << conte << " " << error_viriel << endl ;
	    fiche_linear << conte << " " << error_linear << endl ;
	    fiche_error_m1 << conte << " " << error_m1 << " " << sqrt(hole.area()/16./M_PI) << " " << m1 << endl ;
	    fiche_error_m2 << conte << " " << error_m2 << " " << star.mass_b() << " " << m2 << endl ; 
	    }
	
	// The axis position
	double scaling_axe = pow((2+error_linear)/(2+2*error_linear), 0.1) ;
	axe_pos *= scaling_axe ;
	star.set_mp().set_ori (axe_pos, 0, 0) ;
        hole.set_mp().set_ori (-distance + axe_pos, 0, 0) ;

	// Value of omega
	double new_ome = star.compute_angul() ;
        if (new_ome !=0)
		set_omega(relax*new_ome + (1-relax)*omega) ;
		//set_omega(new_ome) ;
		
	// Converge to the right masses :
	if (search_dist) {
	        double error_dist = (distance-dist)/dist ;
		double scale_d = pow((2+error_dist)/(2+2*error_dist), 0.2) ;
		distance *= scale_d ;
		cout << "Distance = " << distance << endl ;
		star.mp.homothetie(scale_d) ;
		hole.mp.homothetie(scale_d) ;
	/*
	    double scaling_r = pow((2-error_m1)/(2-2*error_m1), 0.2) ;
	    hole.mp.homothetie_interne(scaling_r) ;
	    hole.set_rayon(hole.get_rayon()*scaling_r) ;
	    
	    double scaling_ent = pow((2-error_m2)/(2-2*error_m2), 0.2) ;
	    ent_c *= scaling_ent ;
	    */
	    }
	fiche_ome << conte << " " << omega << " " << star.compute_angul() << endl ;
	fiche_axe << conte << " " << axe_pos << endl ;
	
	cout << "PAS TOTAL : " << conte << " DIFFERENCE : " << erreur << endl ;
	if (erreur < precis)
	    loop = false ;
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
    
}
