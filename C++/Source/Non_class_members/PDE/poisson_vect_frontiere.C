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


char poisson_vect_frontiere_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:28  e_gourgoulhon
 * Initial revision
 *
 * Revision 2.2  2000/10/26  09:08:06  phil
 * *** empty log message ***
 *
 * Revision 2.1  2000/10/26  09:01:18  phil
 * *** empty log message ***
 *
 * Revision 2.0  2000/10/19  09:36:36  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

// Header C : 
#include <stdlib.h>
#include <math.h>

// Headers Lorene :
#include "proto.h"
#include "valeur.h"
#include "cmp.h"
#include "tenseur.h"

    // USING OOhara
void poisson_vect_frontiere (double lambda, const Tenseur& source, Tenseur& shift, 
	    const Valeur& lim_x, const Valeur& lim_y, const Valeur& lim_z, 
	    int num_front, double precision, int itermax) {
	
    // METTRE TOUT PLEIN D'ASSERT
   
    // Confort
    int nt = lim_x.get_mg()->get_nt(num_front+1) ;
    int np = lim_x.get_mg()->get_np(num_front+1) ;
    int nz = lim_x.get_mg()->get_nzone() ;
    
   if (shift.get_etat() == ETATZERO) {
	shift.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    shift.set(i).annule_hard() ;
	shift.set_std_base() ;
    }
    
    Tenseur so (source) ;
    
    // La source scalaire :
    Tenseur cop_so (so) ;
    cop_so.dec2_dzpuis() ;
    cop_so.dec2_dzpuis() ;
    
    Tenseur scal (*so.get_mp()) ;
    scal.set_etat_qcq() ;
    
    Cmp source_scal (contract(cop_so.gradient(), 0, 1)()/(lambda+1)) ;
    source_scal.inc2_dzpuis() ;
    if (source_scal.get_etat()== ETATZERO) {
	source_scal.annule_hard() ;
	source_scal.std_base_scal() ;
	}

    Tenseur copie_so (so) ;
    copie_so.dec_dzpuis() ;
    
    Tenseur source_vect (*so.get_mp(), 1, CON, *source.get_triad()) ;
    Tenseur auxi (*so.get_mp(), 1, COV, *source.get_triad()) ;
    Cmp grad_shift (source_scal.get_mp()) ;
    
    // La condition sur la derivee du scalaire :
    Valeur lim_scal (lim_x.get_mg()) ;
    Tenseur shift_old (*shift.get_mp(), 1, CON, shift.get_mp()->get_bvect_cart()) ;
    
    int conte = 0 ;
    int indic = 1 ;
    
    while (indic ==1) {
	
	shift_old = shift ;
	
	grad_shift = contract(shift.gradient(), 0, 1)() ;
	grad_shift.dec2_dzpuis() ;
	grad_shift.va.coef_i() ;
	
	lim_scal = 1 ; // Permet d'affecter les trucs qui vont bien !
	for (int k=0 ; k<np ; k++)
	    for (int j=0 ; j<nt ; j++)
		lim_scal.set(num_front, k, j, 0) = 
		    grad_shift.va (num_front+1, k, j, 0) ;
	lim_scal.std_base_scal() ;
		        
	// On resout la scalaire :
	scal.set() = source_scal.poisson_dirichlet (lim_scal, num_front) ;
	
	// La source vectorielle :
	source_vect.set_etat_qcq() ;
	auxi = scal.gradient() ;
	auxi.inc_dzpuis() ;
	for (int i=0 ; i<3 ; i++)
	    source_vect.set(i) = copie_so(i) - lambda * auxi(i) ;
	
	indic = 0;
	for (int i=0 ; i<3 ; i++)
	    if (source_vect(i).get_etat()==ETATQCQ)
		indic = 1 ;
	if (indic==0) {
	    for (int i=0 ; i<3 ; i++)
		source_vect.set(i).annule_hard() ;
	    source_vect.set_std_base() ;
	}
	
	// On resout les equations de poisson sur le shift :
	shift.set(0) = source_vect(0).poisson_dirichlet (lim_x, num_front) ;
	shift.set(1) = source_vect(1).poisson_dirichlet (lim_y, num_front) ;
	shift.set(2) = source_vect(2).poisson_dirichlet (lim_z, num_front) ;
	
	double erreur = 0 ;
	for (int i=0 ; i<3 ; i++) 
	    if (max(norme(shift(i))) > precision) {
	    Tbl diff (diffrelmax (shift(i), shift_old(i))) ;
	    for (int j=num_front+1 ; j<nz ; j++)
		if (diff(j)> erreur)
		    erreur = diff(j) ;
	    }
	
	cout << "Pas " << conte << " : Difference " << erreur << endl ;
	conte ++ ;
	
	if ((erreur <precision) || (conte > itermax))
	    indic = -1 ;
	}
}


void poisson_vect_binaire ( double lambda, 
		const Tenseur& source_un, const Tenseur& source_deux, 
		const Valeur& bound_x_un, const Valeur& bound_y_un, 
		const Valeur& bound_z_un, const Valeur& bound_x_deux, 
		const Valeur& bound_y_deux, const Valeur& bound_z_deux, 
		Tenseur& sol_un, Tenseur& sol_deux, int num_front, double precision) {
		    
    // METTRE DES ASSERT
    assert (sol_un.get_etat() != ETATNONDEF) ;
    assert (sol_deux.get_etat() != ETATNONDEF) ;
    
    // Les bases des deux vecteurs doivent etre alignees ou non alignees :
     
    assert (sol_un.get_mp() == source_un.get_mp()) ;
    assert (sol_deux.get_mp() == source_deux.get_mp()) ;
    
    double orientation_un = sol_un.get_mp()->get_rot_phi() ;
    assert ((orientation_un==0) || (orientation_un==M_PI)) ;
    
    double orientation_deux = sol_deux.get_mp()->get_rot_phi() ;
    assert ((orientation_deux==0) || (orientation_deux==M_PI)) ;
    
    int same_orient = (orientation_un == orientation_deux) ? 1 : -1 ;
    
    
    if (sol_un.get_etat() == ETATZERO) {
	sol_un.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    sol_un.set(i).annule_hard() ;
	sol_un.set_std_base() ;
    }
    
     if (sol_deux.get_etat() == ETATZERO) {
	sol_deux.set_etat_qcq() ;
	for (int i=0 ; i<3 ; i++)
	    sol_deux.set(i).annule_hard() ;
	sol_deux.set_std_base() ;
    }
    
    Valeur limite_x_un (bound_x_un.get_mg()) ;
    limite_x_un = bound_x_un ;
    Valeur limite_y_un (bound_y_un.get_mg()) ;
    limite_y_un = bound_y_un ;
    Valeur limite_z_un (bound_z_un.get_mg()) ;
    limite_z_un = bound_z_un ;
      
    Valeur limite_x_deux (bound_x_deux.get_mg()) ;
    limite_x_deux = bound_x_deux ;
    Valeur limite_y_deux (bound_y_deux.get_mg()) ;
    limite_y_deux = bound_y_deux ;
    Valeur limite_z_deux (bound_z_deux.get_mg()) ;
    limite_z_deux = bound_z_deux ;
    
    Valeur limite_chi_un (bound_x_un.get_mg()) ;
    limite_chi_un = 0 ;
    limite_chi_un.std_base_scal() ;
    
    Valeur limite_chi_deux (bound_x_deux.get_mg()) ;
    limite_chi_deux = 0 ;
    limite_chi_deux.std_base_scal() ;
    
    Mtbl xa_mtbl_un (source_un.get_mp()->get_mg()) ;
    xa_mtbl_un.set_etat_qcq() ;
    Mtbl ya_mtbl_un (source_un.get_mp()->get_mg()) ;
    ya_mtbl_un.set_etat_qcq() ;
    Mtbl za_mtbl_un (source_un.get_mp()->get_mg()) ;
    za_mtbl_un.set_etat_qcq() ;
    Mtbl xa_mtbl_deux (source_deux.get_mp()->get_mg()) ;
    xa_mtbl_deux.set_etat_qcq() ;
    Mtbl ya_mtbl_deux (source_deux.get_mp()->get_mg()) ;
    ya_mtbl_deux.set_etat_qcq() ;
    Mtbl za_mtbl_deux (source_deux.get_mp()->get_mg()) ;
    za_mtbl_deux.set_etat_qcq() ;
    
    xa_mtbl_un = source_un.get_mp()->xa ;
    ya_mtbl_un = source_un.get_mp()->ya ;
    za_mtbl_un = source_un.get_mp()->za ;
    
    xa_mtbl_deux = source_deux.get_mp()->xa ;
    ya_mtbl_deux = source_deux.get_mp()->ya ;
    za_mtbl_deux = source_deux.get_mp()->za ;
    
    double xabs, yabs, zabs ;
    double air,  theta,  phi ;
    double valeur ;
    
    int nbrep_un = bound_x_un.get_mg()->get_np(num_front) ;
    int nbret_un = bound_x_un.get_mg()->get_nt(num_front) ;
    int nbrep_deux = bound_x_deux.get_mg()->get_np(num_front) ;
    int nbret_deux = bound_x_deux.get_mg()->get_nt(num_front) ;
    int nz_un = bound_x_un.get_mg()->get_nzone() ;
    int nz_deux = bound_x_deux.get_mg()->get_nzone() ;
    
    // La source de l'equation scalaire sur 1
    Tenseur cop_so_un (source_un) ;
    cop_so_un.dec2_dzpuis() ;
    cop_so_un.dec2_dzpuis() ;
    
    Cmp source_scal_un (contract (cop_so_un.gradient(), 0, 1)()/(lambda+1)) ;
    if (source_scal_un.get_etat() == ETATZERO) {
	source_scal_un.annule_hard() ;
	source_scal_un.std_base_scal() ;
    }
    source_scal_un.inc2_dzpuis() ;
    
    // La source de l'equation scalaire sur 2
    Tenseur cop_so_deux (source_deux) ;
    cop_so_deux.dec2_dzpuis() ;
    cop_so_deux.dec2_dzpuis() ;
    
    Cmp source_scal_deux (contract (cop_so_deux.gradient(), 0, 1)()/(lambda+1)) ;
    if (source_scal_deux.get_etat() == ETATZERO) {
	source_scal_deux.annule_hard() ;
	source_scal_deux.std_base_scal() ;
    }
    source_scal_deux.inc2_dzpuis() ;
    
    // Les copies :
    Tenseur copie_so_un (source_un) ;
    copie_so_un.dec_dzpuis() ;
    
    Tenseur copie_so_deux (source_deux) ;
    copie_so_deux.dec_dzpuis() ;
    
    // ON COMMENCE LA BOUCLE :
    Tenseur sol_un_old (sol_un) ;
    Tenseur sol_deux_old (sol_deux) ;
    
    int indic = 1 ;
    int conte = 0 ;
    
    while (indic == 1) {
	
	// On resout les deux equations scalaires :
	Tenseur chi_un (source_scal_un.poisson_dirichlet (limite_chi_un, num_front)) ;
	Tenseur chi_deux (source_scal_deux.poisson_dirichlet (limite_chi_deux, num_front)) ;
	
	// On calcul les source pour les equation vectorielles :
	Tenseur source_vect_un (copie_so_un) ;
	if (source_vect_un.get_etat() == ETATZERO) {
	    source_vect_un.set_etat_qcq() ;
	    for (int i=0 ; i<3 ; i++) {
		source_vect_un.set(i).annule_hard() ;
		source_vect_un.set(i).set_dzpuis(3) ;
		}
	    source_vect_un.set_std_base() ;
	    }
	Tenseur grad_chi_un (chi_un.gradient()) ;
	grad_chi_un.inc_dzpuis() ;
	for (int i=0 ; i<3 ; i++)
	    source_vect_un.set(i) = source_vect_un(i)-lambda*grad_chi_un(i) ;
	
	Tenseur source_vect_deux (copie_so_deux) ;
	if (source_vect_deux.get_etat() == ETATZERO) {
	    source_vect_deux.set_etat_qcq() ;
	    for (int i=0 ; i<3 ; i++) {
		source_vect_deux.set(i).annule_hard() ;
		source_vect_deux.set(i).set_dzpuis(3) ;
		}
	    source_vect_deux.set_std_base() ;
	    }
	Tenseur grad_chi_deux (chi_deux.gradient()) ;
	grad_chi_deux.inc_dzpuis() ;
	for (int i=0 ; i<3 ; i++)
	    source_vect_deux.set(i) = source_vect_deux(i)-lambda*grad_chi_deux(i) ;
	
	
	sol_un_old = sol_un ;
	sol_deux_old = sol_deux ;
	
	// On resout les equation vectorielles :
	sol_un.set(0) = source_vect_un(0).poisson_dirichlet (limite_x_un, num_front) ;
	sol_un.set(1) = source_vect_un(1).poisson_dirichlet (limite_y_un, num_front) ;
	sol_un.set(2) = source_vect_un(2).poisson_dirichlet (limite_z_un, num_front) ;
	sol_deux.set(0) = source_vect_deux(0).poisson_dirichlet (limite_x_deux, num_front) ;
	sol_deux.set(1) = source_vect_deux(1).poisson_dirichlet (limite_y_deux, num_front) ;
	sol_deux.set(2) = source_vect_deux(2).poisson_dirichlet (limite_z_deux, num_front) ;
	
	
	// On modifie les Cl sur chi :
	Cmp div_shift_un (contract(sol_un.gradient(), 0, 1)()) ;
	div_shift_un.dec2_dzpuis() ;
	div_shift_un.va.coef_i() ;
	
	limite_chi_un = 1 ; // Affectation
	for (int k=0 ; k<nbrep_un ; k++)
	    for (int j=0 ; j<nbret_un ; j++)
		limite_chi_un.set(num_front, k, j, 0) =
		    div_shift_un.va (num_front+1, k, j, 0) ;
	limite_chi_un.std_base_scal() ;
	
	Cmp div_shift_deux (contract(sol_deux.gradient(), 0, 1)()) ;
	div_shift_deux.dec2_dzpuis() ;
	div_shift_deux.va.coef_i() ;
	
	limite_chi_deux = 1 ; // Affectation
	for (int k=0 ; k<nbrep_deux ; k++)
	    for (int j=0 ; j<nbret_deux ; j++)
		limite_chi_deux.set(num_front, k, j, 0) =
		    div_shift_deux.va (num_front+1, k, j, 0) ;
	limite_chi_deux.std_base_scal() ;
	
	
	// On modifie les Cl sur sol_un :
	for (int k=0 ; k<nbrep_un ; k++)
	    for (int j=0 ; j<nbret_un ; j++) {
		xabs = xa_mtbl_un (num_front+1, k, j, 0) ;
		yabs = ya_mtbl_un (num_front+1, k, j, 0) ;
		zabs = za_mtbl_un (num_front+1, k, j, 0) ;
		
		source_deux.get_mp()->convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
		
		valeur = sol_deux(0).val_point(air, theta, phi) ;
		limite_x_un.set(num_front, k, j, 0) = 
			bound_x_un(num_front, k, j, 0) - same_orient*valeur ;
		
		valeur = sol_deux(1).val_point(air, theta, phi) ;
		limite_y_un.set(num_front, k, j, 0) = 
			bound_y_un(num_front, k, j, 0) - same_orient*valeur ;
		
		valeur = sol_deux(2).val_point(air, theta, phi) ;
		limite_z_un.set(num_front, k, j, 0) = 
			bound_z_un(num_front, k, j, 0) - valeur ;
	    }
	    
	// On modifie les Cl sur sol_deux :
	for (int k=0 ; k<nbrep_deux ; k++)
	    for (int j=0 ; j<nbret_deux ; j++) {
		xabs = xa_mtbl_deux (num_front+1, k, j, 0) ;
		yabs = ya_mtbl_deux (num_front+1, k, j, 0) ;
		zabs = za_mtbl_deux (num_front+1, k, j, 0) ;
		
		source_un.get_mp()->convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
		
		valeur = sol_un(0).val_point(air, theta, phi) ;
		limite_x_deux.set(num_front, k, j, 0) = 
			bound_x_deux(num_front, k, j, 0) - same_orient*valeur ;
		
		valeur = sol_un(1).val_point(air, theta, phi) ;
		limite_y_deux.set(num_front, k, j, 0) = 
			bound_y_deux(num_front, k, j, 0) - same_orient*valeur ;
		
		valeur = sol_un(2).val_point(air, theta, phi) ;
		limite_z_deux.set(num_front, k, j, 0) = 
			bound_z_deux(num_front, k, j, 0) - valeur ;
	    }
	    
	double erreur = 0 ;
	
	for (int i=0 ; i<3 ; i++) {
	    Tbl diff_un (diffrelmax (sol_un_old(i), sol_un(i))) ;
	    for (int j=num_front+1 ; j<nz_un ; j++)
		if (erreur<diff_un(j))
		    erreur = diff_un(j) ;
	}
	
	for (int i=0 ; i<3 ; i++) {
	    Tbl diff_deux (diffrelmax (sol_deux_old(i), sol_deux(i))) ;
	    for (int j=num_front+1 ; j<nz_deux ; j++)
		if (erreur<diff_deux(j))
		    erreur = diff_deux(j) ;
	}
	
	cout << "Pas " << conte << " : Difference " << erreur << endl ;
	
	if (erreur < precision)
	    indic = -1 ;
	conte ++ ;
	}
}
