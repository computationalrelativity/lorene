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


char pde_frontiere_bin_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.2  2003/10/03 15:58:50  j_novak
 * Cleaning of some headers
 *
 * Revision 1.1.1.1  2001/11/20 15:19:28  e_gourgoulhon
 * LORENE
 *
 * Revision 2.3  2000/12/04  14:30:16  phil
 * correction CL.
 *
 * Revision 2.2  2000/12/01  15:18:26  phil
 * vire trucs inutiles
 *
 * Revision 2.1  2000/12/01  15:16:49  phil
 * correction version Neumann
 *
 * Revision 2.0  2000/10/19  09:36:07  phil
 * *** empty log message ***
 *
 *
 * $Header$
 *
 */

//standard
#include <stdlib.h>
#include <math.h>

// LORENE
#include "tenseur.h"
#include "proto.h"

// Version avec une fonction de theta, phi.

void dirichlet_binaire (const Cmp& source_un, const Cmp& source_deux, 
			const Valeur& boundary_un, const Valeur& boundary_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) {
    
    // Les verifs sur le mapping :
    assert (source_un.get_mp() == sol_un.get_mp()) ;
    assert (source_deux.get_mp() == sol_deux.get_mp()) ;
    
    Valeur limite_un (boundary_un.get_mg()) ;
    Valeur limite_deux (boundary_deux.get_mg()) ;
    
    Cmp sol_un_old (sol_un.get_mp()) ;
    Cmp sol_deux_old (sol_deux.get_mp()) ;
    
    Mtbl xa_mtbl_un (source_un.get_mp()->xa) ;
    Mtbl ya_mtbl_un (source_un.get_mp()->ya) ;
    Mtbl za_mtbl_un (source_un.get_mp()->za) ;
    Mtbl xa_mtbl_deux (source_deux.get_mp()->xa) ;
    Mtbl ya_mtbl_deux (source_deux.get_mp()->ya) ;
    Mtbl za_mtbl_deux (source_deux.get_mp()->za) ;
    
    double xabs, yabs, zabs ;
    double air,  theta,  phi ;
    double valeur ;
    
    int nbrep_un = boundary_un.get_mg()->get_np(num_front) ;
    int nbret_un = boundary_un.get_mg()->get_nt(num_front) ;
    int nbrep_deux = boundary_deux.get_mg()->get_np(num_front) ;
    int nbret_deux = boundary_deux.get_mg()->get_nt(num_front) ;
    
    int nz_un = boundary_un.get_mg()->get_nzone() ;
    int nz_deux = boundary_deux.get_mg()->get_nzone() ;
    
    // Initialisation valeur limite avant iteration !
   limite_un = 1 ; //Pour initialiser les tableaux
   for (int k=0 ; k<nbrep_un ; k++)
    for (int j=0 ; j<nbret_un ; j++)
	limite_un.set(num_front, k, j, 0) =
		sol_un.va.val_point_jk(num_front+1, -1, j, k) ;
    limite_un.set_base (boundary_un.base) ;

    limite_deux = 1 ;
    for (int k=0 ; k<nbrep_deux ; k++)
	for (int j=0 ; j<nbret_deux ; j++)
	  limite_deux.set(num_front, k, j, 0) =
	    sol_deux.va.val_point_jk(num_front+1, -1, j, k) ;
    limite_deux.set_base (boundary_deux.base) ;


    int conte = 0 ;
    int indic = 1 ;
    
    while (indic==1) {
	
	sol_un_old = sol_un ;
	sol_deux_old = sol_deux ;
	
	sol_un = source_un.poisson_dirichlet(limite_un, num_front) ;
	sol_deux = source_deux.poisson_dirichlet(limite_deux, num_front) ;
	
	xa_mtbl_deux = source_deux.get_mp()->xa ;
	ya_mtbl_deux = source_deux.get_mp()->ya ;
	za_mtbl_deux = source_deux.get_mp()->za ;
	
	
	for (int k=0 ; k<nbrep_deux ; k++)
	    for (int j=0 ; j<nbret_deux ; j++) {
		xabs = xa_mtbl_deux (num_front+1, k, j, 0) ;
		yabs = ya_mtbl_deux (num_front+1, k, j, 0) ;
		zabs = za_mtbl_deux (num_front+1, k, j, 0) ;
		
		source_un.get_mp()->convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
		valeur = sol_un.val_point(air, theta, phi) ;
		
		limite_deux.set(num_front, k, j, 0) = 
			boundary_deux(num_front, k, j, 0) - valeur ;
	    }
	   
	xa_mtbl_un = source_un.get_mp()->xa ;
	ya_mtbl_un = source_un.get_mp()->ya ;
	za_mtbl_un = source_un.get_mp()->za ;
	
	for (int k=0 ; k<nbrep_un ; k++)
	    for (int j=0 ; j<nbret_un ; j++) {
		xabs = xa_mtbl_un (num_front+1, k, j, 0) ;
		yabs = ya_mtbl_un (num_front+1, k, j, 0) ;
		zabs = za_mtbl_un (num_front+1, k, j, 0) ;
		
		source_deux.get_mp()->convert_absolute 
		    (xabs, yabs, zabs, air, theta, phi) ;
		valeur = sol_deux.val_point(air, theta, phi) ;
		
		limite_un.set(num_front, k, j, 0) = 
		    boundary_un(num_front, k, j, 0) - valeur ;
	    }
	
	double erreur = 0 ;
	Tbl diff_un (diffrelmax(sol_un, sol_un_old)) ;
	for (int i=num_front+1 ; i<nz_un ; i++)
	    if (diff_un(i) > erreur)
		erreur = diff_un(i) ;
	
	Tbl diff_deux (diffrelmax(sol_deux, sol_deux_old)) ;
	for (int i=num_front+1 ; i<nz_deux ; i++)
	    if (diff_deux(i) > erreur)
		erreur = diff_deux(i) ;
	
	cout << "Pas " << conte << " : Difference " << erreur << endl ;
	
	conte ++ ;
	if (erreur < precision)
	    indic = -1 ;
    }
					
}

// Version avec des doubles :
void dirichlet_binaire (const Cmp& source_un, const Cmp& source_deux, 
			double bound_un, double bound_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) {
				
    Valeur boundary_un (source_un.get_mp()->get_mg()->get_angu()) ;
    if (bound_un == 0)
	boundary_un.annule_hard() ;
    else
	boundary_un = bound_un ;
    boundary_un.std_base_scal() ;
    
    Valeur boundary_deux (source_deux.get_mp()->get_mg()->get_angu()) ;
    if (bound_deux == 0)
	boundary_deux.annule_hard() ;
    else
	boundary_deux = bound_deux ;
    boundary_deux.std_base_scal() ;
    
    dirichlet_binaire (source_un, source_deux, boundary_un, boundary_deux, 
			sol_un, sol_deux, num_front, precision) ;
}

// Version avec une fonction de theta, phi.

void neumann_binaire (const Cmp& source_un, const Cmp& source_deux, 
			const Valeur& boundary_un, const Valeur& boundary_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) {
    
    // Les verifs sur le mapping :
    assert (source_un.get_mp() == sol_un.get_mp()) ;
    assert (source_deux.get_mp() == sol_deux.get_mp()) ;
    
    // Alignes ou non ?
    double orient_un = source_un.get_mp()->get_rot_phi() ;
    assert ((orient_un==0) || (orient_un==M_PI)) ;
    double orient_deux = source_deux.get_mp()->get_rot_phi() ;
    assert ((orient_deux==0) || (orient_deux==M_PI)) ;
    int same_orient = (orient_un==orient_deux) ? 1 : -1 ;
    
    Valeur limite_un (boundary_un.get_mg()) ;
    Valeur limite_deux (boundary_deux.get_mg()) ;
    
    Cmp sol_un_old (sol_un.get_mp()) ;
    Cmp sol_deux_old (sol_deux.get_mp()) ;
    
    Mtbl xa_mtbl_un (source_un.get_mp()->xa) ;
    Mtbl ya_mtbl_un (source_un.get_mp()->ya) ;
    Mtbl za_mtbl_un (source_un.get_mp()->za) ;
    
    Mtbl cost_mtbl_un (source_un.get_mp()->cost) ;
    Mtbl sint_mtbl_un (source_un.get_mp()->sint) ;
    Mtbl cosp_mtbl_un (source_un.get_mp()->cosp) ;
    Mtbl sinp_mtbl_un (source_un.get_mp()->sinp) ;
    
    
    Mtbl xa_mtbl_deux (source_deux.get_mp()->xa) ;
    Mtbl ya_mtbl_deux (source_deux.get_mp()->ya) ;
    Mtbl za_mtbl_deux (source_deux.get_mp()->za) ;
    
    Mtbl cost_mtbl_deux (source_deux.get_mp()->cost) ;
    Mtbl sint_mtbl_deux (source_deux.get_mp()->sint) ;
    Mtbl cosp_mtbl_deux (source_deux.get_mp()->cosp) ;
    Mtbl sinp_mtbl_deux (source_deux.get_mp()->sinp) ;
    
    double xabs, yabs, zabs ;
    double air,  theta,  phi ;
    double valeur ;
     
    int nbrep_un = boundary_un.get_mg()->get_np(num_front) ;
    int nbret_un = boundary_un.get_mg()->get_nt(num_front) ;
    int nbrep_deux = boundary_deux.get_mg()->get_np(num_front) ;
    int nbret_deux = boundary_deux.get_mg()->get_nt(num_front) ;
    
    int nz_un = boundary_un.get_mg()->get_nzone() ;
    int nz_deux = boundary_deux.get_mg()->get_nzone() ;
    
    // Initialisation des CL :
    limite_un = 1 ;
    limite_deux = 2 ;
    Cmp der_un (sol_un.dsdr()) ;
    Cmp der_deux (sol_deux.dsdr()) ;
    
    for (int k=0 ; k<nbrep_un ; k++)
	for (int j=0 ; j<nbret_un ; j++)
	    limite_un.set(num_front, k, j, 0) =
		der_un.va.val_point_jk(num_front+1, -1, j, k) ;
    limite_un.set_base (boundary_un.base) ;

    for (int k=0 ; k<nbrep_deux ; k++)
	for (int j=0 ; j<nbret_deux ; j++)
	  limite_deux.set(num_front, k, j, 0) =
	    der_deux.va.val_point_jk(num_front+1, -1, j, k) ;
    limite_deux.set_base (boundary_deux.base) ;
    
    int conte = 0 ;
    int indic = 1 ;
    
    while (indic==1) {

	sol_un_old = sol_un ;
	sol_deux_old = sol_deux ;
	
	sol_un = source_un.poisson_neumann(limite_un, num_front) ;
	sol_deux = source_deux.poisson_neumann(limite_deux, num_front) ;
	
	// On veut les derivees de l'un sur l'autre :
	Tenseur copie_un (sol_un) ;
	Tenseur grad_sol_un (copie_un.gradient()) ;
	grad_sol_un.dec2_dzpuis() ;
	grad_sol_un.set(0) = grad_sol_un(0)*same_orient ;
	grad_sol_un.set(1) = grad_sol_un(1)*same_orient ;
	
	for (int k=0 ; k<nbrep_deux ; k++)
	    for (int j=0 ; j<nbret_deux ; j++) {
		xabs = xa_mtbl_deux (num_front+1, k, j, 0) ;
		yabs = ya_mtbl_deux (num_front+1, k, j, 0) ;
		zabs = za_mtbl_deux (num_front+1, k, j, 0) ;
		
		source_un.get_mp()->convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
				
		valeur = sint_mtbl_deux (num_front+1, k, j, 0) * (
cosp_mtbl_deux(num_front+1, k, j, 0)*grad_sol_un(0).val_point(air, theta, phi)+
sinp_mtbl_deux(num_front+1, k, j, 0)*grad_sol_un(1).val_point(air, theta, phi))+
cost_mtbl_deux(num_front+1, k, j, 0)*grad_sol_un(2).val_point(air, theta, phi);

		limite_deux.set(num_front, k, j, 0) = 
			boundary_deux(num_front, k, j, 0) - valeur ;
	    }
	
	Tenseur copie_deux (sol_deux) ;
	Tenseur grad_sol_deux (copie_deux.gradient()) ;
	grad_sol_deux.dec2_dzpuis() ;
	grad_sol_deux.set(0) = grad_sol_deux(0)*same_orient ;
	grad_sol_deux.set(1) = grad_sol_deux(1)*same_orient ;
	
	for (int k=0 ; k<nbrep_un ; k++)
	    for (int j=0 ; j<nbret_un ; j++) {
		xabs = xa_mtbl_un (num_front+1, k, j, 0) ;
		yabs = ya_mtbl_un (num_front+1, k, j, 0) ;
		zabs = za_mtbl_un (num_front+1, k, j, 0) ;
		
		source_deux.get_mp()->convert_absolute 
				(xabs, yabs, zabs, air, theta, phi) ;
				
		valeur = sint_mtbl_un (num_front+1, k, j, 0) * (
cosp_mtbl_un(num_front+1, k, j, 0)*grad_sol_deux(0).val_point(air, theta, phi)+
sinp_mtbl_un(num_front+1, k, j, 0)*grad_sol_deux(1).val_point(air, theta, phi))+
cost_mtbl_un(num_front+1, k, j, 0)*grad_sol_deux(2).val_point(air, theta, phi);

		limite_un.set(num_front, k, j, 0) = 
			boundary_un(num_front, k, j, 0) - valeur ;
	    }
		
	double erreur = 0 ;
	Tbl diff_un (diffrelmax(sol_un, sol_un_old)) ;
	for (int i=num_front+1 ; i<nz_un ; i++)
	    if (diff_un(i) > erreur)
		erreur = diff_un(i) ;
	
	Tbl diff_deux (diffrelmax(sol_deux, sol_deux_old)) ;
	for (int i=num_front+1 ; i<nz_deux ; i++)
	    if (diff_deux(i) > erreur)
		erreur = diff_deux(i) ;
	
	cout << "Pas " << conte << " : Difference " << erreur << endl ;
	conte ++ ;
	
	if (erreur < precision)
	    indic = -1 ;
    }					
}

// Version avec des doubles :
void neumann_binaire (const Cmp& source_un, const Cmp& source_deux, 
			double bound_un, double bound_deux, 
				Cmp& sol_un, Cmp& sol_deux, int num_front, 
				double precision) {
				
    Valeur boundary_un (source_un.get_mp()->get_mg()->get_angu()) ;
    if (bound_un == 0)
	boundary_un.annule_hard () ;
    else 
	boundary_un = bound_un ;
    boundary_un.std_base_scal() ;
    
    Valeur boundary_deux (source_deux.get_mp()->get_mg()->get_angu()) ;
    if (bound_deux == 0)
	boundary_deux.annule_hard() ;
    else
	boundary_deux = bound_deux ;
    boundary_deux.std_base_scal() ;
    
    neumann_binaire (source_un, source_deux, boundary_un, boundary_deux, 
			sol_un, sol_deux, num_front, precision) ;
}			
